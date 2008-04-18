      SUBROUTINE thermcell(ngrid,nlay,ptimestep
     s                  ,pplay,pplev,pphi
     s                  ,pu,pv,pt,po
     s                  ,pduadj,pdvadj,pdtadj,pdoadj
     s                  ,fm0,entr0
c    s                  ,pu_therm,pv_therm
     s                  ,r_aspect,l_mix,w2di,tho)

      use dimens_m
      use dimphy
      use YOMCST
      IMPLICIT NONE

c=======================================================================
c
c   Calcul du transport verticale dans la couche limite en presence
c   de "thermiques" explicitement representes
c
c   Réécriture à partir d'un listing papier à Habas, le 14/02/00
c
c   le thermique est supposé homogène et dissipé par mélange avec
c   son environnement. la longueur l_mix contrôle l'efficacité du
c   mélange
c
c   Le calcul du transport des différentes espèces se fait en prenant
c   en compte:
c     1. un flux de masse montant
c     2. un flux de masse descendant
c     3. un entrainement
c     4. un detrainement
c
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------


c   arguments:
c   ----------

      INTEGER ngrid,nlay,w2di,tho
      real ptimestep,l_mix,r_aspect
      REAL pt(ngrid,nlay),pdtadj(ngrid,nlay)
      REAL pu(ngrid,nlay),pduadj(ngrid,nlay)
      REAL pv(ngrid,nlay),pdvadj(ngrid,nlay)
      REAL po(ngrid,nlay),pdoadj(ngrid,nlay)
      REAL, intent(in):: pplay(ngrid,nlay)
      real, intent(in):: pplev(ngrid,nlay+1)
      real pphi(ngrid,nlay)

      integer idetr
      save idetr
      data idetr/3/

c   local:
c   ------

      INTEGER ig,k,l,lmaxa(klon),lmix(klon)
      real zsortie1d(klon)
c CR: on remplace lmax(klon,klev+1)
      INTEGER lmax(klon),lmin(klon),lentr(klon)
      real linter(klon)
      real zmix(klon), fracazmix(klon) 
c RC 
      real zmax(klon),zw,zz,zw2(klon,klev+1),ztva(klon,klev),zzz

      real zlev(klon,klev+1),zlay(klon,klev)
      REAL zh(klon,klev),zdhadj(klon,klev)
      REAL ztv(klon,klev)
      real zu(klon,klev),zv(klon,klev),zo(klon,klev)
      REAL wh(klon,klev+1)
      real wu(klon,klev+1),wv(klon,klev+1),wo(klon,klev+1)
      real zla(klon,klev+1)
      real zwa(klon,klev+1)
      real zld(klon,klev+1)
      real zwd(klon,klev+1)
      real zsortie(klon,klev)
      real zva(klon,klev)
      real zua(klon,klev)
      real zoa(klon,klev)

      real zha(klon,klev)
      real wa_moy(klon,klev+1)
      real fraca(klon,klev+1)
      real fracc(klon,klev+1)
      real zf,zf2
      real thetath2(klon,klev),wth2(klon,klev)
      common/comtherm/thetath2,wth2

      real count_time
      integer isplit,nsplit,ialt
      parameter (nsplit=10)
      data isplit/0/
      save isplit

      logical sorties
      real rho(klon,klev),rhobarz(klon,klev+1),masse(klon,klev)
      real zpspsk(klon,klev)

c     real wmax(klon,klev),wmaxa(klon)
      real wmax(klon),wmaxa(klon)
      real wa(klon,klev,klev+1)
      real wd(klon,klev+1)
      real larg_part(klon,klev,klev+1)
      real fracd(klon,klev+1)
      real xxx(klon,klev+1)
      real larg_cons(klon,klev+1)
      real larg_detr(klon,klev+1)
      real fm0(klon,klev+1),entr0(klon,klev),detr(klon,klev)
      real pu_therm(klon,klev),pv_therm(klon,klev)
      real fm(klon,klev+1),entr(klon,klev)
      real fmc(klon,klev+1)

cCR:nouvelles variables
      real f_star(klon,klev+1),entr_star(klon,klev)
      real entr_star_tot(klon),entr_star2(klon)
      real f(klon), f0(klon)
      real zlevinter(klon)
      logical first
      data first /.false./
      save first
cRC

      character*2 str2
      character*10 str10

      LOGICAL vtest(klon),down

      EXTERNAL SCOPY

      integer ncorrec,ll
      save ncorrec
      data ncorrec/0/
      
c
c-----------------------------------------------------------------------
c   initialisation:
c   ---------------
c
       sorties=.true.
      IF(ngrid.NE.klon) THEN
         PRINT*
         PRINT*,'STOP dans convadj'
         PRINT*,'ngrid    =',ngrid
         PRINT*,'klon  =',klon
      ENDIF
c
c-----------------------------------------------------------------------
c   incrementation eventuelle de tendances precedentes:
c   ---------------------------------------------------

       print*,'0 OK convect8'

      DO 1010 l=1,nlay
         DO 1015 ig=1,ngrid
            zpspsk(ig,l)=(pplay(ig,l)/pplev(ig,1))**RKAPPA
            zh(ig,l)=pt(ig,l)/zpspsk(ig,l)
            zu(ig,l)=pu(ig,l)
            zv(ig,l)=pv(ig,l)
            zo(ig,l)=po(ig,l)
            ztv(ig,l)=zh(ig,l)*(1.+0.61*zo(ig,l))
1015     CONTINUE
1010  CONTINUE

       print*,'1 OK convect8'
c                       --------------------
c
c
c                       + + + + + + + + + + +
c
c
c  wa, fraca, wd, fracd --------------------   zlev(2), rhobarz
c  wh,wt,wo ...
c
c                       + + + + + + + + + + +  zh,zu,zv,zo,rho
c
c
c                       --------------------   zlev(1)
c                       \\\\\\\\\\\\\\\\\\\\
c
c

c-----------------------------------------------------------------------
c   Calcul des altitudes des couches
c-----------------------------------------------------------------------

      do l=2,nlay
         do ig=1,ngrid
            zlev(ig,l)=0.5*(pphi(ig,l)+pphi(ig,l-1))/RG
         enddo
      enddo
      do ig=1,ngrid
         zlev(ig,1)=0.
         zlev(ig,nlay+1)=(2.*pphi(ig,klev)-pphi(ig,klev-1))/RG
      enddo
      do l=1,nlay
         do ig=1,ngrid
            zlay(ig,l)=pphi(ig,l)/RG
         enddo
      enddo

c      print*,'2 OK convect8'
c-----------------------------------------------------------------------
c   Calcul des densites
c-----------------------------------------------------------------------

      do l=1,nlay
         do ig=1,ngrid
            rho(ig,l)=pplay(ig,l)/(zpspsk(ig,l)*RD*zh(ig,l))
         enddo
      enddo

      do l=2,nlay
         do ig=1,ngrid
            rhobarz(ig,l)=0.5*(rho(ig,l)+rho(ig,l-1))
         enddo
      enddo

      do k=1,nlay
         do l=1,nlay+1
            do ig=1,ngrid
               wa(ig,k,l)=0.
            enddo
         enddo
      enddo

c      print*,'3 OK convect8'
c------------------------------------------------------------------
c   Calcul de w2, quarre de w a partir de la cape
c   a partir de w2, on calcule wa, vitesse de l'ascendance
c
c   ATTENTION: Dans cette version, pour cause d'economie de memoire,
c   w2 est stoke dans wa
c
c   ATTENTION: dans convect8, on n'utilise le calcule des wa
c   independants par couches que pour calculer l'entrainement
c   a la base et la hauteur max de l'ascendance.
c
c   Indicages:
c   l'ascendance provenant du niveau k traverse l'interface l avec
c   une vitesse wa(k,l).
c
c                       --------------------
c
c                       + + + + + + + + + + 
c
c  wa(k,l)   ----       --------------------    l
c             /\
c            /||\       + + + + + + + + + + 
c             ||
c             ||        --------------------
c             ||
c             ||        + + + + + + + + + + 
c             ||
c             ||        --------------------
c             ||__
c             |___      + + + + + + + + + +     k
c
c                       --------------------
c
c
c
c------------------------------------------------------------------

cCR: ponderation entrainement des couches instables
cdef des entr_star tels que entr=f*entr_star      
      do l=1,klev
         do ig=1,ngrid 
            entr_star(ig,l)=0.
         enddo
      enddo
c determination de la longueur de la couche d entrainement
      do ig=1,ngrid
         lentr(ig)=1
      enddo

con ne considere que les premieres couches instables
      do k=nlay-2,1,-1
         do ig=1,ngrid
            if (ztv(ig,k).gt.ztv(ig,k+1).and.
     s          ztv(ig,k+1).le.ztv(ig,k+2)) then
               lentr(ig)=k
            endif
          enddo
      enddo
    
c determination du lmin: couche d ou provient le thermique
      do ig=1,ngrid
         lmin(ig)=1
      enddo
      do ig=1,ngrid
         do l=nlay,2,-1
            if (ztv(ig,l-1).gt.ztv(ig,l)) then
               lmin(ig)=l-1
            endif
         enddo
      enddo
c
c definition de l'entrainement des couches
      do l=1,klev-1
         do ig=1,ngrid 
            if (ztv(ig,l).gt.ztv(ig,l+1).and.
     s          l.ge.lmin(ig).and.l.le.lentr(ig)) then 
                 entr_star(ig,l)=(ztv(ig,l)-ztv(ig,l+1))*
     s                           (zlev(ig,l+1)-zlev(ig,l))
            endif
         enddo
      enddo
c pas de thermique si couches 1->5 stables
      do ig=1,ngrid
         if (lmin(ig).gt.5) then
            do l=1,klev
               entr_star(ig,l)=0.
            enddo
         endif
      enddo 
c calcul de l entrainement total
      do ig=1,ngrid
         entr_star_tot(ig)=0.
      enddo
      do ig=1,ngrid
         do k=1,klev
            entr_star_tot(ig)=entr_star_tot(ig)+entr_star(ig,k)
         enddo
      enddo
c
      print*,'fin calcul entr_star'
      do k=1,klev
         do ig=1,ngrid 
            ztva(ig,k)=ztv(ig,k)
         enddo
      enddo
cRC
c      print*,'7 OK convect8'
      do k=1,klev+1
         do ig=1,ngrid
            zw2(ig,k)=0.
            fmc(ig,k)=0.
cCR
            f_star(ig,k)=0.
cRC
            larg_cons(ig,k)=0.
            larg_detr(ig,k)=0.
            wa_moy(ig,k)=0.
         enddo
      enddo

c      print*,'8 OK convect8'
      do ig=1,ngrid
         linter(ig)=1.
         lmaxa(ig)=1
         lmix(ig)=1
         wmaxa(ig)=0.
      enddo

cCR: 
      do l=1,nlay-2
         do ig=1,ngrid
            if (ztv(ig,l).gt.ztv(ig,l+1)
     s         .and.entr_star(ig,l).gt.1.e-10
     s         .and.zw2(ig,l).lt.1e-10) then
               f_star(ig,l+1)=entr_star(ig,l)
ctest:calcul de dteta
               zw2(ig,l+1)=2.*RG*(ztv(ig,l)-ztv(ig,l+1))/ztv(ig,l+1)
     s                     *(zlev(ig,l+1)-zlev(ig,l))
     s                     *0.4*pphi(ig,l)/(pphi(ig,l+1)-pphi(ig,l))
               larg_detr(ig,l)=0.
            else if ((zw2(ig,l).ge.1e-10).and.
     s               (f_star(ig,l)+entr_star(ig,l).gt.1.e-10)) then
               f_star(ig,l+1)=f_star(ig,l)+entr_star(ig,l)
               ztva(ig,l)=(f_star(ig,l)*ztva(ig,l-1)+entr_star(ig,l)
     s                    *ztv(ig,l))/f_star(ig,l+1)
               zw2(ig,l+1)=zw2(ig,l)*(f_star(ig,l)/f_star(ig,l+1))**2+
     s                     2.*RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)
     s                     *(zlev(ig,l+1)-zlev(ig,l))
            endif
c determination de zmax continu par interpolation lineaire
            if (zw2(ig,l+1).lt.0.) then
ctest
               if (abs(zw2(ig,l+1)-zw2(ig,l)).lt.1e-10) then
                  print*,'pb linter'
               endif
               linter(ig)=(l*(zw2(ig,l+1)-zw2(ig,l))
     s           -zw2(ig,l))/(zw2(ig,l+1)-zw2(ig,l))
               zw2(ig,l+1)=0.
               lmaxa(ig)=l
            else
               if (zw2(ig,l+1).lt.0.) then
                  print*,'pb1 zw2<0'
               endif
               wa_moy(ig,l+1)=sqrt(zw2(ig,l+1))
            endif
            if (wa_moy(ig,l+1).gt.wmaxa(ig)) then
c   lmix est le niveau de la couche ou w (wa_moy) est maximum
               lmix(ig)=l+1
               wmaxa(ig)=wa_moy(ig,l+1)
            endif
         enddo
      enddo
      print*,'fin calcul zw2'
c
c Calcul de la couche correspondant a la hauteur du thermique
      do ig=1,ngrid
         lmax(ig)=lentr(ig)
      enddo
      do ig=1,ngrid
         do l=nlay,lentr(ig)+1,-1
            if (zw2(ig,l).le.1.e-10) then
               lmax(ig)=l-1
            endif
         enddo
      enddo
c pas de thermique si couches 1->5 stables
      do ig=1,ngrid
         if (lmin(ig).gt.5) then
            lmax(ig)=1
            lmin(ig)=1
         endif
      enddo 
c    
c Determination de zw2 max
      do ig=1,ngrid
         wmax(ig)=0.
      enddo

      do l=1,nlay
         do ig=1,ngrid
            if (l.le.lmax(ig)) then
                if (zw2(ig,l).lt.0.)then
                  print*,'pb2 zw2<0'
                endif
                zw2(ig,l)=sqrt(zw2(ig,l))
                wmax(ig)=max(wmax(ig),zw2(ig,l))
            else
                 zw2(ig,l)=0.
            endif
          enddo
      enddo

c   Longueur caracteristique correspondant a la hauteur des thermiques.
      do  ig=1,ngrid
         zmax(ig)=0.
         zlevinter(ig)=zlev(ig,1)
      enddo
      do  ig=1,ngrid
c calcul de zlevinter
          zlevinter(ig)=(zlev(ig,lmax(ig)+1)-zlev(ig,lmax(ig)))*
     s    linter(ig)+zlev(ig,lmax(ig))-lmax(ig)*(zlev(ig,lmax(ig)+1)
     s    -zlev(ig,lmax(ig)))
       zmax(ig)=max(zmax(ig),zlevinter(ig)-zlev(ig,lmin(ig)))
      enddo

      print*,'avant fermeture'
c Fermeture,determination de f
      do ig=1,ngrid
         entr_star2(ig)=0.
      enddo
      do ig=1,ngrid
         if (entr_star_tot(ig).LT.1.e-10) then
            f(ig)=0.
         else
             do k=lmin(ig),lentr(ig)
                entr_star2(ig)=entr_star2(ig)+entr_star(ig,k)**2
     s                    /(rho(ig,k)*(zlev(ig,k+1)-zlev(ig,k)))
             enddo
c Nouvelle fermeture
             f(ig)=wmax(ig)/(max(500.,zmax(ig))*r_aspect
     s             *entr_star2(ig))*entr_star_tot(ig)
ctest
c             if (first) then
c             f(ig)=f(ig)+(f0(ig)-f(ig))*exp(-ptimestep/zmax(ig)
c     s             *wmax(ig))
c             endif
         endif
c         f0(ig)=f(ig)
c         first=.true.
      enddo
      print*,'apres fermeture'

c Calcul de l'entrainement
       do k=1,klev
         do ig=1,ngrid 
            entr(ig,k)=f(ig)*entr_star(ig,k)
         enddo
      enddo
c Calcul des flux
      do ig=1,ngrid
         do l=1,lmax(ig)-1
            fmc(ig,l+1)=fmc(ig,l)+entr(ig,l)
         enddo
      enddo

cRC


c      print*,'9 OK convect8'
c     print*,'WA1 ',wa_moy

c   determination de l'indice du debut de la mixed layer ou w decroit

c   calcul de la largeur de chaque ascendance dans le cas conservatif.
c   dans ce cas simple, on suppose que la largeur de l'ascendance provenant
c   d'une couche est égale à la hauteur de la couche alimentante.
c   La vitesse maximale dans l'ascendance est aussi prise comme estimation
c   de la vitesse d'entrainement horizontal dans la couche alimentante.

      do l=2,nlay
         do ig=1,ngrid
            if (l.le.lmaxa(ig)) then
               zw=max(wa_moy(ig,l),1.e-10)
               larg_cons(ig,l)=zmax(ig)*r_aspect
     s         *fmc(ig,l)/(rhobarz(ig,l)*zw)
            endif
         enddo
      enddo

      do l=2,nlay
         do ig=1,ngrid
            if (l.le.lmaxa(ig)) then
c              if (idetr.eq.0) then
c  cette option est finalement en dur.
                  if ((l_mix*zlev(ig,l)).lt.0.)then
                   print*,'pb l_mix*zlev<0'
                  endif
                  larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
c              else if (idetr.eq.1) then
c                 larg_detr(ig,l)=larg_cons(ig,l)
c    s            *sqrt(l_mix*zlev(ig,l))/larg_cons(ig,lmix(ig))
c              else if (idetr.eq.2) then
c                 larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
c    s            *sqrt(wa_moy(ig,l))
c              else if (idetr.eq.4) then
c                 larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
c    s            *wa_moy(ig,l)
c              endif
            endif
         enddo
       enddo

c      print*,'10 OK convect8'
c     print*,'WA2 ',wa_moy
c   calcul de la fraction de la maille concernée par l'ascendance en tenant
c   compte de l'epluchage du thermique.
c
cCR def de  zmix continu (profil parabolique des vitesses)
      do ig=1,ngrid
           if (lmix(ig).gt.1.) then
c test 
              if (((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))
     s        *((zlev(ig,lmix(ig)))-(zlev(ig,lmix(ig)+1)))
     s        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))
     s        *((zlev(ig,lmix(ig)-1))-(zlev(ig,lmix(ig))))).gt.1e-10)
     s        then
c             
            zmix(ig)=((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))
     s        *((zlev(ig,lmix(ig)))**2-(zlev(ig,lmix(ig)+1))**2)
     s        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))
     s        *((zlev(ig,lmix(ig)-1))**2-(zlev(ig,lmix(ig)))**2))
     s        /(2.*((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))
     s        *((zlev(ig,lmix(ig)))-(zlev(ig,lmix(ig)+1)))
     s        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))
     s        *((zlev(ig,lmix(ig)-1))-(zlev(ig,lmix(ig))))))
            else
            zmix(ig)=zlev(ig,lmix(ig))
            print*,'pb zmix'
            endif
         else 
         zmix(ig)=0.
         endif
ctest
         if ((zmax(ig)-zmix(ig)).lt.0.) then
            zmix(ig)=0.99*zmax(ig)
c            print*,'pb zmix>zmax'
         endif
      enddo
c
c calcul du nouveau lmix correspondant
      do ig=1,ngrid
         do l=1,klev
            if (zmix(ig).ge.zlev(ig,l).and.
     s          zmix(ig).lt.zlev(ig,l+1)) then
              lmix(ig)=l
             endif
          enddo
      enddo
c
      do l=2,nlay
         do ig=1,ngrid
            if(larg_cons(ig,l).gt.1.) then
c     print*,ig,l,lmix(ig),lmaxa(ig),larg_cons(ig,l),'  KKK'
               fraca(ig,l)=(larg_cons(ig,l)-larg_detr(ig,l))
     s            /(r_aspect*zmax(ig))
c test
               fraca(ig,l)=max(fraca(ig,l),0.)
               fraca(ig,l)=min(fraca(ig,l),0.5)
               fracd(ig,l)=1.-fraca(ig,l)
               fracc(ig,l)=larg_cons(ig,l)/(r_aspect*zmax(ig))
            else
c              wa_moy(ig,l)=0.
               fraca(ig,l)=0.
               fracc(ig,l)=0.
               fracd(ig,l)=1.
            endif
         enddo
      enddo                  
cCR: calcul de fracazmix
       do ig=1,ngrid
          fracazmix(ig)=(fraca(ig,lmix(ig)+1)-fraca(ig,lmix(ig)))/
     s     (zlev(ig,lmix(ig)+1)-zlev(ig,lmix(ig)))*zmix(ig)
     s    +fraca(ig,lmix(ig))-zlev(ig,lmix(ig))*(fraca(ig,lmix(ig)+1)
     s    -fraca(ig,lmix(ig)))/(zlev(ig,lmix(ig)+1)-zlev(ig,lmix(ig)))
       enddo
c
       do l=2,nlay
          do ig=1,ngrid
             if(larg_cons(ig,l).gt.1.) then
               if (l.gt.lmix(ig)) then
ctest
                 if (zmax(ig)-zmix(ig).lt.1.e-10) then
c                   print*,'pb xxx'
                   xxx(ig,l)=(lmaxa(ig)+1.-l)/(lmaxa(ig)+1.-lmix(ig))
                 else
                 xxx(ig,l)=(zmax(ig)-zlev(ig,l))/(zmax(ig)-zmix(ig))
                 endif
           if (idetr.eq.0) then
               fraca(ig,l)=fracazmix(ig)
           else if (idetr.eq.1) then
               fraca(ig,l)=fracazmix(ig)*xxx(ig,l)
           else if (idetr.eq.2) then
               fraca(ig,l)=fracazmix(ig)*(1.-(1.-xxx(ig,l))**2)
           else
               fraca(ig,l)=fracazmix(ig)*xxx(ig,l)**2
           endif
c     print*,ig,l,lmix(ig),lmaxa(ig),xxx(ig,l),'LLLLLLL'
               fraca(ig,l)=max(fraca(ig,l),0.)
               fraca(ig,l)=min(fraca(ig,l),0.5)
               fracd(ig,l)=1.-fraca(ig,l)
               fracc(ig,l)=larg_cons(ig,l)/(r_aspect*zmax(ig))
             endif
            endif
         enddo
      enddo
      
      print*,'fin calcul fraca'
c      print*,'11 OK convect8'
c     print*,'Ea3 ',wa_moy
c------------------------------------------------------------------
c   Calcul de fracd, wd
c   somme wa - wd = 0
c------------------------------------------------------------------


      do ig=1,ngrid
         fm(ig,1)=0.
         fm(ig,nlay+1)=0.
      enddo

      do l=2,nlay
           do ig=1,ngrid
              fm(ig,l)=fraca(ig,l)*wa_moy(ig,l)*rhobarz(ig,l)
cCR:test
              if (entr(ig,l-1).lt.1e-10.and.fm(ig,l).gt.fm(ig,l-1)
     s            .and.l.gt.lmix(ig)) then
                 fm(ig,l)=fm(ig,l-1)
c                 write(1,*)'ajustement fm, l',l
              endif
c              write(1,*)'ig,l,fm(ig,l)',ig,l,fm(ig,l)
cRC
           enddo
         do ig=1,ngrid
            if(fracd(ig,l).lt.0.1) then
               stop'fracd trop petit'
            else
c    vitesse descendante "diagnostique"
               wd(ig,l)=fm(ig,l)/(fracd(ig,l)*rhobarz(ig,l))
            endif
         enddo
      enddo

      do l=1,nlay
         do ig=1,ngrid
c           masse(ig,l)=rho(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
            masse(ig,l)=(pplev(ig,l)-pplev(ig,l+1))/RG
         enddo
      enddo

      print*,'12 OK convect8'
c     print*,'WA4 ',wa_moy
cc------------------------------------------------------------------
c   calcul du transport vertical
c------------------------------------------------------------------

      go to 4444
c     print*,'XXXXXXXXXXXXXXX ptimestep= ',ptimestep
      do l=2,nlay-1
         do ig=1,ngrid
            if(fm(ig,l+1)*ptimestep.gt.masse(ig,l)
     s      .and.fm(ig,l+1)*ptimestep.gt.masse(ig,l+1)) then
c     print*,'WARN!!! FM>M ig=',ig,' l=',l,'  FM='
c    s         ,fm(ig,l+1)*ptimestep
c    s         ,'   M=',masse(ig,l),masse(ig,l+1)
             endif
         enddo
      enddo

      do l=1,nlay
         do ig=1,ngrid
            if(entr(ig,l)*ptimestep.gt.masse(ig,l)) then
c     print*,'WARN!!! E>M ig=',ig,' l=',l,'  E=='
c    s         ,entr(ig,l)*ptimestep
c    s         ,'   M=',masse(ig,l)
             endif
         enddo
      enddo

      do l=1,nlay
         do ig=1,ngrid
            if(.not.fm(ig,l).ge.0..or..not.fm(ig,l).le.10.) then
c     print*,'WARN!!! fm exagere ig=',ig,'   l=',l
c    s         ,'   FM=',fm(ig,l)
            endif
            if(.not.masse(ig,l).ge.1.e-10
     s         .or..not.masse(ig,l).le.1.e4) then
            endif
            if(.not.entr(ig,l).ge.0..or..not.entr(ig,l).le.10.) then
c     print*,'WARN!!! entr exagere ig=',ig,'   l=',l
c    s         ,'   E=',entr(ig,l)
            endif
         enddo
      enddo

4444   continue

cCR:redefinition du entr
       do l=1,nlay
         do ig=1,ngrid
            detr(ig,l)=fm(ig,l)+entr(ig,l)-fm(ig,l+1)
            if (detr(ig,l).lt.0.) then
                entr(ig,l)=entr(ig,l)-detr(ig,l)
                detr(ig,l)=0.
c     print*,'WARNING !!! detrainement negatif ',ig,l
            endif
         enddo
      enddo
cRC
      if (w2di.eq.1) then
         fm0=fm0+ptimestep*(fm-fm0)/float(tho)
         entr0=entr0+ptimestep*(entr-entr0)/float(tho)
      else
         fm0=fm
         entr0=entr
      endif

      if (1.eq.1) then
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zh,zdhadj,zha)
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zo,pdoadj,zoa)
      else
         call dqthermcell2(ngrid,nlay,ptimestep,fm0,entr0,masse,fraca
     .    ,zh,zdhadj,zha)
         call dqthermcell2(ngrid,nlay,ptimestep,fm0,entr0,masse,fraca
     .    ,zo,pdoadj,zoa)
      endif

      if (1.eq.0) then
         call dvthermcell2(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,fraca,zmax
     .    ,zu,zv,pduadj,pdvadj,zua,zva)
      else
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zu,pduadj,zua)
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zv,pdvadj,zva)
      endif

      do l=1,nlay
         do ig=1,ngrid
            zf=0.5*(fracc(ig,l)+fracc(ig,l+1))
            zf2=zf/(1.-zf)
            thetath2(ig,l)=zf2*(zha(ig,l)-zh(ig,l))**2
            wth2(ig,l)=zf2*(0.5*(wa_moy(ig,l)+wa_moy(ig,l+1)))**2
         enddo
      enddo



c     print*,'13 OK convect8'
c     print*,'WA5 ',wa_moy
      do l=1,nlay
         do ig=1,ngrid
            pdtadj(ig,l)=zdhadj(ig,l)*zpspsk(ig,l)
         enddo
      enddo


c     do l=1,nlay
c        do ig=1,ngrid
c           if(abs(pdtadj(ig,l))*86400..gt.500.) then
c     print*,'WARN!!! ig=',ig,'  l=',l
c    s         ,'   pdtadj=',pdtadj(ig,l)
c           endif
c           if(abs(pdoadj(ig,l))*86400..gt.1.) then
c     print*,'WARN!!! ig=',ig,'  l=',l
c    s         ,'   pdoadj=',pdoadj(ig,l)
c           endif
c        enddo
c      enddo

      print*,'14 OK convect8'
c------------------------------------------------------------------
c   Calculs pour les sorties
c------------------------------------------------------------------

      if(sorties) then
      do l=1,nlay
         do ig=1,ngrid
            zla(ig,l)=(1.-fracd(ig,l))*zmax(ig)
            zld(ig,l)=fracd(ig,l)*zmax(ig)
            if(1.-fracd(ig,l).gt.1.e-10)
     s      zwa(ig,l)=wd(ig,l)*fracd(ig,l)/(1.-fracd(ig,l))
         enddo
      enddo

cdeja fait
c      do l=1,nlay
c         do ig=1,ngrid
c            detr(ig,l)=fm(ig,l)+entr(ig,l)-fm(ig,l+1)
c            if (detr(ig,l).lt.0.) then
c                entr(ig,l)=entr(ig,l)-detr(ig,l)
c                detr(ig,l)=0.
c     print*,'WARNING !!! detrainement negatif ',ig,l
c            endif
c         enddo
c      enddo

c     print*,'15 OK convect8'

      isplit=isplit+1


	goto 123
123   continue

      endif

c     if(wa_moy(1,4).gt.1.e-10) stop

      print*,'19 OK convect8'
      return
      end

      subroutine dqthermcell(ngrid,nlay,ptimestep,fm,entr,masse
     .    ,q,dq,qa)
      use dimens_m
      use dimphy
      implicit none

c=======================================================================
c
c   Calcul du transport verticale dans la couche limite en presence
c   de "thermiques" explicitement representes
c   calcul du dq/dt une fois qu'on connait les ascendances
c
c=======================================================================


      integer ngrid,nlay

      real ptimestep
      real, intent(in):: masse(ngrid,nlay)
      real fm(ngrid,nlay+1)
      real entr(ngrid,nlay)
      real q(ngrid,nlay)
      real dq(ngrid,nlay)

      real qa(klon,klev),detr(klon,klev),wqd(klon,klev+1)

      integer ig,k

c   calcul du detrainement

      do k=1,nlay
         do ig=1,ngrid
            detr(ig,k)=fm(ig,k)-fm(ig,k+1)+entr(ig,k)
         enddo
      enddo

c   calcul de la valeur dans les ascendances
      do ig=1,ngrid
         qa(ig,1)=q(ig,1)
      enddo

      do k=2,nlay
         do ig=1,ngrid
            if ((fm(ig,k+1)+detr(ig,k))*ptimestep.gt.
     s         1.e-5*masse(ig,k)) then
               qa(ig,k)=(fm(ig,k)*qa(ig,k-1)+entr(ig,k)*q(ig,k))
     s         /(fm(ig,k+1)+detr(ig,k))
            else
               qa(ig,k)=q(ig,k)
            endif
         enddo
      enddo

      do k=2,nlay
         do ig=1,ngrid
c             wqd(ig,k)=fm(ig,k)*0.5*(q(ig,k-1)+q(ig,k))
            wqd(ig,k)=fm(ig,k)*q(ig,k)
         enddo
      enddo
      do ig=1,ngrid
         wqd(ig,1)=0.
         wqd(ig,nlay+1)=0.
      enddo

      do k=1,nlay
         do ig=1,ngrid
            dq(ig,k)=(detr(ig,k)*qa(ig,k)-entr(ig,k)*q(ig,k)
     s               -wqd(ig,k)+wqd(ig,k+1))
     s               /masse(ig,k)
         enddo
      enddo

      return
      end
      subroutine dvthermcell(ngrid,nlay,ptimestep,fm,entr,masse
     .    ,fraca,larga
     .    ,u,v,du,dv,ua,va)
      use dimens_m
      use dimphy
      implicit none

c=======================================================================
c
c   Calcul du transport verticale dans la couche limite en presence
c   de "thermiques" explicitement representes
c   calcul du dq/dt une fois qu'on connait les ascendances
c
c=======================================================================


      integer ngrid,nlay

      real ptimestep
      real masse(ngrid,nlay),fm(ngrid,nlay+1)
      real fraca(ngrid,nlay+1)
      real larga(ngrid)
      real entr(ngrid,nlay)
      real u(ngrid,nlay)
      real ua(ngrid,nlay)
      real du(ngrid,nlay)
      real v(ngrid,nlay)
      real va(ngrid,nlay)
      real dv(ngrid,nlay)

      real qa(klon,klev),detr(klon,klev)
      real wvd(klon,klev+1),wud(klon,klev+1)
      real gamma0,gamma(klon,klev+1)
      real dua,dva
      integer iter

      integer ig,k

c   calcul du detrainement

      do k=1,nlay
         do ig=1,ngrid
            detr(ig,k)=fm(ig,k)-fm(ig,k+1)+entr(ig,k)
         enddo
      enddo

c   calcul de la valeur dans les ascendances
      do ig=1,ngrid
         ua(ig,1)=u(ig,1)
         va(ig,1)=v(ig,1)
      enddo

      do k=2,nlay
         do ig=1,ngrid
            if ((fm(ig,k+1)+detr(ig,k))*ptimestep.gt.
     s         1.e-5*masse(ig,k)) then
c   On itère sur la valeur du coeff de freinage.
c              gamma0=rho(ig,k)*(zlev(ig,k+1)-zlev(ig,k))
               gamma0=masse(ig,k)
     s         *sqrt( 0.5*(fraca(ig,k+1)+fraca(ig,k)) )
     s         *0.5/larga(ig)
c              gamma0=0.
c   la première fois on multiplie le coefficient de freinage
c   par le module du vent dans la couche en dessous.
               dua=ua(ig,k-1)-u(ig,k-1)
               dva=va(ig,k-1)-v(ig,k-1)
               do iter=1,5
                  gamma(ig,k)=gamma0*sqrt(dua**2+dva**2)
                  ua(ig,k)=(fm(ig,k)*ua(ig,k-1)
     s               +(entr(ig,k)+gamma(ig,k))*u(ig,k))
     s               /(fm(ig,k+1)+detr(ig,k)+gamma(ig,k))
                  va(ig,k)=(fm(ig,k)*va(ig,k-1)
     s               +(entr(ig,k)+gamma(ig,k))*v(ig,k))
     s               /(fm(ig,k+1)+detr(ig,k)+gamma(ig,k))
c                 print*,k,ua(ig,k),va(ig,k),u(ig,k),v(ig,k),dua,dva
                  dua=ua(ig,k)-u(ig,k)
                  dva=va(ig,k)-v(ig,k)
               enddo
            else
               ua(ig,k)=u(ig,k)
               va(ig,k)=v(ig,k)
               gamma(ig,k)=0.
            endif
         enddo
      enddo

      do k=2,nlay
         do ig=1,ngrid
            wud(ig,k)=fm(ig,k)*u(ig,k)
            wvd(ig,k)=fm(ig,k)*v(ig,k)
         enddo
      enddo
      do ig=1,ngrid
         wud(ig,1)=0.
         wud(ig,nlay+1)=0.
         wvd(ig,1)=0.
         wvd(ig,nlay+1)=0.
      enddo

      do k=1,nlay
         do ig=1,ngrid
            du(ig,k)=((detr(ig,k)+gamma(ig,k))*ua(ig,k)
     s               -(entr(ig,k)+gamma(ig,k))*u(ig,k)
     s               -wud(ig,k)+wud(ig,k+1))
     s               /masse(ig,k)
            dv(ig,k)=((detr(ig,k)+gamma(ig,k))*va(ig,k)
     s               -(entr(ig,k)+gamma(ig,k))*v(ig,k)
     s               -wvd(ig,k)+wvd(ig,k+1))
     s               /masse(ig,k)
         enddo
      enddo

      return
      end
      subroutine dqthermcell2(ngrid,nlay,ptimestep,fm,entr,masse,frac
     .    ,q,dq,qa)
      use dimens_m
      use dimphy
      implicit none

c=======================================================================
c
c   Calcul du transport verticale dans la couche limite en presence
c   de "thermiques" explicitement representes
c   calcul du dq/dt une fois qu'on connait les ascendances
c
c=======================================================================


      integer ngrid,nlay

      real ptimestep
      real masse(ngrid,nlay),fm(ngrid,nlay+1)
      real entr(ngrid,nlay),frac(ngrid,nlay)
      real q(ngrid,nlay)
      real dq(ngrid,nlay)

      real qa(klon,klev),detr(klon,klev),wqd(klon,klev+1)
      real qe(klon,klev),zf,zf2

      integer ig,k

c   calcul du detrainement

      do k=1,nlay
         do ig=1,ngrid
            detr(ig,k)=fm(ig,k)-fm(ig,k+1)+entr(ig,k)
         enddo
      enddo

c   calcul de la valeur dans les ascendances
      do ig=1,ngrid
         qa(ig,1)=q(ig,1)
         qe(ig,1)=q(ig,1)
      enddo

      do k=2,nlay
         do ig=1,ngrid
            if ((fm(ig,k+1)+detr(ig,k))*ptimestep.gt.
     s         1.e-5*masse(ig,k)) then
               zf=0.5*(frac(ig,k)+frac(ig,k+1))
               zf2=1./(1.-zf)
               qa(ig,k)=(fm(ig,k)*qa(ig,k-1)+zf2*entr(ig,k)*q(ig,k))
     s         /(fm(ig,k+1)+detr(ig,k)+entr(ig,k)*zf*zf2)
               qe(ig,k)=(q(ig,k)-zf*qa(ig,k))*zf2
            else
               qa(ig,k)=q(ig,k)
               qe(ig,k)=q(ig,k)
            endif
         enddo
      enddo

      do k=2,nlay
         do ig=1,ngrid
c             wqd(ig,k)=fm(ig,k)*0.5*(q(ig,k-1)+q(ig,k))
            wqd(ig,k)=fm(ig,k)*qe(ig,k)
         enddo
      enddo
      do ig=1,ngrid
         wqd(ig,1)=0.
         wqd(ig,nlay+1)=0.
      enddo

      do k=1,nlay
         do ig=1,ngrid
            dq(ig,k)=(detr(ig,k)*qa(ig,k)-entr(ig,k)*qe(ig,k)
     s               -wqd(ig,k)+wqd(ig,k+1))
     s               /masse(ig,k)
         enddo
      enddo

      return
      end
      subroutine dvthermcell2(ngrid,nlay,ptimestep,fm,entr,masse
     .    ,fraca,larga
     .    ,u,v,du,dv,ua,va)
      use dimens_m
      use dimphy
      implicit none

c=======================================================================
c
c   Calcul du transport verticale dans la couche limite en presence
c   de "thermiques" explicitement representes
c   calcul du dq/dt une fois qu'on connait les ascendances
c
c=======================================================================


      integer ngrid,nlay

      real ptimestep
      real masse(ngrid,nlay),fm(ngrid,nlay+1)
      real fraca(ngrid,nlay+1)
      real larga(ngrid)
      real entr(ngrid,nlay)
      real u(ngrid,nlay)
      real ua(ngrid,nlay)
      real du(ngrid,nlay)
      real v(ngrid,nlay)
      real va(ngrid,nlay)
      real dv(ngrid,nlay)

      real qa(klon,klev),detr(klon,klev),zf,zf2
      real wvd(klon,klev+1),wud(klon,klev+1)
      real gamma0,gamma(klon,klev+1)
      real ue(klon,klev),ve(klon,klev)
      real dua,dva
      integer iter

      integer ig,k

c   calcul du detrainement

      do k=1,nlay
         do ig=1,ngrid
            detr(ig,k)=fm(ig,k)-fm(ig,k+1)+entr(ig,k)
         enddo
      enddo

c   calcul de la valeur dans les ascendances
      do ig=1,ngrid
         ua(ig,1)=u(ig,1)
         va(ig,1)=v(ig,1)
         ue(ig,1)=u(ig,1)
         ve(ig,1)=v(ig,1)
      enddo

      do k=2,nlay
         do ig=1,ngrid
            if ((fm(ig,k+1)+detr(ig,k))*ptimestep.gt.
     s         1.e-5*masse(ig,k)) then
c   On itère sur la valeur du coeff de freinage.
c              gamma0=rho(ig,k)*(zlev(ig,k+1)-zlev(ig,k))
               gamma0=masse(ig,k)
     s         *sqrt( 0.5*(fraca(ig,k+1)+fraca(ig,k)) )
     s         *0.5/larga(ig)
     s         *1.
c    s         *0.5
c              gamma0=0.
               zf=0.5*(fraca(ig,k)+fraca(ig,k+1))
               zf=0.
               zf2=1./(1.-zf)
c   la première fois on multiplie le coefficient de freinage
c   par le module du vent dans la couche en dessous.
               dua=ua(ig,k-1)-u(ig,k-1)
               dva=va(ig,k-1)-v(ig,k-1)
               do iter=1,5
c   On choisit une relaxation lineaire.
                  gamma(ig,k)=gamma0
c   On choisit une relaxation quadratique.
                  gamma(ig,k)=gamma0*sqrt(dua**2+dva**2)
                  ua(ig,k)=(fm(ig,k)*ua(ig,k-1)
     s               +(zf2*entr(ig,k)+gamma(ig,k))*u(ig,k))
     s               /(fm(ig,k+1)+detr(ig,k)+entr(ig,k)*zf*zf2
     s                 +gamma(ig,k))
                  va(ig,k)=(fm(ig,k)*va(ig,k-1)
     s               +(zf2*entr(ig,k)+gamma(ig,k))*v(ig,k))
     s               /(fm(ig,k+1)+detr(ig,k)+entr(ig,k)*zf*zf2
     s                 +gamma(ig,k))
c                 print*,k,ua(ig,k),va(ig,k),u(ig,k),v(ig,k),dua,dva
                  dua=ua(ig,k)-u(ig,k)
                  dva=va(ig,k)-v(ig,k)
                  ue(ig,k)=(u(ig,k)-zf*ua(ig,k))*zf2
                  ve(ig,k)=(v(ig,k)-zf*va(ig,k))*zf2
               enddo
            else
               ua(ig,k)=u(ig,k)
               va(ig,k)=v(ig,k)
               ue(ig,k)=u(ig,k)
               ve(ig,k)=v(ig,k)
               gamma(ig,k)=0.
            endif
         enddo
      enddo

      do k=2,nlay
         do ig=1,ngrid
            wud(ig,k)=fm(ig,k)*ue(ig,k)
            wvd(ig,k)=fm(ig,k)*ve(ig,k)
         enddo
      enddo
      do ig=1,ngrid
         wud(ig,1)=0.
         wud(ig,nlay+1)=0.
         wvd(ig,1)=0.
         wvd(ig,nlay+1)=0.
      enddo

      do k=1,nlay
         do ig=1,ngrid
            du(ig,k)=((detr(ig,k)+gamma(ig,k))*ua(ig,k)
     s               -(entr(ig,k)+gamma(ig,k))*ue(ig,k)
     s               -wud(ig,k)+wud(ig,k+1))
     s               /masse(ig,k)
            dv(ig,k)=((detr(ig,k)+gamma(ig,k))*va(ig,k)
     s               -(entr(ig,k)+gamma(ig,k))*ve(ig,k)
     s               -wvd(ig,k)+wvd(ig,k+1))
     s               /masse(ig,k)
         enddo
      enddo

      return
      end
