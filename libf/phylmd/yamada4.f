!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/yamada4.F,v 1.1 2004/06/22 11:45:36 lmdzadmin Exp $
!
      SUBROUTINE yamada4(ngrid,dt,g,rconst,plev,temp
     s   ,zlev,zlay,u,v,teta,cd,q2,km,kn,kq,ustar
     s   ,iflag_pbl)
      use dimens_m
      use dimphy
      IMPLICIT NONE
c.......................................................................
c.......................................................................
c
c dt : pas de temps
c g  : g
c zlev : altitude a chaque niveau (interface inferieure de la couche
c        de meme indice)
c zlay : altitude au centre de chaque couche
c u,v : vitesse au centre de chaque couche
c       (en entree : la valeur au debut du pas de temps)
c teta : temperature potentielle au centre de chaque couche
c        (en entree : la valeur au debut du pas de temps)
c cd : cdrag
c      (en entree : la valeur au debut du pas de temps)
c q2 : $q^2$ au bas de chaque couche
c      (en entree : la valeur au debut du pas de temps)
c      (en sortie : la valeur a la fin du pas de temps)
c km : diffusivite turbulente de quantite de mouvement (au bas de chaque
c      couche)
c      (en sortie : la valeur a la fin du pas de temps)
c kn : diffusivite turbulente des scalaires (au bas de chaque couche)
c      (en sortie : la valeur a la fin du pas de temps)
c
c  iflag_pbl doit valoir entre 6 et 9
c      l=6, on prend  systematiquement une longueur d'equilibre
c    iflag_pbl=6 : MY 2.0
c    iflag_pbl=7 : MY 2.0.Fournier
c    iflag_pbl=8 : MY 2.5
c    iflag_pbl=9 : un test ?

c.......................................................................
      REAL, intent(in):: dt
      real g,rconst
      real plev(klon,klev+1),temp(klon,klev)
      real ustar(klon)
      real kmin,qmin,pblhmin(klon),coriol(klon)
      REAL zlev(klon,klev+1)
      REAL zlay(klon,klev)
      REAL u(klon,klev)
      REAL v(klon,klev)
      REAL teta(klon,klev)
      REAL cd(klon)
      REAL q2(klon,klev+1),qpre
      REAL unsdz(klon,klev)
      REAL unsdzdec(klon,klev+1)

      REAL km(klon,klev+1)
      REAL kmpre(klon,klev+1),tmp2
      REAL mpre(klon,klev+1)
      REAL kn(klon,klev+1)
      REAL kq(klon,klev+1)
      real ff(klon,klev+1),delta(klon,klev+1)
      real aa(klon,klev+1),aa0,aa1
      integer iflag_pbl,ngrid


      integer nlay,nlev
      PARAMETER (nlay=klev)
      PARAMETER (nlev=klev+1)

      logical first
      integer ipas
      save first,ipas
      data first,ipas/.true.,0/


      integer ig,k


      real ri,zrif,zalpha,zsm,zsn
      real rif(klon,klev+1),sm(klon,klev+1),alpha(klon,klev)

      real m2(klon,klev+1),dz(klon,klev+1),zq,n2(klon,klev+1)
      real dtetadz(klon,klev+1)
      real m2cstat,mcstat,kmcstat
      real l(klon,klev+1),l0(klon)
      save l0

      real sq(klon),sqz(klon),zz(klon,klev+1)
      integer iter

      real ric,rifc,b1,kap
      save ric,rifc,b1,kap
      data ric,rifc,b1,kap/0.195,0.191,16.6,0.4/

      real frif,falpha,fsm
      real fl,zzz,zl0,zq2,zn2

      real rino(klon,klev+1),smyam(klon,klev),styam(klon,klev)
     s  ,lyam(klon,klev),knyam(klon,klev)
     s  ,w2yam(klon,klev),t2yam(klon,klev)
      common/pbldiag/rino,smyam,styam,lyam,knyam,w2yam,t2yam

      frif(ri)=0.6588*(ri+0.1776-sqrt(ri*ri-0.3221*ri+0.03156))
      falpha(ri)=1.318*(0.2231-ri)/(0.2341-ri)
      fsm(ri)=1.96*(0.1912-ri)*(0.2341-ri)/((1.-ri)*(0.2231-ri))
      fl(zzz,zl0,zq2,zn2)=
     s     max(min(l0(ig)*kap*zlev(ig,k)/(kap*zlev(ig,k)+l0(ig))
     s     ,0.5*sqrt(q2(ig,k))/sqrt(max(n2(ig,k),1.e-10))) ,1.)

      if (.not.(iflag_pbl.ge.6.and.iflag_pbl.le.9)) then
           stop'probleme de coherence dans appel a MY'
      endif

      ipas=ipas+1
      if (0.eq.1.and.first) then
      do ig=1,1000
         ri=(ig-800.)/500.
         if (ri.lt.ric) then
            zrif=frif(ri)
         else
            zrif=rifc
         endif
         if(zrif.lt.0.16) then
            zalpha=falpha(zrif)
            zsm=fsm(zrif)
         else
            zalpha=1.12
            zsm=0.085
         endif
c     print*,ri,rif,zalpha,zsm
      enddo
      endif

c.......................................................................
c  les increments verticaux
c.......................................................................
c
c!!!!! allerte !!!!!c
c!!!!! zlev n'est pas declare a nlev !!!!!c
c!!!!! ---->
                                                      DO ig=1,ngrid
            zlev(ig,nlev)=zlay(ig,nlay)
     &             +( zlay(ig,nlay) - zlev(ig,nlev-1) )
                                                      ENDDO
c!!!!! <----
c!!!!! allerte !!!!!c
c
      DO k=1,nlay
                                                      DO ig=1,ngrid
        unsdz(ig,k)=1.E+0/(zlev(ig,k+1)-zlev(ig,k))
                                                      ENDDO
      ENDDO
                                                      DO ig=1,ngrid
      unsdzdec(ig,1)=1.E+0/(zlay(ig,1)-zlev(ig,1))
                                                      ENDDO
      DO k=2,nlay
                                                      DO ig=1,ngrid
        unsdzdec(ig,k)=1.E+0/(zlay(ig,k)-zlay(ig,k-1))
                                                     ENDDO
      ENDDO
                                                      DO ig=1,ngrid
      unsdzdec(ig,nlay+1)=1.E+0/(zlev(ig,nlay+1)-zlay(ig,nlay))
                                                     ENDDO
c
c.......................................................................

      do k=2,klev
                                                          do ig=1,ngrid
         dz(ig,k)=zlay(ig,k)-zlay(ig,k-1)
         m2(ig,k)=((u(ig,k)-u(ig,k-1))**2+(v(ig,k)-v(ig,k-1))**2)
     s             /(dz(ig,k)*dz(ig,k))
         dtetadz(ig,k)=(teta(ig,k)-teta(ig,k-1))/dz(ig,k)
         n2(ig,k)=g*2.*dtetadz(ig,k)/(teta(ig,k-1)+teta(ig,k))
c        n2(ig,k)=0.
         ri=n2(ig,k)/max(m2(ig,k),1.e-10)
         if (ri.lt.ric) then
            rif(ig,k)=frif(ri)
         else
            rif(ig,k)=rifc
         endif
         if(rif(ig,k).lt.0.16) then
            alpha(ig,k)=falpha(rif(ig,k))
            sm(ig,k)=fsm(rif(ig,k))
         else
            alpha(ig,k)=1.12
            sm(ig,k)=0.085
         endif
         zz(ig,k)=b1*m2(ig,k)*(1.-rif(ig,k))*sm(ig,k)
c     print*,'RIF L=',k,rif(ig,k),ri*alpha(ig,k)


                                                          enddo
      enddo


c====================================================================
c   Au premier appel, on determine l et q2 de facon iterative.
c iterration pour determiner la longueur de melange


      if (first.or.iflag_pbl.eq.6) then
                                                          do ig=1,ngrid
      l0(ig)=10.
                                                          enddo
      do k=2,klev-1
                                                          do ig=1,ngrid
        l(ig,k)=l0(ig)*kap*zlev(ig,k)/(kap*zlev(ig,k)+l0(ig))
                                                          enddo
      enddo

      do iter=1,10
                                                          do ig=1,ngrid
         sq(ig)=1.e-10
         sqz(ig)=1.e-10
                                                          enddo
         do k=2,klev-1
                                                          do ig=1,ngrid
           q2(ig,k)=l(ig,k)**2*zz(ig,k)
           l(ig,k)=fl(zlev(ig,k),l0(ig),q2(ig,k),n2(ig,k))
           zq=sqrt(q2(ig,k))
           sqz(ig)=sqz(ig)+zq*zlev(ig,k)*(zlay(ig,k)-zlay(ig,k-1))
           sq(ig)=sq(ig)+zq*(zlay(ig,k)-zlay(ig,k-1))
                                                          enddo
         enddo
                                                          do ig=1,ngrid
         l0(ig)=0.2*sqz(ig)/sq(ig)
c        l0(ig)=30.
                                                          enddo
c      print*,'ITER=',iter,'  L0=',l0

      enddo

c     print*,'Fin de l initialisation de q2 et l0'

      endif ! first

c====================================================================
c  Calcul de la longueur de melange.
c====================================================================

c   Mise a jour de l0
                                                          do ig=1,ngrid
      sq(ig)=1.e-10
      sqz(ig)=1.e-10
                                                          enddo
      do k=2,klev-1
                                                          do ig=1,ngrid
        zq=sqrt(q2(ig,k))
        sqz(ig)=sqz(ig)+zq*zlev(ig,k)*(zlay(ig,k)-zlay(ig,k-1))
        sq(ig)=sq(ig)+zq*(zlay(ig,k)-zlay(ig,k-1))
                                                          enddo
      enddo
                                                          do ig=1,ngrid
      l0(ig)=0.2*sqz(ig)/sq(ig)
c        l0(ig)=30.
                                                          enddo
c      print*,'ITER=',iter,'  L0=',l0
c   calcul de l(z)
      do k=2,klev
                                                          do ig=1,ngrid
         l(ig,k)=fl(zlev(ig,k),l0(ig),q2(ig,k),n2(ig,k))
         if(first) then
           q2(ig,k)=l(ig,k)**2*zz(ig,k)
         endif
                                                          enddo
      enddo

c====================================================================
c   Yamada 2.0
c====================================================================
      if (iflag_pbl.eq.6) then

      do k=2,klev
                                                          do ig=1,ngrid
         q2(ig,k)=l(ig,k)**2*zz(ig,k)
                                                          enddo
      enddo


      else if (iflag_pbl.eq.7) then
c====================================================================
c   Yamada 2.Fournier
c====================================================================

c  Calcul de l,  km, au pas precedent
      do k=2,klev
                                                          do ig=1,ngrid
c        print*,'SMML=',sm(ig,k),l(ig,k)
         delta(ig,k)=q2(ig,k)/(l(ig,k)**2*sm(ig,k))
         kmpre(ig,k)=l(ig,k)*sqrt(q2(ig,k))*sm(ig,k)
         mpre(ig,k)=sqrt(m2(ig,k))
c        print*,'0L=',k,l(ig,k),delta(ig,k),km(ig,k)
                                                          enddo
      enddo

      do k=2,klev-1
                                                          do ig=1,ngrid
        m2cstat=max(alpha(ig,k)*n2(ig,k)+delta(ig,k)/b1,1.e-12)
        mcstat=sqrt(m2cstat)

c        print*,'M2 L=',k,mpre(ig,k),mcstat
c
c  -----{puis on ecrit la valeur de q qui annule l'equation de m
c        supposee en q3}
c
        IF (k.eq.2) THEN
          kmcstat=1.E+0 / mcstat
     &    *( unsdz(ig,k)*kmpre(ig,k+1)
     &                        *mpre(ig,k+1)
     &      +unsdz(ig,k-1)
     &              *cd(ig)
     &              *( sqrt(u(ig,3)**2+v(ig,3)**2)
     &                -mcstat/unsdzdec(ig,k)
     &                -mpre(ig,k+1)/unsdzdec(ig,k+1) )**2)
     &      /( unsdz(ig,k)+unsdz(ig,k-1) )
        ELSE
          kmcstat=1.E+0 / mcstat
     &    *( unsdz(ig,k)*kmpre(ig,k+1)
     &                        *mpre(ig,k+1)
     &      +unsdz(ig,k-1)*kmpre(ig,k-1)
     &                          *mpre(ig,k-1) )
     &      /( unsdz(ig,k)+unsdz(ig,k-1) )
        ENDIF
c       print*,'T2 L=',k,tmp2
        tmp2=kmcstat
     &      /( sm(ig,k)/q2(ig,k) )
     &      /l(ig,k)
        q2(ig,k)=max(tmp2,1.e-12)**(2./3.)
c       print*,'Q2 L=',k,q2(ig,k)
c
                                                          enddo
      enddo

      else if (iflag_pbl.ge.8) then
c====================================================================
c   Yamada 2.5 a la Didi
c====================================================================


c  Calcul de l,  km, au pas precedent
      do k=2,klev
                                                          do ig=1,ngrid
c        print*,'SMML=',sm(ig,k),l(ig,k)
         delta(ig,k)=q2(ig,k)/(l(ig,k)**2*sm(ig,k))
         if (delta(ig,k).lt.1.e-20) then
c     print*,'ATTENTION   L=',k,'   Delta=',delta(ig,k)
            delta(ig,k)=1.e-20
         endif
         km(ig,k)=l(ig,k)*sqrt(q2(ig,k))*sm(ig,k)
         aa0=
     s   (m2(ig,k)-alpha(ig,k)*n2(ig,k)-delta(ig,k)/b1)
         aa1=
     s   (m2(ig,k)*(1.-rif(ig,k))-delta(ig,k)/b1)
c abder      print*,'AA L=',k,aa0,aa1,aa1/max(m2(ig,k),1.e-20)
         aa(ig,k)=aa1*dt/(delta(ig,k)*l(ig,k))
c     print*,'0L=',k,l(ig,k),delta(ig,k),km(ig,k)
         qpre=sqrt(q2(ig,k))
         if (iflag_pbl.eq.8 ) then
            if (aa(ig,k).gt.0.) then
               q2(ig,k)=(qpre+aa(ig,k)*qpre*qpre)**2
            else
               q2(ig,k)=(qpre/(1.-aa(ig,k)*qpre))**2
            endif
         else ! iflag_pbl=9
            if (aa(ig,k)*qpre.gt.0.9) then
               q2(ig,k)=(qpre*10.)**2
            else
               q2(ig,k)=(qpre/(1.-aa(ig,k)*qpre))**2
            endif
         endif
         q2(ig,k)=min(max(q2(ig,k),1.e-10),1.e4)
c     print*,'Q2 L=',k,q2(ig,k),qpre*qpre
                                                          enddo
      enddo

      endif ! Fin du cas 8

c     print*,'OK8'

c====================================================================
c   Calcul des coefficients de mélange
c====================================================================
      do k=2,klev
c     print*,'k=',k
                                                          do ig=1,ngrid
cabde      print*,'KML=',l(ig,k),q2(ig,k),sm(ig,k)
         zq=sqrt(q2(ig,k))
         km(ig,k)=l(ig,k)*zq*sm(ig,k)
         kn(ig,k)=km(ig,k)*alpha(ig,k)
         kq(ig,k)=l(ig,k)*zq*0.2
c     print*,'KML=',km(ig,k),kn(ig,k)
                                                          enddo
      enddo

c   Traitement des cas noctrunes avec l'introduction d'une longueur
c   minilale.

c====================================================================
c   Traitement particulier pour les cas tres stables.
c   D'apres Holtslag Boville.

      print*,'YAMADA4 0'

                                                          do ig=1,ngrid
      coriol(ig)=1.e-4
      pblhmin(ig)=0.07*ustar(ig)/max(abs(coriol(ig)),2.546e-5)
                                                          enddo

       print*,'pblhmin ',pblhmin
CTest a remettre 21 11 02
c test abd 13 05 02      if(0.eq.1) then
      if(1.eq.1) then
      do k=2,klev
         do ig=1,klon
            if (teta(ig,2).gt.teta(ig,1)) then
               qmin=ustar(ig)*(max(1.-zlev(ig,k)/pblhmin(ig),0.))**2
               kmin=kap*zlev(ig,k)*qmin
            else
               kmin=-1. ! kmin n'est utilise que pour les SL stables.
            endif 
            if (kn(ig,k).lt.kmin.or.km(ig,k).lt.kmin) then
c               print*,'Seuil min Km K=',k,kmin,km(ig,k),kn(ig,k)
c     s           ,sqrt(q2(ig,k)),pblhmin(ig),qmin/sm(ig,k)
               kn(ig,k)=kmin
               km(ig,k)=kmin
               kq(ig,k)=kmin
c   la longueur de melange est suposee etre l= kap z
c   K=l q Sm d'ou q2=(K/l Sm)**2
               q2(ig,k)=(qmin/sm(ig,k))**2
            endif
         enddo
      enddo
      endif

      print*,'YAMADA4 1'
c   Diagnostique pour stokage

      rino=rif
      smyam(:,1:klev)=sm(:,1:klev)
      styam=sm(:,1:klev)*alpha(:,1:klev)
      lyam(1:klon,1:klev)=l(:,1:klev)
      knyam(1:klon,1:klev)=kn(:,1:klev)

c   Estimations de w'2 et T'2 d'apres Abdela et McFarlane

        if(1.eq.0)then
      w2yam=q2(:,1:klev)*0.24
     s    +lyam(:,1:klev)*5.17*kn(:,1:klev)*n2(:,1:klev)
     s   /sqrt(q2(:,1:klev))

      t2yam=9.1*kn(:,1:klev)*dtetadz(:,1:klev)**2/sqrt(q2(:,1:klev))
     s  *lyam(:,1:klev)
	endif

c     print*,'OKFIN'
      first=.false.
      return
      end
