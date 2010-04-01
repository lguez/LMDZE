!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/bilan_dyn.F,v 1.5 2005/03/16 10:12:17 fairhead Exp $
!
      SUBROUTINE bilan_dyn (ntrac,dt_app,dt_cum,
     s  ps,masse,pk,flux_u,flux_v,teta,phi,ucov,vcov,trac)

c   AFAIRE
c   Prevoir en champ nq+1 le diagnostique de l'energie
c   en faisant Qzon=Cv T + L * ...
c             vQ..A=Cp T + L * ...

      USE histcom
      use calendar
      use histwrite_m
      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use comgeom, only: constang_2d, cu_2d, cv_2d, rlatv
      use temps
      use iniprint
      use inigrads_m, only: inigrads

      IMPLICIT NONE


c====================================================================
c
c   Sous-programme consacre à des diagnostics dynamiques de base
c
c 
c   De facon generale, les moyennes des scalaires Q sont ponderees par
c   la masse.
c
c   Les flux de masse sont eux simplement moyennes.
c
c====================================================================

c   Arguments :
c   ===========

      integer ntrac
      real dt_app,dt_cum
      real ps(iip1,jjp1)
      real masse(iip1,jjp1,llm),pk(iip1,jjp1,llm)
      real flux_u(iip1,jjp1,llm)
      real flux_v(iip1,jjm,llm)
      real teta(iip1,jjp1,llm)
      real phi(iip1,jjp1,llm)
      real ucov(iip1,jjp1,llm)
      real vcov(iip1,jjm,llm)
      real trac(iip1,jjp1,llm,ntrac)

c   Local :
c   =======

      integer icum,ncum
      logical first
      real zz,zqy,zfactv(jjm,llm)

      integer nQ
      parameter (nQ=7)


cym      character*6 nom(nQ)
cym      character*6 unites(nQ)
      character*6,save :: nom(nQ)
      character*6,save :: unites(nQ)

      character*10 file
      integer ifile
      parameter (ifile=4)

      integer itemp,igeop,iecin,iang,iu,iovap,iun
      integer i_sortie

      save first,icum,ncum
      save itemp,igeop,iecin,iang,iu,iovap,iun
      save i_sortie

      real time
      integer itau
      save time,itau
      data time,itau/0.,0/

      data first/.true./
      data itemp,igeop,iecin,iang,iu,iovap,iun/1,2,3,4,5,6,7/
      data i_sortie/1/

      real ww

c   variables dynamiques intermédiaires
      REAL vcont(iip1,jjm,llm),ucont(iip1,jjp1,llm)
      REAL ang(iip1,jjp1,llm),unat(iip1,jjp1,llm)
      REAL massebx(iip1,jjp1,llm),masseby(iip1,jjm,llm)
      REAL vorpot(iip1,jjm,llm)
      REAL w(iip1,jjp1,llm),ecin(iip1,jjp1,llm),convm(iip1,jjp1,llm)
      REAL bern(iip1,jjp1,llm)

c   champ contenant les scalaires advectés.
      real Q(iip1,jjp1,llm,nQ)
    
c   champs cumulés
      real ps_cum(iip1,jjp1)
      real masse_cum(iip1,jjp1,llm)
      real flux_u_cum(iip1,jjp1,llm)
      real flux_v_cum(iip1,jjm,llm)
      real Q_cum(iip1,jjp1,llm,nQ)
      real flux_uQ_cum(iip1,jjp1,llm,nQ)
      real flux_vQ_cum(iip1,jjm,llm,nQ)
      real flux_wQ_cum(iip1,jjp1,llm,nQ)
      real dQ(iip1,jjp1,llm,nQ)

      save ps_cum,masse_cum,flux_u_cum,flux_v_cum
      save Q_cum,flux_uQ_cum,flux_vQ_cum

c   champs de tansport en moyenne zonale
      integer ntr,itr
      parameter (ntr=5)

cym      character*10 znom(ntr,nQ)
cym      character*20 znoml(ntr,nQ)
cym      character*10 zunites(ntr,nQ)
      character*10,save :: znom(ntr,nQ)
      character*20,save :: znoml(ntr,nQ)
      character*10,save :: zunites(ntr,nQ)

      integer iave,itot,immc,itrs,istn
      data iave,itot,immc,itrs,istn/1,2,3,4,5/
      character*3 ctrs(ntr)
      data ctrs/'  ','TOT','MMC','TRS','STN'/

      real zvQ(jjm,llm,ntr,nQ),zvQtmp(jjm,llm)
      real zavQ(jjm,ntr,nQ),psiQ(jjm,llm+1,nQ)
      real zmasse(jjm,llm),zamasse(jjm)

      real zv(jjm,llm),psi(jjm,llm+1)

      integer i,j,l,iQ


c   Initialisation du fichier contenant les moyennes zonales.
c   ---------------------------------------------------------

      integer fileid
      integer thoriid, zvertiid
      save fileid

      integer ndex3d(jjm*llm)

C   Variables locales
C
      real zjulian
      character*3 str
      character*10 ctrac
      integer ii,jj
      integer zan, dayref
C
      real rlong(jjm),rlatg(jjm)

      !!print *, "Call sequence information: bilan_dyn"

c=====================================================================
c   Initialisation
c=====================================================================

      time=time+dt_app
      itau=itau+1

      if (first) then


        icum=0
c       initialisation des fichiers
        first=.false.
c   ncum est la frequence de stokage en pas de temps
        ncum=dt_cum/dt_app
        if (abs(ncum*dt_app-dt_cum).gt.1.e-5*dt_app) then
           print *,
     .            'Pb : le pas de cumule doit etre multiple du pas'
           print *,'dt_app=',dt_app
           print *,'dt_cum=',dt_cum
           stop
        endif

        if (i_sortie.eq.1) then
         file='dynzon'
         call inigrads(ifile ,(/0./),180./pi,0.,0.,rlatv,-90.,90.,
     $        180./pi ,presnivs,1. ,dt_cum,file,'dyn_zon ')
        endif

        nom(itemp)='T'
        nom(igeop)='gz'
        nom(iecin)='K'
        nom(iang)='ang'
        nom(iu)='u'
        nom(iovap)='ovap'
        nom(iun)='un'

        unites(itemp)='K'
        unites(igeop)='m2/s2'
        unites(iecin)='m2/s2'
        unites(iang)='ang'
        unites(iu)='m/s'
        unites(iovap)='kg/kg'
        unites(iun)='un'


c   Initialisation du fichier contenant les moyennes zonales.
c   ---------------------------------------------------------

      zan = annee_ref
      dayref = day_ref
      CALL ymds2ju(zan, 1, dayref, 0.0, zjulian)
      
      rlong=0.
      rlatg=rlatv*180./pi
       
      call histbeg_totreg('dynzon', rlong(:1), rlatg,
     .             1, 1, 1, jjm,
     .             itau_dyn, zjulian, dt_cum, thoriid, fileid)

C
C  Appel a histvert pour la grille verticale
C
      call histvert(fileid, 'presnivs', 'Niveaux sigma','mb',
     .              llm, presnivs, zvertiid)
C
C  Appels a histdef pour la definition des variables a sauvegarder
      do iQ=1,nQ
         do itr=1,ntr
            if(itr.eq.1) then
               znom(itr,iQ)=nom(iQ)
               znoml(itr,iQ)=nom(iQ)
               zunites(itr,iQ)=unites(iQ)
            else
               znom(itr,iQ)=ctrs(itr)//'v'//nom(iQ)
               znoml(itr,iQ)='transport : v * '//nom(iQ)//' '//ctrs(itr)
               zunites(itr,iQ)='m/s * '//unites(iQ)
            endif
         enddo
      enddo

c   Declarations des champs avec dimension verticale
c      print*,'1HISTDEF'
      do iQ=1,nQ
         do itr=1,ntr
      IF (prt_level > 5)
     . print *,'var ',itr,iQ
     .      ,znom(itr,iQ),znoml(itr,iQ),zunites(itr,iQ)
            call histdef(fileid,znom(itr,iQ),znoml(itr,iQ),
     .        zunites(itr,iQ),1,jjm,thoriid,llm,1,llm,zvertiid,
     .        'ave(X)',dt_cum,dt_cum)
         enddo
c   Declarations pour les fonctions de courant
c      print*,'2HISTDEF'
          call histdef(fileid,'psi'//nom(iQ)
     .      ,'stream fn. '//znoml(itot,iQ),
     .      zunites(itot,iQ),1,jjm,thoriid,llm,1,llm,zvertiid,
     .      'ave(X)',dt_cum,dt_cum)
      enddo


c   Declarations pour les champs de transport d'air
c      print*,'3HISTDEF'
      call histdef(fileid, 'masse', 'masse',
     .             'kg', 1, jjm, thoriid, llm, 1, llm, zvertiid,
     .             'ave(X)', dt_cum, dt_cum)
      call histdef(fileid, 'v', 'v',
     .             'm/s', 1, jjm, thoriid, llm, 1, llm, zvertiid,
     .             'ave(X)', dt_cum, dt_cum)
c   Declarations pour les fonctions de courant
c      print*,'4HISTDEF'
          call histdef(fileid,'psi','stream fn. MMC ','mega t/s',
     .      1,jjm,thoriid,llm,1,llm,zvertiid,
     .      'ave(X)',dt_cum,dt_cum)


c   Declaration des champs 1D de transport en latitude
c      print*,'5HISTDEF'
      do iQ=1,nQ
         do itr=2,ntr
            call histdef(fileid,'a'//znom(itr,iQ),znoml(itr,iQ),
     .        zunites(itr,iQ),1,jjm,thoriid,1,1,1,-99,
     .        'ave(X)',dt_cum,dt_cum)
         enddo
      enddo


c      print*,'8HISTDEF'
               CALL histend(fileid)


      endif


c=====================================================================
c   Calcul des champs dynamiques
c   ----------------------------

c   énergie cinétique
      ucont(:,:,:)=0
      CALL covcont(llm,ucov,vcov,ucont,vcont)
      CALL enercin(vcov,ucov,vcont,ucont,ecin)

c   moment cinétique
      do l=1,llm
         ang(:,:,l)=ucov(:,:,l)+constang_2d(:,:)
         unat(:,:,l)=ucont(:,:,l)*cu_2d(:,:)
      enddo

      Q(:,:,:,itemp)=teta(:,:,:)*pk(:,:,:)/cpp
      Q(:,:,:,igeop)=phi(:,:,:)
      Q(:,:,:,iecin)=ecin(:,:,:)
      Q(:,:,:,iang)=ang(:,:,:)
      Q(:,:,:,iu)=unat(:,:,:)
      Q(:,:,:,iovap)=trac(:,:,:,1)
      Q(:,:,:,iun)=1.


c=====================================================================
c   Cumul
c=====================================================================
c
      if(icum.EQ.0) then
         ps_cum=0.
         masse_cum=0.
         flux_u_cum=0.
         flux_v_cum=0.
         Q_cum=0.
         flux_vQ_cum=0.
         flux_uQ_cum=0.
      endif

      IF (prt_level > 5)
     . print *,'dans bilan_dyn ',icum,'->',icum+1
      icum=icum+1

c   accumulation des flux de masse horizontaux
      ps_cum=ps_cum+ps
      masse_cum=masse_cum+masse
      flux_u_cum=flux_u_cum+flux_u
      flux_v_cum=flux_v_cum+flux_v
      do iQ=1,nQ
      Q_cum(:,:,:,iQ)=Q_cum(:,:,:,iQ)+Q(:,:,:,iQ)*masse(:,:,:)
      enddo

c=====================================================================
c  FLUX ET TENDANCES
c=====================================================================

c   Flux longitudinal
c   -----------------
      do iQ=1,nQ
         do l=1,llm
            do j=1,jjp1
               do i=1,iim
                  flux_uQ_cum(i,j,l,iQ)=flux_uQ_cum(i,j,l,iQ)
     s            +flux_u(i,j,l)*0.5*(Q(i,j,l,iQ)+Q(i+1,j,l,iQ))
               enddo
               flux_uQ_cum(iip1,j,l,iQ)=flux_uQ_cum(1,j,l,iQ)
            enddo
         enddo
      enddo

c    flux méridien
c    -------------
      do iQ=1,nQ
         do l=1,llm
            do j=1,jjm
               do i=1,iip1
                  flux_vQ_cum(i,j,l,iQ)=flux_vQ_cum(i,j,l,iQ)
     s            +flux_v(i,j,l)*0.5*(Q(i,j,l,iQ)+Q(i,j+1,l,iQ))
               enddo
            enddo
         enddo
      enddo


c    tendances
c    ---------

c   convergence horizontale
      call  convflu(flux_uQ_cum,flux_vQ_cum,llm*nQ,dQ)

c   calcul de la vitesse verticale
      call convmas(flux_u_cum,flux_v_cum,convm)
      CALL vitvert(convm,w)

      do iQ=1,nQ
         do l=1,llm-1
            do j=1,jjp1
               do i=1,iip1
                  ww=-0.5*w(i,j,l+1)*(Q(i,j,l,iQ)+Q(i,j,l+1,iQ))
                  dQ(i,j,l  ,iQ)=dQ(i,j,l  ,iQ)-ww
                  dQ(i,j,l+1,iQ)=dQ(i,j,l+1,iQ)+ww
               enddo
            enddo
         enddo
      enddo
      IF (prt_level > 5)
     . print *,'Apres les calculs fait a chaque pas'
c=====================================================================
c   PAS DE TEMPS D'ECRITURE
c=====================================================================
      if (icum.eq.ncum) then
c=====================================================================

      IF (prt_level > 5)
     . print *,'Pas d ecriture'

c   Normalisation
      do iQ=1,nQ
         Q_cum(:,:,:,iQ)=Q_cum(:,:,:,iQ)/masse_cum(:,:,:)
      enddo
      zz=1./float(ncum)
      ps_cum=ps_cum*zz
      masse_cum=masse_cum*zz
      flux_u_cum=flux_u_cum*zz
      flux_v_cum=flux_v_cum*zz
      flux_uQ_cum=flux_uQ_cum*zz
      flux_vQ_cum=flux_vQ_cum*zz
      dQ=dQ*zz


c   A retravailler eventuellement
c   division de dQ par la masse pour revenir aux bonnes grandeurs
      do iQ=1,nQ
         dQ(:,:,:,iQ)=dQ(:,:,:,iQ)/masse_cum(:,:,:)
      enddo
 
c=====================================================================
c   Transport méridien
c=====================================================================

c   cumul zonal des masses des mailles
c   ----------------------------------
      zv=0.
      zmasse=0.
      call massbar(masse_cum,massebx,masseby)
      do l=1,llm
         do j=1,jjm
            do i=1,iim
               zmasse(j,l)=zmasse(j,l)+masseby(i,j,l)
               zv(j,l)=zv(j,l)+flux_v_cum(i,j,l)
            enddo
            zfactv(j,l)=cv_2d(1,j)/zmasse(j,l)
         enddo
      enddo

c     print*,'3OK'
c   --------------------------------------------------------------
c   calcul de la moyenne zonale du transport :
c   ------------------------------------------
c
c                                     --
c TOT : la circulation totale       [ vq ]
c
c                                      -     -
c MMC : mean meridional circulation [ v ] [ q ]
c
c                                     ----      --       - -
c TRS : transitoires                [ v'q'] = [ vq ] - [ v q ]
c
c                                     - * - *       - -       -     -
c STT : stationaires                [ v   q   ] = [ v q ] - [ v ] [ q ]
c
c                                              - -
c    on utilise aussi l'intermediaire TMP :  [ v q ]
c
c    la variable zfactv transforme un transport meridien cumule
c    en kg/s * unte-du-champ-transporte en m/s * unite-du-champ-transporte
c
c   --------------------------------------------------------------


c   ----------------------------------------
c   Transport dans le plan latitude-altitude
c   ----------------------------------------

      zvQ=0.
      psiQ=0.
      do iQ=1,nQ
         zvQtmp=0.
         do l=1,llm
            do j=1,jjm
c              print*,'j,l,iQ=',j,l,iQ
c   Calcul des moyennes zonales du transort total et de zvQtmp
               do i=1,iim
                  zvQ(j,l,itot,iQ)=zvQ(j,l,itot,iQ)
     s                            +flux_vQ_cum(i,j,l,iQ)
                  zqy=      0.5*(Q_cum(i,j,l,iQ)*masse_cum(i,j,l)+
     s                           Q_cum(i,j+1,l,iQ)*masse_cum(i,j+1,l))
                  zvQtmp(j,l)=zvQtmp(j,l)+flux_v_cum(i,j,l)*zqy
     s             /(0.5*(masse_cum(i,j,l)+masse_cum(i,j+1,l)))
                  zvQ(j,l,iave,iQ)=zvQ(j,l,iave,iQ)+zqy
               enddo
c              print*,'aOK'
c   Decomposition
               zvQ(j,l,iave,iQ)=zvQ(j,l,iave,iQ)/zmasse(j,l)
               zvQ(j,l,itot,iQ)=zvQ(j,l,itot,iQ)*zfactv(j,l)
               zvQtmp(j,l)=zvQtmp(j,l)*zfactv(j,l)
               zvQ(j,l,immc,iQ)=zv(j,l)*zvQ(j,l,iave,iQ)*zfactv(j,l)
               zvQ(j,l,itrs,iQ)=zvQ(j,l,itot,iQ)-zvQtmp(j,l)
               zvQ(j,l,istn,iQ)=zvQtmp(j,l)-zvQ(j,l,immc,iQ)
            enddo
         enddo
c   fonction de courant meridienne pour la quantite Q
         do l=llm,1,-1
            do j=1,jjm
               psiQ(j,l,iQ)=psiQ(j,l+1,iQ)+zvQ(j,l,itot,iQ)
            enddo
         enddo
      enddo

c   fonction de courant pour la circulation meridienne moyenne
      psi=0.
      do l=llm,1,-1
         do j=1,jjm
            psi(j,l)=psi(j,l+1)+zv(j,l)
            zv(j,l)=zv(j,l)*zfactv(j,l)
         enddo
      enddo

c     print*,'4OK'
c   sorties proprement dites
      if (i_sortie.eq.1) then
      do iQ=1,nQ
         do itr=1,ntr
            call histwrite(fileid,znom(itr,iQ),itau,zvQ(:,:,itr,iQ))
         enddo
         call histwrite(fileid,'psi'//nom(iQ),itau,psiQ(:,1:llm,iQ))
      enddo

      call histwrite(fileid,'masse',itau,zmasse)
      call histwrite(fileid,'v',itau,zv)
      psi=psi*1.e-9
      call histwrite(fileid,'psi',itau,psi(:,1:llm))

      endif


c   -----------------
c   Moyenne verticale
c   -----------------

      zamasse=0.
      do l=1,llm
         zamasse(:)=zamasse(:)+zmasse(:,l)
      enddo
      zavQ=0.
      do iQ=1,nQ
         do itr=2,ntr
            do l=1,llm
               zavQ(:,itr,iQ)=zavQ(:,itr,iQ)+zvQ(:,l,itr,iQ)*zmasse(:,l)
            enddo
            zavQ(:,itr,iQ)=zavQ(:,itr,iQ)/zamasse(:)
            call histwrite(fileid,'a'//znom(itr,iQ),itau,zavQ(:,itr,iQ))
         enddo
      enddo

c     on doit pouvoir tracer systematiquement la fonction de courant.

c=====================================================================
c/////////////////////////////////////////////////////////////////////
      icum=0                  !///////////////////////////////////////
      endif ! icum.eq.ncum    !///////////////////////////////////////
c/////////////////////////////////////////////////////////////////////
c=====================================================================

      return
      end
