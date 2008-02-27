!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/read_reanalyse.F,v 1.3 2005/04/15 12:31:21 lmdzadmin Exp $
!
c
c
      subroutine read_reanalyse(timestep,psi
     s   ,u,v,t,q,masse,ps,mode,nlevnc)

c   mode=0 variables naturelles
c   mode=1 variabels GCM

c -----------------------------------------------------------------
c   Declarations
c -----------------------------------------------------------------
      use dimens_m
      use paramet_m
      use comvert
      use comgeom
      use guide_m
      IMPLICIT NONE

c common
c ------

      include "netcdf.inc"


c arguments
c ---------
      integer nlevnc
      integer timestep,mode,l

      real psi(iip1,jjp1)
      real u(iip1,jjp1,llm),v(iip1,jjm,llm)
      real t(iip1,jjp1,llm),ps(iip1,jjp1),q(iip1,jjp1,llm)
      real masse(iip1,jjp1,llm),pk(iip1,jjp1,llm)


c local
c -----
      integer ncidu,varidu,ncidv,varidv,ncidt,varidt,ncidps,varidps
      integer ncidpl
      integer varidpl,ncidQ,varidQ
      save ncidu,varidu,ncidv,varidv,ncidt,varidt,ncidps,varidps
      save ncidpl
      save varidpl,ncidQ,varidQ

      real*4 unc(iip1,jjp1,nlevnc),vnc(iip1,jjm,nlevnc)
      real*4 tnc(iip1,jjp1,nlevnc),psnc(iip1,jjp1)
      real*4 Qnc(iip1,jjp1,nlevnc)
      real*4 pl(nlevnc)

      integer start(4),count(4),status

      real rcode
      logical first
      save first

      data first/.true./



c -----------------------------------------------------------------
c   Initialisation de la lecture des fichiers
c -----------------------------------------------------------------
      if (first) then
           ncidpl=-99
           print*,'Intitialisation de read reanalsye'

c Vent zonal
            if (guide_u) then
            ncidu=NCOPN('u.nc',NCNOWRIT,rcode)
            varidu=NCVID(ncidu,'UWND',rcode)
            print*,'ncidu,varidu',ncidu,varidu
            if (ncidpl.eq.-99) ncidpl=ncidu
            endif

c Vent meridien
            if (guide_v) then
            ncidv=NCOPN('v.nc',NCNOWRIT,rcode)
            varidv=NCVID(ncidv,'VWND',rcode)
            print*,'ncidv,varidv',ncidv,varidv
            if (ncidpl.eq.-99) ncidpl=ncidv
            endif

c Temperature
            if (guide_T) then
            ncidt=NCOPN('T.nc',NCNOWRIT,rcode)
            varidt=NCVID(ncidt,'AIR',rcode)
            print*,'ncidt,varidt',ncidt,varidt
            if (ncidpl.eq.-99) ncidpl=ncidt
            endif

c Humidite
            if (guide_Q) then
            ncidQ=NCOPN('hur.nc',NCNOWRIT,rcode)
            varidQ=NCVID(ncidQ,'RH',rcode)
            print*,'ncidQ,varidQ',ncidQ,varidQ
            if (ncidpl.eq.-99) ncidpl=ncidQ
            endif

c Pression de surface
            if (guide_P) then
            ncidps=NCOPN('ps.nc',NCNOWRIT,rcode)
            varidps=NCVID(ncidps,'SP',rcode)
            print*,'ncidps,varidps',ncidps,varidps
            endif

c Coordonnee verticale
            if (ncep) then
               print*,'Vous etes entrain de lire des donnees NCEP'
               varidpl=NCVID(ncidpl,'LEVEL',rcode)
            else
               print*,'Vous etes entrain de lire des donnees ECMWF'
               varidpl=NCVID(ncidpl,'PRESSURE',rcode)
            endif
            print*,'ncidu,varidpl',ncidu,varidpl
      endif
      print*,'ok1'

c Niveaux de pression
      print*,'WARNING!!! Il n y a pas de test de coherence'
      print*,'sur le nombre de niveaux verticaux dans le fichier nc'
      status=NF_GET_VARA_REAL(ncidpl,varidpl,1,nlevnc,pl)
c  passage en pascal
      pl(:)=100.*pl(:)
      if (first) then
       do l=1,nlevnc
          print*,'PL(',l,')=',pl(l)
       enddo
      endif

c -----------------------------------------------------------------
c   lecture des champs u, v, T, ps
c -----------------------------------------------------------------

c  dimensions pour les champs scalaires et le vent zonal
c  -----------------------------------------------------

      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=timestep

      count(1)=iip1
      count(2)=jjp1
      count(3)=nlevnc
      count(4)=1

c mise a zero des tableaux 
c ------------------------
       unc(:,:,:)=0.
       vnc(:,:,:)=0.
       tnc(:,:,:)=0.
       Qnc(:,:,:)=0.

c  Vent zonal
c  ----------

      if (guide_u) then
      print*,'avant la lecture de UNCEP nd de niv:',nlevnc
      status=NF_GET_VARA_REAL(ncidu,varidu,start,count,unc)
c     call dump2d(iip1,jjp1,unc,'VENT NCEP   ')
c     call dump2d(iip1,40,unc(1,1,nlevnc),'VENT NCEP   ')
      print*,'WARNING!!! Correction bidon pour palier a un '
      print*,'probleme dans la creation des fichiers nc'
      call correctbid(iim,jjp1*nlevnc,unc)
      call dump2d(iip1,jjp1,unc,'UNC COUCHE 1 ')
      endif

c  Temperature
c  -----------

      print*,'ncidt=',ncidt,'varidt=',varidt,'start=',start
      print*,'count=',count
      if (guide_T) then
      status=NF_GET_VARA_REAL(ncidt,varidt,start,count,tnc)
      call dump2d(iip1,jjp1,tnc,'TNC COUCHE 1 AAA ')
      call correctbid(iim,jjp1*nlevnc,tnc)
      call dump2d(iip1,jjp1,tnc,'TNC COUCHE 1 BBB ')
      endif

c  Humidite
c  --------

      if (guide_Q) then
      status=NF_GET_VARA_REAL(ncidQ,varidQ,start,count,Qnc)
      call correctbid(iim,jjp1*nlevnc,Qnc)
      call dump2d(iip1,jjp1,Qnc,'QNC COUCHE 1 ')
      endif

      count(2)=jjm
c  Vent meridien
c  -------------

      if (guide_v) then
      status=NF_GET_VARA_REAL(ncidv,varidv,start,count,vnc)
      call correctbid(iim,jjm*nlevnc,vnc)
      call dump2d(iip1,jjm,vnc,'VNC COUCHE 1 ')
      endif

      start(3)=timestep
      start(4)=0
      count(2)=jjp1
      count(3)=1
      count(4)=0

c  Pression de surface
c  -------------------

      if (guide_P) then
      status=NF_GET_VARA_REAL(ncidps,varidps,start,count,psnc)
      call dump2d(iip1,jjp1,psnc,'PSNC COUCHE 1 ')
      call correctbid(iim,jjp1,psnc)
      endif



c -----------------------------------------------------------------
c  Interpollation verticale sur les niveaux modele
c -----------------------------------------------------------------
      call reanalyse2nat(nlevnc,psi,unc,vnc,tnc,Qnc,psnc,pl,u,v,t,Q
     s    ,ps,masse,pk)

      call dump2d(iip1,jjm,v,'V COUCHE APRES ')


c -----------------------------------------------------------------
c  Passage aux variables du modele (vents covariants, temperature
c  potentielle et humidite specifique)
c -----------------------------------------------------------------
      call nat2gcm(u,v,t,Q,pk,u,v,t,Q)
      print*,'TIMESTEP ',timestep
      if(mode.ne.1) stop'mode pas egal 0'
c     call dump2d(iip1,jjm,v,'VCOV COUCHE 1 ')

c   Lignes introduites a une epoque pour un probleme oublie...
c     do l=1,llm
c        do i=1,iip1
c           v(i,1,l)=0.
c           v(i,jjm,l)=0.
c        enddo
c     enddo
      first=.false.

      return
      end


c===========================================================================
      subroutine reanalyse2nat(nlevnc,psi
     s   ,unc,vnc,tnc,qnc,psnc,pl,u,v,t,q
     s   ,ps,masse,pk)
c===========================================================================

c -----------------------------------------------------------------
c   Inversion Nord/sud de la grille + interpollation sur les niveaux
c   verticaux du modele.
c -----------------------------------------------------------------

      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use comgeom
      use exner_hyb_m, only: exner_hyb
      use guide_m
      use pression_m, only: pression

      implicit none


      integer nlevnc
      real psi(iip1,jjp1)
      real u(iip1,jjp1,llm),v(iip1,jjm,llm)
      real t(iip1,jjp1,llm),ps(iip1,jjp1),q(iip1,jjp1,llm)

      real pl(nlevnc)
      real unc(iip1,jjp1,nlevnc),vnc(iip1,jjm,nlevnc)
      real tnc(iip1,jjp1,nlevnc),psnc(iip1,jjp1)
      real qnc(iip1,jjp1,nlevnc)

      real zu(iip1,jjp1,llm),zv(iip1,jjm,llm)
      real zt(iip1,jjp1,llm),zq(iip1,jjp1,llm)

      real pext(iip1,jjp1,llm)
      real pbarx(iip1,jjp1,llm),pbary(iip1,jjm,llm)
      real plunc(iip1,jjp1,llm),plvnc(iip1,jjm,llm)
      real plsnc(iip1,jjp1,llm)

      real p(iip1,jjp1,llmp1),pk(iip1,jjp1,llm),pks(iip1,jjp1)
      real pkf(iip1,jjp1,llm)
      real masse(iip1,jjp1,llm),pls(iip1,jjp1,llm)
      real prefkap,unskap


      integer i,j,l


c -----------------------------------------------------------------
c   calcul de la pression au milieu des couches.
c -----------------------------------------------------------------

      CALL pression( ip1jmp1, ap, bp, psi, p )
      call massdair(p,masse)
      CALL exner_hyb(psi,p,pks,pk,pkf)

c    ....  Calcul de pls , pression au milieu des couches ,en Pascals
      unskap=1./kappa
      prefkap =  preff  ** kappa
c     PRINT *,' Pref kappa unskap  ',preff,kappa,unskap
      DO l = 1, llm
       DO j=1,jjp1
        DO i =1, iip1
        pls(i,j,l) = preff * ( pk(i,j,l)/cpp) ** unskap
        ENDDO
       ENDDO
       ENDDO


c -----------------------------------------------------------------
c   calcul des pressions pour les grilles u et v
c -----------------------------------------------------------------

      do l=1,llm
      do j=1,jjp1
         do i=1,iip1
            pext(i,j,l)=pls(i,j,l)*aire_2d(i,j)
         enddo
      enddo
      enddo
      call massbar(pext, pbarx, pbary )
      do l=1,llm
      do j=1,jjp1
         do i=1,iip1
            plunc(i,jjp1+1-j,l)=pbarx(i,j,l)/aireu_2d(i,j)
            plsnc(i,jjp1+1-j,l)=pls(i,j,l)
         enddo
      enddo
      enddo
      do l=1,llm
      do j=1,jjm
         do i=1,iip1
            plvnc(i,jjm+1-j,l)=pbary(i,j,l)/airev_2d(i,j)
         enddo
      enddo
      enddo

c -----------------------------------------------------------------

      if (guide_P) then
      do j=1,jjp1
         do i=1,iim
            ps(i,j)=psnc(i,jjp1+1-j)
         enddo
         ps(iip1,j)=ps(1,j)
      enddo
      endif


c -----------------------------------------------------------------
      call pres2lev(unc,zu,nlevnc,llm,pl,plunc,iip1,jjp1)
      call pres2lev(vnc,zv,nlevnc,llm,pl,plvnc,iip1,jjm )
      call pres2lev(tnc,zt,nlevnc,llm,pl,plsnc,iip1,jjp1)
      call pres2lev(qnc,zq,nlevnc,llm,pl,plsnc,iip1,jjp1)

c     call dump2d(iip1,jjp1,ps,'PS    ')
c     call dump2d(iip1,jjp1,psu,'PS    ')
c     call dump2d(iip1,jjm,psv,'PS    ')
c  Inversion Nord/Sud
      do l=1,llm
         do j=1,jjp1
            do i=1,iim
               u(i,j,l)=zu(i,jjp1+1-j,l)
               t(i,j,l)=zt(i,jjp1+1-j,l)
               q(i,j,l)=zq(i,jjp1+1-j,l)
            enddo
            u(iip1,j,l)=u(1,j,l)
            t(iip1,j,l)=t(1,j,l)
            q(iip1,j,l)=q(1,j,l)
         enddo
      enddo

      do l=1,llm
         do j=1,jjm
            do i=1,iim
               v(i,j,l)=zv(i,jjm+1-j,l)
            enddo
            v(iip1,j,l)=v(1,j,l)
         enddo
      enddo

      return
      end

c===========================================================================
      subroutine nat2gcm(u,v,t,rh,pk,ucov,vcov,teta,q)
c===========================================================================

      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use comgeom
      use q_sat_m, only: q_sat
      use guide_m
      implicit none


      real u(iip1,jjp1,llm),v(iip1,jjm,llm)
      real t(iip1,jjp1,llm),pk(iip1,jjp1,llm),rh(iip1,jjp1,llm)
      real ps(iip1,jjp1)

      real ucov(iip1,jjp1,llm),vcov(iip1,jjm,llm)
      real teta(iip1,jjp1,llm),q(iip1,jjp1,llm)

      real pres(iip1,jjp1,llm),qsat(iip1,jjp1,llm)

      real unskap

      integer i,j,l


      print*,'Entree dans nat2gcm'
c    ucov(:,:,:)=0.
c    do l=1,llm
c       ucov(:,2:jjm,l)=u(:,2:jjm,l)*cu_2d(:,2:jjm)
c    enddo
c    ucov(iip1,:,:)=ucov(1,:,:)

c    teta(:,:,:)=t(:,:,:)*cpp/pk(:,:,:)
c    teta(iip1,:,:)=teta(1,:,:)
     
c   calcul de ucov et de la temperature potentielle
      do l=1,llm
         do j=1,jjp1
            do i=1,iim
               ucov(i,j,l)=u(i,j,l)*cu_2d(i,j)
               teta(i,j,l)=t(i,j,l)*cpp/pk(i,j,l)
            enddo
            ucov(iip1,j,l)=ucov(1,j,l)
            teta(iip1,j,l)=teta(1,j,l)
         enddo
         do i=1,iip1
            ucov(i,1,l)=0.
            ucov(i,jjp1,l)=0.
            teta(i,1,l)=teta(1,1,l)
            teta(i,jjp1,l)=teta(1,jjp1,l)
         enddo
      enddo

c   calcul de ucov
      do l=1,llm
         do j=1,jjm
            do i=1,iim
               vcov(i,j,l)=v(i,j,l)*cv_2d(i,j)
            enddo
            vcov(iip1,j,l)=vcov(1,j,l)
         enddo
      enddo

c     call dump2d(iip1,jjp1,teta,'TETA EN BAS   ')
c     call dump2d(iip1,jjp1,teta(1,1,llm),'TETA EN HAUT   ')

c  Humidite relative -> specifique
c  -------------------------------
      if (1.eq.0) then
c   FINALEMENT ON GUIDE EN HUMIDITE RELATIVE
      print*,'calcul de unskap'
      unskap   = 1./ kappa
      print*,'calcul de pres'
      pres(:,:,:)=preff*(pk(:,:,:)/cpp)**unskap
      print*,'calcul de qsat'
      qsat = q_sat(t, pres)
      print*,'calcul de q'
c   ATTENTION : humidites relatives en %
      rh(:,:,:)=max(rh(:,:,:)*0.01,1.e-6)
      q(:,:,:)=qsat(:,:,:)*rh(:,:,:)
      print*,'calcul de q OK'

      call dump2d(iip1,jjp1,pres,'PRESSION PREMIERE COUCHE   ')
      call dump2d(iip1,jjp1,q,'HUMIDITE SPECIFIQUE COUCHE 1   ') 
      endif


      return
      end



c===========================================================================
      subroutine correctbid(iim,nl,x)
c===========================================================================
      integer iim,nl
      real x(iim+1,nl)
      integer i,l
      real zz

      do l=1,nl
         do i=2,iim-1
            if(abs(x(i,l)).gt.1.e10) then
               zz=0.5*(x(i-1,l)+x(i+1,l))
c              print*,'correction ',i,l,x(i,l),zz
               x(i,l)=zz
            endif
         enddo
      enddo
      return
      end
