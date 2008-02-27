!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/iniphysiq.F,v 1.1.1.1 2004/05/19 12:53:08 lmdzadmin Exp $
!
c
c
      SUBROUTINE iniphysiq(ngrid,nlayer,
     $           punjours,
     $           pdayref,ptimestep,
     $           plat,plon,parea,pcu,pcv,
     $           prad,pg,pr,pcpp)
      use dimens_m
      use dimphy
      use comgeomphy
      use suphec_m, only: suphec

      IMPLICIT NONE
c
c=======================================================================
c
c   subject:
c   --------
c
c   Initialisation for the physical parametrisations of the LMD 
c   martian atmospheric general circulation modele.
c
c   author: Frederic Hourdin 15 / 10 /93
c   -------
c
c   arguments:
c   ----------
c
c   input:
c   ------
c
c    ngrid                 Size of the horizontal grid.
c                          All internal loops are performed on that grid.
c    nlayer                Number of vertical layers.
c    pdayref               Day of reference for the simulation
c    firstcall             True at the first call
c    lastcall              True at the last call
c    pday                  Number of days counted from the North. Spring
c                          equinoxe.
c
c=======================================================================
c
c-----------------------------------------------------------------------
c   declarations:
c   -------------
 

      REAL prad,pg,pr,pcpp,punjours
 
      INTEGER ngrid,nlayer
      REAL plat(ngrid),plon(ngrid),parea(klon),pcu(klon),pcv(klon)
      INTEGER pdayref
 
      REAL ptimestep
 
      IF (nlayer.NE.klev) THEN
         PRINT*,'STOP in inifis'
         PRINT*,'Probleme de dimensions :'
         PRINT*,'nlayer     = ',nlayer
         PRINT*,'klev   = ',klev
         STOP
      ENDIF

      IF (ngrid.NE.klon) THEN
         PRINT*,'STOP in inifis'
         PRINT*,'Probleme de dimensions :'
         PRINT*,'ngrid     = ',ngrid
         PRINT*,'klon   = ',klon
         STOP
      ENDIF

      airephy=parea
      cuphy=pcu
      cvphy=pcv
      rlond = plon
      rlatd = plat

      call suphec
      print*,'ATTENTION !!! TRAVAILLER SUR INIPHYSIQ'
      print*,'CONTROLE DES LATITUDES, LONGITUDES, PARAMETRES ...'

      END
