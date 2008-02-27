!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/YOEGWD.h,v 1.1.1.1 2004/05/19 12:53:08 lmdzadmin Exp $
!
C     -----------------------------------------------------------------
C*    *COMMON* *YOEGWD* - PARAMETERS FOR GRAVITY WAVE DRAG CALCULATIONS
C     -----------------------------------------------------------------
C
      integer NKTOPG,NSTRA
      real GFRCRIT,GKWAKE,GRCRIT,GVCRIT,GKDRAG,GKLIFT
      real GHMAX,GRAHILO,GSIGCR,GSSEC,GTSEC,GVSEC
      COMMON/YOEGWD/ GFRCRIT,GKWAKE,GRCRIT,GVCRIT,GKDRAG,GKLIFT
     *        ,GHMAX,GRAHILO,GSIGCR,NKTOPG,NSTRA,GSSEC,GTSEC,GVSEC
C


