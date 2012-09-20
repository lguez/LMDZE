!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/laplacien_rot.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
      SUBROUTINE laplacien_rot ( klevel, rotin, rotout,ghx,ghy )
c
c    P. Le Van
c
c   ************************************************************
c    ...  calcul de  ( rotat x nxgrad )  du rotationnel rotin  .
c   ************************************************************
c
c     klevel et rotin  sont des arguments  d'entree pour le s-prog
c      rotout           est  un argument  de sortie pour le s-prog
c
      use dimens_m
      use paramet_m
      use comgeom
      use filtreg_m, only: filtreg
      IMPLICIT NONE
c

c 
c   ..........    variables  en  arguments     .............
c
      INTEGER, intent(in):: klevel
      REAL rotin( ip1jm,klevel ), rotout( ip1jm,klevel )
c
c   ..........    variables   locales       ................
c
      REAL ghy(ip1jm,klevel), ghx(ip1jmp1,klevel)
c   ........................................................
c
c
      CALL  filtreg ( rotin ,   jjm, klevel,   2, 1, .FALSE.)

      CALL   nxgrad ( klevel, rotin,   ghx ,  ghy               )
      CALL   rotatf  ( klevel, ghx  ,   ghy , rotout             )
c
      RETURN
      END
