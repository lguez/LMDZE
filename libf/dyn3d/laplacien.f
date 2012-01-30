!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/laplacien.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE laplacien ( klevel, teta, divgra )
c
c     P. Le Van
c
c   ************************************************************
c    ....     calcul de  (div( grad ))   de   teta  .....
c   ************************************************************
c     klevel et teta  sont des arguments  d'entree pour le s-prog
c      divgra     est  un argument  de sortie pour le s-prog
c
      use grad_m, only: grad
      use dimens_m
      use paramet_m
      use comgeom
      use filtreg_m, only: filtreg
      use divergf_m, only: divergf

      IMPLICIT NONE
c

c
c    .........      variables  en arguments   ..............
c
      INTEGER, intent(in):: klevel
      REAL teta( ip1jmp1,klevel ), divgra( ip1jmp1,klevel )
c
c    ............     variables  locales      ..............
c
      REAL ghy(ip1jm,llm), ghx(ip1jmp1,llm)
c    .......................................................


c
      CALL SCOPY ( ip1jmp1 * klevel, teta, 1, divgra, 1 )

      CALL filtreg( divgra,  jjp1, klevel,  2, 1, .TRUE., 1 )
      CALL   grad ( klevel,divgra,   ghx , ghy              )
      CALL  divergf ( klevel, ghx , ghy  , divgra           )

      RETURN
      END
