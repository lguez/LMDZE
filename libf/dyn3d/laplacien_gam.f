!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/laplacien_gam.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE laplacien_gam ( klevel, cuvsga, cvusga, unsaigam ,
     *                        unsapolnga, unsapolsga, teta, divgra )

c  P. Le Van
c
c   ************************************************************
c
c      ....   calcul de  (div( grad ))   de   teta  .....
c   ************************************************************
c    klevel et teta  sont des arguments  d'entree pour le s-prog
c      divgra     est  un argument  de sortie pour le s-prog
c
      use grad_m, only: grad
      use dimens_m
      use paramet_m
      use comgeom
      IMPLICIT NONE
c

c
c    ............     variables  en arguments    ..........
c
      INTEGER, intent(in):: klevel
      REAL teta( ip1jmp1,klevel ), divgra( ip1jmp1,klevel )
      REAL cuvsga(ip1jm) , cvusga( ip1jmp1 ),unsaigam(ip1jmp1),
     *     unsapolnga, unsapolsga
c
c    ...........    variables  locales    .................
c
      REAL ghy(ip1jm,llm), ghx(ip1jmp1,llm)
c    ......................................................

      CALL SCOPY ( ip1jmp1 * klevel, teta, 1, divgra, 1 )
c
      CALL   grad ( klevel, divgra, ghx, ghy )
c
      CALL  diverg_gam ( klevel, cuvsga, cvusga,  unsaigam  ,
     *                 unsapolnga, unsapolsga, ghx , ghy , divgra )

c


      RETURN
      END
