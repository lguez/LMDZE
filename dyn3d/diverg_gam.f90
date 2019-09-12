module diverg_gam_m

  IMPLICIT NONE

contains

  SUBROUTINE diverg_gam(klevel, cuvscvgam, cvuscugam, unsairegam, unsapolnga, &
       unsapolsga, x, y, div)

    ! From LMDZ4/libf/dyn3d/diverg_gam.F, version 1.1.1.1 2004/05/19
    ! 12:53:05

    ! Author: P. Le Van

    ! Calcule la divergence \`a tous les niveaux d'un vecteur de
    ! composantes covariantes x et y.

    USE dimensions, only: iim
    USE paramet_m, only: ip1jmp1, ip1jm, iip1, iip2, ip1jmi1

    ! div      est  un argument  de sortie pour le s-prog

    ! ATTENTION : pendant ce s-pg , ne pas toucher au COMMON/scratch/  .

    ! ---------------------------------------------------------------------

    ! ..........          variables en arguments    ...................

    INTEGER, INTENT (IN) :: klevel
    REAL cuvscvgam(ip1jm), cvuscugam(ip1jmp1), unsairegam(ip1jmp1)
    REAL unsapolnga, unsapolsga
    REAL, intent(in):: x(ip1jmp1, klevel), y(ip1jm, klevel)
    real div(ip1jmp1, klevel)

    ! ...............     variables  locales   .........................

    REAL aiy1(iip1), aiy2(iip1)
    REAL sumypn, sumyps
    INTEGER l, ij

    ! ...................................................................

    DO l = 1, klevel

       DO ij = iip2, ip1jm - 1
          div(ij+1, l) = (cvuscugam(ij+1)*x(ij+1,l)-cvuscugam(ij)*x(ij,l)+ &
               cuvscvgam(ij-iim)*y(ij-iim,l)-cuvscvgam(ij+1)*y(ij+1,l))* &
               unsairegam(ij+1)
       END DO

       ! ....  correction pour  div( 1,j,l)  ......
       ! ....   div(1,j,l)= div(iip1,j,l) ....

       ! DIR$ IVDEP
       DO ij = iip2, ip1jm, iip1
          div(ij, l) = div(ij+iim, l)
       END DO

       ! ....  calcul  aux poles  .....

       DO ij = 1, iim
          aiy1(ij) = cuvscvgam(ij)*y(ij, l)
          aiy2(ij) = cuvscvgam(ij+ip1jmi1)*y(ij+ip1jmi1, l)
       END DO
       sumypn = sum(aiy1(:iim)) * unsapolnga
       sumyps = sum(aiy2(:iim)) * unsapolsga

       DO ij = 1, iip1
          div(ij, l) = -sumypn
          div(ij+ip1jm, l) = sumyps
       END DO
    END DO

  END SUBROUTINE diverg_gam

end module diverg_gam_m
