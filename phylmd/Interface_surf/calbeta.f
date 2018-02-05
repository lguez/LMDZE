module calbeta_m

  IMPLICIT NONE

contains

  SUBROUTINE calbeta(indice, snow, qsol, beta, vcal, vdif)

    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: April 14th, 1994

    ! Calcul de quelques paramètres pour appliquer la couche limite.

    USE indicesol, ONLY: is_lic, is_oce, is_sic, is_ter

    INTEGER, intent(in):: indice
    REAL, intent(in):: snow(:)

    REAL, intent(in):: qsol(:) ! (knon)
    ! column-density of water in soil, in kg m-2

    REAL, intent(out):: beta(:), vcal(:), vdif(:) ! (knon)

    ! Local:

    REAL, PARAMETER:: tau_gl = 86400. * 5. 
    ! temps de relaxation pour la glace de mer

    REAL, PARAMETER:: max_eau_sol = 150. ! in kg m-2

    ! Pour une épaisseur du sol de 15 cm :
    REAL, PARAMETER:: calsol = 1. / (2.5578E6 * 0.15)
    REAL, PARAMETER:: calsno = 1. / (2.3867E6 * 0.15)
    REAL, PARAMETER:: calice = 1. / (5.1444E6 * 0.15)

    !------------------------------------------------------------

    select case (indice)
    case(is_oce)
       beta = 1.
       vcal = 0.
       vdif = 0.

    case (is_sic)
       beta = 1.
       vcal = merge(calsno, calice, snow > 0.)
       vdif = 1. / tau_gl

    case (is_ter)
       beta = min(2. * qsol / max_eau_sol, 1.)
       vcal = merge(calsno, calsol, snow > 0.)
       vdif = 0.

    case (is_lic)
       beta = 1.
       vcal = merge(calsno, calice, snow > 0.)
       vdif = 0.
    END select

  END SUBROUTINE calbeta

end module calbeta_m
