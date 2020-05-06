module albsno_m

  IMPLICIT none

contains

  SUBROUTINE albsno(agesno, alb_neig, snow_fall)

    ! From phylmd/interface_surf.F90, version 1.8 2005/05/25 13:10:09

    ! Calcul \^age de la neige.

    use comconst, only: dtphys

    REAL, intent(inout):: agesno(:) ! (knon) age of snow, in days
    real, intent(out):: alb_neig(:) ! (knon)

    real, intent(in):: snow_fall(:) !(knon)
    ! precipitation, solid water mass flux (kg / m2 / s), positive down

    !------------------------------------------------------------------------

    ! D\'esert partout:
    alb_neig = 0.55 + 0.3 * EXP(- agesno / 5.)

    agesno = (agesno + (1. - agesno / 50.) * dtphys / 86400.) &
         * EXP(- snow_fall * dtphys / 0.3)

  END SUBROUTINE albsno

end module albsno_m
