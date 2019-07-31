module albsno_m

  ! From phylmd/interface_surf.F90, version 1.8 2005/05/25 13:10:09

  IMPLICIT none

contains

  SUBROUTINE albsno(agesno, alb_neig, snow_fall)

    use comconst, only: dtphys

    REAL, intent(inout):: agesno(:) ! (knon) age of snow, in days
    real, intent(out):: alb_neig(:) ! (knon)
    real, intent(in):: snow_fall(:) !(knon)

    !------------------------------------------------------------------------

    ! D\'esert partout:
    alb_neig = 0.55 + 0.3 * EXP(- agesno / 5.)

    ! Modulation en fonction de l'\^age de la neige :
    agesno = (agesno + (1. - agesno / 50.) * dtphys / 86400.) &
         * EXP(- snow_fall * dtphys / 0.3)

  END SUBROUTINE albsno

end module albsno_m
