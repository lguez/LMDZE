module comconst

  use nr_util, only: twopi

  implicit none

  real, parameter:: daysec = 86400. ! number of seconds per day

  REAL, protected:: dtvr ! time step for dynamics, in s
  REAL, protected:: dtphys ! time step for physics, in s

  real, parameter:: rad = 6371229. ! radius of the Earth (in m)

  real, parameter:: cpp = 1004.70885
  ! specific heat capacity at constant pressure of dry air, in J K-1 kg-1

  real, parameter:: kappa = 0.2857143 ! r / cpp

  real, parameter:: r = cpp * kappa
  ! specific ideal gas constant for dry air, in J K-1 kg-1

  real, parameter:: g = 9.8 ! acceleration of gravity (in m s-2)

  real, parameter:: omeg = twopi / daysec
  ! angular speed of rotation of the Earth (in rad s-1)

  private twopi

contains

  SUBROUTINE iniconst

    ! From dyn3d/iniconst.F,v 1.1.1.1 2004/05/19 12:53:05
    ! P. Le Van

    USE conf_gcm_m, ONLY: day_step, iphysiq

    IMPLICIT NONE

    !-----------------------------------------------------------------------

    dtvr = daysec / real(day_step)
    dtphys  = iphysiq * dtvr

    print *, 'dtvr = ', dtvr
    print *, 'dtphys = ', dtphys
    PRINT *, 'cpp = ', cpp
    PRINT *, 'R = ', r
    PRINT *, 'kappa = ', kappa

  END SUBROUTINE iniconst

end module comconst
