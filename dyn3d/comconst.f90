module comconst

  use nr_util, only: twopi

  implicit none

  real, parameter:: daysec = 86400. ! number of seconds per day

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

end module comconst
