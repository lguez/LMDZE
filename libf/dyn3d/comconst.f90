module comconst

  use dimens_m, only: jjm

  implicit none

  INTEGER im, jm, lllm, imp1
  integer, parameter:: jmp1 = jjm + 1
  integer lllmm1, lllmp1, lcl
  REAL dtvr ! time step for dynamics (in s)
  real, parameter:: daysec = 86400. ! number of seconds per day
  REAL pi, dtphys, dtdiss
  real, parameter:: rad = 6371229. ! radius of the Earth (in m)
  real r
  real, parameter:: cpp = 1004.70885, kappa = 0.2857143
  REAL cotot, unsim
  real, parameter:: g = 9.8 ! acceleration of gravity (in m s-2)
  real omeg ! angular speed of rotation of the Earth (in rad s-1)

  private jjm

contains

  subroutine initialize

    print *, "Call sequence information: initialize"
    pi     =  Acos(-1.)
    omeg   = 2 * pi / daysec

  end subroutine initialize

end module comconst
