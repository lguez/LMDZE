module radepsi

  ! From phylmd/radepsi.h, version 1.1.1.1 2004/05/19 12:53:08

  implicit none

  double precision, parameter:: ZEELOG = 1.D-07
  ! (1.d-10 (not good for 32-bit machines))

  double precision, parameter:: ZEPSC  = 1.D-20
  double precision, parameter:: ZEPSCO = 1.D-10
  double precision, parameter:: ZEPSCQ = 1.D-10
  double precision, parameter:: ZEPSCT = 1.D-20
  double precision, parameter:: ZEPSCW = 1.D-20
  double precision, parameter:: ZEPSEC = 1.0d-12
  double precision, parameter:: ZEPSCR = 1.0d-10

  double precision, parameter:: REPSCT=1.0E-10

end module radepsi
