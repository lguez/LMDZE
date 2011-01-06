program test_disvert

  ! This is a program in Fortran 95.
  ! This program writes the vertical coordinates at level interfaces.

  use dimens_m, only: llm
  use comvert, only: pa, disvert, ap, bp, preff
  use new_unit_m, only: new_unit

  implicit none

  integer unit, l

  REAL s(llm+1) ! hybrid sigma-pressure coordinate
  ! "s(l)" : à l'interface des couches "l" et "l-1"

  real p(llm + 1) ! pressure at level, in Pa
  real z(llm) ! pressure-altitude (km)

  !---------------------

  print *, "llm = ", llm
  pa = 5e4
  call disvert
  s = ap / pa + bp
  p = ap + bp * preff
  z = - 7. * log(ap(:llm) / preff + bp(:llm))

  call new_unit(unit)
  open(unit, file="test_disvert.csv", status="replace", action="write")
  ! Title line:
  write(unit, fmt=*) '"ap (hPa)" "bp" "s" "p (hPa)" "z (km)"'
  do l = 1, llm
     write(unit, fmt=*) ap(l) / 100., bp(l), s(l), p(l) / 100., z(l)
  end do
  close(unit)
  print *, 'The file "test_disvert.csv" has been created.'

  open(unit, file="delta_z.csv", status="replace", action="write")
  ! Title line:
  write(unit, fmt=*) '"z middle (km)" "delta z (km)"'
  do l = 1, llm - 1
     write(unit, fmt=*) (z(l) + z(l+1)) / 2, z(l+1) - z(l)
  end do
  close(unit)
  print *, 'The file "delta_z.csv" has been created.'

end program test_disvert
