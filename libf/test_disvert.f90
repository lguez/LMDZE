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

  !---------------------

  print *, "llm = ", llm
  pa = 5e4
  call disvert
  s = ap / pa + bp
  p = ap + bp * preff

  call new_unit(unit)
  open(unit, file="test_disvert.csv", status="replace", action="write")
  write(unit, fmt=*) '"ap (hPa)" "bp" "s" "p (hPa)"' ! title line
  do l = 1, llm
     write(unit, fmt=*) ap(l) / 100., bp(l), s(l), p(l) / 100.
  end do
  close(unit)
  print *, 'The file "test_disvert.csv" has been created.'

end program test_disvert
