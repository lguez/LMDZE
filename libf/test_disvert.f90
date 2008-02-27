program test_disvert

  use dimens_m, only: llm
  use comvert, only: pa, disvert, ap, bp, preff
  use comconst, only: initialize
  use in_out, only: new_unit

  implicit none

  integer unit, l
  REAL sigma(llm+1) ! "sigma(l)" : à l'interface des couches "l" et "l-1"
  real p(llm + 1) ! pressure at level, in Pa

  !---------------------

  print *, "llm = ", llm
  call initialize
  pa = 5e4
  call disvert
  sigma = ap / pa + bp
  p = ap + bp * preff

  unit = new_unit()
  open(unit, file="test_disvert.csv", status="replace", action="write")
  write(unit, fmt=*) '"ap (Pa)" "bp" "sigma" "p (Pa)"' ! title line
  do l = 1, llm
     write(unit, fmt=*) ap(l), bp(l), sigma(l), p(l)
  end do
  close(unit)
  print *, 'The file "test_disvert.csv" has been created.'

end program test_disvert
