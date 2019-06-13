program test_orbite

  use orbite_m, only: orbite

  implicit none

  REAL xjour ! jour de l'ann\'ee \`a compter du premier janvier

  REAL longi
  ! longitude vraie de la Terre dans son orbite solaire, par rapport
  ! au point vernal (21 mars), en degr\'es

  REAL dist ! distance terre-soleil, en ua
  integer i
  integer, parameter:: n = 100
  real delta_x

  !------------------------------------------------------------------------

  xjour = 1.
  delta_x = 359. / (n - 1)
  print *, '"" degrees ua'
  print *, "xjour longi dist"

  do i = 1, n
     call orbite(xjour, longi, dist)
     print *, xjour, longi, dist
     xjour = xjour + delta_x
  end do

end program test_orbite
