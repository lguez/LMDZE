SUBROUTINE minmax(imax, xi, zmin, zmax)

  ! From dyn3d/minmax.F, version 1.1.1.1 2004/05/19 12:53:07
  ! P. Le Van

  INTEGER, intent(in):: imax
  REAL, intent(in):: xi(imax)
  REAL, intent(out):: zmin, zmax

  !------------------------------

  zmin = minval(xi)
  zmax = maxval(xi)

END SUBROUTINE minmax
