module pressure_var

  use dimens_m, only: iim, jjm, llm

  implicit none

  real pls(iim + 1, jjm + 1, llm)
  ! (pressure at mid-layer of LMDZ grid, in Pa)
  ! "pls(i, j, l)" is at longitude "rlonv(i)", latitude "rlatu(j)",
  ! for layer "l")

  REAL p3d(iim + 1, jjm + 1, llm+1) ! pressure at layer interfaces, in Pa
  ! ("p3d(i, j, l)" is at longitude "rlonv(i)", latitude "rlatu(j)",
  ! for interface "l")

  private iim, jjm, llm

end module pressure_var
