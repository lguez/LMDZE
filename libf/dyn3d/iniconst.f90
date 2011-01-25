SUBROUTINE iniconst

  ! From dyn3d/iniconst.F,v 1.1.1.1 2004/05/19 12:53:05
  ! P. Le Van

  USE dimens_m, ONLY : iim, jjm, llm
  USE comconst, ONLY : cpp, dtphys, dtvr, im, imp1, jm, kappa, lllm, lllmm1, &
       lllmp1, r, unsim
  USE comvert, ONLY : disvert
  USE conf_gcm_m, ONLY : iphysiq

  IMPLICIT NONE

  !-----------------------------------------------------------------------

  !   dimension des boucles:
  im      = iim
  jm      = jjm
  lllm    = llm
  imp1    = iim 
  lllmm1  = llm - 1
  lllmp1  = llm + 1

  dtphys  = iphysiq * dtvr
  unsim   = 1./iim
  r       = cpp * kappa
  PRINT *, 'cpp = ', cpp
  PRINT *, 'R = ', r
  PRINT *, 'kappa = ', kappa
  CALL disvert

END SUBROUTINE iniconst
