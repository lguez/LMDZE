SUBROUTINE iniconst

  ! From dyn3d/iniconst.F,v 1.1.1.1 2004/05/19 12:53:05
  ! P. Le Van

  use dimens_m
  use paramet_m
  use comconst, only: cpp, dtdiss, dtphys, dtvr, im, imp1, jm, kappa, lllm, &
       lllmm1, lllmp1, pi, r, unsim
  use comvert, only: disvert
  use conf_gcm_m, only: idissip, iphysiq
  use temps

  IMPLICIT NONE

  !-----------------------------------------------------------------------

  !   dimension des boucles:
  im      = iim
  jm      = jjm
  lllm    = llm
  imp1    = iim 
  lllmm1  = llm - 1
  lllmp1  = llm + 1

  dtdiss  = idissip * dtvr
  dtphys  = iphysiq * dtvr
  unsim   = 1./iim
  pi      = 2.*ASIN( 1. )
  r       = cpp * kappa
  PRINT *, 'cpp = ', cpp
  PRINT *, 'R = ', r
  PRINT *, 'kappa = ', kappa
  CALL disvert

END SUBROUTINE iniconst
