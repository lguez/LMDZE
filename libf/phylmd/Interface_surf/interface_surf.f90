MODULE interface_surf

  ! From phylmd/interface_surf.F90, version 1.8 2005/05/25 13:10:09
  ! L. Fairhead, LMD, february 2000

  IMPLICIT none

  ! run_off ruissellement total
  REAL, ALLOCATABLE, DIMENSION(:), SAVE :: run_off, run_off_lic
  real, allocatable, dimension(:), save :: coastalflow, riverflow

  REAL, ALLOCATABLE, DIMENSION(:, :), SAVE :: tmp_rriv, tmp_rcoa, tmp_rlic
  ! pour simuler la fonte des glaciers antarctiques
  REAL, save :: tau_calv

END MODULE interface_surf
