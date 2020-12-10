target_sources(LMDZE_objects PRIVATE conf_guide.f90 coordij.f90
  guide.f90 init_tau2alpha.f90 tau2alpha.f90)

add_subdirectory(Read_reanalyse)
