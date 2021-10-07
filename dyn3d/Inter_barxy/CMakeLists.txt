target_sources(LMDZE_objects PRIVATE inter_barx.f90 inter_barxy.f90
  inter_bary.f90 ord_coord.f90 ord_coordm.f90)

add_executable(test_inter_barxy test_inter_barxy.f90)
target_link_libraries(test_inter_barxy PRIVATE LMDZE_objects NR_util::nr_util)
