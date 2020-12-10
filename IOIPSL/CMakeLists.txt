target_sources(LMDZE_objects PRIVATE errioipsl.f90 find_str.f90
  histcom_var.f90 histvar_seq.f90 histwrite.f90 histwrite_real.f90
  ioget_calendar.f90 mathop2.f90 mathop.f90 strlowercase.f90)

foreach(my_dir IN ITEMS Mathelp Histcom Calendar)
  add_subdirectory(${my_dir})
endforeach(my_dir)
