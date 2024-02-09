module init_iophy_m

  implicit none

contains

  subroutine init_iophy

    ! Libraries:
    use xios, only: xios_duration, xios_set_timestep, &
         xios_close_context_definition, xios_set_time_origin, xios_date, &
         xios_set_start_date, xios_define_calendar

    use conf_gcm_m, only: dtphys
    use dynetat0_m, only: day_ini
    use dynetat0_chosen_m, only: annee_ref, day_ref

    ! Local:
    TYPE(xios_duration) dtime

    !----------------------------------------------------------------------

    CALL xios_define_calendar(type = "D360")
    dtime%second = dtphys
    call xios_set_timestep(dtime)
    call xios_set_time_origin(xios_date(annee_ref, 1,day_ref, 0, 0, 0))
    CALL xios_set_start_date(xios_date(annee_ref, 1, day_ini, 0, 0, 0))
    call xios_close_context_definition

  end subroutine init_iophy

end module init_iophy_m
