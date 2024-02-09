module init_iophy_m

  implicit none

contains

  subroutine init_iophy

    ! Libraries:
    use jumble, only: rad_to_deg
    use xios, only: xios_set_axis_attr, xios_set_domain_attr, xios_duration, &
         xios_set_timestep, xios_close_context_definition, &
         xios_set_time_origin, xios_date, xios_set_start_date, &
         xios_define_calendar

    use conf_gcm_m, only: dtphys
    use dimensions, only: iim, jjm, llm
    use disvert_m, only: presnivs
    use dynetat0_m, only: rlatu, rlonv, day_ini
    use dynetat0_chosen_m, only: annee_ref, day_ref

    ! Local:
    TYPE(xios_duration) dtime

    !----------------------------------------------------------------------

    CALL xios_define_calendar(type = "D360")
    dtime%second = dtphys
    call xios_set_timestep(dtime)
    call xios_set_time_origin(xios_date(annee_ref, 1,day_ref, 0, 0, 0))
    CALL xios_set_start_date(xios_date(annee_ref, 1, day_ini, 0, 0, 0))
    CALL xios_set_domain_attr("dom_glo", type = "rectilinear", &
         ni_glo = iim, nj_glo = jjm + 1, &
         lonvalue_1d = dble(rlonv(:iim) * rad_to_deg), &
         latvalue_1d = dble(rlatu * rad_to_deg))
    CALL xios_set_axis_attr("presnivs", n_glo = llm, &
         value = dble(presnivs / 100.))
    call xios_close_context_definition

  end subroutine init_iophy

end module init_iophy_m
