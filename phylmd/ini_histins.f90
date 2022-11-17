module ini_histins_m

  implicit none

  integer, save:: nid_ins

contains

  subroutine ini_histins

    ! From phylmd/ini_histins.h, version 1.2, 2005/05/25 13:10:09

    ! Libraries:
    use jumble, only: pi

    use clesphys, only: ecrit_ins, ok_instan
    use clesphys2, only: conv_emanuel
    use conf_gcm_m, only: dtphys
    use dimensions, only: iim, jjm, llm, nqmx
    use disvert_m, only: presnivs
    use dynetat0_chosen_m, only: day_ref, annee_ref
    use dynetat0_m, only: rlatu, rlonv
    USE histbeg_totreg_m, ONLY : histbeg_totreg
    USE histdef_m, ONLY : histdef
    USE histend_m, ONLY : histend
    USE histvert_m, ONLY : histvert
    use indicesol, only: nbsrf, clnsurf
    use infotrac_init_m, only: tname, ttext
    use phyetat0_m, only: itau_phy
    USE ymds2ju_m, only: ymds2ju

    ! Local:
    double precision julian
    real zsto, zout
    integer nhori, nvert, nsrf, iq, it

    !-------------------------------------------------------------------

    print *, 'Call sequence information: ini_histins'

    test_ok_instan: IF (ok_instan) THEN
       zsto = dtphys * ecrit_ins
       zout = dtphys * ecrit_ins
       CALL ymds2ju(annee_ref, 1, day_ref, 0.0, julian)
       CALL histbeg_totreg("histins", rlonv(:iim) / pi * 180., &
            rlatu / pi * 180., 1, iim, &
            1, jjm + 1, itau_phy, julian, dtphys, nhori, nid_ins)
       print *, 'itau_phy = ', itau_phy
       print *, "julian = ", julian
       CALL histvert(nid_ins, "presnivs", "Vertical levels", "mb", &
            presnivs/100., nvert)

       ! Once:
       CALL histdef(nid_ins, "phis", "surface geopotential", "m2 s-2", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, &
            "once", zsto, zout)
       CALL histdef(nid_ins, "aire", "Grid area", "-", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, &
            "once", zsto, zout)

       ! Champs 2D:

       CALL histdef(nid_ins, "tsol", "Surface Temperature", "K", &
            iim, jjm + 1, nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "t2m", "Temperature 2m", "K", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "q2m", "Specific humidity 2m", "Kg/Kg", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "u10m", "Vent zonal 10m", "m/s", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "v10m", "Vent meridien 10m", "m/s", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "psol", "Surface Pressure", "Pa", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "plul", "Large-scale Precip.", "mm/day", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "pluc", "Convective Precip.", "mm/day", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "cdrm", "Momentum drag coef.", "-", &
            iim, jjm + 1, nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "cdrh", "Heat drag coef.", "-", &
            iim, jjm + 1, nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "precip", "Precipitation Totale liq+sol",  &
            "kg/(s*m2)", &
            iim, jjm + 1, nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "snow", "Snow fall", "kg/(s*m2)", &
            iim, jjm + 1, nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "topl", "OLR", "W/m2", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "evap", "evaporation", "kg/(s*m2)", &
            iim, jjm + 1, nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "sols", "Solar rad. at surf.", "W/m2", &
            iim, jjm + 1, nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "rls", "surface net downward longwave flux", &
            "W/m2", iim, jjm + 1, nhori, 1, 1, 1, -99, "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "solldown", "Down. IR rad. at surface",  &
            "W/m2", iim, jjm + 1, nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "bils", "Surf. total heat flux", "W/m2", &
            iim, jjm + 1, nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "sens", "Sensible heat flux", "W/m2", &
            iim, jjm + 1, nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "fder", "Heat flux derivative", "W/m2", &
            iim, jjm + 1, nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "dtsvdfo", "Boundary-layer dTs(o)", "K/s", &
            iim, jjm + 1, nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "dtsvdft", "Boundary-layer dTs(t)", "K/s", &
            iim, jjm + 1, nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "dtsvdfg", "Boundary-layer dTs(g)", "K/s", &
            iim, jjm + 1, nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "dtsvdfi", "Boundary-layer dTs(g)", "K/s", &
            iim, jjm + 1, nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "msnow", "surface snow amount", "kg/m2", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "zxfqcalving", "ice calving", "kg m-2 s-1", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "run_off_lic", "land ice melt to ocean", &
            "kg m-2 s-1", iim, jjm + 1, nhori, 1, 1, 1, -99, "inst(X)", zsto, &
            zout)
       CALL histdef(nid_ins, "flat", "latent heat flux", "W m-2", iim, &
            jjm + 1, nhori, 1, 1, 1, -99, "inst(X)", zsto, zout)

       DO nsrf = 1, nbsrf
          call histdef(nid_ins, "fract_"//clnsurf(nsrf),  &
               "Fraction "//clnsurf(nsrf), "1",   &
               iim, jjm + 1, nhori, 1, 1, 1, -99, &
               "inst(X)", zsto, zout)
          call histdef(nid_ins, "sens_"//clnsurf(nsrf),  &
               "Sensible heat flux "//clnsurf(nsrf), "W/m2",   &
               iim, jjm + 1, nhori, 1, 1, 1, -99, &
               "inst(X)", zsto, zout)
          call histdef(nid_ins, "tsol_"//clnsurf(nsrf),  &
               "Surface Temperature"//clnsurf(nsrf), "W/m2",   &
               iim, jjm + 1, nhori, 1, 1, 1, -99, &
               "inst(X)", zsto, zout)
          call histdef(nid_ins, "lat_"//clnsurf(nsrf),  &
               "Latent heat flux "//clnsurf(nsrf), "W/m2",   &
               iim, jjm + 1, nhori, 1, 1, 1, -99, &
               "inst(X)", zsto, zout)
          call histdef(nid_ins, "taux_"//clnsurf(nsrf),  &
               "Zonal wind stress"//clnsurf(nsrf), "Pa", &
               iim, jjm + 1, nhori, 1, 1, 1, -99, &
               "inst(X)", zsto, zout)
          call histdef(nid_ins, "tauy_"//clnsurf(nsrf),  &
               "Meridional xind stress "//clnsurf(nsrf), "Pa",   &
               iim, jjm + 1, nhori, 1, 1, 1, -99, &
               "inst(X)", zsto, zout)
          call histdef(nid_ins, "albe_"//clnsurf(nsrf),  &
               "Albedo "//clnsurf(nsrf), "-",   &
               iim, jjm + 1, nhori, 1, 1, 1, -99, &
               "inst(X)", zsto, zout)
          call histdef(nid_ins, "rugs_"//clnsurf(nsrf),  &
               "rugosite "//clnsurf(nsrf), "-",   &
               iim, jjm + 1, nhori, 1, 1, 1, -99, &
               "inst(X)", zsto, zout)
          call histdef(nid_ins, "u10m_"//clnsurf(nsrf),  &
               "zonal wind 10 m "//clnsurf(nsrf), "m s-1",   &
               iim, jjm + 1, nhori, 1, 1, 1, -99, &
               "inst(X)", zsto, zout)
          call histdef(nid_ins, "v10m_"//clnsurf(nsrf),  &
               "meridional wind 10 m "//clnsurf(nsrf), "m s-1",   &
               iim, jjm + 1, nhori, 1, 1, 1, -99, &
               "inst(X)", zsto, zout)
       END DO

       CALL histdef(nid_ins, "rugs", "rugosity", "-", &
            iim, jjm + 1, nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "albs", "Surface albedo", "-", &
            iim, jjm + 1, nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "s_pblh", "Boundary Layer Height", "m", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "s_pblt", "T at Boundary Layer Height",  &
            "K", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "s_lcl", "Condensation level", "m", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "s_capCL", &
            "Convective available potential energy for atmospheric boundary " &
            // "layer", "J/m2", iim, jjm + 1, nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "s_oliqCL", "Liq Water in BL", "kg/m2", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "s_cteiCL", "Instability criteria (ABL)", "K", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "s_therm", "Exces du thermique", "K", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "qsurf", "Surface Air humidity", "", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "ffonte", "Thermal flux for snow melting", &
            "W m-2", iim, jjm + 1, nhori, 1, 1, 1, -99, "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "ue", "zonal dry static energy transport", &
            "W m-1", iim, jjm + 1, nhori, 1, 1, 1, -99, "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "ve", "meridional dry static energy transport", &
            "W m-1", iim, jjm + 1, nhori, 1, 1, 1, -99, "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "uq", "zonal humidity transport", "kg s-1 m-1", &
            iim, jjm + 1, nhori, 1, 1, 1, -99, "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "vq", "meridional humidity transport", &
            "kg s-1 m-1", iim, jjm + 1, nhori, 1, 1, 1, -99, "inst(X)", zsto, &
            zout)

       if (conv_emanuel) then
          CALL histdef(nid_ins, "ptop", "cloud top pressure", &
               "Pa", iim, jjm + 1, nhori, 1, 1, 1, -99, "inst(X)", zsto, zout)
          CALL histdef(nid_ins, "dnwd0", "unsaturated downdraft", &
               "kg/m2/s", iim, jjm + 1, nhori, llm, 1, llm, nvert, "inst(X)", &
               zsto, zout)
       end if

       ! Champs 3D:

       CALL histdef(nid_ins, "tro3", "ozone mole fraction", "-", &
            iim, jjm + 1, nhori, llm, 1, llm, nvert, "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "temp", "Temperature", "K", &
            iim, jjm + 1, nhori, llm, 1, llm, nvert, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "vitu", "Zonal wind", "m/s", &
            iim, jjm + 1, nhori, llm, 1, llm, nvert, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "vitv", "Merid wind", "m/s", &
            iim, jjm + 1, nhori, llm, 1, llm, nvert, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "geop", "Geopotential height", "m", &
            iim, jjm + 1, nhori, llm, 1, llm, nvert, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "pres", "Air pressure", "Pa", &
            iim, jjm + 1, nhori, llm, 1, llm, nvert, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "dtvdf", "Boundary-layer dT", "K/s", &
            iim, jjm + 1, nhori, llm, 1, llm, nvert, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "dqvdf", "Boundary-layer dQ", "Kg/Kg/s", &
            iim, jjm + 1, nhori, llm, 1, llm, nvert, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "zmasse", "column density of air in cell", &
            "kg m-2", iim, jjm + 1, nhori, llm, 1, llm, nvert, "inst(X)", &
            zsto, zout)
       CALL histdef(nid_ins, "rhum", "Relative humidity", &
            "", iim, jjm + 1, nhori, llm, 1, llm, nvert, "inst(X)", &
            zsto, zout)
       CALL histdef(nid_ins, "d_t_ec", "kinetic dissipation dT", &
            "K/s", iim, jjm + 1, nhori, llm, 1, llm, nvert, "inst(X)", &
            zsto, zout)
       CALL histdef(nid_ins, "dtsw0", "CS SW radiation dT", &
            "K/s", iim, jjm + 1, nhori, llm, 1, llm, nvert, "inst(X)", &
            zsto, zout)
       CALL histdef(nid_ins, "dtlw0", "CS LW radiation dT", &
            "K/s", iim, jjm + 1, nhori, llm, 1, llm, nvert, "inst(X)", &
            zsto, zout)
       CALL histdef(nid_ins, "pmflxr", "convective precipitation liquid", "", &
            iim, jjm + 1, nhori, llm, 1, llm, nvert, "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "re", "cloud droplet effective radius", &
            "micrometer", iim, jjm + 1, nhori, llm, 1, llm, nvert, &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "fl", &
            "denominator of cloud droplet effective radius", "", iim, &
            jjm + 1, nhori, llm, 1, llm, nvert, "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "rld", "LW downward radiation", "W m-2", iim, &
            jjm + 1, nhori, llm, 1, llm, nvert, "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "rldcs", "LW CS downward radiation", "W m-2", &
            iim, jjm + 1, nhori, llm, 1, llm, nvert, "inst(X)", zsto, zout)

       DO it = 1, nqmx - 2
          ! champ 2D
          iq=it+2
          CALL histdef(nid_ins, tname(iq), ttext(iq), "U/kga", iim, jjm+1, &
               nhori, llm, 1, llm, nvert, "inst(X)", zsto, zout)
          CALL histdef(nid_ins, "fl"//tname(iq), "Flux "//ttext(iq), &
               "U/m2/s", iim, jjm+1, nhori, llm, 1, llm, nvert, &
               "inst(X)", zsto, zout)
          CALL histdef(nid_ins, "d_tr_th_"//tname(iq), &
               "tendance thermique"// ttext(iq), "?", &
               iim, jjm+1, nhori, llm, 1, llm, nvert, &
               "inst(X)", zsto, zout)
          CALL histdef(nid_ins, "d_tr_cv_"//tname(iq), &
               "tendance convection"// ttext(iq), "?", &
               iim, jjm+1, nhori, llm, 1, llm, nvert, &
               "inst(X)", zsto, zout)
          CALL histdef(nid_ins, "d_tr_cl_"//tname(iq), &
               "tendance couche limite"// ttext(iq), "?", &
               iim, jjm+1, nhori, llm, 1, llm, nvert, &
               "inst(X)", zsto, zout)
       ENDDO

       CALL histend(nid_ins)
    ENDIF test_ok_instan

  end subroutine ini_histins

end module ini_histins_m
