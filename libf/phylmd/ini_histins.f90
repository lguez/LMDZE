module ini_histins_m

  implicit none

contains

  subroutine ini_histins(dtime, ok_instan, nid_ins)

    ! From phylmd/ini_histins.h, v 1.2 2005/05/25 13:10:09

    use dimens_m, only: iim, jjm, llm
    use dimphy, only: klon
    use temps, only: itau_phy, day_ref, annee_ref
    use clesphys, only: ecrit_ins
    use indicesol, only: nbsrf, clnsurf
    USE calendar, only: ymds2ju
    USE histbeg_totreg_m, ONLY : histbeg_totreg
    USE histdef_m, ONLY : histdef
    USE histend_m, ONLY : histend
    USE histvert_m, ONLY : histvert
    use phyetat0_m, only: rlon, rlat
    use comvert, only: presnivs

    REAL, intent(in):: dtime ! pas temporel de la physique (s)
    logical, intent(in):: ok_instan
    integer, intent(out):: nid_ins

    REAL zx_lon(iim, jjm + 1), zx_lat(iim, jjm + 1)
    real zjulian, zsto, zout
    integer i, nhori, nvert, idayref, nsrf

    !-------------------------------------------------------------------

    IF (ok_instan) THEN

       zsto = dtime * ecrit_ins
       zout = dtime * ecrit_ins

       idayref = day_ref
       CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)

       CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), rlon, zx_lon)
       DO i = 1, iim
          zx_lon(i, 1) = rlon(i+1)
          zx_lon(i, (jjm + 1)) = rlon(i+1)
       ENDDO
       CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), rlat, zx_lat)
       CALL histbeg_totreg("histins", zx_lon(:, 1), zx_lat(1, :), 1, iim, 1, &
            jjm + 1, itau_phy, zjulian, dtime, nhori, nid_ins)
       write(*, *)'Inst ', itau_phy, zjulian
       CALL histvert(nid_ins, "presnivs", "Vertical levels", "mb", &
            llm, presnivs/100., nvert)

       CALL histdef(nid_ins, "phis", "Surface geop. height", "-", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "once", zsto, zout)

       CALL histdef(nid_ins, "aire", "Grid area", "-", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "once", zsto, zout)

       ! Champs 2D:

       CALL histdef(nid_ins, "tsol", "Surface Temperature", "K", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "t2m", "Temperature 2m", "K", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "q2m", "Specific humidity 2m", "Kg/Kg", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "u10m", "Vent zonal 10m", "m/s", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "v10m", "Vent meridien 10m", "m/s", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "psol", "Surface Pressure", "Pa", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "plul", "Large-scale Precip.", "mm/day", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "pluc", "Convective Precip.", "mm/day", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "cdrm", "Momentum drag coef.", "-", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "cdrh", "Heat drag coef.", "-", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "precip", "Precipitation Totale liq+sol",  &
            "kg/(s*m2)", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "snow", "Snow fall", "kg/(s*m2)", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)

       !        CALL histdef(nid_ins, "snow_mass", "Snow Mass", "kg/m2",
       !    .                iim, (jjm + 1), nhori, 1, 1, 1, -99,
       !    .                "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "topl", "OLR", "W/m2", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "evap", "Evaporation", "kg/(s*m2)", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "sols", "Solar rad. at surf.", "W/m2", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "soll", "IR rad. at surface", "W/m2", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "solldown", "Down. IR rad. at surface",  &
            "W/m2", iim, (jjm + 1), nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "bils", "Surf. total heat flux", "W/m2", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "sens", "Sensible heat flux", "W/m2", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "fder", "Heat flux derivation", "W/m2", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "dtsvdfo", "Boundary-layer dTs(o)", "K/s", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "dtsvdft", "Boundary-layer dTs(t)", "K/s", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "dtsvdfg", "Boundary-layer dTs(g)", "K/s", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "dtsvdfi", "Boundary-layer dTs(g)", "K/s", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)

       DO nsrf = 1, nbsrf

          call histdef(nid_ins, "pourc_"//clnsurf(nsrf),  &
               "% "//clnsurf(nsrf), "%",   &
               iim, (jjm + 1), nhori, 1, 1, 1, -99, &
               "inst(X)", zsto, zout)

          call histdef(nid_ins, "fract_"//clnsurf(nsrf),  &
               "Fraction "//clnsurf(nsrf), "1",   &
               iim, (jjm + 1), nhori, 1, 1, 1, -99, &
               "inst(X)", zsto, zout)

          call histdef(nid_ins, "sens_"//clnsurf(nsrf),  &
               "Sensible heat flux "//clnsurf(nsrf), "W/m2",   &
               iim, (jjm + 1), nhori, 1, 1, 1, -99, &
               "inst(X)", zsto, zout)

          call histdef(nid_ins, "tsol_"//clnsurf(nsrf),  &
               "Surface Temperature"//clnsurf(nsrf), "W/m2",   &
               iim, (jjm + 1), nhori, 1, 1, 1, -99, &
               "inst(X)", zsto, zout)

          call histdef(nid_ins, "lat_"//clnsurf(nsrf),  &
               "Latent heat flux "//clnsurf(nsrf), "W/m2",   &
               iim, (jjm + 1), nhori, 1, 1, 1, -99, &
               "inst(X)", zsto, zout)

          call histdef(nid_ins, "taux_"//clnsurf(nsrf),  &
               "Zonal wind stress"//clnsurf(nsrf), "Pa", &
               iim, (jjm + 1), nhori, 1, 1, 1, -99, &
               "inst(X)", zsto, zout)

          call histdef(nid_ins, "tauy_"//clnsurf(nsrf),  &
               "Meridional xind stress "//clnsurf(nsrf), "Pa",   &
               iim, (jjm + 1), nhori, 1, 1, 1, -99, &
               "inst(X)", zsto, zout)

          call histdef(nid_ins, "albe_"//clnsurf(nsrf),  &
               "Albedo "//clnsurf(nsrf), "-",   &
               iim, (jjm + 1), nhori, 1, 1, 1, -99, &
               "inst(X)", zsto, zout)

          call histdef(nid_ins, "rugs_"//clnsurf(nsrf),  &
               "rugosite "//clnsurf(nsrf), "-",   &
               iim, (jjm + 1), nhori, 1, 1, 1, -99, &
               "inst(X)", zsto, zout)
          !XXX
       END DO
       CALL histdef(nid_ins, "rugs", "rugosity", "-", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "albs", "Surface albedo", "-", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)
       CALL histdef(nid_ins, "albslw", "Surface albedo LW", "-", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99,  &
            "inst(X)", zsto, zout)

       !IM cf. AM 081204 BEG
       ! HBTM2
       CALL histdef(nid_ins, "s_pblh", "Boundary Layer Height", "m", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "s_pblt", "T at Boundary Layer Height",  &
            "K", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "s_lcl", "Condensation level", "m", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "s_capCL", "Conv avlbl pot ener for ABL", "J/m2", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "s_oliqCL", "Liq Water in BL", "kg/m2", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "s_cteiCL", "Instability criteria (ABL)", "K", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "s_therm", "Exces du thermique", "K", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "s_trmb1", "deep_cape(HBTM2)", "J/m2", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "s_trmb2", "inhibition (HBTM2)", "J/m2", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "s_trmb3", "Point Omega (HBTM2)", "m", &
            iim, (jjm + 1), nhori, 1, 1, 1, -99, &
            "inst(X)", zsto, zout)

       !IM cf. AM 081204 END

       ! Champs 3D:

       CALL histdef(nid_ins, "temp", "Temperature", "K", &
            iim, (jjm + 1), nhori, llm, 1, llm, nvert, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "vitu", "Zonal wind", "m/s", &
            iim, (jjm + 1), nhori, llm, 1, llm, nvert, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "vitv", "Merid wind", "m/s", &
            iim, (jjm + 1), nhori, llm, 1, llm, nvert, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "geop", "Geopotential height", "m", &
            iim, (jjm + 1), nhori, llm, 1, llm, nvert, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "pres", "Air pressure", "Pa", &
            iim, (jjm + 1), nhori, llm, 1, llm, nvert, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "dtvdf", "Boundary-layer dT", "K/s", &
            iim, (jjm + 1), nhori, llm, 1, llm, nvert, &
            "inst(X)", zsto, zout)

       CALL histdef(nid_ins, "dqvdf", "Boundary-layer dQ", "Kg/Kg/s", &
            iim, (jjm + 1), nhori, llm, 1, llm, nvert, &
            "inst(X)", zsto, zout)

       CALL histend(nid_ins)
    ENDIF

  end subroutine ini_histins

end module ini_histins_m
