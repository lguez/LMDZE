module ini_hist

  ! This module is clean: no C preprocessor directive, no include line.

  IMPLICIT none

contains

  subroutine ini_histhf(dtime, presnivs, nid_hf, nid_hf3d)

    ! From phylmd/ini_histhf.h, version 1.3 2005/05/25 13:10:09

    use dimens_m, only: iim, jjm, llm
    use temps, only: day_ref, annee_ref, itau_phy
    use dimphy, only: klon
    USE ioipsl, only: ymds2ju, histbeg_totreg, histvert, histend
    use phyetat0_m, only: rlon, rlat

    REAL, intent(in):: dtime ! pas temporel de la physique (s)
    real, intent(in):: presnivs(:)
    integer, intent(out):: nid_hf, nid_hf3d

    REAL zx_lon(iim, jjm + 1), zx_lat(iim, jjm + 1)
    integer idayref
    real zjulian
    integer i, nhori, nvert

    !-----------------------------------------------

    idayref = day_ref
    CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)

    CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), rlon, zx_lon)
    DO i = 1, iim
       zx_lon(i, 1) = rlon(i+1)
       zx_lon(i, (jjm + 1)) = rlon(i+1)
    ENDDO

    CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), rlat, zx_lat)
    CALL histbeg_totreg("histhf", zx_lon(:, 1), zx_lat(1, :), 1, iim, 1, &
         (jjm + 1), itau_phy, zjulian, dtime, nhori, nid_hf)

    CALL histvert(nid_hf, "presnivs", "Vertical levels", "mb", &
         llm, presnivs/100., nvert)

    call ini_histhf3d(dtime, presnivs, nid_hf3d)
    CALL histend(nid_hf)

  end subroutine ini_histhf

  !******************************************************************

  subroutine ini_histhf3d(dtime, presnivs, nid_hf3d)

    ! From phylmd/ini_histhf3d.h, v 1.2 2005/05/25 13:10:09

    ! sorties hf 3d

    use dimens_m, only: iim, jjm, llm
    use dimphy, only: klon, nbtr
    use temps, only: itau_phy, day_ref, annee_ref
    use clesphys, only: ecrit_hf
    use phyetat0_m, only: rlon, rlat
    USE ioipsl, only: ymds2ju, histbeg_totreg, histvert, histend, histdef

    REAL, intent(in):: dtime ! pas temporel de la physique (s)
    real, intent(in):: presnivs(:)
    integer, intent(out):: nid_hf3d

    real zstohf, zout
    REAL zx_lon(iim, jjm + 1), zx_lat(iim, jjm + 1)
    real zjulian
    integer i, nhori, nvert, idayref

    !------------------------------------------

    zstohf = dtime * REAL(ecrit_hf)
    zout = dtime * REAL(ecrit_hf)

    idayref = day_ref
    CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)

    CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), rlon, zx_lon)
    DO i = 1, iim
       zx_lon(i, 1) = rlon(i+1)
       zx_lon(i, (jjm + 1)) = rlon(i+1)
    ENDDO

    CALL gr_fi_ecrit(1, klon, iim, (jjm + 1), rlat, zx_lat)
    CALL histbeg_totreg("histhf3d", zx_lon(:, 1), zx_lat(1, :), 1, iim, 1, &
         (jjm + 1), itau_phy, zjulian, dtime, nhori, nid_hf3d)

    CALL histvert(nid_hf3d, "presnivs", "Vertical levels", "mb", &
         llm, presnivs/100., nvert)

    ! Champs 3D:

    CALL histdef(nid_hf3d, "temp", "Air temperature", "K", &
         iim, (jjm + 1), nhori, llm, 1, llm, nvert, &
         "ave(X)", zstohf, zout)

    CALL histdef(nid_hf3d, "ovap", "Specific humidity", "kg/kg", &
         iim, (jjm + 1), nhori, llm, 1, llm, nvert, &
         "ave(X)", zstohf, zout)

    CALL histdef(nid_hf3d, "vitu", "Zonal wind", "m/s", &
         iim, (jjm + 1), nhori, llm, 1, llm, nvert, &
         "ave(X)", zstohf, zout)

    CALL histdef(nid_hf3d, "vitv", "Meridional wind", "m/s", &
         iim, (jjm + 1), nhori, llm, 1, llm, nvert, &
         "ave(X)", zstohf, zout)

    if (nbtr >= 3) then
       CALL histdef(nid_hf3d, "O3", "Ozone mass fraction", "?", iim, &
            (jjm + 1), nhori, llm, 1, llm, nvert, "ave(X)", zstohf, &
            zout)
    end if

    CALL histend(nid_hf3d)

  end subroutine ini_histhf3d

  !******************************************************************

  subroutine ini_histday(dtime, presnivs, ok_journe, nid_day, nq)

    ! From phylmd/ini_histday.h, v 1.3 2005/05/25 13:10:09

    use dimens_m, only: iim, jjm, llm
    use temps, only: itau_phy, day_ref, annee_ref
    USE ioipsl, only: ymds2ju, histbeg_totreg, histvert, histend, histdef
    use phyetat0_m, only: rlon, rlat
    use clesphys, only: ecrit_day
    use grid_change, only: gr_phy_write_2d

    REAL, intent(in):: dtime ! pas temporel de la physique (s)
    real, intent(in):: presnivs(:)
    logical, intent(in):: ok_journe
    integer, intent(out):: nid_day
    INTEGER, intent(in):: nq ! nombre de traceurs (y compris vapeur d'eau)

    ! Variables local to the procedure:
    REAL zx_lat(iim, jjm + 1)
    integer nhori, nvert
    real zjulian

    !--------------------------------

    IF (ok_journe) THEN
       CALL ymds2ju(annee_ref, 1, day_ref, 0., zjulian)
       zx_lat = gr_phy_write_2d(rlat)
       CALL histbeg_totreg("histday", rlon(2: iim+1), zx_lat(1, :), 1, iim, &
            1, jjm + 1, itau_phy, zjulian, dtime, nhori, nid_day)
       CALL histvert(nid_day, "presnivs", "Vertical levels", "mb", &
            llm, presnivs/100., nvert)
       if (nq <= 4) then
          call histdef(nid_day, "Sigma_O3_Royer", &
               "column-density of ozone, in a cell, from Royer", "DU", &
               pxsize=iim, pysize=jjm+1, phoriid=nhori, pzsize=llm, &
               par_oriz=1, par_szz=llm, pzid=nvert, popp="ave(X)", &
               pfreq_opp=dtime, pfreq_wrt=real(ecrit_day))
       end if
       CALL histend(nid_day)
    ENDIF

  end subroutine ini_histday

  !****************************************************

  subroutine ini_histins(dtime, presnivs, ok_instan, nid_ins)

    ! From phylmd/ini_histins.h, v 1.2 2005/05/25 13:10:09

    use dimens_m, only: iim, jjm, llm
    use dimphy, only: klon
    use temps, only: itau_phy, day_ref, annee_ref
    use clesphys, only: ecrit_ins
    use indicesol, only: nbsrf, clnsurf
    USE ioipsl, only: ymds2ju, histbeg_totreg, histvert, histend, histdef
    use phyetat0_m, only: rlon, rlat

    REAL, intent(in):: dtime ! pas temporel de la physique (s)
    real, intent(in):: presnivs(:)
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

  !*************************************************

  subroutine ini_histrac(nid_tra, pdtphys, presnivs, nq_phys, lessivage)

    ! From phylmd/ini_histrac.h, version 1.10 2006/02/21 08:08:30

    use dimens_m, only: iim, jjm, llm
    use ioipsl, only: ymds2ju, histbeg_totreg, histvert, histdef, histend
    use temps, only: annee_ref, day_ref, itau_phy
    use iniadvtrac_m, only: niadv, tnom, ttext
    use dimphy, only: klon
    use clesphys, only: ecrit_tra
    use grid_change, only: gr_phy_write_2d
    use phyetat0_m, only: rlon, rlat

    INTEGER, intent(out):: nid_tra
    real, intent(in):: pdtphys  ! pas d'integration pour la physique (s)
    REAL, intent(in):: presnivs(:)

    integer, intent(in):: nq_phys
    ! (nombre de traceurs auxquels on applique la physique)

    logical, intent(in):: lessivage

    ! Variables local to the procedure:

    REAL zjulian
    REAL zx_lat(iim, jjm+1)
    INTEGER nhori, nvert
    REAL zsto, zout
    integer it, iq, iiq

    !---------------------------------------------------------

    CALL ymds2ju(annee_ref, month=1, day=day_ref, sec=0.0, julian=zjulian)
    zx_lat(:, :) = gr_phy_write_2d(rlat)
    CALL histbeg_totreg("histrac", rlon(2:iim+1), zx_lat(1, :), &
         1, iim, 1, jjm+1, itau_phy, zjulian, pdtphys, nhori, nid_tra)
    CALL histvert(nid_tra, "presnivs", "Vertical levels", "mb", llm, &
         presnivs, nvert)

    zsto = pdtphys
    zout = pdtphys * REAL(ecrit_tra)

    CALL histdef(nid_tra, "phis", "Surface geop. height", "-", &
         iim, jjm+1, nhori, 1, 1, 1, -99, &
         "once",  zsto, zout)
    CALL histdef(nid_tra, "aire", "Grid area", "-", &
         iim, jjm+1, nhori, 1, 1, 1, -99, &
         "once",  zsto, zout)
    CALL histdef(nid_tra, "zmasse", "column density of air in cell", &
         "kg m-2", iim, jjm + 1, nhori, llm, 1, llm, nvert, "ave(X)", &
         zsto, zout)

    DO it = 1, nq_phys
       ! champ 2D
       iq=it+2
       iiq=niadv(iq)
       CALL histdef(nid_tra, tnom(iq), ttext(iiq), "U/kga", iim, jjm+1, &
            nhori, llm, 1, llm, nvert, "ave(X)", zsto, zout)
       if (lessivage) THEN
          CALL histdef(nid_tra, "fl"//tnom(iq), "Flux "//ttext(iiq), &
               "U/m2/s", iim, jjm+1, nhori, llm, 1, llm, nvert, &
               "ave(X)", zsto, zout)
       endif

       !---Ajout Olivia
       CALL histdef(nid_tra, "d_tr_th_"//tnom(iq), &
            "tendance thermique"// ttext(iiq), "?", &
            iim, jjm+1, nhori, llm, 1, llm, nvert, &
            "ave(X)", zsto, zout)
       CALL histdef(nid_tra, "d_tr_cv_"//tnom(iq), &
            "tendance convection"// ttext(iiq), "?", &
            iim, jjm+1, nhori, llm, 1, llm, nvert, &
            "ave(X)", zsto, zout)
       CALL histdef(nid_tra, "d_tr_cl_"//tnom(iq), &
            "tendance couche limite"// ttext(iiq), "?", &
            iim, jjm+1, nhori, llm, 1, llm, nvert, &
            "ave(X)", zsto, zout)
       !---fin Olivia	 

    ENDDO

    CALL histdef(nid_tra, "pplay", "", "-", &
         iim, jjm+1, nhori, llm, 1, llm, nvert, &
         "inst(X)", zout, zout)
    CALL histdef(nid_tra, "t", "", "-", &
         iim, jjm+1, nhori, llm, 1, llm, nvert, &
         "inst(X)", zout, zout)

    CALL histend(nid_tra)

  end subroutine ini_histrac

end module ini_hist
