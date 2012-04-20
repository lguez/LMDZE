module ini_histhf3d_m

  implicit none

contains

  subroutine ini_histhf3d(dtime, nid_hf3d)

    ! From phylmd/ini_histhf3d.h, v 1.2 2005/05/25 13:10:09

    ! sorties hf 3d

    use dimens_m, only: iim, jjm, llm
    use dimphy, only: klon, nbtr
    use temps, only: itau_phy, day_ref, annee_ref
    use clesphys, only: ecrit_hf
    use phyetat0_m, only: rlon, rlat
    USE calendar, only: ymds2ju
    USE histbeg_totreg_m, ONLY : histbeg_totreg
    USE histdef_m, ONLY : histdef
    USE histend_m, ONLY : histend
    USE histvert_m, ONLY : histvert
    use comvert, only: presnivs

    REAL, intent(in):: dtime ! pas temporel de la physique (s)
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

end module ini_histhf3d_m
