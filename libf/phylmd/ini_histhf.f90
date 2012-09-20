module ini_histhf_m

  implicit none

contains

  subroutine ini_histhf(dtime, nid_hf, nid_hf3d)

    ! From phylmd/ini_histhf.h, version 1.3 2005/05/25 13:10:09

    use dimens_m, only: iim, jjm, llm
    use temps, only: day_ref, annee_ref, itau_phy
    use dimphy, only: klon
    USE calendar, only: ymds2ju
    USE histbeg_totreg_m, ONLY : histbeg_totreg
    USE histend_m, ONLY : histend
    USE histvert_m, ONLY : histvert
    use phyetat0_m, only: rlon, rlat
    use disvert_m, only: presnivs
    use ini_histhf3d_m, only: ini_histhf3d

    REAL, intent(in):: dtime ! pas temporel de la physique (s)
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

    call ini_histhf3d(dtime, nid_hf3d)
    CALL histend(nid_hf)

  end subroutine ini_histhf

end module ini_histhf_m
