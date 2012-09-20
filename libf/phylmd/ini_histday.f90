module ini_histday_m

  implicit none

contains

  subroutine ini_histday(dtime, ok_journe, nid_day, nq)

    ! From phylmd/ini_histday.h, v 1.3 2005/05/25 13:10:09

    use dimens_m, only: iim, jjm, llm
    use temps, only: itau_phy, day_ref, annee_ref
    USE calendar, only: ymds2ju
    USE histbeg_totreg_m, ONLY : histbeg_totreg
    USE histdef_m, ONLY : histdef
    USE histend_m, ONLY : histend
    USE histvert_m, ONLY : histvert
    use phyetat0_m, only: rlon, rlat
    use clesphys, only: ecrit_day
    use grid_change, only: gr_phy_write_2d
    use disvert_m, only: presnivs

    REAL, intent(in):: dtime ! pas temporel de la physique (s)
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
               xsize=iim, ysize=jjm+1, horiid=nhori, pzsize=llm, par_oriz=1, &
               par_szz=llm, pzid=nvert, popp="ave(X)", pfreq_opp=dtime, &
               pfreq_wrt=real(ecrit_day))
       end if
       CALL histend(nid_day)
    ENDIF

  end subroutine ini_histday

end module ini_histday_m
