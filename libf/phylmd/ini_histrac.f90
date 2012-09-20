module ini_histrac_m

  implicit none

contains


  subroutine ini_histrac(nid_tra, pdtphys, nq_phys, lessivage)

    ! From phylmd/ini_histrac.h, version 1.10 2006/02/21 08:08:30

    use dimens_m, only: iim, jjm, llm
    USE calendar, only: ymds2ju
    USE histbeg_totreg_m, ONLY : histbeg_totreg
    USE histdef_m, ONLY : histdef
    USE histend_m, ONLY : histend
    USE histvert_m, ONLY : histvert
    use temps, only: annee_ref, day_ref, itau_phy
    use iniadvtrac_m, only: niadv, tnom, ttext
    use dimphy, only: klon
    use clesphys, only: ecrit_tra
    use grid_change, only: gr_phy_write_2d
    use phyetat0_m, only: rlon, rlat
    use disvert_m, only: presnivs

    INTEGER, intent(out):: nid_tra
    real, intent(in):: pdtphys  ! pas d'integration pour la physique (s)

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
    CALL histdef(nid_tra, "T", "temperature", "K", iim, jjm+1, nhori, llm, &
         1, llm, nvert, "inst(X)", zout, zout)

    CALL histend(nid_tra)

  end subroutine ini_histrac

end module ini_histrac_m
