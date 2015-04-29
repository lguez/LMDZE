module ini_histrac_m

  implicit none

contains

  subroutine ini_histrac(nid_tra, pdtphys, nq_phys, lessivage)

    ! From phylmd/ini_histrac.h, version 1.10 2006/02/21 08:08:30

    use clesphys, only: ecrit_tra
    use dimens_m, only: iim, jjm, llm
    use disvert_m, only: presnivs
    use dimphy, only: klon
    use dynetat0_m, only: day_ref, annee_ref
    use grid_change, only: gr_phy_write_2d
    USE histbeg_totreg_m, ONLY: histbeg_totreg
    USE histdef_m, ONLY : histdef
    USE histend_m, ONLY : histend
    USE histvert_m, ONLY : histvert
    use iniadvtrac_m, only: tname, ttext
    use phyetat0_m, only: rlon, rlat
    use temps, only: itau_phy
    USE ymds2ju_m, only: ymds2ju

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
    integer it, iq

    !---------------------------------------------------------

    CALL ymds2ju(annee_ref, month=1, day=day_ref, sec=0.0, julian=zjulian)
    zx_lat(:, :) = gr_phy_write_2d(rlat)
    CALL histbeg_totreg("histrac", rlon(2:iim+1), zx_lat(1, :), &
         1, iim, 1, jjm+1, itau_phy, zjulian, pdtphys, nhori, nid_tra)
    CALL histvert(nid_tra, "presnivs", "Vertical levels", "mb", presnivs, nvert)

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
       CALL histdef(nid_tra, tname(iq), ttext(iq), "U/kga", iim, jjm+1, &
            nhori, llm, 1, llm, nvert, "ave(X)", zsto, zout)
       if (lessivage) THEN
          CALL histdef(nid_tra, "fl"//tname(iq), "Flux "//ttext(iq), &
               "U/m2/s", iim, jjm+1, nhori, llm, 1, llm, nvert, &
               "ave(X)", zsto, zout)
       endif

       !---Ajout Olivia
       CALL histdef(nid_tra, "d_tr_th_"//tname(iq), &
            "tendance thermique"// ttext(iq), "?", &
            iim, jjm+1, nhori, llm, 1, llm, nvert, &
            "ave(X)", zsto, zout)
       CALL histdef(nid_tra, "d_tr_cv_"//tname(iq), &
            "tendance convection"// ttext(iq), "?", &
            iim, jjm+1, nhori, llm, 1, llm, nvert, &
            "ave(X)", zsto, zout)
       CALL histdef(nid_tra, "d_tr_cl_"//tname(iq), &
            "tendance couche limite"// ttext(iq), "?", &
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
