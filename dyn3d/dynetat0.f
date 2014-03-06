module dynetat0_m

  IMPLICIT NONE

  INTEGER day_ini

contains

  SUBROUTINE dynetat0(vcov, ucov, teta, q, masse, ps, phis, time_0)

    ! From dynetat0.F, version 1.2, 2004/06/22 11:45:30
    ! Authors: P. Le Van, L. Fairhead
    ! This procedure reads the initial state of the atmosphere.

    use comconst, only: dtvr
    use comgeom, only: rlonu, rlatu, rlonv, rlatv, cu_2d, cv_2d, aire_2d
    use conf_gcm_m, only: fxyhypb, ysinus
    use dimens_m, only: iim, jjm, llm, nqmx
    use disvert_m, only: pa
    use ener, only: etot0, ang0, ptot0, stot0, ztot0
    use iniadvtrac_m, only: tname
    use netcdf, only: NF90_NOWRITE, NF90_NOERR
    use netcdf95, only: NF95_GET_VAR, nf95_open, nf95_inq_varid, NF95_CLOSE, &
         NF95_Gw_VAR
    use nr_util, only: assert
    use serre, only: clon, clat, grossismy, grossismx
    use temps, only: day_ref, itau_dyn, annee_ref

    REAL, intent(out):: vcov(: , :, :) ! (iim + 1, jjm, llm)
    REAL, intent(out):: ucov(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(out):: teta(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(out):: q(:, :, :, :) ! (iim + 1, jjm + 1, llm, nqmx)
    REAL, intent(out):: masse(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(out):: ps(:, :) ! (iim + 1, jjm + 1) in Pa
    REAL, intent(out):: phis(:, :) ! (iim + 1, jjm + 1)
    REAL, intent(out):: time_0

    ! Local variables: 
    INTEGER iq
    REAL, pointer:: tab_cntrl(:) ! tableau des paramètres du run
    INTEGER ierr, ncid, varid

    !-----------------------------------------------------------------------

    print *, "Call sequence information: dynetat0"

    call assert((/size(ucov, 1), size(vcov, 1), size(masse, 1), size(ps, 1), &
         size(phis, 1), size(q, 1), size(teta, 1)/) == iim + 1, "dynetat0 iim")
    call assert((/size(ucov, 2), size(vcov, 2) + 1, size(masse, 2), &
         size(ps, 2), size(phis, 2), size(q, 2), size(teta, 2)/) == jjm + 1, &
         "dynetat0 jjm")
    call assert((/size(vcov, 3), size(ucov, 3), size(teta, 3), size(q, 3), &
         size(masse, 3)/) == llm, "dynetat0 llm")
    call assert(size(q, 4) == nqmx, "dynetat0 q nqmx")

    ! Fichier état initial :
    call nf95_open("start.nc", NF90_NOWRITE, ncid)

    call nf95_inq_varid(ncid, "controle", varid)
    call NF95_Gw_VAR(ncid, varid, tab_cntrl)

    call assert(int(tab_cntrl(1)) == iim, "dynetat0 tab_cntrl iim") 
    call assert(int(tab_cntrl(2)) == jjm, "dynetat0 tab_cntrl jjm") 
    call assert(int(tab_cntrl(3)) == llm, "dynetat0 tab_cntrl llm") 

    day_ref = int(tab_cntrl(4))
    annee_ref = int(tab_cntrl(5))

    IF (dtvr /= tab_cntrl(12)) THEN
       print *, 'Warning: the time steps from day_step and "start.nc" ' // &
          'are different.'
       print *, 'dtvr from day_step: ', dtvr
       print *, 'dtvr from "start.nc": ', tab_cntrl(12)
       print *, 'Using the value from day_step.'
    ENDIF

    etot0 = tab_cntrl(13)
    ptot0 = tab_cntrl(14)
    ztot0 = tab_cntrl(15)
    stot0 = tab_cntrl(16)
    ang0 = tab_cntrl(17)
    pa = tab_cntrl(18)
    clon = tab_cntrl(20)
    clat = tab_cntrl(21)
    grossismx = tab_cntrl(22)
    grossismy = tab_cntrl(23)
    fxyhypb = tab_cntrl(24) == 1.
    if (.not. fxyhypb) ysinus = tab_cntrl(27) == 1.
    itau_dyn = tab_cntrl(31)

    call NF95_INQ_VARID (ncid, "rlonu", varid)
    call NF95_GET_VAR(ncid, varid, rlonu)

    call NF95_INQ_VARID (ncid, "rlatu", varid)
    call NF95_GET_VAR(ncid, varid, rlatu)

    call NF95_INQ_VARID (ncid, "rlonv", varid)
    call NF95_GET_VAR(ncid, varid, rlonv)

    call NF95_INQ_VARID (ncid, "rlatv", varid)
    call NF95_GET_VAR(ncid, varid, rlatv)

    call NF95_INQ_VARID (ncid, "cu", varid)
    call NF95_GET_VAR(ncid, varid, cu_2d)

    call NF95_INQ_VARID (ncid, "cv", varid)
    call NF95_GET_VAR(ncid, varid, cv_2d)

    call NF95_INQ_VARID (ncid, "aire", varid)
    call NF95_GET_VAR(ncid, varid, aire_2d)

    call NF95_INQ_VARID (ncid, "phisinit", varid)
    call NF95_GET_VAR(ncid, varid, phis)

    call NF95_INQ_VARID (ncid, "temps", varid)
    call NF95_GET_VAR(ncid, varid, time_0)

    day_ini = tab_cntrl(30) + INT(time_0)
    time_0 = time_0 - INT(time_0)
    ! {0 <= time0 < 1}

    deallocate(tab_cntrl) ! pointer

    call NF95_INQ_VARID (ncid, "ucov", varid)
    call NF95_GET_VAR(ncid, varid, ucov)

    call NF95_INQ_VARID (ncid, "vcov", varid)
    call NF95_GET_VAR(ncid, varid, vcov)

    call NF95_INQ_VARID (ncid, "teta", varid)
    call NF95_GET_VAR(ncid, varid, teta)

    DO iq = 1, nqmx
       call NF95_INQ_VARID(ncid, tname(iq), varid, ierr)
       IF (ierr /= NF90_NOERR) THEN
          PRINT *, 'dynetat0: "' // tname(iq) // '" not found, ' // &
               "setting it to zero..."
          q(:, :, :, iq) = 0.
       ELSE
          call NF95_GET_VAR(ncid, varid, q(:, :, :, iq))
       ENDIF
    ENDDO

    call NF95_INQ_VARID (ncid, "masse", varid)
    call NF95_GET_VAR(ncid, varid, masse)

    call NF95_INQ_VARID (ncid, "ps", varid)
    call NF95_GET_VAR(ncid, varid, ps)
    ! Check that there is a single value at each pole:
    call assert(ps(1, 1) == ps(2:, 1), "dynetat0 ps north pole")
    call assert(ps(1, jjm + 1) == ps(2:, jjm + 1), "dynetat0 ps south pole")

    call NF95_CLOSE(ncid)

  END SUBROUTINE dynetat0

end module dynetat0_m
