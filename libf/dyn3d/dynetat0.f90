module dynetat0_m

  IMPLICIT NONE

  INTEGER day_ini

contains

  SUBROUTINE dynetat0(vcov, ucov, teta, q, masse, ps, phis, time_0)

    ! From dynetat0.F, version 1.2, 2004/06/22 11:45:30
    ! Authors: P. Le Van, L. Fairhead
    ! Objet : lecture de l'�tat initial

    use comconst, only: im, dtvr, jm, lllm, omeg
    use comvert, only: pa
    use comgeom, only: rlonu, rlatu, rlonv, rlatv, cu_2d, cv_2d, aire_2d
    use dimens_m, only: iim, jjm, llm, nqmx
    use ener, only: etot0, ang0, ptot0, stot0, ztot0
    use iniadvtrac_m, only: tname
    use logic, only: fxyhypb, ysinus
    use serre, only: clon, clat, grossismy, grossismx
    use netcdf95, only: NF95_GET_VAR, nf95_open, nf95_inq_varid, NF95_CLOSE
    use netcdf, only: NF90_NOWRITE, NF90_NOERR
    use nr_util, only: assert
    use temps, only: day_ref, itau_dyn, annee_ref

    ! Arguments:
    REAL, intent(out):: vcov(: , :), ucov(:, :), teta(:, :)
    REAL, intent(out):: q(:, :, :), masse(:, :)
    REAL, intent(out):: ps(:) ! in Pa
    REAL, intent(out):: phis(:, :)
    REAL, intent(out):: time_0

    ! Variables 
    INTEGER iq
    INTEGER, PARAMETER:: length = 100
    REAL tab_cntrl(length) ! tableau des param�tres du run
    INTEGER ierr, ncid, varid

    !-----------------------------------------------------------------------

    print *, "Call sequence information: dynetat0"

    call assert(size(vcov, 1) == (iim + 1) * jjm, "dynetat0 vcov 1")
    call assert((/size(ucov, 1), size(teta, 1), size(q, 1), size(masse, 1), &
         size(ps)/) == (iim + 1) * (jjm + 1), "dynetat0 (iim + 1) * (jjm + 1)")
    call assert(shape(phis) == (/iim + 1, jjm + 1/), "dynetat0 phis")
    call assert((/size(vcov, 2), size(ucov, 2), size(teta, 2), size(q, 2), &
         size(masse, 2)/) == llm, "dynetat0 llm")
    call assert(size(q, 3) == nqmx, "dynetat0 q 3")

    ! Fichier �tat initial :
    call nf95_open("start.nc", NF90_NOWRITE, ncid)

    call nf95_inq_varid(ncid, "controle", varid)
    call NF95_GET_VAR(ncid, varid, tab_cntrl)

    im = int(tab_cntrl(1))
    jm = int(tab_cntrl(2))
    lllm = int(tab_cntrl(3))
    call assert(im == iim, "dynetat0 im iim") 
    call assert(jm == jjm, "dynetat0 jm jjm") 
    call assert(lllm == llm, "dynetat0 lllm llm") 

    day_ref = int(tab_cntrl(4))
    annee_ref = int(tab_cntrl(5))

    omeg = tab_cntrl(7)
    PRINT *, 'omeg = ', omeg

    dtvr = tab_cntrl(12)
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

    call NF95_INQ_VARID (ncid, "ucov", varid)
    call NF95_GET_VAR(ncid, varid, ucov, count_nc=(/iim + 1, jjm + 1, llm/))

    call NF95_INQ_VARID (ncid, "vcov", varid)
    call NF95_GET_VAR(ncid, varid, vcov, count_nc=(/iim + 1, jjm, llm/))

    call NF95_INQ_VARID (ncid, "teta", varid)
    call NF95_GET_VAR(ncid, varid, teta, count_nc=(/iim + 1, jjm + 1, llm/))

    DO iq = 1, nqmx
       call NF95_INQ_VARID(ncid, tname(iq), varid, ierr)
       IF (ierr /= NF90_NOERR) THEN
          PRINT *, 'dynetat0: "' // tname(iq) // '" not found, ' // &
               "setting it to zero..."
          q(:, :, iq) = 0.
       ELSE
          call NF95_GET_VAR(ncid, varid, q(:, :, iq), &
               count_nc=(/iim + 1, jjm + 1, llm/))
       ENDIF
    ENDDO

    call NF95_INQ_VARID (ncid, "masse", varid)
    call NF95_GET_VAR(ncid, varid, masse, count_nc=(/iim + 1, jjm + 1, llm/))

    call NF95_INQ_VARID (ncid, "ps", varid)
    call NF95_GET_VAR(ncid, varid, ps, count_nc=(/iim + 1, jjm + 1/))

    call NF95_CLOSE(ncid)

  END SUBROUTINE dynetat0

end module dynetat0_m
