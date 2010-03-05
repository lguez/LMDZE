module dynetat0_m

  IMPLICIT NONE

  INTEGER day_ini

contains

  SUBROUTINE dynetat0(vcov, ucov, teta, q, masse, ps, phis, time_0)

    ! From dynetat0.F, version 1.2 2004/06/22 11:45:30
    ! Authors:  P. Le Van, L. Fairhead
    ! Objet : lecture de l'état initial

    use dimens_m, only: iim, jjm, llm, nqmx
    use comconst, only: im, cpp, dtvr, g, kappa, jm, lllm, omeg, rad
    use comvert, only: pa
    use logic, only: fxyhypb, ysinus
    use comgeom, only: rlonu, rlatu, rlonv, rlatv, cu_2d, cv_2d, aire_2d
    use serre, only: clon, clat, grossismy, grossismx
    use temps, only: day_ref, itau_dyn, annee_ref
    use ener, only: etot0, ang0, ptot0, stot0, ztot0
    use iniadvtrac_m, only: tname
    use netcdf95, only: nf95_open, nf95_inq_varid, handle_err, NF95_CLOSE
    use netcdf, only: NF90_NOWRITE, NF90_GET_VAR, NF90_NOERR
    use numer_rec, only: assert

    !   Arguments:
    REAL, intent(out):: vcov(: , :), ucov(:, :), teta(:, :)
    REAL, intent(out):: q(:, :, :), masse(:, :)
    REAL, intent(out):: ps(:) ! in Pa
    REAL, intent(out):: phis(:, :)
    REAL, intent(out):: time_0

    !   Variables 
    INTEGER length, iq
    PARAMETER (length = 100)
    REAL tab_cntrl(length) ! tableau des parametres du run
    INTEGER ierr, ncid, nvarid

    !-----------------------------------------------------------------------

    print *, "Call sequence information: dynetat0"

    call assert(size(vcov, 1) == (iim + 1) * jjm, "dynetat0 vcov 1")
    call assert((/size(ucov, 1), size(teta, 1), size(q, 1), size(masse, 1), &
         size(ps)/) == (iim + 1) * (jjm + 1), "dynetat0 (iim + 1) * (jjm + 1)")
    call assert(shape(phis) == (/iim + 1, jjm + 1/), "dynetat0 phis")
    call assert((/size(vcov, 2), size(ucov, 2), size(teta, 2), size(q, 2), &
         size(masse, 2)/) == llm, "dynetat0 llm")
    call assert(size(q, 3) == nqmx, "dynetat0 q 3")

    ! Fichier état initial :
    call nf95_open("start.nc", NF90_NOWRITE, ncid)

    call nf95_inq_varid(ncid, "controle", nvarid)
    ierr = NF90_GET_VAR(ncid, nvarid, tab_cntrl)
    call handle_err("dynetat0, controle", ierr, ncid)

    im         = int(tab_cntrl(1))
    jm         = int(tab_cntrl(2))
    lllm       = int(tab_cntrl(3))
    day_ref    = int(tab_cntrl(4))
    annee_ref  = int(tab_cntrl(5))
    omeg       = tab_cntrl(7)
    dtvr       = tab_cntrl(12)
    etot0      = tab_cntrl(13)
    ptot0      = tab_cntrl(14)
    ztot0      = tab_cntrl(15)
    stot0      = tab_cntrl(16)
    ang0       = tab_cntrl(17)
    pa         = tab_cntrl(18)
    clon       = tab_cntrl(20)
    clat       = tab_cntrl(21)
    grossismx  = tab_cntrl(22)
    grossismy  = tab_cntrl(23)

    IF (tab_cntrl(24) == 1.)  THEN
       fxyhypb  = .TRUE.
    ELSE
       fxyhypb = .FALSE.
       ysinus  = .FALSE.
       IF (tab_cntrl(27) == 1.) ysinus = .TRUE. 
    ENDIF

    day_ini = tab_cntrl(30)
    itau_dyn = tab_cntrl(31)

    PRINT *, 'rad = ', rad
    PRINT *, 'omeg = ', omeg
    PRINT *, 'g = ', g
    PRINT *, 'cpp = ', cpp
    PRINT *, 'kappa = ', kappa

    IF (im /= iim)  THEN
       PRINT 1, im, iim
       STOP 1
    ELSE  IF (jm /= jjm)  THEN
       PRINT 2, jm, jjm
       STOP 1
    ELSE  IF (lllm /= llm)  THEN
       PRINT 3, lllm, llm
       STOP 1
    ENDIF

    call NF95_INQ_VARID (ncid, "rlonu", nvarid)
    ierr = NF90_GET_VAR(ncid, nvarid, rlonu)
    call handle_err("dynetat0, rlonu", ierr, ncid)

    call NF95_INQ_VARID (ncid, "rlatu", nvarid)
    ierr = NF90_GET_VAR(ncid, nvarid, rlatu)
    call handle_err("dynetat0, rlatu", ierr, ncid)

    call NF95_INQ_VARID (ncid, "rlonv", nvarid)
    ierr = NF90_GET_VAR(ncid, nvarid, rlonv)
    call handle_err("dynetat0, rlonv", ierr, ncid)

    call NF95_INQ_VARID (ncid, "rlatv", nvarid)
    ierr = NF90_GET_VAR(ncid, nvarid, rlatv)
    call handle_err("dynetat0, rlatv", ierr, ncid)

    call NF95_INQ_VARID (ncid, "cu", nvarid)
    ierr = NF90_GET_VAR(ncid, nvarid, cu_2d)
    call handle_err("dynetat0, cu", ierr, ncid)

    call NF95_INQ_VARID (ncid, "cv", nvarid)
    ierr = NF90_GET_VAR(ncid, nvarid, cv_2d)
    call handle_err("dynetat0, cv", ierr, ncid)

    call NF95_INQ_VARID (ncid, "aire", nvarid)
    ierr = NF90_GET_VAR(ncid, nvarid, aire_2d)
    call handle_err("dynetat0, aire", ierr, ncid)

    call NF95_INQ_VARID (ncid, "phisinit", nvarid)
    ierr = NF90_GET_VAR(ncid, nvarid, phis)
    call handle_err("dynetat0, phisinit", ierr, ncid)

    call NF95_INQ_VARID (ncid, "temps", nvarid)
    ierr = NF90_GET_VAR(ncid, nvarid, time_0)
    call handle_err("dynetat0, temps", ierr, ncid)

    call NF95_INQ_VARID (ncid, "ucov", nvarid)
    ierr = NF90_GET_VAR(ncid, nvarid, ucov, count=(/iim + 1, jjm + 1, llm/))
    call handle_err("dynetat0, ucov", ierr, ncid)

    call NF95_INQ_VARID (ncid, "vcov", nvarid)
    ierr = NF90_GET_VAR(ncid, nvarid, vcov, count=(/iim + 1, jjm, llm/))
    call handle_err("dynetat0, vcov", ierr, ncid)

    call NF95_INQ_VARID (ncid, "teta", nvarid)
    ierr = NF90_GET_VAR(ncid, nvarid, teta, count=(/iim + 1, jjm + 1, llm/))
    call handle_err("dynetat0, teta", ierr, ncid)

    DO iq = 1, nqmx
       call NF95_INQ_VARID(ncid, tname(iq), nvarid, ierr)
       IF (ierr  /=  NF90_NOERR) THEN
          PRINT *, 'dynetat0: le champ "' // tname(iq) // '" est absent, ' // &
               "il est donc initialisé à zéro."
          q(:, :, iq) = 0.
       ELSE
          ierr = NF90_GET_VAR(ncid, nvarid, q(:, :, iq), &
               count=(/iim + 1, jjm + 1, llm/))
          call handle_err("dynetat0, " // tname(iq), ierr, ncid)
       ENDIF
    ENDDO

    call NF95_INQ_VARID (ncid, "masse", nvarid)
    ierr = NF90_GET_VAR(ncid, nvarid, masse, count=(/iim + 1, jjm + 1, llm/))
    call handle_err("dynetat0, masse", ierr, ncid)

    call NF95_INQ_VARID (ncid, "ps", nvarid)
    ierr = NF90_GET_VAR(ncid, nvarid, ps, count=(/iim + 1, jjm + 1/))
    call handle_err("dynetat0, ps", ierr, ncid)

    call NF95_CLOSE(ncid)

    day_ini=day_ini+INT(time_0)
    time_0=time_0-INT(time_0)
    ! {0 <= time0 < 1}

1   FORMAT(//10x, 'la valeur de im =', i4, 2x, &
         'lue sur le fichier de demarrage est differente de la valeur ' &
         // 'parametree iim =', i4//)
2   FORMAT(//10x, 'la valeur de jm =', i4, 2x, &
         'lue sur le fichier de demarrage est differente de la valeur ' &
         // 'parametree jjm =', i4//)
3   FORMAT(//10x, 'la valeur de lmax =', i4, 2x, &
         'lue sur le fichier demarrage est differente de la valeur ' &
         // 'parametree llm =', i4//)

  END SUBROUTINE dynetat0

end module dynetat0_m
