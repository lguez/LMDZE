module dynetat0_m

  ! This module is clean: no C preprocessor directive, no include line.

  IMPLICIT NONE

contains

  SUBROUTINE dynetat0(vcov, ucov, teta, q, masse, ps, phis, time)

    ! From dynetat0.F, version 1.2 2004/06/22 11:45:30

    ! Authors:  P. Le Van, L. Fairhead
    ! Objet : lecture de l'état initial

    use dimens_m, only: iim, jjm, llm, nqmx
    use comconst, only: im, cpp, dtvr, g, kappa, jm, lllm, omeg, rad
    use comvert, only: pa
    use logic, only: fxyhypb, ysinus
    use comgeom, only: rlonu, rlatu, rlonv, rlatv, cu_2d, cv_2d, aire_2d
    use serre, only: clon, clat, grossismy, grossismx
    use temps, only: day_ref, day_ini, itau_dyn, annee_ref
    use ener, only: etot0, ang0, ptot0, stot0, ztot0
    use advtrac_m, only: tname
    use netcdf95, only: nf95_open, nf95_inq_varid, handle_err, NF95_CLOSE
    use netcdf, only: NF90_NOWRITE, NF90_GET_VAR, NF90_NOERR
    use numer_rec, only: assert

    !   Arguments:
    REAL, intent(out):: vcov(: , :), ucov(:, :), teta(:, :)
    REAL, intent(out):: q(:, :, :), masse(:, :)
    REAL, intent(out):: ps(:) ! in Pa
    REAL, intent(out):: phis(:, :)
    REAL, intent(out):: time

    !   Variables 
    INTEGER length, iq
    PARAMETER (length = 100)
    REAL tab_cntrl(length) ! tableau des parametres du run
    INTEGER ierr, nid, nvarid

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
    call nf95_open("start.nc", NF90_NOWRITE, nid)

    call nf95_inq_varid(nid, "controle", nvarid)
    ierr = NF90_GET_VAR(nid, nvarid, tab_cntrl)
    call handle_err("dynetat0, controle", ierr, nid)

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

    call NF95_INQ_VARID (nid, "rlonu", nvarid)
    ierr = NF90_GET_VAR(nid, nvarid, rlonu)
    call handle_err("dynetat0, rlonu", ierr, nid)

    call NF95_INQ_VARID (nid, "rlatu", nvarid)
    ierr = NF90_GET_VAR(nid, nvarid, rlatu)
    call handle_err("dynetat0, rlatu", ierr, nid)

    call NF95_INQ_VARID (nid, "rlonv", nvarid)
    ierr = NF90_GET_VAR(nid, nvarid, rlonv)
    call handle_err("dynetat0, rlonv", ierr, nid)

    call NF95_INQ_VARID (nid, "rlatv", nvarid)
    ierr = NF90_GET_VAR(nid, nvarid, rlatv)
    call handle_err("dynetat0, rlatv", ierr, nid)

    call NF95_INQ_VARID (nid, "cu", nvarid)
    ierr = NF90_GET_VAR(nid, nvarid, cu_2d)
    call handle_err("dynetat0, cu", ierr, nid)

    call NF95_INQ_VARID (nid, "cv", nvarid)
    ierr = NF90_GET_VAR(nid, nvarid, cv_2d)
    call handle_err("dynetat0, cv", ierr, nid)

    call NF95_INQ_VARID (nid, "aire", nvarid)
    ierr = NF90_GET_VAR(nid, nvarid, aire_2d)
    call handle_err("dynetat0, aire", ierr, nid)

    call NF95_INQ_VARID (nid, "phisinit", nvarid)
    ierr = NF90_GET_VAR(nid, nvarid, phis)
    call handle_err("dynetat0, phisinit", ierr, nid)

    call NF95_INQ_VARID (nid, "temps", nvarid)
    ierr = NF90_GET_VAR(nid, nvarid, time)
    call handle_err("dynetat0, temps", ierr, nid)

    call NF95_INQ_VARID (nid, "ucov", nvarid)
    ierr = NF90_GET_VAR(nid, nvarid, ucov, count=(/iim + 1, jjm + 1, llm/))
    call handle_err("dynetat0, ucov", ierr, nid)

    call NF95_INQ_VARID (nid, "vcov", nvarid)
    ierr = NF90_GET_VAR(nid, nvarid, vcov, count=(/iim + 1, jjm, llm/))
    call handle_err("dynetat0, vcov", ierr, nid)

    call NF95_INQ_VARID (nid, "teta", nvarid)
    ierr = NF90_GET_VAR(nid, nvarid, teta, count=(/iim + 1, jjm + 1, llm/))
    call handle_err("dynetat0, teta", ierr, nid)

    DO iq = 1, nqmx
       call NF95_INQ_VARID(nid, tname(iq), nvarid, ierr)
       IF (ierr  /=  NF90_NOERR) THEN
          PRINT *, 'dynetat0: le champ "' // tname(iq) // '" est absent, ' // &
               "il est donc initialisé à zéro."
          q(:, :, iq) = 0.
       ELSE
          ierr = NF90_GET_VAR(nid, nvarid, q(:, :, iq), &
               count=(/iim + 1, jjm + 1, llm/))
          call handle_err("dynetat0, " // tname(iq), ierr, nid)
       ENDIF
    ENDDO

    call NF95_INQ_VARID (nid, "masse", nvarid)
    ierr = NF90_GET_VAR(nid, nvarid, masse, count=(/iim + 1, jjm + 1, llm/))
    call handle_err("dynetat0, masse", ierr, nid)

    call NF95_INQ_VARID (nid, "ps", nvarid)
    ierr = NF90_GET_VAR(nid, nvarid, ps, count=(/iim + 1, jjm + 1/))
    call handle_err("dynetat0, ps", ierr, nid)

    call NF95_CLOSE(nid)

    day_ini=day_ini+INT(time)
    time=time-INT(time)

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
