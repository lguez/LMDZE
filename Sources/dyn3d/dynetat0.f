module dynetat0_m

  use dimens_m, only: iim, jjm

  IMPLICIT NONE

  private iim, jjm

  INTEGER day_ini 
  ! day number at the beginning of the run, based at value 1 on
  ! January 1st of annee_ref

  integer:: day_ref = 1 ! jour de l'ann\'ee de l'\'etat initial
  ! (= 350 si 20 d\'ecembre par exemple)

  integer:: annee_ref = 1998 ! Annee de l'etat initial (avec 4 chiffres)

  REAL clon ! longitude of the center of the zoom, in rad
  real clat ! latitude of the center of the zoom, in rad

  real grossismx, grossismy
  ! facteurs de grossissement du zoom, selon la longitude et la latitude
  ! = 2 si 2 fois, = 3 si 3 fois, etc.

  real dzoomx, dzoomy
  ! extensions en longitude et latitude de la zone du zoom (fractions
  ! de la zone totale)

  real taux, tauy
  ! raideur de la transition de l'int\'erieur \`a l'ext\'erieur du zoom
  
  real rlatu(jjm + 1)
  ! (latitudes of points of the "scalar" and "u" grid, in rad)

  real rlatv(jjm) 
  ! (latitudes of points of the "v" grid, in rad, in decreasing order)

  real rlonu(iim + 1) ! longitudes of points of the "u" grid, in rad

  real rlonv(iim + 1)
  ! (longitudes of points of the "scalar" and "v" grid, in rad)

  real xprimu(iim + 1), xprimv(iim + 1)
  ! xprim[uv] = 2 pi / iim * (derivative of the longitudinal zoom
  ! function)(rlon[uv])

  REAL xprimm025(iim + 1), xprimp025(iim + 1)
  REAL rlatu1(jjm), rlatu2(jjm), yprimu1(jjm), yprimu2(jjm)

  save

contains

  SUBROUTINE dynetat0(vcov, ucov, teta, q, masse, ps, phis)

    ! From dynetat0.F, version 1.2, 2004/06/22 11:45:30
    ! Authors: P. Le Van, L. Fairhead
    ! This procedure reads the initial state of the atmosphere.

    use comconst, only: dtvr
    use conf_gcm_m, only: raz_date
    use dimens_m, only: iim, jjm, llm, nqmx
    use disvert_m, only: pa
    use ener, only: etot0, ang0, ptot0, stot0, ztot0
    use iniadvtrac_m, only: tname
    use netcdf, only: NF90_NOWRITE, NF90_NOERR
    use netcdf95, only: NF95_GET_VAR, nf95_open, nf95_inq_varid, NF95_CLOSE, &
         NF95_Gw_VAR
    use nr_util, only: assert
    use temps, only: itau_dyn
    use unit_nml_m, only: unit_nml

    REAL, intent(out):: vcov(: , :, :) ! (iim + 1, jjm, llm)
    REAL, intent(out):: ucov(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(out):: teta(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(out):: q(:, :, :, :) ! (iim + 1, jjm + 1, llm, nqmx)
    REAL, intent(out):: masse(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(out):: ps(:, :) ! (iim + 1, jjm + 1) in Pa
    REAL, intent(out):: phis(:, :) ! (iim + 1, jjm + 1)

    ! Local variables: 
    INTEGER iq
    REAL, pointer:: tab_cntrl(:) ! tableau des param\`etres du run
    INTEGER ierr, ncid, varid

    namelist /dynetat0_nml/ day_ref, annee_ref

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

    ! Fichier \'etat initial :
    call nf95_open("start.nc", NF90_NOWRITE, ncid)

    call nf95_inq_varid(ncid, "controle", varid)
    call NF95_Gw_VAR(ncid, varid, tab_cntrl)

    call assert(int(tab_cntrl(1)) == iim, "dynetat0 tab_cntrl iim") 
    call assert(int(tab_cntrl(2)) == jjm, "dynetat0 tab_cntrl jjm") 
    call assert(int(tab_cntrl(3)) == llm, "dynetat0 tab_cntrl llm") 

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
    dzoomx = tab_cntrl(25)
    dzoomy = tab_cntrl(26)
    taux = tab_cntrl(28)
    tauy = tab_cntrl(29)

    print *, "Enter namelist 'dynetat0_nml'."
    read(unit=*, nml=dynetat0_nml)
    write(unit_nml, nml=dynetat0_nml)

    if (raz_date) then
       print *, 'Resetting the date, using the namelist.'
       day_ini = day_ref
       itau_dyn = 0
    else
       day_ref = tab_cntrl(4)
       annee_ref = tab_cntrl(5)
       itau_dyn = tab_cntrl(31)
       day_ini = tab_cntrl(30)
    end if

    print *, "day_ini = ", day_ini

    deallocate(tab_cntrl) ! pointer

    call NF95_INQ_VARID (ncid, "rlonu", varid)
    call NF95_GET_VAR(ncid, varid, rlonu)

    call NF95_INQ_VARID (ncid, "rlatu", varid)
    call NF95_GET_VAR(ncid, varid, rlatu)

    call NF95_INQ_VARID (ncid, "rlonv", varid)
    call NF95_GET_VAR(ncid, varid, rlonv)

    call NF95_INQ_VARID (ncid, "rlatv", varid)
    call NF95_GET_VAR(ncid, varid, rlatv)

    CALL nf95_inq_varid(ncid, 'xprimu', varid)
    CALL nf95_get_var(ncid, varid, xprimu)

    CALL nf95_inq_varid(ncid, 'xprimv', varid)
    CALL nf95_get_var(ncid, varid, xprimv)

    CALL nf95_inq_varid(ncid, 'xprimm025', varid)
    CALL nf95_get_var(ncid, varid, xprimm025)

    CALL nf95_inq_varid(ncid, 'xprimp025', varid)
    CALL nf95_get_var(ncid, varid, xprimp025)

    call NF95_INQ_VARID (ncid, "rlatu1", varid)
    call NF95_GET_VAR(ncid, varid, rlatu1)

    call NF95_INQ_VARID (ncid, "rlatu2", varid)
    call NF95_GET_VAR(ncid, varid, rlatu2)

    CALL nf95_inq_varid(ncid, 'yprimu1', varid)
    CALL nf95_get_var(ncid, varid, yprimu1)

    CALL nf95_inq_varid(ncid, 'yprimu2', varid)
    CALL nf95_get_var(ncid, varid, yprimu2)

    call NF95_INQ_VARID (ncid, "phisinit", varid)
    call NF95_GET_VAR(ncid, varid, phis)

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
