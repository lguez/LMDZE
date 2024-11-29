MODULE dynredem_m

  IMPLICIT NONE

  INTEGER, protected:: ncid

CONTAINS

  SUBROUTINE dynredem(vcov, ucov, teta, q, masse, ps, iday_end, itau)

    ! From dyn3d/dynredem.F, version 1.2, 2004/06/22 11:45:30
    ! \'Ecriture du fichier de red\'emarrage au format NetCDF

    use jumble, only: assert
    USE netcdf95, ONLY: nf95_clobber, nf95_float, nf95_global, nf95_unlimited, &
         nf95_create, nf95_def_dim, nf95_def_var, nf95_enddef, nf95_put_att, &
         nf95_put_var, nf95_close, nf95_inq_varid

    use caldyn0_m, only: ang0, etot0, ptot0, stot0, ztot0
    USE comconst, ONLY: daysec, ra
    USE dimensions, ONLY: iim, jjm, llm, nqmx
    USE disvert_m, ONLY: ap, bp, preff, presnivs
    use dynetat0_m, only: rlatu, rlatv, rlonu, rlonv, rlatu1, rlatu2, yprimu1, &
         yprimu2, xprimp025, xprimm025, xprimu, xprimv
    use dynetat0_chosen_m, only: day_ref, annee_ref, clat, clon, dzoomx, &
         dzoomy, grossismx, grossismy, taux, tauy
    use grid_noro_m, only: phis
    USE infotrac_init_m, ONLY: tname, ttext
    USE ju2ymds_m, ONLY: ju2ymds
    USE paramet_m, ONLY: iip1, jjp1, llmp1
    use suphec_m, only: rg, rcpd, rkappa, romega
    use ymds2ju_m, only: ymds2ju

    REAL, INTENT(IN):: vcov(:, :, :) ! (iim + 1, jjm, llm)
    REAL, INTENT(IN):: ucov(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, INTENT(IN):: teta(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, INTENT(IN):: q(:, :, :, :) ! (iim + 1, jjm + 1, llm, nqmx)
    REAL, INTENT(IN):: masse(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, INTENT(IN):: ps(:, :) ! (iim + 1, jjm + 1)
    INTEGER, INTENT(IN):: iday_end
    INTEGER, INTENT(IN):: itau

    ! Local:

    INTEGER iq
    INTEGER, PARAMETER:: length = 100
    REAL tab_cntrl(length) ! tableau des param\`etres du run

    ! Pour NetCDF :
    INTEGER idim_index
    INTEGER idim_rlonu, idim_rlonv, idim_rlatu, idim_rlatv
    INTEGER idim_s, idim_sig
    INTEGER dimid_temps
    INTEGER varid
    integer varid_controle, varid_rlonu, varid_rlatu, varid_rlonv, varid_rlatv
    integer varid_xprimu, varid_xprimv, varid_xprimm025, varid_xprimp025
    integer varid_rlatu1, varid_rlatu2, varid_yprimu1, varid_yprimu2, varid_ap
    integer varid_bp, varid_presnivs, varid_phis

    double precision julian
    INTEGER year, day, month
    CHARACTER(len=30) unites

    !-----------------------------------------------------------------------

    PRINT *, 'Call sequence information: dynredem'
    call assert((/size(vcov, 1), size(ucov, 1), size(teta, 1), size(q, 1), &
         size(masse, 1), size(ps, 1)/) == iim + 1, "dynredem1 iim")
    call assert((/size(vcov, 2) + 1, size(ucov, 2), size(teta, 2), size(q, 2), &
         size(masse, 2), size(ps, 2)/) == jjm + 1, "dynredem1 jjm")
    call assert((/size(vcov, 3), size(ucov, 3), size(teta, 3), size(q, 3), &
         size(masse, 3)/) == llm, "dynredem1 llm")
    call assert(size(q, 4) == nqmx, "dynredem1 nqmx")
    CALL ymds2ju(annee_ref, 1, iday_end, 0., julian)
    CALL ju2ymds(julian, year, month, day)
    tab_cntrl(1) = iim
    tab_cntrl(2) = jjm
    tab_cntrl(3) = llm
    tab_cntrl(4) = day_ref
    tab_cntrl(5) = annee_ref
    tab_cntrl(6) = ra
    tab_cntrl(7) = romega
    tab_cntrl(8) = rg
    tab_cntrl(9) = rcpd
    tab_cntrl(10) = rkappa
    tab_cntrl(11) = daysec
    tab_cntrl(12) = 0.
    tab_cntrl(13) = etot0
    tab_cntrl(14) = ptot0
    tab_cntrl(15) = ztot0
    tab_cntrl(16) = stot0
    tab_cntrl(17) = ang0
    tab_cntrl(18) = 0.
    tab_cntrl(19) = preff

    ! Param\`etres pour le zoom :
    tab_cntrl(20) = clon
    tab_cntrl(21) = clat
    tab_cntrl(22) = grossismx
    tab_cntrl(23) = grossismy
    tab_cntrl(24) = 1.
    tab_cntrl(25) = dzoomx
    tab_cntrl(26) = dzoomy
    tab_cntrl(27) = 0.
    tab_cntrl(28) = taux
    tab_cntrl(29) = tauy

    tab_cntrl(30) = iday_end
    tab_cntrl(31) = itau
    tab_cntrl(32:) = 0.

    CALL nf95_create("restart.nc", nf95_clobber, ncid)
    CALL nf95_put_att(ncid, nf95_global, 'title', &
         'start file for the dynamics code')

    ! Definir les dimensions du fichiers:

    CALL nf95_def_dim(ncid, 'index', length, idim_index)
    CALL nf95_def_dim(ncid, 'rlonu', iip1, idim_rlonu)
    CALL nf95_def_dim(ncid, 'rlatu', jjp1, idim_rlatu)
    CALL nf95_def_dim(ncid, 'rlonv', iip1, idim_rlonv)
    CALL nf95_def_dim(ncid, 'rlatv', jjm, idim_rlatv)
    CALL nf95_def_dim(ncid, 'sigs', llm, idim_s)
    CALL nf95_def_dim(ncid, 'sig', llmp1, idim_sig)
    CALL nf95_def_dim(ncid, 'temps', nf95_unlimited, dimid_temps)

    ! Definir et enregistrer certains champs invariants:

    CALL nf95_def_var(ncid, 'controle', nf95_float, idim_index, varid_controle)
    CALL nf95_put_att(ncid, varid_controle, 'title', 'Parametres de controle')

    CALL nf95_def_var(ncid, 'rlonu', nf95_float, idim_rlonu, varid_rlonu)
    CALL nf95_put_att(ncid, varid_rlonu, 'title', 'Longitudes des points U')

    CALL nf95_def_var(ncid, 'rlatu', nf95_float, idim_rlatu, varid_rlatu)
    CALL nf95_put_att(ncid, varid_rlatu, 'title', 'Latitudes des points U')

    CALL nf95_def_var(ncid, 'rlonv', nf95_float, idim_rlonv, varid_rlonv)
    CALL nf95_put_att(ncid, varid_rlonv, 'title', 'Longitudes des points V')

    CALL nf95_def_var(ncid, 'rlatv', nf95_float, idim_rlatv, varid_rlatv)
    CALL nf95_put_att(ncid, varid_rlatv, 'title', 'Latitudes des points V')

    CALL nf95_def_var(ncid, 'xprimu', nf95_float, idim_rlonu, varid_xprimu)
    CALL nf95_put_att(ncid, varid_xprimu, 'title', 'dx / dX aux points u')

    CALL nf95_def_var(ncid, 'xprimv', nf95_float, idim_rlonv, varid_xprimv)
    CALL nf95_put_att(ncid, varid_xprimv, 'title', 'dx / dX aux points v')

    CALL nf95_def_var(ncid, 'xprimm025', nf95_float, idim_rlonu, &
         varid_xprimm025)
    CALL nf95_def_var(ncid, 'xprimp025', nf95_float, idim_rlonu, &
         varid_xprimp025)

    CALL nf95_def_var(ncid, 'rlatu1', nf95_float, idim_rlatv, varid_rlatu1)
    CALL nf95_def_var(ncid, 'rlatu2', nf95_float, idim_rlatv, varid_rlatu2)
    CALL nf95_def_var(ncid, 'yprimu1', nf95_float, idim_rlatv, varid_yprimu1)
    CALL nf95_def_var(ncid, 'yprimu2', nf95_float, idim_rlatv, varid_yprimu2)

    CALL nf95_def_var(ncid, 'ap', nf95_float, idim_sig, varid_ap)
    CALL nf95_put_att(ncid, varid_ap, 'title', 'Coefficient A pour hybride')

    CALL nf95_def_var(ncid, 'bp', nf95_float, idim_sig, varid_bp)
    CALL nf95_put_att(ncid, varid_bp, 'title', 'Coefficient B pour hybride')

    CALL nf95_def_var(ncid, 'presnivs', nf95_float, idim_s, varid_presnivs)

    ! Geopentiel au sol:

    CALL nf95_def_var(ncid, 'phis', nf95_float, (/idim_rlonv, idim_rlatu/), &
         varid_phis)
    CALL nf95_put_att(ncid, varid_phis, 'standard_name', 'surface_geopotential')
    CALL nf95_put_att(ncid, varid_phis, 'units', 'm2 s-2')

    ! Definir les variables pour pouvoir les enregistrer plus tard:

    CALL nf95_def_var(ncid, 'temps', nf95_float, dimid_temps, varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Temps de simulation')
    WRITE(unites, fmt = 200) year, month, day
200 FORMAT ('days since ', I4, '-', I2.2, '-', I2.2, ' 00:00:00')
    CALL nf95_put_att(ncid, varid, 'units', unites)

    CALL nf95_def_var(ncid, 'ucov', nf95_float, &
         (/idim_rlonu, idim_rlatu, idim_s, dimid_temps/), varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Vitesse U')

    CALL nf95_def_var(ncid, 'vcov', nf95_float, &
         (/idim_rlonv, idim_rlatv, idim_s, dimid_temps/), varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Vitesse V')

    CALL nf95_def_var(ncid, 'teta', nf95_float, &
         (/idim_rlonv, idim_rlatu, idim_s, dimid_temps/), varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Temperature')

    DO iq = 1, nqmx
       CALL nf95_def_var(ncid, tname(iq), nf95_float, &
            (/idim_rlonv, idim_rlatu, idim_s, dimid_temps/), varid)
       CALL nf95_put_att(ncid, varid, 'title', ttext(iq))
    END DO

    CALL nf95_def_var(ncid, 'masse', nf95_float, &
         (/idim_rlonv, idim_rlatu, idim_s, dimid_temps/), varid)
    CALL nf95_put_att(ncid, varid, 'title', 'C est quoi ?')

    CALL nf95_def_var(ncid, 'ps', nf95_float, &
         (/idim_rlonv, idim_rlatu, dimid_temps/), varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Pression au sol')

    CALL nf95_enddef(ncid)

    CALL nf95_put_var(ncid, varid_controle, tab_cntrl)
    CALL nf95_put_var(ncid, varid_rlonu, rlonu)
    CALL nf95_put_var(ncid, varid_rlatu, rlatu)
    CALL nf95_put_var(ncid, varid_rlonv, rlonv)
    CALL nf95_put_var(ncid, varid_rlatv, rlatv)
    CALL nf95_put_var(ncid, varid_xprimu, xprimu)
    CALL nf95_put_var(ncid, varid_xprimv, xprimv)
    CALL nf95_put_var(ncid, varid_xprimm025, xprimm025)
    CALL nf95_put_var(ncid, varid_xprimp025, xprimp025)
    call NF95_PUT_VAR(ncid, varid_rlatu1, rlatu1)
    call NF95_PUT_VAR(ncid, varid_rlatu2, rlatu2)
    CALL nf95_put_var(ncid, varid_yprimu1, yprimu1)
    CALL nf95_put_var(ncid, varid_yprimu2, yprimu2)
    CALL nf95_put_var(ncid, varid_ap, ap)
    CALL nf95_put_var(ncid, varid_bp, bp)
    CALL nf95_put_var(ncid, varid_presnivs, presnivs)
    CALL nf95_put_var(ncid, varid_phis, phis)

    ! \'Ecriture de la coordonn\'ee temps:
    call nf95_inq_varid(ncid, 'temps', varid)
    call nf95_put_var(ncid, varid, values = 0.)

    ! \'Ecriture des champs

    call nf95_inq_varid(ncid, 'ucov', varid)
    call nf95_put_var(ncid, varid, ucov)
    call nf95_inq_varid(ncid, 'vcov', varid)
    call nf95_put_var(ncid, varid, vcov)
    call nf95_inq_varid(ncid, 'teta', varid)
    call nf95_put_var(ncid, varid, teta)

    DO iq = 1, nqmx
       call nf95_inq_varid(ncid, tname(iq), varid)
       call nf95_put_var(ncid, varid, q(:, :, :, iq))
    END DO

    call nf95_inq_varid(ncid, 'masse', varid)
    call nf95_put_var(ncid, varid, masse)
    call nf95_inq_varid(ncid, 'ps', varid)
    call nf95_put_var(ncid, varid, ps)
    call nf95_close(ncid)

  END SUBROUTINE dynredem

END MODULE dynredem_m
