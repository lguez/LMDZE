module etat0dyn_m

  IMPLICIT NONE

contains

  SUBROUTINE etat0dyn(tsol_2d)

    ! From "etat0_netcdf.F", version 1.3, 2005/05/25 13:10:09

    use jumble, only: assert
    use netcdf95, only: nf95_close, nf95_get_var, nf95_gw_var, nf95_put_var, &
         nf95_inq_varid, nf95_open

    use caldyn0_m, only: caldyn0
    use comgeom, only: aire_2d, apoln, apols, cu_2d, cv_2d
    use dimensions, only: iim, jjm, llm, nqmx
    use disvert_m, only: ap, bp, preff
    use dynetat0_m, only: rlatu, rlatv, rlonu, rlonv
    use dynetat0_chosen_m, only: day_ref
    use dynredem0_m, only: dynredem0
    use dynredem1_m, only: dynredem1
    use exner_hyb_m, only: exner_hyb
    use geopot_m, only: geopot
    use massdair_m, only: massdair
    use q_sat_m, only: q_sat
    use regr_lat_time_coefoz_m, only: regr_lat_time_coefoz
    use regr_pr_o3_m, only: regr_pr_o3
    use start_init_dyn_m, only: start_init_dyn
    use start_inter_3d_m, only: start_inter_3d
    use suphec_m, only: rcpd, rkappa

    real, intent(in):: tsol_2d(:, :) ! (iim + 1, jjm + 1)
    ! both soil temperature and surface temperature, in K

    ! Local:

    REAL, dimension(iim + 1, jjm + 1, llm):: ucov, t3d, teta
    REAL vcov(iim + 1, jjm, llm)

    REAL q(iim + 1, jjm + 1, llm, nqmx)
    ! Mass fractions of trace species. "q(i, j, l)" is at longitude
    ! "rlonv(i)", latitude "rlatu(j)" and pressure level "pls(i, j,
    ! l)".

    real ps(iim + 1, jjm + 1)
    INTEGER l
    REAL pk(iim + 1, jjm + 1, llm) ! fonction d'Exner aux milieux des couches 
    real pks(iim + 1, jjm + 1)
    REAL masse(iim + 1, jjm + 1, llm)
    REAL phi(iim + 1, jjm + 1, llm)

    real pls(iim + 1, jjm + 1, llm)
    ! (pressure at mid-layer of LMDZ grid, in Pa)
    ! "pls(i, j, l)" is at longitude "rlonv(i)", latitude "rlatu(j)",
    ! for layer "l")

    REAL p3d(iim + 1, jjm + 1, llm+1) ! pressure at layer interfaces, in Pa
    ! ("p3d(i, j, l)" is at longitude "rlonv(i)", latitude "rlatu(j)",
    ! for interface "l")

    !---------------------------------

    print *, "Call sequence information: etat0dyn"
    CALL start_init_dyn(tsol_2d, ps)

    ! Compute pressure on intermediate levels:
    forall(l = 1: llm + 1) p3d(:, :, l) = ap(l) + bp(l) * ps
    CALL exner_hyb(ps, p3d, pks, pk)
    call assert(MINVAL(pk) /= MAXVAL(pk), '"pk" should not be constant')

    pls = preff * (pk / rcpd)**(1. / rkappa)
    PRINT *, "minval(pls) = ", minval(pls)
    print *, "maxval(pls) = ", maxval(pls)

    call start_inter_3d('U', rlonv, rlatv, pls, ucov)
    forall (l = 1: llm) ucov(:iim, :, l) = ucov(:iim, :, l) * cu_2d(:iim, :)
    ucov(iim+1, :, :) = ucov(1, :, :)

    call start_inter_3d('V', rlonu, rlatu(:jjm), pls(:, :jjm, :), vcov)
    forall (l = 1: llm) vcov(:iim, :, l) = vcov(:iim, :, l) * cv_2d(:iim, :)
    vcov(iim + 1, :, :) = vcov(1, :, :)

    call start_inter_3d('TEMP', rlonu, rlatv, pls, t3d)
    PRINT *, 'minval(t3d) = ', minval(t3d)
    print *, "maxval(t3d) = ", maxval(t3d)

    teta(:iim, :, :) = t3d(:iim, :, :) * rcpd / pk(:iim, :, :)
    teta(iim + 1, :, :) = teta(1, :, :)
    DO l = 1, llm
       teta(:, 1, l) = SUM(aire_2d(:, 1) * teta(:, 1, l)) / apoln
       teta(:, jjm + 1, l) = SUM(aire_2d(:, jjm + 1) * teta(:, jjm + 1, l)) &
            / apols
    ENDDO

    ! Water vapor:
    call start_inter_3d('R', rlonu, rlatv, pls, q(:, :, :, 1))
    q(:, :, :, 1) = 0.01 * q(:, :, :, 1) * q_sat(t3d, pls)
    WHERE (q(:, :, :, 1) < 0.) q(:, :, :, 1) = 1E-10
    
    forall (l = 1:llm)
       q(:, 1, l, 1) = SUM(aire_2d(:, 1) * q(:, 1, l, 1)) / apoln
       q(:, jjm + 1, l, 1) &
            = SUM(aire_2d(:, jjm + 1) * q(:, jjm + 1, l, 1)) / apols
    END forall

    q(:, :, :, 2:4) = 0. ! liquid water, radon and lead

    if (nqmx >= 5) then
       ! Ozone:
       call regr_lat_time_coefoz
       call regr_pr_o3(p3d, q(:, :, :, 5))
       ! Convert from mole fraction to mass fraction:
       q(:, :, :, 5) = q(:, :, :, 5) * 48. / 29.
    end if

    ! Calcul interm\'ediaire :
    CALL massdair(p3d, masse)

    forall (l = 1:llm)
       masse(:, 1, l) = SUM(aire_2d(:iim, 1) * masse(:iim, 1, l)) / apoln
       masse(:, jjm + 1, l) = &
            SUM(aire_2d(:iim, jjm + 1) * masse(:iim, jjm + 1, l)) / apols
    END forall

    CALL geopot(teta, pk , pks, phi)
    CALL caldyn0(ucov, vcov, teta, ps, pk, phi)
    CALL dynredem0(day_ref)
    CALL dynredem1(vcov, ucov, teta, q, masse, ps, itau = 0)

  END SUBROUTINE etat0dyn

end module etat0dyn_m
