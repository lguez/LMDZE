module caldyn0_m

  IMPLICIT NONE

  REAL, save, protected:: ang0, etot0, ptot0, ztot0, stot0

contains


  SUBROUTINE caldyn0(ucov, vcov, teta, ps, pk, phi)

    ! From dyn3d/caldyn0.F, version 1.1.1.1, 2004/05/19 12:53:07
    ! Authors:  P. Le Van, F. Forget
    ! Objet : calcul des tendances dynamiques

    use bernoui_m, only: bernoui
    USE comgeom, ONLY: airesurg
    use convmas_m, only: convmas
    use covcont_m, only: covcont
    USE dimensions, ONLY: iim, jjm, llm
    USE disvert_m, ONLY: ap, bp
    use enercin_m, only: enercin
    use flumass_m, only: flumass
    use massbar_m, only: massbar
    use massbarxy_m, only: massbarxy
    use massdair_m, only: massdair
    USE paramet_m, ONLY: iip1, ip1jmp1, jjp1, llmp1
    use sortvarc_m, only: sortvarc
    use tourpot_m, only: tourpot
    use vitvert_m, only: vitvert

    REAL, INTENT(IN):: ucov(:, :, :) ! (iim + 1, jjm + 1, llm) vent covariant
    REAL, INTENT(IN):: vcov(:, :, :) ! (iim + 1, jjm, llm) ! vent covariant
    REAL, INTENT(IN):: teta(ip1jmp1, llm)
    REAL, INTENT (IN):: ps(ip1jmp1)
    REAL, INTENT (IN):: pk(iip1, jjp1, llm)
    REAL, INTENT (IN):: phi(iim + 1, jjm + 1, llm)

    ! Local:
    real masse(iim + 1, jjm + 1, llm)
    REAL w(iim + 1, jjm + 1, llm)
    REAL pbaru(ip1jmp1, llm), pbarv((iim + 1) * jjm, llm)
    REAL vcont(iim + 1, jjm, llm), ucont(iim + 1, jjm + 1, llm)
    REAL p(ip1jmp1, llmp1)
    REAL massebx(ip1jmp1, llm), masseby((iim + 1) * jjm, llm)
    REAL vorpot(iim + 1, jjm, llm)
    real ecin(iim + 1, jjm + 1, llm), convm(iim + 1, jjm + 1, llm)
    REAL massebxy(iim + 1, jjm, llm), dp(iim + 1, jjm + 1)
    INTEGER l
    real ang, etot, ptot, ztot, stot, rmsdpdt, rmsv

    !-----------------------------------------------------------------------

    PRINT *, 'Call sequence information: caldyn0'

    CALL covcont(ucov, vcov, ucont, vcont)
    forall (l = 1: llm + 1) p(:, l) = ap(l) + bp(l) * ps
    CALL massdair(p, masse)
    CALL massbar(masse, massebx, masseby)
    CALL massbarxy(masse, massebxy)
    CALL flumass(massebx, masseby, vcont, ucont, pbaru, pbarv)
    CALL convmas(pbaru, pbarv, convm)
    dp = convm(:, :, 1) / airesurg
    w = vitvert(convm)
    CALL tourpot(vcov, ucov, massebxy, vorpot)
    CALL enercin(vcov, ucov, vcont, ucont, ecin)
    CALL sortvarc(ucov, teta, ps, masse, pk, vorpot, phi, bernoui(phi, ecin), &
         dp, ang, etot, ptot, ztot, stot, rmsdpdt, rmsv)
    PRINT *, 'ptot = ', ptot
    PRINT *, 'etot = ', etot
    PRINT *, 'ztot = ', ztot
    PRINT *, 'stot = ', stot
    PRINT *, 'ang = ', ang
    etot0 = etot
    ptot0 = ptot
    ztot0 = ztot
    stot0 = stot
    ang0  = ang

  END SUBROUTINE caldyn0

  !********************************************************************

  subroutine read_caldyn0(ncid_start)

    use netcdf95, only: NF95_GET_VAR, nf95_inq_varid

    integer, intent(in):: ncid_start

    ! Local:
    REAL tab_cntrl(5) ! partie du tableau des param\`etres du run
    integer varid

    !------------------------------------------------------------------

    call nf95_inq_varid(ncid_start, "controle", varid)
    call NF95_Get_VAR(ncid_start, varid, tab_cntrl, start = [13])

    etot0 = tab_cntrl(1)
    ptot0 = tab_cntrl(2)
    ztot0 = tab_cntrl(3)
    stot0 = tab_cntrl(4)
    ang0 = tab_cntrl(5)

  end subroutine read_caldyn0

end module caldyn0_m
