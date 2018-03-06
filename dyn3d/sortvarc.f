module sortvarc_m

  IMPLICIT NONE

  real, save:: ang, etot, ptot, ztot, stot, rmsdpdt, rmsv

contains

  SUBROUTINE sortvarc(ucov, teta, ps, masse, pk, phis, vorpot, phi, bern, dp, &
       resetvarc)

    ! From dyn3d/sortvarc.F, version 1.1.1.1, 2004/05/19 12:53:07
    ! Author: P. Le Van
    ! Objet : sortie des variables de contr\^ole

    USE comconst, ONLY: daysec, g, omeg, rad
    USE comgeom, ONLY: aire_2d, cu_2d
    USE dimens_m, ONLY: iim, jjm, llm
    use dynetat0_m, ONLY: rlatu
    USE ener, ONLY: ang0, etot0, ptot0, stot0, ztot0
    use filtreg_scal_m, only: filtreg_scal
    use massbarxy_m, only: massbarxy
    USE paramet_m, ONLY: jjp1

    REAL, INTENT(IN):: ucov(iim + 1, jjm + 1, llm)
    REAL, INTENT(IN):: teta(iim + 1, jjm + 1, llm)
    REAL, INTENT(IN):: ps(iim + 1, jjm + 1)
    REAL, INTENT(IN):: masse(iim + 1, jjm + 1, llm)
    REAL, INTENT(IN):: pk(iim + 1, jjm + 1, llm)
    REAL, INTENT(IN):: phis(iim + 1, jjm + 1)
    REAL, INTENT(IN):: vorpot(:, :, :) ! (iim + 1, jjm, llm)
    REAL, intent(in):: phi(iim + 1, jjm + 1, llm)
    real, intent(in):: bern(iim + 1, jjm + 1, llm)
    REAL, intent(in):: dp(iim + 1, jjm + 1)
    logical, intent(in):: resetvarc

    ! Local:
    REAL bernf(iim + 1, jjm + 1, llm)
    REAL etotl(llm), angl(llm), ge(iim, 2:jjm)
    REAL cosphi(2:jjm)
    REAL radsg, radomeg
    REAL massebxy(iim + 1, jjm, llm)
    INTEGER j, l

    !-----------------------------------------------------------------------

    PRINT *, "Call sequence information: sortvarc"

    rmsdpdt = daysec * 0.01 * sqrt(sum(dp(:iim, :)**2) / (iim * jjp1))

    ! Calcul du moment  angulaire :
    
    radsg = rad / g
    radomeg = rad * omeg
    cosphi = cos(rlatu(2:jjm))

    DO l = 1, llm
       forall (j = 2:jjm) ge(:, j) = (ucov(:iim, j, l) / cu_2d(:iim, j) &
            + radomeg * cosphi(j)) * masse(:iim, j, l) * cosphi(j)
       angl(l) = radsg * sum(ge)
    END DO

    ang = sum(angl)

    ! Calcul  de l'energie, de l'enstrophie, de l'entropie et de rmsv :

    bernf = bern
    CALL filtreg_scal(bernf, direct = .false., intensive = .false.)

    ptot = sum(ps(:iim, :) * aire_2d(:iim, :))

    forall (l = 1:llm) etotl(l) = sum(masse(:iim, :, l) * (phis(:iim, :) &
         + teta(:iim, :, l) * pk(:iim, :, l) + bernf(:iim, :, l) &
         - phi(:iim, :, l)))
    etot = sum(etotl)

    CALL massbarxy(masse, massebxy)
    ztot = sum(vorpot(:iim, :, :)**2 * massebxy(:iim, :, :))

    stot = sum(masse(:iim, :, :) * teta(:iim, :, :))
    rmsv = 2. &
         * sum(masse(:iim, :, :) * max(bernf(:iim, :, :) - phi(:iim, :, :), 0.))

    IF (resetvarc .or. ptot0 == 0.) then
       print *, 'sortvarc: recomputed initial values.'
       etot0 = etot
       ptot0 = ptot
       ztot0 = ztot
       stot0 = stot
       ang0  = ang
       PRINT *, 'ptot0 = ', ptot0
       PRINT *, 'etot0 = ', etot0
       PRINT *, 'ztot0 = ', ztot0
       PRINT *, 'stot0 = ', stot0
       PRINT *, 'ang0 = ', ang0
    END IF

    IF (.not. resetvarc) then
       etot = etot / etot0
       rmsv = sqrt(rmsv / ptot)
       ptot = ptot / ptot0
       ztot = ztot / ztot0
       stot = stot / stot0
       ang = ang / ang0
    end IF

  END SUBROUTINE sortvarc

end module sortvarc_m
