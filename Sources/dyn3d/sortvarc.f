module sortvarc_m

  IMPLICIT NONE

  real, save:: ang, etot, ptot, ztot, stot, rmsdpdt, rmsv

contains

  SUBROUTINE sortvarc(ucov, teta, ps, masse, pk, phis, vorpot, phi, &
       bern, dp, resetvarc)

    ! From dyn3d/sortvarc.F, version 1.1.1.1 2004/05/19 12:53:07
    ! Author: P. Le Van
    ! Objet : sortie des variables de contr\^ole

    USE comconst, ONLY: daysec, g, omeg, rad
    USE comgeom, ONLY: aire_2d, cu_2d
    USE dimens_m, ONLY: iim, jjm, llm
    use dynetat0_m, ONLY: rlatu
    USE ener, ONLY: ang0, etot0, ptot0, stot0, ztot0
    use filtreg_scal_m, only: filtreg_scal
    use massbarxy_m, only: massbarxy
    USE paramet_m, ONLY: iip1, ip1jm, jjp1

    REAL, INTENT(IN):: ucov(iim + 1, jjm + 1, llm)
    REAL, INTENT(IN):: teta(iim + 1, jjm + 1, llm)
    REAL, INTENT(IN):: ps(iim + 1, jjm + 1)
    REAL, INTENT(IN):: masse(iim + 1, jjm + 1, llm)
    REAL, INTENT(IN):: pk(iim + 1, jjm + 1, llm)
    REAL, INTENT(IN):: phis(iim + 1, jjm + 1)
    REAL, INTENT(IN):: vorpot(ip1jm, llm)
    REAL, intent(in):: phi(iim + 1, jjm + 1, llm)
    real, intent(in):: bern(iim + 1, jjm + 1, llm)
    REAL, intent(in):: dp(iim + 1, jjm + 1)
    logical, intent(in):: resetvarc

    ! Local:
    REAL vor(ip1jm), bernf(iim + 1, jjm + 1, llm), ztotl(llm)
    REAL etotl(llm), stotl(llm), rmsvl(llm), angl(llm), ge(iim + 1, jjm + 1)
    REAL cosphi(2:jjm)
    REAL radsg, radomeg
    REAL massebxy(ip1jm, llm)
    INTEGER j, l, ij
    REAL ssum

    !-----------------------------------------------------------------------

    PRINT *, "Call sequence information: sortvarc"

    CALL massbarxy(masse, massebxy)

    ! Calcul  de  rmsdpdt
    ge = dp*dp
    rmsdpdt = sum(ge) - sum(ge(1, :))
    rmsdpdt = daysec*1.E-2*sqrt(rmsdpdt / (iim * jjp1))
    bernf = bern
    CALL filtreg_scal(bernf, direct = .false., intensive = .false.)

    ! Calcul du moment  angulaire
    radsg = rad/g
    radomeg = rad*omeg
    cosphi = cos(rlatu(2:jjm))

    ! Calcul  de l'energie, de l'enstrophie, de l'entropie et de rmsv

    DO l = 1, llm
       DO ij = 1, ip1jm
          vor(ij) = vorpot(ij, l)*vorpot(ij, l)*massebxy(ij, l)
       END DO
       ztotl(l) = (ssum(ip1jm, vor, 1)-ssum(jjm, vor, iip1))

       ge = masse(:, :, l) * (phis + teta(:, :, l) * pk(:, :, l) &
            + bernf(:, :, l) - phi(:, :, l))
       etotl(l) = sum(ge) - sum(ge(1, :))

       ge = masse(:, :, l)*teta(:, :, l)
       stotl(l) = sum(ge) - sum(ge(1, :))

       ge = masse(:, :, l) * max(bernf(:, :, l) - phi(:, :, l), 0.)
       rmsvl(l) = 2.*(sum(ge)-sum(ge(1, :)))

       forall (j = 2:jjm) ge(:, j) = (ucov(:, j, l) / cu_2d(:, j) &
            + radomeg * cosphi(j)) * masse(:, j, l) * cosphi(j)
       angl(l) = radsg * (sum(ge(:, 2:jjm)) - sum(ge(1, 2:jjm)))
    END DO

    ge = ps * aire_2d
    ptot = sum(ge) - sum(ge(1, :))
    etot = sum(etotl)
    ztot = sum(ztotl)
    stot = sum(stotl)
    rmsv = sum(rmsvl)
    ang = sum(angl)

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
       etot = etot/etot0
       rmsv = sqrt(rmsv/ptot)
       ptot = ptot/ptot0
       ztot = ztot/ztot0
       stot = stot/stot0
       ang = ang/ang0
    end IF

  END SUBROUTINE sortvarc

end module sortvarc_m
