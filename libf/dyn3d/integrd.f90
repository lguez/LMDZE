module integrd_m

  IMPLICIT NONE

contains

  SUBROUTINE integrd(vcovm1, ucovm1, tetam1, psm1, massem1, dv, dudyn, &
       dteta, dp, vcov, ucov, teta, q, ps, masse, finvmaold, dt, leapf)

    ! From dyn3d/integrd.F, version 1.1.1.1 2004/05/19 12:53:05
    ! Author: P. Le Van 
    ! Objet: incrémentation des tendances dynamiques

    USE comgeom, ONLY : aire, apoln, apols
    USE dimens_m, ONLY : iim, jjm, llm
    USE disvert_m, ONLY : ap, bp
    USE filtreg_m, ONLY : filtreg
    use massdair_m, only: massdair
    use nr_util, only: assert
    USE paramet_m, ONLY : iip1, iip2, ip1jm, ip1jmp1, jjp1, llmp1
    use qminimum_m, only: qminimum

    ! Arguments: 

    REAL vcov(ip1jm, llm), ucov((iim + 1) * (jjm + 1), llm)
    real, intent(inout):: teta((iim + 1) * (jjm + 1), llm)
    REAL q(:, :, :, :) ! (iim + 1, jjm + 1, llm, nq)
    REAL, intent(inout):: ps((iim + 1) * (jjm + 1))
    REAL masse((iim + 1) * (jjm + 1), llm)

    REAL vcovm1(ip1jm, llm), ucovm1((iim + 1) * (jjm + 1), llm)
    REAL, intent(inout):: tetam1((iim + 1) * (jjm + 1), llm)
    REAL, intent(inout):: psm1((iim + 1) * (jjm + 1))
    real massem1((iim + 1) * (jjm + 1), llm)

    REAL dv(ip1jm, llm), dudyn((iim + 1) * (jjm + 1), llm)
    REAL dteta((iim + 1) * (jjm + 1), llm), dp((iim + 1) * (jjm + 1))
    REAL finvmaold((iim + 1) * (jjm + 1), llm)
    LOGICAL, INTENT (IN) :: leapf
    real, intent(in):: dt

    ! Local variables: 

    INTEGER nq
    REAL vscr(ip1jm), uscr((iim + 1) * (jjm + 1)), hscr((iim + 1) * (jjm + 1))
    real pscr((iim + 1) * (jjm + 1))
    REAL massescr((iim + 1) * (jjm + 1), llm)
    real finvmasse((iim + 1) * (jjm + 1), llm)
    REAL p((iim + 1) * (jjm + 1), llmp1)
    REAL tpn, tps, tppn(iim), tpps(iim)
    REAL qpn, qps, qppn(iim), qpps(iim)
    REAL deltap((iim + 1) * (jjm + 1), llm)

    INTEGER l, ij, iq

    REAL ssum

    !-----------------------------------------------------------------------

    call assert(size(q, 1) == iim + 1, size(q, 2) == jjm + 1, &
         size(q, 3) == llm, "integrd")
    nq = size(q, 4)

    DO l = 1, llm
       DO ij = 1, iip1
          ucov(ij, l) = 0.
          ucov(ij+ip1jm, l) = 0.
          uscr(ij) = 0.
          uscr(ij+ip1jm) = 0.
       END DO
    END DO

    massescr = masse

    ! Integration de ps :

    pscr = ps
    ps = psm1 + dt * dp

    DO ij = 1, (iim + 1) * (jjm + 1)
       IF (ps(ij) < 0.) THEN
          PRINT *, 'integrd: au point ij = ', ij, &
               ', negative surface pressure ', ps(ij)
          STOP 1
       END IF
    END DO

    DO ij = 1, iim
       tppn(ij) = aire(ij)*ps(ij)
       tpps(ij) = aire(ij+ip1jm) * ps(ij+ip1jm)
    END DO
    tpn = ssum(iim, tppn, 1)/apoln
    tps = ssum(iim, tpps, 1)/apols
    DO ij = 1, iip1
       ps(ij) = tpn
       ps(ij+ip1jm) = tps
    END DO

    ! Calcul de la nouvelle masse d'air au dernier temps integre t+1

    forall (l = 1: llm + 1) p(:, l) = ap(l) + bp(l) * ps
    CALL massdair(p, masse)

    finvmasse = masse
    CALL filtreg(finvmasse, jjp1, llm, -2, 2, .TRUE.)

    ! integration de ucov, vcov, h

    DO l = 1, llm
       DO ij = iip2, ip1jm
          uscr(ij) = ucov(ij, l)
          ucov(ij, l) = ucovm1(ij, l) + dt*dudyn(ij, l)
       END DO

       DO ij = 1, ip1jm
          vscr(ij) = vcov(ij, l)
          vcov(ij, l) = vcovm1(ij, l) + dt*dv(ij, l)
       END DO

       hscr = teta(:, l)
       teta(:, l) = tetam1(:, l) * massem1(:, l) / masse(:, l) &
            + dt * dteta(:, l) / masse(:, l)

       ! Calcul de la valeur moyenne, unique aux poles pour teta

       DO ij = 1, iim
          tppn(ij) = aire(ij)*teta(ij, l)
          tpps(ij) = aire(ij+ip1jm)*teta(ij+ip1jm, l)
       END DO
       tpn = ssum(iim, tppn, 1)/apoln
       tps = ssum(iim, tpps, 1)/apols

       DO ij = 1, iip1
          teta(ij, l) = tpn
          teta(ij+ip1jm, l) = tps
       END DO

       IF (leapf) THEN
          ucovm1(:, l)  =uscr
          vcovm1(:, l) = vscr
          tetam1(:, l) = hscr
       END IF
    END DO

    DO l = 1, llm
       DO ij = 1, (iim + 1) * (jjm + 1)
          deltap(ij, l) = p(ij, l) - p(ij, l+1)
       END DO
    END DO

    CALL qminimum(q, nq, deltap)

    ! Calcul de la valeur moyenne, unique aux poles pour q

    DO iq = 1, nq
       DO l = 1, llm
          DO ij = 1, iim
             qppn(ij) = aire(ij)*q(ij, 1, l, iq)
             qpps(ij) = aire(ij+ip1jm)*q(ij, jjm + 1, l, iq)
          END DO
          qpn = ssum(iim, qppn, 1)/apoln
          qps = ssum(iim, qpps, 1)/apols

          DO ij = 1, iip1
             q(ij, 1, l, iq) = qpn
             q(ij, jjm + 1, l, iq) = qps
          END DO
       END DO
    END DO

    finvmaold = finvmasse

    ! Fin de l'integration de q

    IF (leapf) THEN
       psm1 = pscr
       massem1 = massescr
    END IF

  END SUBROUTINE integrd

end module integrd_m
