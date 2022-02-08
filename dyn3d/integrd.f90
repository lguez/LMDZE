module integrd_m

  IMPLICIT NONE

contains

  SUBROUTINE integrd(vcovm1, ucovm1, tetam1, psm1, massem1, dv, du, dteta, &
       dp, vcov, ucov, teta, q, ps, masse, dt, leapf)

    ! From dyn3d/integrd.F, version 1.1.1.1, 2004/05/19 12:53:05
    ! Author: P. Le Van 
    ! Objet : incrémentation des tendances dynamiques

    ! Libraries:
    use jumble, only: assert

    USE comgeom, ONLY : aire_2d, apoln, apols
    USE dimensions, ONLY : iim, jjm, llm
    USE disvert_m, ONLY : ap, bp
    use massdair_m, only: massdair
    USE paramet_m, ONLY : ip1jm
    use qminimum_m, only: qminimum

    REAL, intent(inout):: vcovm1(ip1jm, llm), ucovm1((iim + 1) * (jjm + 1), llm)
    REAL, intent(inout):: tetam1(iim + 1, jjm + 1, llm)
    REAL, intent(inout):: psm1(:, :) ! (iim + 1, jjm + 1)
    real, intent(inout):: massem1(iim + 1, jjm + 1, llm)
    REAL, intent(in):: dv(ip1jm, llm), du((iim + 1) * (jjm + 1), llm)
    REAL, intent(in):: dteta(iim + 1, jjm + 1, llm)
    REAL, intent(in):: dp(:, :) ! (iim + 1, jjm + 1)
    REAL, intent(inout):: vcov(ip1jm, llm), ucov((iim + 1) * (jjm + 1), llm)
    real, intent(inout):: teta(iim + 1, jjm + 1, llm)
    REAL, intent(inout):: q(:, :, :, :) ! (iim + 1, jjm + 1, llm, nq)
    REAL, intent(inout):: ps(:, :) ! (iim + 1, jjm + 1) pression au sol, en Pa
    REAL, intent(inout):: masse(iim + 1, jjm + 1, llm)
    real, intent(in):: dt ! time step, in s
    LOGICAL, INTENT (IN) :: leapf

    ! Local: 
    REAL finvmaold(iim + 1, jjm + 1, llm)
    INTEGER nq
    REAL vscr(ip1jm), uscr((iim + 1) * (jjm + 1)), hscr(iim + 1, jjm + 1)
    real pscr(iim + 1, jjm + 1)
    REAL p(iim + 1, jjm + 1, llm + 1)
    REAL deltap(iim + 1, jjm + 1, llm)
    INTEGER l, ij, iq, i, j

    !-----------------------------------------------------------------------

    call assert(size(q, 1) == iim + 1, size(q, 2) == jjm + 1, &
         size(q, 3) == llm, "integrd")
    nq = size(q, 4)

    DO l = 1, llm
       DO ij = 1, iim + 1
          ucov(ij, l) = 0.
          ucov(ij+ip1jm, l) = 0.
          uscr(ij) = 0.
          uscr(ij+ip1jm) = 0.
       END DO
    END DO

    ! Integration de ps :

    pscr = ps
    ps = psm1 + dt * dp

    DO j = 1, jjm + 1
       do i = 1, iim + 1
          IF (ps(i, j) < 0.) THEN
             PRINT *, 'integrd: au point i, j = ', i, j, &
                  ', negative surface pressure ', ps(i, j)
             STOP 1
          END IF
       END DO
    end DO

    ps(:, 1) = sum(aire_2d(:iim, 1) * ps(:iim, 1)) / apoln
    ps(:, jjm + 1) = sum(aire_2d(:iim, jjm + 1) * ps(:iim, jjm + 1)) / apols
    

    ! Calcul de la nouvelle masse d'air au dernier temps integre t+1

    forall (l = 1: llm + 1) p(:, :, l) = ap(l) + bp(l) * ps
    CALL massdair(p, finvmaold)

    ! integration de ucov, vcov, h

    DO l = 1, llm
       DO ij = iim + 2, ip1jm
          uscr(ij) = ucov(ij, l)
          ucov(ij, l) = ucovm1(ij, l) + dt * du(ij, l)
       END DO

       DO ij = 1, ip1jm
          vscr(ij) = vcov(ij, l)
          vcov(ij, l) = vcovm1(ij, l) + dt * dv(ij, l)
       END DO

       hscr = teta(:, :, l)
       teta(:, :, l) = tetam1(:, :, l) * massem1(:, :, l) / finvmaold(:, :, l) &
            + dt * dteta(:, :, l) / finvmaold(:, :, l)

       ! Calcul de la valeur moyenne, unique aux poles pour teta
       teta(:, 1, l) = sum(aire_2d(:iim, 1) * teta(:iim, 1, l)) / apoln
       teta(:, jjm + 1, l) = sum(aire_2d(:iim, jjm + 1) &
            * teta(:iim, jjm + 1, l)) / apols

       IF (leapf) THEN
          ucovm1(:, l) = uscr
          vcovm1(:, l) = vscr
          tetam1(:, :, l) = hscr
       END IF
    END DO

    forall (l = 1:llm) deltap(:, :, l) = p(:, :, l) - p(:, :, l + 1)
    CALL qminimum(q, nq, deltap)

    ! Calcul de la valeur moyenne, unique aux poles pour q
    DO iq = 1, nq
       DO l = 1, llm
          q(:, 1, l, iq) = sum(aire_2d(:iim, 1) * q(:iim, 1, l, iq)) / apoln
          q(:, jjm + 1, l, iq) = sum(aire_2d(:iim, jjm + 1) &
               * q(:iim, jjm + 1, l, iq)) / apols
       END DO
    END DO

    IF (leapf) THEN
       psm1 = pscr
       massem1 = masse
    END IF

    masse = finvmaold

  END SUBROUTINE integrd

end module integrd_m
