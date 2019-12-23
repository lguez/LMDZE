SUBROUTINE advnqx(q, qg, qd)

  ! Auteurs:   Calcul des valeurs de q aux point u.

  ! --------------------------------------------------------------------
  USE dimensions
  USE paramet_m
  USE conf_gcm_m
  IMPLICIT NONE



  ! Arguments:
  ! ----------
  REAL q(ip1jmp1, llm), qg(ip1jmp1, llm), qd(ip1jmp1, llm)

  ! Local
  ! ---------

  INTEGER ij, l

  REAL dxqu(ip1jmp1), zqu(ip1jmp1)
  REAL zqmax(ip1jmp1), zqmin(ip1jmp1)
  LOGICAL extremum(ip1jmp1)

  INTEGER mode
  SAVE mode
  DATA mode/1/

  ! calcul des pentes en u:
  ! -----------------------
  IF (mode==0) THEN
    DO l = 1, llm
      DO ij = 1, ip1jm
        qd(ij, l) = q(ij, l)
        qg(ij, l) = q(ij, l)
      END DO
    END DO
  ELSE
    DO l = 1, llm
      DO ij = iip2, ip1jm - 1
        dxqu(ij) = q(ij+1, l) - q(ij, l)
        zqu(ij) = 0.5*(q(ij+1,l)+q(ij,l))
      END DO
      DO ij = iip1 + iip1, ip1jm, iip1
        dxqu(ij) = dxqu(ij-iim)
        zqu(ij) = zqu(ij-iim)
      END DO
      DO ij = iip2, ip1jm - 1
        zqu(ij) = zqu(ij) - dxqu(ij+1)/12.
      END DO
      DO ij = iip1 + iip1, ip1jm, iip1
        zqu(ij) = zqu(ij-iim)
      END DO
      DO ij = iip2 + 1, ip1jm
        zqu(ij) = zqu(ij) + dxqu(ij-1)/12.
      END DO
      DO ij = iip1 + iip1, ip1jm, iip1
        zqu(ij-iim) = zqu(ij)
      END DO

      ! calcul des valeurs max et min acceptees aux interfaces

      DO ij = iip2, ip1jm - 1
        zqmax(ij) = max(q(ij+1,l), q(ij,l))
        zqmin(ij) = min(q(ij+1,l), q(ij,l))
      END DO
      DO ij = iip1 + iip1, ip1jm, iip1
        zqmax(ij) = zqmax(ij-iim)
        zqmin(ij) = zqmin(ij-iim)
      END DO
      DO ij = iip2 + 1, ip1jm
        extremum(ij) = dxqu(ij)*dxqu(ij-1) <= 0.
      END DO
      DO ij = iip1 + iip1, ip1jm, iip1
        extremum(ij-iim) = extremum(ij)
      END DO
      DO ij = iip2, ip1jm
        zqu(ij) = min(max(zqmin(ij),zqu(ij)), zqmax(ij))
      END DO
      DO ij = iip2 + 1, ip1jm
        IF (extremum(ij)) THEN
          qg(ij, l) = q(ij, l)
          qd(ij, l) = q(ij, l)
        ELSE
          qd(ij, l) = zqu(ij)
          qg(ij, l) = zqu(ij-1)
        END IF
      END DO
      DO ij = iip1 + iip1, ip1jm, iip1
        qd(ij-iim, l) = qd(ij, l)
        qg(ij-iim, l) = qg(ij, l)
      END DO

      GO TO 8888

      DO ij = iip2 + 1, ip1jm
        IF (extremum(ij) .AND. .NOT. extremum(ij-1)) qd(ij-1, l) = q(ij, l)
      END DO

      DO ij = iip1 + iip1, ip1jm, iip1
        qd(ij-iim, l) = qd(ij, l)
      END DO
      DO ij = iip2, ip1jm - 1
        IF (extremum(ij) .AND. .NOT. extremum(ij+1)) qg(ij+1, l) = q(ij, l)
      END DO

      DO ij = iip1 + iip1, ip1jm, iip1
        qg(ij, l) = qg(ij-iim, l)
      END DO
8888  CONTINUE
    END DO
  END IF
  RETURN
END SUBROUTINE advnqx
