SUBROUTINE advnqy(q, qs, qn)

  ! Auteurs:   Calcul des valeurs de q aux point v.

  ! --------------------------------------------------------------------
  USE dimensions
  USE paramet_m
  USE conf_gcm_m
  IMPLICIT NONE



  ! Arguments:
  ! ----------
  REAL q(ip1jmp1, llm), qs(ip1jmp1, llm), qn(ip1jmp1, llm)

  ! Local
  ! ---------

  INTEGER ij, l

  REAL dyqv(ip1jm), zqv(ip1jm, llm)
  REAL zqmax(ip1jm), zqmin(ip1jm)
  LOGICAL extremum(ip1jmp1)

  INTEGER mode
  SAVE mode
  DATA mode/1/

  IF (mode==0) THEN
    DO l = 1, llm
      DO ij = 1, ip1jmp1
        qn(ij, l) = q(ij, l)
        qs(ij, l) = q(ij, l)
      END DO
    END DO
  ELSE

    ! calcul des pentes en u:
    ! -----------------------
    DO l = 1, llm
      DO ij = 1, ip1jm
        dyqv(ij) = q(ij, l) - q(ij+iip1, l)
      END DO

      DO ij = iip2, ip1jm - iip1
        zqv(ij, l) = 0.5*(q(ij+iip1,l)+q(ij,l))
        zqv(ij, l) = zqv(ij, l) + (dyqv(ij+iip1)-dyqv(ij-iip1))/12.
      END DO

      DO ij = iip2, ip1jm
        extremum(ij) = dyqv(ij)*dyqv(ij-iip1) <= 0.
      END DO

      ! Pas de pentes aux poles
      DO ij = 1, iip1
        zqv(ij, l) = q(ij, l)
        zqv(ip1jm-iip1+ij, l) = q(ip1jm+ij, l)
        extremum(ij) = .TRUE.
        extremum(ip1jmp1-iip1+ij) = .TRUE.
      END DO

      ! calcul des valeurs max et min acceptees aux interfaces
      DO ij = 1, ip1jm
        zqmax(ij) = max(q(ij+iip1,l), q(ij,l))
        zqmin(ij) = min(q(ij+iip1,l), q(ij,l))
      END DO

      DO ij = 1, ip1jm
        zqv(ij, l) = min(max(zqmin(ij),zqv(ij,l)), zqmax(ij))
      END DO

      DO ij = iip2, ip1jm
        IF (extremum(ij)) THEN
          qs(ij, l) = q(ij, l)
          qn(ij, l) = q(ij, l)
          ! if (.not.extremum(ij-iip1)) qs(ij-iip1,l)=q(ij,l)
          ! if (.not.extremum(ij+iip1)) qn(ij+iip1,l)=q(ij,l)
        ELSE
          qs(ij, l) = zqv(ij, l)
          qn(ij, l) = zqv(ij-iip1, l)
        END IF
      END DO

      DO ij = 1, iip1
        qs(ij, l) = q(ij, l)
        qn(ij, l) = q(ij, l)
        qs(ip1jm+ij, l) = q(ip1jm+ij, l)
        qn(ip1jm+ij, l) = q(ip1jm+ij, l)
      END DO

    END DO
  END IF
  RETURN
END SUBROUTINE advnqy

