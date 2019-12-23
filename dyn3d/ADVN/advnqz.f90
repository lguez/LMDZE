SUBROUTINE advnqz(q, qh, qb)

  ! Auteurs:   Calcul des valeurs de q aux point v.

  ! --------------------------------------------------------------------
  USE dimensions
  USE paramet_m
  USE conf_gcm_m
  IMPLICIT NONE



  ! Arguments:
  ! ----------
  REAL q(ip1jmp1, llm), qh(ip1jmp1, llm), qb(ip1jmp1, llm)

  ! Local
  ! ---------

  INTEGER ij, l

  REAL dzqw(ip1jmp1, llm+1), zqw(ip1jmp1, llm+1)
  REAL zqmax(ip1jmp1, llm), zqmin(ip1jmp1, llm)
  LOGICAL extremum(ip1jmp1, llm)

  INTEGER mode
  SAVE mode

  DATA mode/1/

  ! calcul des pentes en u:
  ! -----------------------

  IF (mode==0) THEN
    DO l = 1, llm
      DO ij = 1, ip1jmp1
        qb(ij, l) = q(ij, l)
        qh(ij, l) = q(ij, l)
      END DO
    END DO
  ELSE
    DO l = 2, llm
      DO ij = 1, ip1jmp1
        dzqw(ij, l) = q(ij, l-1) - q(ij, l)
        zqw(ij, l) = 0.5*(q(ij,l-1)+q(ij,l))
      END DO
    END DO
    DO ij = 1, ip1jmp1
      dzqw(ij, 1) = 0.
      dzqw(ij, llm+1) = 0.
    END DO
    DO l = 2, llm
      DO ij = 1, ip1jmp1
        zqw(ij, l) = zqw(ij, l) + (dzqw(ij,l+1)-dzqw(ij,l-1))/12.
      END DO
    END DO
    DO l = 2, llm - 1
      DO ij = 1, ip1jmp1
        extremum(ij, l) = dzqw(ij, l)*dzqw(ij, l+1) <= 0.
      END DO
    END DO

    ! Pas de pentes en bas et en haut
    DO ij = 1, ip1jmp1
      zqw(ij, 2) = q(ij, 1)
      zqw(ij, llm) = q(ij, llm)
      extremum(ij, 1) = .TRUE.
      extremum(ij, llm) = .TRUE.
    END DO

    ! calcul des valeurs max et min acceptees aux interfaces
    DO l = 2, llm
      DO ij = 1, ip1jmp1
        zqmax(ij, l) = max(q(ij,l-1), q(ij,l))
        zqmin(ij, l) = min(q(ij,l-1), q(ij,l))
      END DO
    END DO

    DO l = 2, llm
      DO ij = 1, ip1jmp1
        zqw(ij, l) = min(max(zqmin(ij,l),zqw(ij,l)), zqmax(ij,l))
      END DO
    END DO

    DO l = 2, llm - 1
      DO ij = 1, ip1jmp1
        IF (extremum(ij,l)) THEN
          qh(ij, l) = q(ij, l)
          qb(ij, l) = q(ij, l)
        ELSE
          qh(ij, l) = zqw(ij, l+1)
          qb(ij, l) = zqw(ij, l)
        END IF
      END DO
    END DO
    ! do l=2,llm-1
    ! do ij=1,ip1jmp1
    ! if(extremum(ij,l)) then
    ! if (.not.extremum(ij,l-1)) qh(ij,l-1)=q(ij,l)
    ! if (.not.extremum(ij,l+1)) qb(ij,l+1)=q(ij,l)
    ! endif
    ! enddo
    ! enddo

    DO ij = 1, ip1jmp1
      qb(ij, 1) = q(ij, 1)
      qh(ij, 1) = q(ij, 1)
      qb(ij, llm) = q(ij, llm)
      qh(ij, llm) = q(ij, llm)
    END DO

  END IF

  RETURN
END SUBROUTINE advnqz

