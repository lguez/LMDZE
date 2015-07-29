SUBROUTINE xtp(imr, jnp, iml, j1, j2, jn, js, pu, dq, q, uc, fx1, xmass, &
    iord)
  DIMENSION uc(imr, *), dc(-iml:imr+iml+1), xmass(imr, jnp), fx1(imr+1), &
    dq(imr, jnp), qtmp(-iml:imr+1+iml)
  DIMENSION pu(imr, jnp), q(imr, jnp), isave(imr)

  imp = imr + 1

  ! van Leer at high latitudes
  jvan = max(1, jnp/18)
  j1vl = j1 + jvan
  j2vl = j2 - jvan

  DO j = j1, j2

    DO i = 1, imr
      qtmp(i) = q(i, j)
    END DO

    IF (j>=jn .OR. j<=js) GO TO 2222
    ! ************* Eulerian **********

    qtmp(0) = q(imr, j)
    qtmp(-1) = q(imr-1, j)
    qtmp(imp) = q(1, j)
    qtmp(imp+1) = q(2, j)

    IF (iord==1 .OR. j==j1 .OR. j==j2) THEN
      DO i = 1, imr
        iu = float(i) - uc(i, j)
        fx1(i) = qtmp(iu)
      END DO
    ELSE
      CALL xmist(imr, iml, qtmp, dc)
      dc(0) = dc(imr)

      IF (iord==2 .OR. j<=j1vl .OR. j>=j2vl) THEN
        DO i = 1, imr
          iu = float(i) - uc(i, j)
          fx1(i) = qtmp(iu) + dc(iu)*(sign(1.,uc(i,j))-uc(i,j))
        END DO
      ELSE
        CALL fxppm(imr, iml, uc(1,j), qtmp, dc, fx1, iord)
      END IF

    END IF

    DO i = 1, imr
      fx1(i) = fx1(i)*xmass(i, j)
    END DO

    GO TO 1309

    ! ***** Conservative (flux-form) Semi-Lagrangian transport *****

2222 CONTINUE

    DO i = -iml, 0
      qtmp(i) = q(imr+i, j)
      qtmp(imp-i) = q(1-i, j)
    END DO

    IF (iord==1 .OR. j==j1 .OR. j==j2) THEN
      DO i = 1, imr
        itmp = int(uc(i,j))
        isave(i) = i - itmp
        iu = i - uc(i, j)
        fx1(i) = (uc(i,j)-itmp)*qtmp(iu)
      END DO
    ELSE
      CALL xmist(imr, iml, qtmp, dc)

      DO i = -iml, 0
        dc(i) = dc(imr+i)
        dc(imp-i) = dc(1-i)
      END DO

      DO i = 1, imr
        itmp = int(uc(i,j))
        rut = uc(i, j) - itmp
        isave(i) = i - itmp
        iu = i - uc(i, j)
        fx1(i) = rut*(qtmp(iu)+dc(iu)*(sign(1.,rut)-rut))
      END DO
    END IF

    DO i = 1, imr
      IF (uc(i,j)>1.) THEN
        ! DIR$ NOVECTOR
        DO ist = isave(i), i - 1
          fx1(i) = fx1(i) + qtmp(ist)
        END DO
      ELSE IF (uc(i,j)<-1.) THEN
        DO ist = i, isave(i) - 1
          fx1(i) = fx1(i) - qtmp(ist)
        END DO
        ! DIR$ VECTOR
      END IF
    END DO
    DO i = 1, imr
      fx1(i) = pu(i, j)*fx1(i)
    END DO

    ! ***************************************

1309 fx1(imp) = fx1(1)
    DO i = 1, imr
      dq(i, j) = dq(i, j) + fx1(i) - fx1(i+1)
    END DO

    ! ***************************************

  END DO
  RETURN
END SUBROUTINE xtp

