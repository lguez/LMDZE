
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/extrapol.F,v 1.1.1.1 2004/05/19
! 12:53:07 lmdzadmin Exp $



SUBROUTINE extrapol(pfild, kxlon, kylat, pmask, norsud, ldper, knbor, pwork)
  IMPLICIT NONE

  ! OASIS routine (Adaptation: Laurent Li, le 14 mars 1997)
  ! Fill up missed values by using the neighbor points

  INTEGER kxlon, kylat ! longitude and latitude dimensions (Input)
  INTEGER knbor ! minimum neighbor number (Input)
  LOGICAL norsud ! True if field is from North to South (Input)
  LOGICAL ldper ! True if take into account the periodicity (Input)
  REAL pmask ! mask value (Input)
  REAL pfild(kxlon, kylat) ! field to be extrapolated (Input/Output)
  REAL pwork(kxlon, kylat) ! working space

  REAL zwmsk
  INTEGER incre, idoit, i, j, k, inbor, ideb, ifin, ilon, jlat
  INTEGER ix(9), jy(9) ! index arrays for the neighbors coordinates
  REAL zmask(9)

  ! We search over the eight closest neighbors

  ! j+1  7  8  9
  ! j  4  5  6    Current point 5 --> (i,j)
  ! j-1  1  2  3
  ! i-1 i i+1


  IF (norsud) THEN
    DO j = 1, kylat
      DO i = 1, kxlon
        pwork(i, j) = pfild(i, kylat-j+1)
      END DO
    END DO
    DO j = 1, kylat
      DO i = 1, kxlon
        pfild(i, j) = pwork(i, j)
      END DO
    END DO
  END IF

  incre = 0

  DO j = 1, kylat
    DO i = 1, kxlon
      pwork(i, j) = pfild(i, j)
    END DO
  END DO

  ! * To avoid problems in floating point tests
  zwmsk = pmask - 1.0

200 CONTINUE
  incre = incre + 1
  DO j = 1, kylat
    DO i = 1, kxlon
      IF (pfild(i,j)>zwmsk) THEN
        pwork(i, j) = pfild(i, j)
        inbor = 0
        ideb = 1
        ifin = 9

        ! * Fill up ix array
        ix(1) = max(1, i-1)
        ix(2) = i
        ix(3) = min(kxlon, i+1)
        ix(4) = max(1, i-1)
        ix(5) = i
        ix(6) = min(kxlon, i+1)
        ix(7) = max(1, i-1)
        ix(8) = i
        ix(9) = min(kxlon, i+1)

        ! * Fill up iy array
        jy(1) = max(1, j-1)
        jy(2) = max(1, j-1)
        jy(3) = max(1, j-1)
        jy(4) = j
        jy(5) = j
        jy(6) = j
        jy(7) = min(kylat, j+1)
        jy(8) = min(kylat, j+1)
        jy(9) = min(kylat, j+1)

        ! * Correct latitude bounds if southernmost or northernmost points
        IF (j==1) ideb = 4
        IF (j==kylat) ifin = 6

        ! * Account for periodicity in longitude

        IF (ldper) THEN
          IF (i==kxlon) THEN
            ix(3) = 1
            ix(6) = 1
            ix(9) = 1
          ELSE IF (i==1) THEN
            ix(1) = kxlon
            ix(4) = kxlon
            ix(7) = kxlon
          END IF
        ELSE
          IF (i==1) THEN
            ix(1) = i
            ix(2) = i + 1
            ix(3) = i
            ix(4) = i + 1
            ix(5) = i
            ix(6) = i + 1
          END IF
          IF (i==kxlon) THEN
            ix(1) = i - 1
            ix(2) = i
            ix(3) = i - 1
            ix(4) = i
            ix(5) = i - 1
            ix(6) = i
          END IF

          IF (i==1 .OR. i==kxlon) THEN
            jy(1) = max(1, j-1)
            jy(2) = max(1, j-1)
            jy(3) = j
            jy(4) = j
            jy(5) = min(kylat, j+1)
            jy(6) = min(kylat, j+1)

            ideb = 1
            ifin = 6
            IF (j==1) ideb = 3
            IF (j==kylat) ifin = 4
          END IF
        END IF ! end for ldper test

        ! * Find unmasked neighbors

        DO k = ideb, ifin
          zmask(k) = 0.
          ilon = ix(k)
          jlat = jy(k)
          IF (pfild(ilon,jlat)<zwmsk) THEN
            zmask(k) = 1.
            inbor = inbor + 1
          END IF
        END DO

        ! * Not enough points around point P are unmasked; interpolation on P
        ! will be done in a future call to extrap.

        IF (inbor>=knbor) THEN
          pwork(i, j) = 0.
          DO k = ideb, ifin
            ilon = ix(k)
            jlat = jy(k)
            pwork(i, j) = pwork(i, j) + pfild(ilon, jlat)*zmask(k)/float( &
              inbor)
          END DO
        END IF

      END IF
    END DO
  END DO

  ! *    3. Writing back unmasked field in pfild
  ! ------------------------------------

  ! * pfild then contains:
  ! - Values which were not masked
  ! - Interpolated values from the inbor neighbors
  ! - Values which are not yet interpolated

  idoit = 0
  DO j = 1, kylat
    DO i = 1, kxlon
      IF (pwork(i,j)>zwmsk) idoit = idoit + 1
      pfild(i, j) = pwork(i, j)
    END DO
  END DO

  IF (idoit/=0) GO TO 200
  ! cc      PRINT*, "Number of extrapolation steps incre =", incre

  IF (norsud) THEN
    DO j = 1, kylat
      DO i = 1, kxlon
        pwork(i, j) = pfild(i, kylat-j+1)
      END DO
    END DO
    DO j = 1, kylat
      DO i = 1, kxlon
        pfild(i, j) = pwork(i, j)
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE extrapol
