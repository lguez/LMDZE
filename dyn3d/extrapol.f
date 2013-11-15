!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/extrapol.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
C
C
      SUBROUTINE extrapol (pfild, kxlon, kylat, pmask,
     .                   norsud, ldper, knbor, pwork)
      IMPLICIT none
c
c OASIS routine (Adaptation: Laurent Li, le 14 mars 1997)
c Fill up missed values by using the neighbor points
c
      INTEGER kxlon, kylat ! longitude and latitude dimensions (Input)
      INTEGER knbor ! minimum neighbor number (Input)
      LOGICAL norsud ! True if field is from North to South (Input)
      LOGICAL ldper ! True if take into account the periodicity (Input)
      REAL pmask ! mask value (Input)
      REAL pfild(kxlon,kylat) ! field to be extrapolated (Input/Output)
      REAL pwork(kxlon,kylat) ! working space
c
      REAL zwmsk
      INTEGER incre, idoit, i, j, k, inbor, ideb, ifin, ilon, jlat
      INTEGER ix(9), jy(9) ! index arrays for the neighbors coordinates
      REAL zmask(9)
C
C  We search over the eight closest neighbors
C
C            j+1  7  8  9
C              j  4  5  6    Current point 5 --> (i,j)
C            j-1  1  2  3
C                i-1 i i+1
c
c
      IF (norsud) THEN
         DO j = 1, kylat
         DO i = 1, kxlon
            pwork(i,j) = pfild(i,kylat-j+1)
         ENDDO
         ENDDO
         DO j = 1, kylat
         DO i = 1, kxlon
            pfild(i,j) = pwork(i,j)
         ENDDO
         ENDDO
      ENDIF
c
      incre = 0
c
      DO j = 1, kylat
      DO i = 1, kxlon
         pwork(i,j) = pfild(i,j)
      ENDDO
      ENDDO
c
C* To avoid problems in floating point tests
      zwmsk = pmask - 1.0
c
200   CONTINUE
      incre = incre + 1
      DO 99999 j = 1, kylat
      DO 99999 i = 1, kxlon
      IF (pfild(i,j).GT. zwmsk) THEN
         pwork(i,j) = pfild(i,j)
         inbor = 0
         ideb = 1
         ifin = 9
C
C* Fill up ix array
         ix(1) = MAX (1,i-1)
         ix(2) = i
         ix(3) = MIN (kxlon,i+1)
         ix(4) = MAX (1,i-1)
         ix(5) = i
         ix(6) = MIN (kxlon,i+1)
         ix(7) = MAX (1,i-1)
         ix(8) = i
         ix(9) = MIN (kxlon,i+1)
C
C* Fill up iy array
         jy(1) = MAX (1,j-1)
         jy(2) = MAX (1,j-1)
         jy(3) = MAX (1,j-1)
         jy(4) = j
         jy(5) = j
         jy(6) = j
         jy(7) = MIN (kylat,j+1)
         jy(8) = MIN (kylat,j+1)
         jy(9) = MIN (kylat,j+1)
C
C* Correct latitude bounds if southernmost or northernmost points
         IF (j .EQ. 1) ideb = 4
         IF (j .EQ. kylat) ifin = 6
C
C* Account for periodicity in longitude
C
         IF (ldper) THEN 
            IF (i .EQ. kxlon) THEN
               ix(3) = 1
               ix(6) = 1
               ix(9) = 1
            ELSE IF (i .EQ. 1) THEN
               ix(1) = kxlon
               ix(4) = kxlon
               ix(7) = kxlon
            ENDIF
         ELSE
            IF (i .EQ. 1) THEN
               ix(1) = i
               ix(2) = i + 1
               ix(3) = i
               ix(4) = i + 1
               ix(5) = i
               ix(6) = i + 1
            ENDIF 
            IF (i .EQ. kxlon) THEN
               ix(1) = i -1
               ix(2) = i
               ix(3) = i - 1
               ix(4) = i
               ix(5) = i - 1
               ix(6) = i
            ENDIF
C
            IF (i .EQ. 1 .OR. i .EQ. kxlon) THEN 
               jy(1) = MAX (1,j-1)
               jy(2) = MAX (1,j-1)
               jy(3) = j
               jy(4) = j
               jy(5) = MIN (kylat,j+1)
               jy(6) = MIN (kylat,j+1)
C
               ideb = 1
               ifin = 6
               IF (j .EQ. 1) ideb = 3
               IF (j .EQ. kylat) ifin = 4
            ENDIF
         ENDIF ! end for ldper test
C
C* Find unmasked neighbors
C
         DO 230 k = ideb, ifin
            zmask(k) = 0.
            ilon = ix(k)
            jlat = jy(k)
            IF (pfild(ilon,jlat) .LT. zwmsk) THEN
               zmask(k) = 1.
               inbor = inbor + 1
            ENDIF
 230     CONTINUE
C
C* Not enough points around point P are unmasked; interpolation on P 
C  will be done in a future call to extrap.
C
         IF (inbor .GE. knbor) THEN
            pwork(i,j) = 0.
            DO k = ideb, ifin
               ilon = ix(k)
               jlat = jy(k)
               pwork(i,j) = pwork(i,j)
     $                      + pfild(ilon,jlat) * zmask(k)/FLOAT(inbor)
            ENDDO
         ENDIF
C
      ENDIF
99999 CONTINUE
C
C*    3. Writing back unmasked field in pfild
C        ------------------------------------
C
C* pfild then contains:
C     - Values which were not masked
C     - Interpolated values from the inbor neighbors
C     - Values which are not yet interpolated
C
      idoit = 0
      DO j = 1, kylat
      DO i = 1, kxlon
         IF (pwork(i,j) .GT. zwmsk) idoit = idoit + 1
         pfild(i,j) = pwork(i,j)
      ENDDO
      ENDDO
c
      IF (idoit .ne. 0) GOTO 200
ccc      PRINT*, "Number of extrapolation steps incre =", incre
c
      IF (norsud) THEN
         DO j = 1, kylat
         DO i = 1, kxlon
            pwork(i,j) = pfild(i,kylat-j+1)
         ENDDO
         ENDDO
         DO j = 1, kylat
         DO i = 1, kxlon
            pfild(i,j) = pwork(i,j)
         ENDDO
         ENDDO
      ENDIF
c
      RETURN
      END
