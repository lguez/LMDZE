module readsulfate_m

  IMPLICIT none

contains

  SUBROUTINE readsulfate(dayvrai, time, first, sulfate)

    ! From LMDZ4/libf/phylmd/readsulfate.F, version 1.2 2005/05/19
    ! 08:27:15 fairhead

    ! This routine reads monthly mean values of sulfate aerosols and
    ! returns a linearly interpolated daily-mean field.

    ! Author: Johannes Quaas (quaas@lmd.jussieu.fr)
    ! April 26th, 2001

    ! ATTENTION!! runs are supposed to start with Jan, 1. 1930
    ! (rday = 1)

    ! The model always has 360 days per year.
    ! SO4 concentration rather than mixing ratio
    ! 10yr-mean-values to interpolate
    ! Introduce flag to read in just one decade

    use dimens_m, only: iim, jjm
    use dimphy, only: klon, klev
    use dynetat0_m, only: annee_ref
    use getso4fromfile_m, only: getso4fromfile

    integer, intent(in):: dayvrai
    ! current day number, based at value 1 on January 1st of annee_ref

    REAL, intent(in):: time ! heure de la journ\'ee en fraction de jour

    LOGICAL, intent(in):: first ! First timestep
    ! (and therefore initialization necessary)

    real, intent(out):: sulfate(klon, klev)
    ! mass of sulfate (monthly mean data, from file) (micro g SO4 / m3)

    ! Local:
    INTEGER i, ig, k, it
    INTEGER j, iday, iyr, iyr1, iyr2
    CHARACTER(len = 4) cyear
    INTEGER im, day1, day2, im2
    real so4_1(iim, jjm + 1, klev, 12)
    real so4_2(iim, jjm + 1, klev, 12) ! sulfate distributions
    double precision, save:: so4(klon, klev, 12) ! SO4 in right dimension
    double precision, save:: so4_out(klon, klev)
    LOGICAL lnewday

    !----------------------------------------------------------------

    iday = dayvrai

    ! Get the year of the run
    iyr = iday/360

    ! Get the day of the actual year:
    iday = iday - iyr*360

    ! 0.02 is about 0.5/24, namly less than half an hour
    lnewday = time < 0.02

    ! All has to be done only, if a new day begins

    test_newday: IF (lnewday .OR. first) THEN
       im = iday/30 + 1 ! the actual month

       ! annee_ref is the initial year (defined in temps.h)
       iyr = iyr + annee_ref

       ! Do I have to read new data? (Is this the first day of a year?)
       IF (first .OR. iday == 1.) THEN
          ! Initialize field
          DO it = 1, 12
             DO k = 1, klev
                DO i = 1, klon
                   so4(i, k, it) = 0.
                ENDDO
             ENDDO
          ENDDO

          IF (iyr < 1850) THEN
             cyear = '.nat'
             print *, 'getso4 iyr = ', iyr, ' ', cyear
             CALL getso4fromfile(cyear, so4_1)
          ELSE IF (iyr >= 2100) THEN
             cyear = '2100'
             print *, 'getso4 iyr = ', iyr, ' ', cyear
             CALL getso4fromfile(cyear, so4_1)
          ELSE
             ! Read in data:
             ! a) from actual 10-yr-period

             IF (iyr < 1900) THEN
                iyr1 = 1850
                iyr2 = 1900
             ELSE IF (iyr >= 1900.and.iyr < 1920) THEN
                iyr1 = 1900
                iyr2 = 1920
             ELSE
                iyr1 = INT(iyr/10)*10
                iyr2 = INT(1 + iyr/10)*10
             ENDIF
             WRITE(cyear, '(I4)') iyr1
             print *, 'getso4 iyr = ', iyr, ' ', cyear
             CALL getso4fromfile(cyear, so4_1)

             ! Read two decades:
             ! b) from the next following one
             WRITE(cyear, '(I4)') iyr2
             print *, 'getso4 iyr = ', iyr, ' ', cyear
             CALL getso4fromfile(cyear, so4_2)

             ! Interpolate linarily to the actual year:
             DO it = 1, 12
                DO k = 1, klev
                   DO j = 1, jjm
                      DO i = 1, iim
                         so4_1(i, j, k, it) = so4_1(i, j, k, it) &
                              - REAL(iyr - iyr1)/REAL(iyr2 - iyr1) &
                              * (so4_1(i, j, k, it) - so4_2(i, j, k, it))
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          ! Transform the horizontal 2D-field into the physics-field
          ! (Also the levels and the latitudes have to be inversed)

          DO it = 1, 12
             DO k = 1, klev
                ! a) at the poles, use the zonal mean:
                DO i = 1, iim
                   ! North pole
                   so4(1, k, it) = so4(1, k, it) &
                        + so4_1(i, jjm + 1, klev + 1 - k, it)
                   ! South pole
                   so4(klon, k, it) = so4(klon, k, it) &
                        + so4_1(i, 1, klev + 1 - k, it)
                ENDDO
                so4(1, k, it) = so4(1, k, it)/REAL(iim)
                so4(klon, k, it) = so4(klon, k, it)/REAL(iim)

                ! b) the values between the poles:
                ig = 1
                DO j = 2, jjm
                   DO i = 1, iim
                      ig = ig + 1
                      if (ig > klon) stop 1
                      so4(ig, k, it) = so4_1(i, jjm + 1 - j, klev + 1 - k, it)
                   ENDDO
                ENDDO
                IF (ig /= klon - 1) then
                   print *, 'Error in readsulfate (var conversion)'
                   STOP 1
                end IF
             ENDDO ! Loop over k (vertical)
          ENDDO ! Loop over it (months)
       ENDIF ! Had to read new data?

       ! Interpolate to actual day:
       IF (iday < im*30 - 15) THEN
          ! in the first half of the month use month before and actual month
          im2 = im - 1
          day1 = im2*30 + 15
          day2 = im2*30 - 15
          IF (im2 <= 0) THEN
             ! the month is january, thus the month before december
             im2 = 12
          ENDIF
          DO k = 1, klev
             DO i = 1, klon
                sulfate(i, k) = so4(i, k, im2) &
                     - REAL(iday - day2)/REAL(day1 - day2) &
                     * (so4(i, k, im2) - so4(i, k, im))
                IF (sulfate(i, k) < 0.) THEN
                   IF (iday - day2 < 0.) write(*, *) 'iday - day2', iday - day2
                   IF (so4(i, k, im2) - so4(i, k, im) < 0.) &
                        write(*, *) 'so4(i, k, im2) - so4(i, k, im)', &
                        so4(i, k, im2) - so4(i, k, im)
                   IF (day1 - day2 < 0.) write(*, *) 'day1 - day2', day1 - day2
                   stop 1
                endif
             ENDDO
          ENDDO
       ELSE
          ! the second half of the month
          im2 = im + 1
          IF (im2 > 12) THEN
             ! the month is december, the following thus january
             im2 = 1
          ENDIF
          day2 = im*30 - 15
          day1 = im*30 + 15
          DO k = 1, klev
             DO i = 1, klon
                sulfate(i, k) = so4(i, k, im2) &
                     - REAL(iday - day2)/REAL(day1 - day2) &
                     * (so4(i, k, im2) - so4(i, k, im))
                IF (sulfate(i, k) < 0.) THEN
                   IF (iday - day2 < 0.) write(*, *) 'iday - day2', iday - day2
                   IF (so4(i, k, im2) - so4(i, k, im) < 0.) &
                        write(*, *) 'so4(i, k, im2) - so4(i, k, im)', &
                        so4(i, k, im2) - so4(i, k, im)
                   IF (day1 - day2 < 0.) write(*, *) 'day1 - day2', day1 - day2
                   stop 1
                endif
             ENDDO
          ENDDO
       ENDIF

       DO k = 1, klev
          DO i = 1, klon
             so4_out(i, k) = sulfate(i, k)
          ENDDO
       ENDDO
    ELSE
       ! If no new day, use old data:
       DO k = 1, klev
          DO i = 1, klon
             sulfate(i, k) = so4_out(i, k)
          ENDDO
       ENDDO
    ENDIF test_newday

  END SUBROUTINE readsulfate

end module readsulfate_m
