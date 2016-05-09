module readsulfate_preind_m

  IMPLICIT none

contains

  SUBROUTINE readsulfate_preind(dayvrai, time, first, sulfate)

    ! Read and calculate pre-industrial values of sulfate. This
    ! routine reads monthly mean values of sulfate aerosols and
    ! returns a linearly interpolated daily-mean field. It does so for
    ! the preindustriel values of the sulfate, to a large part
    ! analogous to the routine readsulfate.

    ! Author:
    ! Johannes Quaas (quaas@lmd.jussieu.fr)
    ! June 26th, 2001

    use dimens_m, only: iim, jjm
    use dimphy, only: klon, klev
    use dynetat0_m, only: annee_ref
    use getso4fromfile_m, only: getso4fromfile

    integer, intent(in):: dayvrai
    ! current day number, based at value 1 on January 1st of annee_ref

    REAL, intent(in):: time ! heure de la journ\'ee en fraction de jour

    LOGICAL, intent(in):: first ! First timestep
    ! (and therefore initialization necessary)

    real, intent(out):: sulfate (klon, klev)
    ! number concentration sulfate (monthly mean data, from file)

    ! Local:
    INTEGER i, ig, k, it
    INTEGER j, iday, iyr
    INTEGER im, day1, day2, im2
    real so4_1(iim, jjm + 1, klev, 12)
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

       IF (first) THEN
          ! Initialize field
          DO it = 1, 12
             DO k = 1, klev
                DO i = 1, klon
                   so4(i, k, it) = 0.
                ENDDO
             ENDDO
          ENDDO

          CALL getso4fromfile('.nat', so4_1)

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

  END SUBROUTINE readsulfate_preind

end module readsulfate_preind_m
