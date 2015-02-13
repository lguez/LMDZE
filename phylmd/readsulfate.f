module readsulfate_m

  IMPLICIT none

contains

  SUBROUTINE readsulfate(r_day, first, sulfate)

    ! From LMDZ4/libf/phylmd/readsulfate.F, version 1.2 2005/05/19
    ! 08:27:15 fairhead

    ! This routine reads in monthly mean values of sulfate aerosols and
    ! returns a linearly interpolated daily-mean field.

    ! Author: Johannes Quaas (quaas@lmd.jussieu.fr)
    ! 26/04/01

    ! Modifications:
    ! 21/06/01: Make integrations of more than one year possible ;-)
    ! ATTENTION!! runs are supposed to start with Jan, 1. 1930
    ! (rday=1)

    ! 27/06/01: Correction: The model always has 360 days per year!
    ! 27/06/01: SO4 concentration rather than mixing ratio
    ! 27/06/01: 10yr-mean-values to interpolate
    ! 20/08/01: Correct the error through integer-values in interpolations
    ! 21/08/01: Introduce flag to read in just one decade

    USE dimens_m, ONLY: iim, jjm
    USE dimphy, ONLY: klev, klon
    USE dynetat0_m, ONLY: annee_ref
    use getso4fromfile_m, only: getso4fromfile

    ! Input:

    real, intent(in):: r_day                   ! Day of integration
    LOGICAL, intent(in):: first                 ! First timestep
    ! (and therefore initialization necessary)

    ! Output:

    real sulfate (klon, klev)  ! Mass of sulfate (monthly mean data,
    !  from file) [ug SO4/m3]

    ! Local Variables:

    INTEGER i, ig, k, it
    INTEGER j, iday, ny, iyr, iyr1, iyr2
    parameter (ny=jjm+1)

    CHARACTER(len=4) cyear

    INTEGER im, day1, day2, im2
    double precision so4_1(iim, jjm+1, klev, 12)
    double precision so4_2(iim, jjm+1, klev, 12)   ! The sulfate distributions

    double precision so4(klon, klev, 12)  ! SO4 in right dimension
    SAVE so4
    double precision so4_out(klon, klev)
    SAVE so4_out

    LOGICAL lnewday
    LOGICAL lonlyone
    PARAMETER (lonlyone=.FALSE.)

    !--------------------------------------------------------------------

    iday = INT(r_day)

    ! Get the year of the run
    iyr  = iday/360

    ! Get the day of the actual year:
    iday = iday-iyr*360

    ! 0.02 is about 0.5/24, namly less than half an hour
    lnewday = (r_day-FLOAT(iday).LT.0.02)

    ! All has to be done only, if a new day begins!

    IF (lnewday.OR.first) THEN
       im = iday/30 +1 ! the actual month
       ! annee_ref is the initial year (defined in temps.h)
       iyr = iyr + annee_ref

       ! Do I have to read new data? (Is this the first day of a year?)
       IF (first.OR.iday.EQ.1.) THEN
          ! Initialize values
          DO it=1,12
             DO k=1,klev
                DO i=1,klon
                   so4(i,k,it)=0.
                ENDDO
             ENDDO
          ENDDO

          IF (iyr .lt. 1850) THEN
             cyear='.nat'
             WRITE(*,*) 'getso4  iyr=', iyr,'   ',cyear
             CALL getso4fromfile(cyear, so4_1)
          ELSE IF (iyr .ge. 2100) THEN
             cyear='2100'
             WRITE(*,*) 'getso4  iyr=', iyr,'   ',cyear
             CALL getso4fromfile(cyear, so4_1)
          ELSE

             ! Read in data:
             ! a) from actual 10-yr-period

             IF (iyr.LT.1900) THEN
                iyr1 = 1850
                iyr2 = 1900
             ELSE IF (iyr.ge.1900.and.iyr.lt.1920) THEN
                iyr1 = 1900
                iyr2 = 1920
             ELSE
                iyr1 = INT(iyr/10)*10
                iyr2 = INT(1+iyr/10)*10
             ENDIF
             WRITE(cyear,'(I4)') iyr1
             WRITE(*,*) 'getso4  iyr=', iyr,'   ',cyear
             CALL getso4fromfile(cyear, so4_1)

             ! If to read two decades:
             IF (.NOT.lonlyone) THEN

                ! b) from the next following one
                WRITE(cyear,'(I4)') iyr2
                WRITE(*,*) 'getso4  iyr=', iyr,'   ',cyear
                CALL getso4fromfile(cyear, so4_2)

             ENDIF

             ! Interpolate linarily to the actual year:
             DO it=1,12
                DO k=1,klev
                   DO j=1,jjm
                      DO i=1,iim
                         so4_1(i,j,k,it)=so4_1(i,j,k,it) &
                              - FLOAT(iyr-iyr1)/FLOAT(iyr2-iyr1) &
                              * (so4_1(i,j,k,it) - so4_2(i,j,k,it))
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO

          ENDIF !lonlyone

          ! Transform the horizontal 2D-field into the physics-field
          ! (Also the levels and the latitudes have to be inversed)

          DO it=1,12
             DO k=1, klev
                ! a) at the poles, use the zonal mean:
                DO i=1,iim
                   ! North pole
                   so4(1,k,it)=so4(1,k,it)+so4_1(i,jjm+1,klev+1-k,it)
                   ! South pole
                   so4(klon,k,it)=so4(klon,k,it)+so4_1(i,1,klev+1-k,it)
                ENDDO
                so4(1,k,it)=so4(1,k,it)/FLOAT(iim)
                so4(klon,k,it)=so4(klon,k,it)/FLOAT(iim)

                ! b) the values between the poles:
                ig=1
                DO j=2,jjm
                   DO i=1,iim
                      ig=ig+1
                      if (ig.gt.klon) write (*,*) 'shit'
                      so4(ig,k,it) = so4_1(i,jjm+1-j,klev+1-k,it)
                   ENDDO
                ENDDO
                IF (ig.NE.klon-1) STOP 'Error in readsulfate (var conversion)'
             ENDDO ! Loop over k (vertical)
          ENDDO ! Loop over it (months)

       ENDIF ! Had to read new data?

       ! Interpolate to actual day:
       IF (iday.LT.im*30-15) THEN
          ! in the first half of the month use month before and actual month
          im2=im-1
          day2 = im2*30-15
          day1 = im2*30+15
          IF (im2.LE.0) THEN
             ! the month is january, thus the month before december
             im2=12
          ENDIF
          DO k=1,klev
             DO i=1,klon
                sulfate(i,k) = so4(i,k,im2)   &
                     - FLOAT(iday-day2)/FLOAT(day1-day2) &
                     * (so4(i,k,im2) - so4(i,k,im))
                IF (sulfate(i,k).LT.0.) THEN
                   IF (iday-day2.LT.0.) write(*,*) 'iday-day2',iday-day2
                   IF (so4(i,k,im2) - so4(i,k,im).LT.0.) &
                        write(*,*) 'so4(i,k,im2) - so4(i,k,im)', &
                        so4(i,k,im2) - so4(i,k,im)
                   IF (day1-day2.LT.0.) write(*,*) 'day1-day2',day1-day2
                   stop 'sulfate'
                endif
             ENDDO
          ENDDO
       ELSE
          ! the second half of the month
          im2=im+1
          IF (im2.GT.12) THEN
             ! the month is december, the following thus january
             im2=1
          ENDIF
          day2 = im*30-15
          day1 = im*30+15
          DO k=1,klev
             DO i=1,klon
                sulfate(i,k) = so4(i,k,im2)   &
                     - FLOAT(iday-day2)/FLOAT(day1-day2) &
                     * (so4(i,k,im2) - so4(i,k,im))
                IF (sulfate(i,k).LT.0.) THEN
                   IF (iday-day2.LT.0.) write(*,*) 'iday-day2',iday-day2
                   IF (so4(i,k,im2) - so4(i,k,im).LT.0.) &
                        write(*,*) 'so4(i,k,im2) - so4(i,k,im)', &
                        so4(i,k,im2) - so4(i,k,im)
                   IF (day1-day2.LT.0.) write(*,*) 'day1-day2',day1-day2
                   stop 'sulfate'
                endif
             ENDDO
          ENDDO
       ENDIF

       !JLD      ! The sulfate concentration [molec cm-3] is read in.
       !JLD      ! Convert it into mass [ug SO4/m3]
       !JLD      ! masse_so4 in [g/mol], n_avogadro in [molec/mol]
       ! The sulfate mass [ug SO4/m3] is read in.
       DO k=1,klev
          DO i=1,klon
             !JLD            sulfate(i,k) = sulfate(i,k)*masse_so4
             !JLD     .           /n_avogadro*1.e+12
             so4_out(i,k) = sulfate(i,k)
             IF (so4_out(i,k).LT.0)  &
                  stop 'WAS SOLL DER SCHEISS ? '
          ENDDO
       ENDDO
    ELSE ! if no new day, use old data:
       DO k=1,klev
          DO i=1,klon
             sulfate(i,k) = so4_out(i,k)
             IF (so4_out(i,k).LT.0)  &
                  stop 'WAS SOLL DER SCHEISS ? '
          ENDDO
       ENDDO
    ENDIF ! Did I have to do anything (was it a new day?)

  END SUBROUTINE readsulfate

end module readsulfate_m
