module readsulfate_preind_m

  IMPLICIT none

contains

  SUBROUTINE readsulfate_preind(r_day, first, pi_sulfate)

    ! Read in /calculate pre-industrial values of sulfate

    use dimens_m
    use dimphy
    use temps
    use SUPHEC_M
    use chem
    use getso4fromfile_m, only: getso4fromfile

    ! Content:
    ! --------
    ! This routine reads in monthly mean values of sulfate aerosols and
    ! returns a linearly interpolated daily-mean field.
    !
    ! It does so for the preindustriel values of the sulfate, to a large part
    ! analogous to the routine readsulfate.
    !
    ! Only Pb: Variables must be saved and don t have to be overwritten!
    !
    ! Author:
    ! -------
    ! Johannes Quaas (quaas@lmd.jussieu.fr)
    ! 26/06/01
    !
    ! Input:
    ! ------
    real, intent(in)::  r_day                   ! Day of integration
    LOGICAL, intent(in):: first                 ! First timestep
    ! (and therefore initialization necessary)
    !
    ! Output:
    ! -------
    real pi_sulfate (klon, klev)  ! Number conc. sulfate (monthly mean data,
    !  from file)
    !
    ! Local Variables:
    ! ----------------
    INTEGER i, ig, k, it
    INTEGER j, iday, ny, iyr
    parameter (ny=jjm+1)

    INTEGER im, day1, day2, im2
    double precision pi_so4_1(iim, jjm+1, klev, 12)

    double precision pi_so4(klon, klev, 12)  ! SO4 in right dimension
    SAVE pi_so4
    double precision pi_so4_out(klon, klev)
    SAVE pi_so4_out

    CHARACTER(len=4) cyear
    LOGICAL lnewday



    iday = INT(r_day)

    ! Get the year of the run
    iyr  = iday/360

    ! Get the day of the actual year:
    iday = iday-iyr*360

    ! 0.02 is about 0.5/24, namly less than half an hour
    lnewday = (r_day-FLOAT(iday).LT.0.02)

    ! ---------------------------------------------
    ! All has to be done only, if a new day begins!
    ! ---------------------------------------------

    IF (lnewday.OR.first) THEN
       im = iday/30 +1 ! the actual month

       ! annee_ref is the initial year (defined in temps.h)
       iyr = iyr + annee_ref


       IF (first) THEN
          cyear='.nat'
          CALL getso4fromfile(cyear,pi_so4_1)

          ! Transform the horizontal 2D-field into the physics-field
          ! (Also the levels and the latitudes have to be inversed)

          ! Initialize field
          DO it=1,12
             DO k=1,klev
                DO i=1,klon
                   pi_so4(i,k,it)=0.
                ENDDO
             ENDDO
          ENDDO

          write (*,*) 'preind: finished reading', FLOAT(iim)
          DO it=1,12
             DO k=1, klev
                ! a) at the poles, use the zonal mean:
                DO i=1,iim
                   ! North pole
                   pi_so4(1,k,it)=pi_so4(1,k,it)+pi_so4_1(i,jjm+1,klev+1-k,it)
                   ! South pole
                   pi_so4(klon,k,it)=pi_so4(klon,k,it)+pi_so4_1(i,1,klev+1-k,it)
                ENDDO
                pi_so4(1,k,it)=pi_so4(1,k,it)/FLOAT(iim)
                pi_so4(klon,k,it)=pi_so4(klon,k,it)/FLOAT(iim)

                ! b) the values between the poles:
                ig=1
                DO j=2,jjm
                   DO i=1,iim
                      ig=ig+1
                      if (ig.gt.klon) write (*,*) 'shit'
                      pi_so4(ig,k,it) = pi_so4_1(i,jjm+1-j,klev+1-k,it)
                   ENDDO
                ENDDO
                IF (ig.NE.klon-1) STOP 'Error in readsulfate (var conversion)'
             ENDDO ! Loop over k (vertical)
          ENDDO ! Loop over it (months)

       ENDIF                     ! Had to read new data?


       ! Interpolate to actual day:
       IF (iday.LT.im*30-15) THEN
          ! in the first half of the month use month before and actual month
          im2=im-1
          day1 = im2*30+15
          day2 = im2*30-15
          IF (im2.LE.0) THEN
             ! the month is january, thus the month before december
             im2=12
          ENDIF
          DO k=1,klev
             DO i=1,klon
                pi_sulfate(i,k) = pi_so4(i,k,im2)   &
                     - FLOAT(iday-day2)/FLOAT(day1-day2) &
                     * (pi_so4(i,k,im2) - pi_so4(i,k,im))
                IF (pi_sulfate(i,k).LT.0.) THEN
                   IF (iday-day2.LT.0.) write(*,*) 'iday-day2',iday-day2
                   IF (pi_so4(i,k,im2) - pi_so4(i,k,im).LT.0.) &
                        write(*,*) 'pi_so4(i,k,im2) - pi_so4(i,k,im)', &
                        pi_so4(i,k,im2) - pi_so4(i,k,im)
                   IF (day1-day2.LT.0.) write(*,*) 'day1-day2',day1-day2
                   stop 'pi_sulfate'
                endif
             ENDDO
          ENDDO
       ELSE
          ! the second half of the month
          im2=im+1
          day1 = im*30+15
          IF (im2.GT.12) THEN
             ! the month is december, the following thus january
             im2=1
          ENDIF
          day2 = im*30-15

          DO k=1,klev
             DO i=1,klon
                pi_sulfate(i,k) = pi_so4(i,k,im2)   &
                     - FLOAT(iday-day2)/FLOAT(day1-day2) &
                     * (pi_so4(i,k,im2) - pi_so4(i,k,im))
                IF (pi_sulfate(i,k).LT.0.) THEN
                   IF (iday-day2.LT.0.) write(*,*) 'iday-day2',iday-day2
                   IF (pi_so4(i,k,im2) - pi_so4(i,k,im).LT.0.) &
                        write(*,*) 'pi_so4(i,k,im2) - pi_so4(i,k,im)', &
                        pi_so4(i,k,im2) - pi_so4(i,k,im)
                   IF (day1-day2.LT.0.) write(*,*) 'day1-day2',day1-day2
                   stop 'pi_sulfate'
                endif
             ENDDO
          ENDDO
       ENDIF


       !JLD      ! The sulfate concentration [molec cm-3] is read in.
       !JLD      ! Convert it into mass [ug SO4/m3]
       !JLD      ! masse_so4 in [g/mol], n_avogadro in [molec/mol]
       DO k=1,klev
          DO i=1,klon
             !JLD            pi_sulfate(i,k) = pi_sulfate(i,k)*masse_so4
             !JLD     .           /n_avogadro*1.e+12
             pi_so4_out(i,k) = pi_sulfate(i,k)
          ENDDO
       ENDDO

    ELSE ! If no new day, use old data:
       DO k=1,klev
          DO i=1,klon
             pi_sulfate(i,k) = pi_so4_out(i,k)
          ENDDO
       ENDDO
    ENDIF ! Was this the beginning of a new day?

  END SUBROUTINE readsulfate_preind

end module readsulfate_preind_m
