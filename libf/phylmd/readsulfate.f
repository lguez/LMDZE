!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/readsulfate.F,v 1.2 2005/05/19 08:27:15 fairhead Exp $
!
      SUBROUTINE readsulfate(r_day, first, sulfate)
      
      use dimens_m      
      use dimphy      
      use temps      
      use SUPHEC_M
      use chem      
      IMPLICIT none
      
c Content: 
c --------
c This routine reads in monthly mean values of sulfate aerosols and 
c returns a linearly interpolated dayly-mean field.      
c 
c
c Author:
c -------
c Johannes Quaas (quaas@lmd.jussieu.fr) 
c 26/04/01
c
c Modifications:
c --------------
c 21/06/01: Make integrations of more than one year possible ;-)     
c           ATTENTION!! runs are supposed to start with Jan, 1. 1930
c                       (rday=1)      
c
c 27/06/01: Correction: The model always has 360 days per year!
c 27/06/01: SO4 concentration rather than mixing ratio      
c 27/06/01: 10yr-mean-values to interpolate     
c 20/08/01: Correct the error through integer-values in interpolations      
c 21/08/01: Introduce flag to read in just one decade
c      
c 
c Input:
c ------
      REAL*8, intent(in):: r_day                   ! Day of integration
      LOGICAL, intent(in):: first                 ! First timestep 
                                    ! (and therefore initialization necessary)
c      
c Output:      
c -------     
      REAL*8  sulfate (klon, klev)  ! Mass of sulfate (monthly mean data, 
                                  !  from file) [ug SO4/m3]
c      
c Local Variables:
c ----------------      
      INTEGER i, ig, k, it
      INTEGER j, iday, ny, iyr, iyr1, iyr2
      parameter (ny=jjm+1)
      
      INTEGER ismaller
CJLD      INTEGER idec1, idec2 ! The two decadal data read ini
      CHARACTER*4 cyear
      
      INTEGER im, day1, day2, im2
      REAL*8 so4_1(iim, jjm+1, klev, 12)
      REAL*8 so4_2(iim, jjm+1, klev, 12)   ! The sulfate distributions
      
      REAL*8 so4(klon, klev, 12)  ! SO4 in right dimension
      SAVE so4
      REAL*8 so4_out(klon, klev)
      SAVE so4_out
      
      LOGICAL lnewday 
      LOGICAL lonlyone
      PARAMETER (lonlyone=.FALSE.)

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
                  so4_1(i,j,k,it)=so4_1(i,j,k,it)
     .                 - FLOAT(iyr-iyr1)/FLOAT(iyr2-iyr1)
     .                 * (so4_1(i,j,k,it) - so4_2(i,j,k,it))
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
               sulfate(i,k) = so4(i,k,im2)  
     .              - FLOAT(iday-day2)/FLOAT(day1-day2)
     .              * (so4(i,k,im2) - so4(i,k,im))
               IF (sulfate(i,k).LT.0.) THEN
                  IF (iday-day2.LT.0.) write(*,*) 'iday-day2',iday-day2
                  IF (so4(i,k,im2) - so4(i,k,im).LT.0.)
     . write(*,*) 'so4(i,k,im2) - so4(i,k,im)',
     . so4(i,k,im2) - so4(i,k,im)
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
               sulfate(i,k) = so4(i,k,im2)  
     .              - FLOAT(iday-day2)/FLOAT(day1-day2)
     .              * (so4(i,k,im2) - so4(i,k,im))
               IF (sulfate(i,k).LT.0.) THEN
                  IF (iday-day2.LT.0.) write(*,*) 'iday-day2',iday-day2
                  IF (so4(i,k,im2) - so4(i,k,im).LT.0.)
     . write(*,*) 'so4(i,k,im2) - so4(i,k,im)',
     . so4(i,k,im2) - so4(i,k,im)
                  IF (day1-day2.LT.0.) write(*,*) 'day1-day2',day1-day2
                  stop 'sulfate'
               endif
            ENDDO
         ENDDO
      ENDIF

      
CJLD      ! The sulfate concentration [molec cm-3] is read in. 
CJLD      ! Convert it into mass [ug SO4/m3]
CJLD      ! masse_so4 in [g/mol], n_avogadro in [molec/mol]
      ! The sulfate mass [ug SO4/m3] is read in. 
      DO k=1,klev
         DO i=1,klon
CJLD            sulfate(i,k) = sulfate(i,k)*masse_so4
CJLD     .           /n_avogadro*1.e+12
            so4_out(i,k) = sulfate(i,k)
            IF (so4_out(i,k).LT.0) 
     .          stop 'WAS SOLL DER SCHEISS ? '
         ENDDO
      ENDDO
      ELSE ! if no new day, use old data:
      DO k=1,klev
         DO i=1,klon
            sulfate(i,k) = so4_out(i,k)
            IF (so4_out(i,k).LT.0) 
     .          stop 'WAS SOLL DER SCHEISS ? '
         ENDDO
      ENDDO
         

      ENDIF ! Did I have to do anything (was it a new day?)
      
      RETURN
      END

      
      
      
      
c-----------------------------------------------------------------------------
c Read in /calculate pre-industrial values of sulfate      
c-----------------------------------------------------------------------------
      
      SUBROUTINE readsulfate_preind (r_day, first, pi_sulfate)
      
      use dimens_m      
      use dimphy      
      use temps      
      use SUPHEC_M
      use chem      
      IMPLICIT none
      
c Content: 
c --------
c This routine reads in monthly mean values of sulfate aerosols and 
c returns a linearly interpolated dayly-mean field.      
c 
c It does so for the preindustriel values of the sulfate, to a large part
c analogous to the routine readsulfate above.      
c
c Only Pb: Variables must be saved and don t have to be overwritten!
c      
c Author:
c -------
c Johannes Quaas (quaas@lmd.jussieu.fr) 
c 26/06/01
c
c Modifications:
c --------------
c see above 
c      
c 
c Input:
c ------
      REAL*8, intent(in)::  r_day                   ! Day of integration
      LOGICAL, intent(in):: first                 ! First timestep 
                                    ! (and therefore initialization necessary)
c      
c Output:      
c -------     
      REAL*8  pi_sulfate (klon, klev)  ! Number conc. sulfate (monthly mean data, 
                                  !  from file)
c      
c Local Variables:
c ----------------      
      INTEGER i, ig, k, it
      INTEGER j, iday, ny, iyr
      parameter (ny=jjm+1)
      
      INTEGER im, day1, day2, im2, ismaller
      REAL*8 pi_so4_1(iim, jjm+1, klev, 12)
      
      REAL*8 pi_so4(klon, klev, 12)  ! SO4 in right dimension
      SAVE pi_so4
      REAL*8 pi_so4_out(klon, klev)
      SAVE pi_so4_out
      
      CHARACTER*4 cyear
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
               pi_sulfate(i,k) = pi_so4(i,k,im2)  
     .              - FLOAT(iday-day2)/FLOAT(day1-day2)
     .              * (pi_so4(i,k,im2) - pi_so4(i,k,im))
               IF (pi_sulfate(i,k).LT.0.) THEN
                  IF (iday-day2.LT.0.) write(*,*) 'iday-day2',iday-day2
                  IF (pi_so4(i,k,im2) - pi_so4(i,k,im).LT.0.)
     . write(*,*) 'pi_so4(i,k,im2) - pi_so4(i,k,im)',
     . pi_so4(i,k,im2) - pi_so4(i,k,im)
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
               pi_sulfate(i,k) = pi_so4(i,k,im2)  
     .              - FLOAT(iday-day2)/FLOAT(day1-day2)
     .              * (pi_so4(i,k,im2) - pi_so4(i,k,im))
               IF (pi_sulfate(i,k).LT.0.) THEN
                  IF (iday-day2.LT.0.) write(*,*) 'iday-day2',iday-day2
                  IF (pi_so4(i,k,im2) - pi_so4(i,k,im).LT.0.)
     . write(*,*) 'pi_so4(i,k,im2) - pi_so4(i,k,im)',
     . pi_so4(i,k,im2) - pi_so4(i,k,im)
                  IF (day1-day2.LT.0.) write(*,*) 'day1-day2',day1-day2
                  stop 'pi_sulfate'
               endif
            ENDDO
         ENDDO
      ENDIF

      
CJLD      ! The sulfate concentration [molec cm-3] is read in. 
CJLD      ! Convert it into mass [ug SO4/m3]
CJLD      ! masse_so4 in [g/mol], n_avogadro in [molec/mol]
      DO k=1,klev
         DO i=1,klon
CJLD            pi_sulfate(i,k) = pi_sulfate(i,k)*masse_so4
CJLD     .           /n_avogadro*1.e+12
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
      RETURN
      END

      
      
      
      
      
      
      
      
      
c-----------------------------------------------------------------------------
c Routine for reading SO4 data from files
c-----------------------------------------------------------------------------
            

      SUBROUTINE getso4fromfile (cyr, so4)
      use dimens_m      
      use dimphy
      include "netcdf.inc"
      CHARACTER*15 fname
      CHARACTER*4 cyr
      
      CHARACTER*6 cvar
      INTEGER START(3), COUNT(3)
      INTEGER  STATUS, NCID, VARID
      INTEGER imth, i, j, k, ny
      PARAMETER (ny=jjm+1)
      
            
      REAL*8 so4mth(iim, ny, klev)
c      REAL*8 so4mth(klev, ny, iim)
      REAL*8 so4(iim, ny, klev, 12)

 
      fname = 'so4.run'//cyr//'.cdf'

      write (*,*) 'reading ', fname
      STATUS = NF_OPEN (fname, NF_NOWRITE, NCID)
      IF (STATUS .NE. NF_NOERR) write (*,*) 'err in open ',status
            
      DO imth=1, 12
         IF (imth.eq.1) THEN
            cvar='SO4JAN'
         ELSEIF (imth.eq.2) THEN
            cvar='SO4FEB'
         ELSEIF (imth.eq.3) THEN
            cvar='SO4MAR'
         ELSEIF (imth.eq.4) THEN
            cvar='SO4APR'
         ELSEIF (imth.eq.5) THEN
            cvar='SO4MAY'
         ELSEIF (imth.eq.6) THEN
            cvar='SO4JUN'
         ELSEIF (imth.eq.7) THEN
            cvar='SO4JUL'
         ELSEIF (imth.eq.8) THEN
            cvar='SO4AUG'
         ELSEIF (imth.eq.9) THEN
            cvar='SO4SEP'
         ELSEIF (imth.eq.10) THEN
            cvar='SO4OCT'
         ELSEIF (imth.eq.11) THEN
            cvar='SO4NOV'
         ELSEIF (imth.eq.12) THEN
            cvar='SO4DEC'
         ENDIF
         start(1)=1
         start(2)=1
         start(3)=1
         count(1)=iim
         count(2)=ny
         count(3)=klev
c         write(*,*) 'here i am'
         STATUS = NF_INQ_VARID (NCID, cvar, VARID)
         write (*,*) ncid,imth,cvar, varid
c         STATUS = NF_INQ_VARID (NCID, VARMONTHS(i), VARID(i))
         IF (STATUS .NE. NF_NOERR) write (*,*) 'err in read ',status      
         STATUS = NF_GET_VARA_DOUBLE
     .    (NCID, VARID, START,COUNT, so4mth)
         IF (STATUS .NE. NF_NOERR) write (*,*) 'err in read data',status
         
         DO k=1,klev
            DO j=1,jjm+1
               DO i=1,iim
                  IF (so4mth(i,j,k).LT.0.) then
                     write(*,*) 'this is shit'
                     write(*,*) 'so4(',i,j,k,') =',so4mth(i,j,k)
                  endif
                  so4(i,j,k,imth)=so4mth(i,j,k)
c                  so4(i,j,k,imth)=so4mth(k,j,i)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
      STATUS = NF_CLOSE(NCID)
      END ! subroutine getso4fromfile
      















