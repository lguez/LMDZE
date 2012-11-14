SUBROUTINE getso4fromfile (cyr, so4)

  ! Routine for reading SO4 data from files

  use dimens_m
  use dimphy
  use netcdf
  CHARACTER*15 fname
  CHARACTER*4 cyr

  CHARACTER*6 cvar
  INTEGER START(3), COUNT(3)
  INTEGER  STATUS, NCID, VARID
  INTEGER imth, i, j, k, ny
  PARAMETER (ny=jjm+1)


  double precision so4mth(iim, ny, klev)
  double precision so4(iim, ny, klev, 12)


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
     !         write(*,*) 'here i am'
     STATUS = NF_INQ_VARID (NCID, cvar, VARID)
     write (*,*) ncid,imth,cvar, varid
     !         STATUS = NF_INQ_VARID (NCID, VARMONTHS(i), VARID(i))
     IF (STATUS .NE. NF_NOERR) write (*,*) 'err in read ',status
     STATUS = NF_GET_VARA_DOUBLE &
          (NCID, VARID, START,COUNT, so4mth)
     IF (STATUS .NE. NF_NOERR) write (*,*) 'err in read data',status

     DO k=1,klev
        DO j=1,jjm+1
           DO i=1,iim
              IF (so4mth(i,j,k).LT.0.) then
                 write(*,*) 'this is shit'
                 write(*,*) 'so4(',i,j,k,') =',so4mth(i,j,k)
              endif
              so4(i,j,k,imth)=so4mth(i,j,k)
              !                  so4(i,j,k,imth)=so4mth(k,j,i)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  STATUS = NF_CLOSE(NCID)
END SUBROUTINE getso4fromfile
