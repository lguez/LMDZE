module getso4fromfile_m

  implicit none

contains

  SUBROUTINE getso4fromfile(cyr, so4)

    ! Routine for reading SO4 data from files

    USE dimens_m, ONLY: iim, jjm
    USE dimphy, ONLY: klev
    USE netcdf, ONLY: nf90_nowrite
    USE netcdf95, ONLY: nf95_close, nf95_get_var, nf95_inq_varid, nf95_open

    CHARACTER(len=*), intent(in):: cyr
    double precision so4(iim, jjm + 1, klev, 12)

    ! Local:

    CHARACTER(len=15) fname

    CHARACTER*6 cvar
    INTEGER START(3), COUNT(3)
    INTEGER  NCID, VARID
    INTEGER imth, i, j, k

    double precision so4mth(iim, jjm + 1, klev)

    !---------------------------------------------------------------------

    fname = 'so4.run'//cyr//'.cdf'

    write(*,*) 'reading ', fname
    call NF95_OPEN(fname, NF90_NOWRITE, NCID)

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
       count(2)=jjm + 1
       count(3)=klev
       call NF95_INQ_VARID(NCID, cvar, VARID)
       write(*,*) ncid,imth,cvar, varid

       call NF95_GET_VAR(NCID, VARID, so4mth, START,COUNT)

       DO k=1,klev
          DO j=1,jjm+1
             DO i=1,iim
                IF (so4mth(i,j,k).LT.0.) then
                   write(*,*) 'this is shit'
                   write(*,*) 'so4(',i,j,k,') =',so4mth(i,j,k)
                endif
                so4(i,j,k,imth)=so4mth(i,j,k)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    call NF95_CLOSE(NCID)

  END SUBROUTINE getso4fromfile

end module getso4fromfile_m
