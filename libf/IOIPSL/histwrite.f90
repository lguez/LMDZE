MODULE histwrite_m

  ! From histcom.f90, v 2.1 2004/04/21 09:27:10

  use histcom_var

  implicit none

  PRIVATE
  PUBLIC histwrite

  INTERFACE histwrite
     !- The "histwrite" procedures give the data to the input-output system.
     !- They trigger the operations to be performed
     !- and the writing to the file if needed.

     !- We test the work to be done at this time here so that at a
     !- later stage we can call different operations and write subroutines
     !- for the REAL and INTEGER interfaces.

     ! INTEGER, INTENT(IN):: pfileid
     ! The ID of the file on which this variable is to be written.
     ! The variable should have been defined in this file before.

     ! CHARACTER(LEN=*), INTENT(IN):: pvarname
     ! short name of the variable

     ! INTEGER, INTENT(IN):: pitau
     ! current timestep

     ! REAL, INTENT(IN):: pdata(:) or (:, :) or (:, :, :)
     ! values of the variable

     ! INTEGER, INTENT(IN):: nbindex
     ! number of indices provided
     ! If it is equal to the size of the full field as provided in histdef
     ! then nothing is done.

     ! INTEGER, INTENT(IN):: nindex(nbindex)
     ! The indices used to expand the variable (pdata) onto the full field

     ! The difference between the procedures is the rank of "pdata".

     MODULE PROCEDURE histwrite_r1d, histwrite_r2d, histwrite_r3d
  END INTERFACE

CONTAINS

  SUBROUTINE histwrite_r1d(pfileid, pvarname, pitau, pdata)

    USE errioipsl, ONLY : histerr
    use calendar, only: isittime
    USE mathelp, ONLY : mathop

    INTEGER, INTENT(IN) :: pfileid, pitau
    REAL, INTENT(IN) :: pdata(:)
    CHARACTER(LEN=*), INTENT(IN) :: pvarname

    ! Variables local to the procedure:
    integer nbindex, nindex(size(pdata))
    LOGICAL :: do_oper, do_write, largebuf
    INTEGER :: varid, io, nbpt_in, nbpt_out
    REAL, ALLOCATABLE, SAVE :: buff_tmp(:)
    INTEGER, SAVE :: buff_tmp_sz
    CHARACTER(LEN=7) :: tmp_opp

    !---------------------------------------------------------------------

    nbindex = size(nindex)
    nindex = 0

    ! 1.0 Try to catch errors like specifying the wrong file ID.
    !     Thanks Marine for showing us what errors users can make !

    IF ( (pfileid < 1).OR.(pfileid > nb_files) ) THEN
       CALL histerr (3, "histwrite", &
            &    'Illegal file ID in the histwrite of variable', pvarname, ' ')
    ENDIF

    ! 1.1 Find the id of the variable to be written and the real time

    CALL histvar_seq (pfileid, pvarname, varid)

    ! 2.0 do nothing for never operation

    tmp_opp = topp(pfileid, varid)

    IF (TRIM(tmp_opp) == "never") THEN
       last_opp_chk(pfileid, varid) = -99
       last_wrt_chk(pfileid, varid) = -99
    ENDIF

    ! 3.0 We check if we need to do an operation

    IF (last_opp_chk(pfileid, varid) == pitau) THEN
       CALL histerr (3, "histwrite", &
            &    'This variable as already been analysed at the present', &
            &    'time step', ' ')
    ENDIF

    CALL isittime &
         &  (pitau, date0(pfileid), deltat(pfileid), freq_opp(pfileid, varid), &
         &   last_opp(pfileid, varid), last_opp_chk(pfileid, varid), do_oper)

    ! 4.0 We check if we need to write the data

    IF (last_wrt_chk(pfileid, varid) == pitau) THEN
       CALL histerr (3, "histwrite", &
            &    'This variable as already been written for the present', &
            &    'time step', ' ')
    ENDIF

    CALL isittime &
         &  (pitau, date0(pfileid), deltat(pfileid), freq_wrt(pfileid, varid), &
         &   last_wrt(pfileid, varid), last_wrt_chk(pfileid, varid), do_write)

    ! 5.0 histwrite called

    IF (do_oper.OR.do_write) THEN

       !-- 5.1 Get the sizes of the data we will handle

       IF (datasz_in(pfileid, varid, 1) <= 0) THEN
          !---- There is the risk here that the user has over-sized the array.
          !---- But how can we catch this ?
          !---- In the worst case we will do impossible operations
          !---- on part of the data !
          datasz_in(pfileid, varid, 1) = SIZE(pdata)
          datasz_in(pfileid, varid, 2) = -1
          datasz_in(pfileid, varid, 3) = -1
       ENDIF

       !-- 5.2 The maximum size of the data will give the size of the buffer

       IF (datasz_max(pfileid, varid) <= 0) THEN
          largebuf = .FALSE.
          DO io=1, nbopp(pfileid, varid)
             IF (INDEX(fuchnbout, sopps(pfileid, varid, io)) > 0) THEN
                largebuf = .TRUE.
             ENDIF
          ENDDO
          IF (largebuf) THEN
             datasz_max(pfileid, varid) = &
                  &        scsize(pfileid, varid, 1) &
                  &       *scsize(pfileid, varid, 2) &
                  &       *scsize(pfileid, varid, 3)
          ELSE
             datasz_max(pfileid, varid) = &
                  &        datasz_in(pfileid, varid, 1)
          ENDIF
       ENDIF

       IF (.NOT.ALLOCATED(buff_tmp)) THEN
          ALLOCATE (buff_tmp(datasz_max(pfileid, varid)))
          buff_tmp_sz = datasz_max(pfileid, varid)
       ELSE IF (datasz_max(pfileid, varid) > buff_tmp_sz) THEN
          DEALLOCATE (buff_tmp)
          ALLOCATE (buff_tmp(datasz_max(pfileid, varid)))
          buff_tmp_sz = datasz_max(pfileid, varid)
       ENDIF

       !-- We have to do the first operation anyway.
       !-- Thus we do it here and change the ranke
       !-- of the data at the same time. This should speed up things.

       nbpt_in = datasz_in(pfileid, varid, 1)
       nbpt_out = datasz_max(pfileid, varid)
       CALL mathop (sopps(pfileid, varid, 1), nbpt_in, pdata, &
            &               missing_val, nbindex, nindex, &
            &               scal(pfileid, varid, 1), nbpt_out, buff_tmp)
       CALL histwrite_real (pfileid, varid, pitau, nbpt_out, &
            &            buff_tmp, nbindex, nindex, do_oper, do_write)
    ENDIF

    ! 6.0 Manage time steps

    IF ((TRIM(tmp_opp) /= "once").AND.(TRIM(tmp_opp) /= "never")) THEN
       last_opp_chk(pfileid, varid) = pitau
       last_wrt_chk(pfileid, varid) = pitau
    ELSE
       last_opp_chk(pfileid, varid) = -99
       last_wrt_chk(pfileid, varid) = -99
    ENDIF
    !---------------------------
  END SUBROUTINE histwrite_r1d

  !===

  SUBROUTINE histwrite_r2d (pfileid, pvarname, pitau, pdata)
    !---------------------------------------------------------------------

    use calendar, only: isittime
    USE errioipsl, ONLY : histerr
    USE mathelp, ONLY : mathop

    INTEGER, INTENT(IN) :: pfileid, pitau
    REAL, DIMENSION(:, :), INTENT(IN) :: pdata
    CHARACTER(LEN=*), INTENT(IN) :: pvarname

    integer nbindex, nindex(size(pdata))
    LOGICAL :: do_oper, do_write, largebuf
    INTEGER :: varid, io, nbpt_in(1:2), nbpt_out
    REAL, ALLOCATABLE, SAVE :: buff_tmp(:)
    INTEGER, SAVE :: buff_tmp_sz
    CHARACTER(LEN=7) :: tmp_opp

    !---------------------------------------------------------------------

    nbindex = size(nindex)
    nindex = 0

    ! 1.0 Try to catch errors like specifying the wrong file ID.
    !     Thanks Marine for showing us what errors users can make !

    IF ( (pfileid < 1).OR.(pfileid > nb_files) ) THEN
       CALL histerr (3, "histwrite", &
            &    'Illegal file ID in the histwrite of variable', pvarname, ' ')
    ENDIF

    ! 1.1 Find the id of the variable to be written and the real time

    CALL histvar_seq (pfileid, pvarname, varid)

    ! 2.0 do nothing for never operation

    tmp_opp = topp(pfileid, varid)

    IF (TRIM(tmp_opp) == "never") THEN
       last_opp_chk(pfileid, varid) = -99
       last_wrt_chk(pfileid, varid) = -99
    ENDIF

    ! 3.0 We check if we need to do an operation

    IF (last_opp_chk(pfileid, varid) == pitau) THEN
       CALL histerr (3, "histwrite", &
            &    'This variable as already been analysed at the present', &
            &    'time step', ' ')
    ENDIF

    CALL isittime &
         &  (pitau, date0(pfileid), deltat(pfileid), freq_opp(pfileid, varid), &
         &   last_opp(pfileid, varid), last_opp_chk(pfileid, varid), do_oper)

    ! 4.0 We check if we need to write the data

    IF (last_wrt_chk(pfileid, varid) == pitau) THEN
       CALL histerr (3, "histwrite", &
            &    'This variable as already been written for the present', &
            &    'time step', ' ')
    ENDIF

    CALL isittime &
         &  (pitau, date0(pfileid), deltat(pfileid), freq_wrt(pfileid, varid), &
         &   last_wrt(pfileid, varid), last_wrt_chk(pfileid, varid), do_write)

    ! 5.0 histwrite called

    IF (do_oper.OR.do_write) THEN

       !-- 5.1 Get the sizes of the data we will handle

       IF (datasz_in(pfileid, varid, 1) <= 0) THEN
          !---- There is the risk here that the user has over-sized the array.
          !---- But how can we catch this ?
          !---- In the worst case we will do impossible operations
          !---- on part of the data !
          datasz_in(pfileid, varid, 1) = SIZE(pdata, DIM=1)
          datasz_in(pfileid, varid, 2) = SIZE(pdata, DIM=2)
          datasz_in(pfileid, varid, 3) = -1
       ENDIF

       !-- 5.2 The maximum size of the data will give the size of the buffer

       IF (datasz_max(pfileid, varid) <= 0) THEN
          largebuf = .FALSE.
          DO io=1, nbopp(pfileid, varid)
             IF (INDEX(fuchnbout, sopps(pfileid, varid, io)) > 0) THEN
                largebuf = .TRUE.
             ENDIF
          ENDDO
          IF (largebuf) THEN
             datasz_max(pfileid, varid) = &
                  &        scsize(pfileid, varid, 1) &
                  &       *scsize(pfileid, varid, 2) &
                  &       *scsize(pfileid, varid, 3)
          ELSE
             datasz_max(pfileid, varid) = &
                  &        datasz_in(pfileid, varid, 1) &
                  &       *datasz_in(pfileid, varid, 2)
          ENDIF
       ENDIF

       IF (.NOT.ALLOCATED(buff_tmp)) THEN
          ALLOCATE (buff_tmp(datasz_max(pfileid, varid)))
          buff_tmp_sz = datasz_max(pfileid, varid)
       ELSE IF (datasz_max(pfileid, varid) > buff_tmp_sz) THEN
          DEALLOCATE (buff_tmp)
          ALLOCATE (buff_tmp(datasz_max(pfileid, varid)))
          buff_tmp_sz = datasz_max(pfileid, varid)
       ENDIF

       !-- We have to do the first operation anyway.
       !-- Thus we do it here and change the ranke
       !-- of the data at the same time. This should speed up things.

       nbpt_in(1:2) = datasz_in(pfileid, varid, 1:2)
       nbpt_out = datasz_max(pfileid, varid)
       CALL mathop (sopps(pfileid, varid, 1), nbpt_in, pdata, &
            &               missing_val, nbindex, nindex, &
            &               scal(pfileid, varid, 1), nbpt_out, buff_tmp)
       CALL histwrite_real (pfileid, varid, pitau, nbpt_out, &
            &            buff_tmp, nbindex, nindex, do_oper, do_write)
    ENDIF

    ! 6.0 Manage time steps

    IF ((TRIM(tmp_opp) /= "once").AND.(TRIM(tmp_opp) /= "never")) THEN
       last_opp_chk(pfileid, varid) = pitau
       last_wrt_chk(pfileid, varid) = pitau
    ELSE
       last_opp_chk(pfileid, varid) = -99
       last_wrt_chk(pfileid, varid) = -99
    ENDIF
    !---------------------------
  END SUBROUTINE histwrite_r2d

  !===

  SUBROUTINE histwrite_r3d (pfileid, pvarname, pitau, pdata)
    !---------------------------------------------------------------------

    use calendar, only: isittime
    USE errioipsl, ONLY : histerr
    USE mathelp, ONLY : mathop

    INTEGER, INTENT(IN) :: pfileid, pitau
    REAL, DIMENSION(:, :, :), INTENT(IN) :: pdata
    CHARACTER(LEN=*), INTENT(IN) :: pvarname

    integer nbindex, nindex(size(pdata))
    LOGICAL :: do_oper, do_write, largebuf
    INTEGER :: varid, io, nbpt_in(1:3), nbpt_out
    REAL, ALLOCATABLE, SAVE :: buff_tmp(:)
    INTEGER, SAVE :: buff_tmp_sz
    CHARACTER(LEN=7) :: tmp_opp

    !---------------------------------------------------------------------

    nbindex = size(nindex)
    nindex = 0

    ! 1.0 Try to catch errors like specifying the wrong file ID.
    !     Thanks Marine for showing us what errors users can make !

    IF ( (pfileid < 1).OR.(pfileid > nb_files) ) THEN
       CALL histerr (3, "histwrite", &
            &    'Illegal file ID in the histwrite of variable', pvarname, ' ')
    ENDIF

    ! 1.1 Find the id of the variable to be written and the real time

    CALL histvar_seq (pfileid, pvarname, varid)

    ! 2.0 do nothing for never operation

    tmp_opp = topp(pfileid, varid)

    IF (TRIM(tmp_opp) == "never") THEN
       last_opp_chk(pfileid, varid) = -99
       last_wrt_chk(pfileid, varid) = -99
    ENDIF

    ! 3.0 We check if we need to do an operation

    IF (last_opp_chk(pfileid, varid) == pitau) THEN
       CALL histerr (3, "histwrite", &
            &    'This variable as already been analysed at the present', &
            &    'time step', ' ')
    ENDIF

    CALL isittime &
         &  (pitau, date0(pfileid), deltat(pfileid), freq_opp(pfileid, varid), &
         &   last_opp(pfileid, varid), last_opp_chk(pfileid, varid), do_oper)

    ! 4.0 We check if we need to write the data

    IF (last_wrt_chk(pfileid, varid) == pitau) THEN
       CALL histerr (3, "histwrite", &
            &    'This variable as already been written for the present', &
            &    'time step', ' ')
    ENDIF

    CALL isittime &
         &  (pitau, date0(pfileid), deltat(pfileid), freq_wrt(pfileid, varid), &
         &   last_wrt(pfileid, varid), last_wrt_chk(pfileid, varid), do_write)

    ! 5.0 histwrite called

    IF (do_oper.OR.do_write) THEN

       !-- 5.1 Get the sizes of the data we will handle

       IF (datasz_in(pfileid, varid, 1) <= 0) THEN
          !---- There is the risk here that the user has over-sized the array.
          !---- But how can we catch this ?
          !---- In the worst case we will do impossible operations
          !---- on part of the data !
          datasz_in(pfileid, varid, 1) = SIZE(pdata, DIM=1)
          datasz_in(pfileid, varid, 2) = SIZE(pdata, DIM=2)
          datasz_in(pfileid, varid, 3) = SIZE(pdata, DIM=3)
       ENDIF

       !-- 5.2 The maximum size of the data will give the size of the buffer

       IF (datasz_max(pfileid, varid) <= 0) THEN
          largebuf = .FALSE.
          DO io =1, nbopp(pfileid, varid)
             IF (INDEX(fuchnbout, sopps(pfileid, varid, io)) > 0) THEN
                largebuf = .TRUE.
             ENDIF
          ENDDO
          IF (largebuf) THEN
             datasz_max(pfileid, varid) = &
                  &        scsize(pfileid, varid, 1) &
                  &       *scsize(pfileid, varid, 2) &
                  &       *scsize(pfileid, varid, 3)
          ELSE
             datasz_max(pfileid, varid) = &
                  &        datasz_in(pfileid, varid, 1) &
                  &       *datasz_in(pfileid, varid, 2) &
                  &       *datasz_in(pfileid, varid, 3)
          ENDIF
       ENDIF

       IF (.NOT.ALLOCATED(buff_tmp)) THEN
          ALLOCATE (buff_tmp(datasz_max(pfileid, varid)))
          buff_tmp_sz = datasz_max(pfileid, varid)
       ELSE IF (datasz_max(pfileid, varid) > buff_tmp_sz) THEN
          DEALLOCATE (buff_tmp)
          ALLOCATE (buff_tmp(datasz_max(pfileid, varid)))
          buff_tmp_sz = datasz_max(pfileid, varid)
       ENDIF

       !-- We have to do the first operation anyway.
       !-- Thus we do it here and change the ranke
       !-- of the data at the same time. This should speed up things.

       nbpt_in(1:3) = datasz_in(pfileid, varid, 1:3)
       nbpt_out = datasz_max(pfileid, varid)
       CALL mathop (sopps(pfileid, varid, 1), nbpt_in, pdata, &
            &               missing_val, nbindex, nindex, &
            &               scal(pfileid, varid, 1), nbpt_out, buff_tmp)
       CALL histwrite_real (pfileid, varid, pitau, nbpt_out, &
            &            buff_tmp, nbindex, nindex, do_oper, do_write)
    ENDIF

    ! 6.0 Manage time steps

    IF ((TRIM(tmp_opp) /= "once").AND.(TRIM(tmp_opp) /= "never")) THEN
       last_opp_chk(pfileid, varid) = pitau
       last_wrt_chk(pfileid, varid) = pitau
    ELSE
       last_opp_chk(pfileid, varid) = -99
       last_wrt_chk(pfileid, varid) = -99
    ENDIF
    !---------------------------
  END SUBROUTINE histwrite_r3d

  !===

  SUBROUTINE histwrite_real(pfileid, varid, pitau, nbdpt, buff_tmp, nbindex, &
       nindex, do_oper, do_write)

    !- This subroutine is internal and does the calculations and writing
    !- if needed. At a later stage it should be split into an operation
    !- and writing subroutines.
    !---------------------------------------------------------------------

    USE mathelp, ONLY : mathop, trans_buff, moycum
    use netcdf, only: NF90_PUT_VAR

    INTEGER, INTENT(IN) :: pfileid, pitau, varid, &
         &                      nbindex, nindex(nbindex), nbdpt
    REAL, DIMENSION(:)  :: buff_tmp
    LOGICAL, INTENT(IN) :: do_oper, do_write

    INTEGER :: tsz, ncid, ncvarid
    INTEGER :: i, iret, ipt, itax
    INTEGER :: io, nbin, nbout
    INTEGER, DIMENSION(4) :: corner, edges
    INTEGER :: itime

    REAL :: rtime
    CHARACTER(LEN=7) :: tmp_opp

    REAL, ALLOCATABLE, SAVE :: buff_tmp2(:)
    INTEGER, SAVE          :: buff_tmp2_sz
    REAL, ALLOCATABLE, SAVE :: buffer_used(:)
    INTEGER, SAVE          :: buffer_sz

    !---------------------------------------------------------------------

    ! The sizes which can be encoutered

    tsz = zsize(pfileid, varid, 1)*zsize(pfileid, varid, 2)*zsize(pfileid, varid, 3)

    ! 1.0 We allocate the memory needed to store the data between write
    !     and the temporary space needed for operations.
    !     We have to keep precedent buffer if needed

    IF (.NOT. ALLOCATED(buffer)) THEN
       ALLOCATE(buffer(buff_pos))
       buffer_sz = buff_pos
       buffer(:)=0.0
    ELSE IF (buffer_sz < buff_pos) THEN
       IF (SUM(buffer)/=0.0) THEN
          ALLOCATE (buffer_used(buffer_sz))
          buffer_used(:)=buffer(:)
          DEALLOCATE (buffer)
          ALLOCATE (buffer(buff_pos))
          buffer_sz = buff_pos
          buffer(:SIZE(buffer_used))=buffer_used
          DEALLOCATE (buffer_used)
       ELSE
          DEALLOCATE (buffer)
          ALLOCATE (buffer(buff_pos))
          buffer_sz = buff_pos
          buffer(:)=0.0
       ENDIF
    ENDIF

    ! The buffers are only deallocated when more space is needed. This
    ! reduces the umber of allocates but increases memory needs.

    IF (.NOT.ALLOCATED(buff_tmp2)) THEN
       ALLOCATE (buff_tmp2(datasz_max(pfileid, varid)))
       buff_tmp2_sz = datasz_max(pfileid, varid)
    ELSE IF ( datasz_max(pfileid, varid) > buff_tmp2_sz) THEN
       DEALLOCATE (buff_tmp2)
       ALLOCATE (buff_tmp2(datasz_max(pfileid, varid)))
       buff_tmp2_sz = datasz_max(pfileid, varid)
    ENDIF

    rtime = pitau * deltat(pfileid)
    tmp_opp = topp(pfileid, varid)

    ! 3.0 Do the operations or transfer the slab of data into buff_tmp

    ! 3.1 DO the Operations only if needed

    IF ( do_oper ) THEN
       i = pfileid
       nbout = nbdpt

       !-- 3.4 We continue the sequence of operations
       !--     we started in the interface routine

       DO io = 2, nbopp(i, varid), 2
          nbin = nbout
          nbout = datasz_max(i, varid)
          CALL mathop(sopps(i, varid, io), nbin, buff_tmp, missing_val, &
               &      nbindex, nindex, scal(i, varid, io), nbout, buff_tmp2)

          nbin = nbout
          nbout = datasz_max(i, varid)
          CALL mathop(sopps(i, varid, io+1), nbin, buff_tmp2, missing_val, &
               &      nbindex, nindex, scal(i, varid, io+1), nbout, buff_tmp)
       ENDDO

       !   3.5 Zoom into the data

       CALL trans_buff &
            &      (zorig(i, varid, 1), zsize(i, varid, 1), &
            &       zorig(i, varid, 2), zsize(i, varid, 2), &
            &       zorig(i, varid, 3), zsize(i, varid, 3), &
            &       scsize(i, varid, 1), scsize(i, varid, 2), scsize(i, varid, 3), &
            &       buff_tmp, buff_tmp2_sz, buff_tmp2)

       !-- 5.0 Do the operations if needed. In the case of instantaneous
       !--     output we do not transfer to the buffer.

       ipt = point(pfileid, varid)

       IF (     (TRIM(tmp_opp) /= "inst") &
            &    .AND.(TRIM(tmp_opp) /= "once") ) THEN
          CALL moycum(tmp_opp, tsz, buffer(ipt:), &
               &       buff_tmp2, nb_opp(pfileid, varid))
       ENDIF

       last_opp(pfileid, varid) = pitau
       nb_opp(pfileid, varid) = nb_opp(pfileid, varid)+1

    ENDIF

    ! 6.0 Write to file if needed

    IF ( do_write ) THEN

       ncvarid = ncvar_ids(pfileid, varid)
       ncid = ncdf_ids(pfileid)

       !-- 6.1 Do the operations that are needed before writting

       IF (     (TRIM(tmp_opp) /= "inst") &
            &    .AND.(TRIM(tmp_opp) /= "once") ) THEN
          rtime = (rtime+last_wrt(pfileid, varid)*deltat(pfileid))/2.0
       ENDIF

       !-- 6.2 Add a value to the time axis of this variable if needed

       IF (     (TRIM(tmp_opp) /= "l_max") &
            &    .AND.(TRIM(tmp_opp) /= "l_min") &
            &    .AND.(TRIM(tmp_opp) /= "once") ) THEN

          itax = var_axid(pfileid, varid)
          itime = nb_wrt(pfileid, varid)+1

          IF (tax_last(pfileid, itax) < itime) THEN
             iret = NF90_PUT_VAR (ncid, tdimid(pfileid, itax), (/ rtime /), &
                  &                            start=(/ itime /), count=(/ 1 /))
             tax_last(pfileid, itax) = itime
          ENDIF
       ELSE
          itime=1
       ENDIF

       !-- 6.3 Write the data. Only in the case of instantaneous output
       !       we do not write the buffer.

       IF (scsize(pfileid, varid, 3) == 1) THEN
          IF (regular(pfileid)) THEN
             corner(1:4) = (/ 1, 1, itime, 0 /)
             edges(1:4) = (/ zsize(pfileid, varid, 1), &
                  &                      zsize(pfileid, varid, 2), &
                  &                       1, 0 /)
          ELSE
             corner(1:4) = (/ 1, itime, 0, 0 /)
             edges(1:4) = (/ zsize(pfileid, varid, 1), 1, 0, 0 /)
          ENDIF
       ELSE
          IF ( regular(pfileid) ) THEN
             corner(1:4) = (/ 1, 1, 1, itime /)
             edges(1:4) = (/ zsize(pfileid, varid, 1), &
                  &                      zsize(pfileid, varid, 2), &
                  &                      zsize(pfileid, varid, 3), 1 /)
          ELSE
             corner(1:4) = (/ 1, 1, itime, 0 /)
             edges(1:4) = (/ zsize(pfileid, varid, 1), &
                  &                      zsize(pfileid, varid, 3), 1, 0 /)
          ENDIF
       ENDIF

       ipt = point(pfileid, varid)

       IF (     (TRIM(tmp_opp) /= "inst") &
            &      .AND.(TRIM(tmp_opp) /= "once") ) THEN
          iret = NF90_PUT_VAR (ncid, ncvarid, buffer(ipt:), &
               &                       start=corner(1:4), count=edges(1:4))
       ELSE
          iret = NF90_PUT_VAR (ncid, ncvarid, buff_tmp2, &
               &                       start=corner(1:4), count=edges(1:4))
       ENDIF

       last_wrt(pfileid, varid) = pitau
       nb_wrt(pfileid, varid) = nb_wrt(pfileid, varid)+1
       nb_opp(pfileid, varid) = 0
       !---
       !   After the write the file can be synchronized so that no data is
       !   lost in case of a crash. This feature gives up on the benefits of
       !   buffering and should only be used in debuging mode. A flag is
       !   needed here to switch to this mode.
       !---
       !   iret = NF90_SYNC (ncid)

    ENDIF
    !----------------------------
  END SUBROUTINE histwrite_real

  !*************************************************************

  SUBROUTINE histvar_seq (pfid, pvarname, pvid)

    !- This subroutine optimized the search for the variable in the table.
    !- In a first phase it will learn the succession of the variables
    !- called and then it will use the table to guess what comes next.
    !- It is the best solution to avoid lengthy searches through array
    !- vectors.

    !- ARGUMENTS :

    !- pfid  : id of the file on which we work
    !- pvarname : The name of the variable we are looking for
    !- pvid     : The var id we found

    USE stringop, ONLY: find_str
    USE errioipsl, ONLY : histerr

    INTEGER, INTENT(in)  :: pfid
    CHARACTER(LEN=*), INTENT(IN) :: pvarname
    INTEGER, INTENT(out) :: pvid

    LOGICAL, SAVE :: learning(nb_files_max)=.TRUE.
    INTEGER, SAVE :: overlap(nb_files_max) = -1
    INTEGER, SAVE :: varseq(nb_files_max, nb_var_max*3)
    INTEGER, SAVE :: varseq_len(nb_files_max) = 0
    INTEGER, SAVE :: varseq_pos(nb_files_max)
    INTEGER, SAVE :: varseq_err(nb_files_max) = 0
    INTEGER      :: nb, sp, nx, pos, ib
    CHARACTER(LEN=20), DIMENSION(nb_var_max) :: tab_str20
    CHARACTER(LEN=20) :: str20
    CHARACTER(LEN=70) :: str70
    INTEGER      :: tab_str20_length(nb_var_max)

    !---------------------------------------------------------------------
    nb = nb_var(pfid)

    IF (learning(pfid)) THEN

       !-- 1.0 We compute the length over which we are going
       !--     to check the overlap

       IF (overlap(pfid) <= 0) THEN
          IF (nb_var(pfid) > 6) THEN
             overlap(pfid) = nb_var(pfid)/3*2
          ELSE
             overlap(pfid) = nb_var(pfid)
          ENDIF
       ENDIF

       !-- 1.1 Find the position of this string

       str20 = pvarname
       tab_str20(1:nb) = name(pfid, 1:nb)
       tab_str20_length(1:nb) = name_length(pfid, 1:nb)

       CALL find_str (nb, tab_str20, tab_str20_length, str20, pos)

       IF (pos > 0) THEN
          pvid = pos
       ELSE
          CALL histerr (3, "histvar_seq", &
               &      'The name of the variable you gave has not been declared', &
               &      'You should use subroutine histdef for declaring variable', &
               &      TRIM(str20))
       ENDIF

       !-- 1.2 If we have not given up we store the position
       !--     in the sequence of calls

       IF ( varseq_err(pfid) .GE. 0 ) THEN
          sp = varseq_len(pfid)+1
          IF (sp <= nb_var_max*3) THEN
             varseq(pfid, sp) = pvid
             varseq_len(pfid) = sp
          ELSE
             CALL histerr (2, "histvar_seq", &
                  &       'The learning process has failed and we give up. '// &
                  &       'Either you sequence is', &
                  &       'too complex or I am too dumb. '// &
                  &       'This will only affect the efficiency', &
                  &       'of your code. Thus if you wish to save time'// &
                  &       ' contact the IOIPSL team. ')
             WRITE(*, *) 'The sequence we have found up to now :'
             WRITE(*, *) varseq(pfid, 1:sp-1)
             varseq_err(pfid) = -1
          ENDIF

          !---- 1.3 Check if we have found the right overlap

          IF (varseq_len(pfid) .GE. overlap(pfid)*2) THEN

             !------ We skip a few variables if needed as they could come
             !------ from the initialisation of the model.

             DO ib = 0, sp-overlap(pfid)*2
                IF ( learning(pfid) .AND.&
                     & SUM(ABS(varseq(pfid, ib+1:ib+overlap(pfid)) -&
                     & varseq(pfid, sp-overlap(pfid)+1:sp))) == 0 ) THEN
                   learning(pfid) = .FALSE.
                   varseq_len(pfid) = sp-overlap(pfid)-ib
                   varseq_pos(pfid) = overlap(pfid)+ib
                   varseq(pfid, 1:varseq_len(pfid)) = &
                        &            varseq(pfid, ib+1:ib+varseq_len(pfid))
                ENDIF
             ENDDO
          ENDIF
       ENDIF
    ELSE

       !-- 2.0 Now we know how the calls to histwrite are sequenced
       !--     and we can get a guess at the var ID

       nx = varseq_pos(pfid)+1
       IF (nx > varseq_len(pfid)) nx = 1

       pvid = varseq(pfid, nx)

       IF (    (INDEX(name(pfid, pvid), pvarname) <= 0)         &
            &    .OR.(name_length(pfid, pvid) /= len_trim(pvarname)) ) THEN
          str20 = pvarname
          tab_str20(1:nb) = name(pfid, 1:nb)
          tab_str20_length(1:nb) = name_length(pfid, 1:nb)
          CALL find_str (nb, tab_str20, tab_str20_length, str20, pos)
          IF (pos > 0) THEN
             pvid = pos
          ELSE
             CALL histerr(3, "histvar_seq", &
                  &  'The name of the variable you gave has not been declared', &
                  &  'You should use subroutine histdef for declaring variable', str20)
          ENDIF
          varseq_err(pfid) = varseq_err(pfid)+1
       ELSE

          !---- We only keep the new position if we have found the variable
          !---- this way. This way an out of sequence call to histwrite does
          !---- not defeat the process.

          varseq_pos(pfid) = nx
       ENDIF

       IF (varseq_err(pfid) .GE. 10) THEN
          WRITE(str70, '("for file ", I3)') pfid
          CALL histerr(2, "histvar_seq", &
               &  'There were 10 errors in the learned sequence of variables', &
               &  str70, 'This looks like a bug, please report it.')
          varseq_err(pfid) = 0
       ENDIF
    ENDIF

  END SUBROUTINE histvar_seq

END MODULE histwrite_m
