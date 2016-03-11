MODULE histwrite_m

  ! From histcom.f90, version 2.1 2004/04/21 09:27:10

  USE errioipsl, ONLY: histerr
  use histbeg_totreg_m, ONLY: nb_files, date0, deltat
  USE histcom_var, ONLY: datasz_in, freq_opp, freq_wrt, fuchnbout, last_opp, &
       last_opp_chk, last_wrt, last_wrt_chk, missing_val, nb_files_max, &
       nb_var_max, nbopp, scal, scsize, sopps, topp
  use histvar_seq_m, only: histvar_seq
  use histwrite_real_m, only: histwrite_real
  use isittime_m, only: isittime
  USE mathop_m, ONLY: mathop

  implicit none

  INTEGER, SAVE:: datasz_max(nb_files_max, nb_var_max) = -1

  INTERFACE histwrite
     ! The "histwrite" procedures give the data to the input-output system.
     ! They trigger the operations to be performed and the writing to
     ! the file if needed.

     ! We test the work to be done at this time here so that at a
     ! later stage we can call different operations and write subroutines
     ! for the REAL and INTEGER interfaces.

     ! INTEGER, INTENT(IN):: fileid
     ! The ID of the file on which this variable is to be written.
     ! The variable should have been defined in this file before.

     ! CHARACTER(LEN=*), INTENT(IN):: varname
     ! short name of the variable

     ! INTEGER, INTENT(IN):: itau
     ! current timestep

     ! REAL, INTENT(IN):: pdata(:) or (:, :) or (:, :, :)
     ! values of the variable

     ! The difference between the procedures is the rank of "pdata".

     MODULE PROCEDURE histwrite_r1d, histwrite_r2d, histwrite_r3d
  END INTERFACE histwrite

  PRIVATE
  public histwrite

CONTAINS

  SUBROUTINE histwrite_r1d(fileid, varname, itau, pdata)

    INTEGER, INTENT(IN):: fileid, itau
    CHARACTER(LEN=*), INTENT(IN):: varname
    REAL, INTENT(IN):: pdata(:)

    ! Variables local to the procedure:
    integer nbindex, nindex(size(pdata))
    LOGICAL:: do_oper, do_write, largebuf
    INTEGER:: varid, io, nbpt_in, nbpt_out
    REAL, ALLOCATABLE, SAVE:: buff_tmp(:)
    INTEGER, SAVE:: buff_tmp_sz
    CHARACTER(LEN=7):: tmp_opp

    !--------------------------------------------------------------------

    nbindex = size(nindex)
    nindex = 0

    ! 1.0 Try to catch errors like specifying the wrong file ID.

    IF ((fileid < 1) .OR. (fileid > nb_files)) THEN
       CALL histerr(3, "histwrite", &
            'Illegal file ID in the histwrite of variable', varname, ' ')
    ENDIF

    ! 1.1 Find the id of the variable to be written and the real time

    CALL histvar_seq(fileid, varname, varid)

    ! 2.0 do nothing for never operation

    tmp_opp = topp(fileid, varid)

    IF (TRIM(tmp_opp) == "never") THEN
       last_opp_chk(fileid, varid) = -99
       last_wrt_chk(fileid, varid) = -99
    ENDIF

    ! 3.0 We check if we need to do an operation

    IF (last_opp_chk(fileid, varid) == itau) THEN
       CALL histerr(3, "histwrite", &
            'This variable as already been analysed at the present', &
            'time step', ' ')
    ENDIF

    CALL isittime(itau, date0(fileid), deltat(fileid), &
         freq_opp(fileid, varid), last_opp(fileid, varid), &
         last_opp_chk(fileid, varid), do_oper)

    ! 4.0 We check if we need to write the data

    IF (last_wrt_chk(fileid, varid) == itau) THEN
       CALL histerr(3, "histwrite", &
            'This variable as already been written for the present', &
            'time step', ' ')
    ENDIF

    CALL isittime(itau, date0(fileid), deltat(fileid), &
         freq_wrt(fileid, varid), last_wrt(fileid, varid), &
         last_wrt_chk(fileid, varid), do_write)

    ! 5.0 histwrite called

    IF (do_oper .OR. do_write) THEN
       ! 5.1 Get the sizes of the data we will handle

       IF (datasz_in(fileid, varid, 1) <= 0) THEN
          ! There is the risk here that the user has over-sized the array.
          ! But how can we catch this ?
          ! In the worst case we will do impossible operations
          ! on part of the data !
          datasz_in(fileid, varid, 1) = SIZE(pdata)
          datasz_in(fileid, varid, 2) = -1
          datasz_in(fileid, varid, 3) = -1
       ENDIF

       ! 5.2 The maximum size of the data will give the size of the buffer

       IF (datasz_max(fileid, varid) <= 0) THEN
          largebuf = .FALSE.
          DO io=1, nbopp(fileid, varid)
             IF (INDEX(fuchnbout, sopps(fileid, varid, io)) > 0) THEN
                largebuf = .TRUE.
             ENDIF
          ENDDO
          IF (largebuf) THEN
             datasz_max(fileid, varid) = scsize(fileid, varid, 1) &
                  * scsize(fileid, varid, 2) *scsize(fileid, varid, 3)
          ELSE
             datasz_max(fileid, varid) = datasz_in(fileid, varid, 1)
          ENDIF
       ENDIF

       IF (.NOT.ALLOCATED(buff_tmp)) THEN
          ALLOCATE(buff_tmp(datasz_max(fileid, varid)))
          buff_tmp_sz = datasz_max(fileid, varid)
       ELSE IF (datasz_max(fileid, varid) > buff_tmp_sz) THEN
          DEALLOCATE(buff_tmp)
          ALLOCATE(buff_tmp(datasz_max(fileid, varid)))
          buff_tmp_sz = datasz_max(fileid, varid)
       ENDIF

       ! We have to do the first operation anyway. Thus we do it here
       ! and change the ranke of the data at the same time. This
       ! should speed up things.

       nbpt_in = datasz_in(fileid, varid, 1)
       nbpt_out = datasz_max(fileid, varid)
       CALL mathop(sopps(fileid, varid, 1), nbpt_in, pdata, missing_val, &
            nbindex, nindex, scal(fileid, varid, 1), nbpt_out, buff_tmp)
       CALL histwrite_real(datasz_max, fileid, varid, itau, nbpt_out, &
            buff_tmp, nbindex, nindex, do_oper, do_write)
    ENDIF

    ! 6.0 Manage time steps

    IF ((TRIM(tmp_opp) /= "once") .AND. (TRIM(tmp_opp) /= "never")) THEN
       last_opp_chk(fileid, varid) = itau
       last_wrt_chk(fileid, varid) = itau
    ELSE
       last_opp_chk(fileid, varid) = -99
       last_wrt_chk(fileid, varid) = -99
    ENDIF

  END SUBROUTINE histwrite_r1d

  !************************************************************************

  SUBROUTINE histwrite_r2d (fileid, varname, itau, pdata)

    INTEGER, INTENT(IN):: fileid, itau
    REAL, INTENT(IN):: pdata(:, :)
    CHARACTER(LEN=*), INTENT(IN):: varname

    integer nbindex, nindex(size(pdata))
    LOGICAL:: do_oper, do_write, largebuf
    INTEGER:: varid, io, nbpt_in(1:2), nbpt_out
    REAL, ALLOCATABLE, SAVE:: buff_tmp(:)
    INTEGER, SAVE:: buff_tmp_sz
    CHARACTER(LEN=7):: tmp_opp

    !--------------------------------------------------------------------

    nbindex = size(nindex)
    nindex = 0

    ! 1.0 Try to catch errors like specifying the wrong file ID.

    IF ( (fileid < 1).OR.(fileid > nb_files) ) THEN
       CALL histerr (3, "histwrite", &
            'Illegal file ID in the histwrite of variable', varname, ' ')
    ENDIF

    ! 1.1 Find the id of the variable to be written and the real time

    CALL histvar_seq (fileid, varname, varid)

    ! 2.0 do nothing for never operation

    tmp_opp = topp(fileid, varid)

    IF (TRIM(tmp_opp) == "never") THEN
       last_opp_chk(fileid, varid) = -99
       last_wrt_chk(fileid, varid) = -99
    ENDIF

    ! 3.0 We check if we need to do an operation

    IF (last_opp_chk(fileid, varid) == itau) CALL histerr (3, "histwrite", &
         'This variable as already been analysed at the present', &
         'time step', ' ')

    CALL isittime(itau, date0(fileid), deltat(fileid), &
         freq_opp(fileid, varid), last_opp(fileid, varid), &
         last_opp_chk(fileid, varid), do_oper)

    ! 4.0 We check if we need to write the data

    IF (last_wrt_chk(fileid, varid) == itau) THEN
       CALL histerr (3, "histwrite", &
            'This variable as already been written for the present', &
            'time step', ' ')
    ENDIF

    CALL isittime &
         (itau, date0(fileid), deltat(fileid), freq_wrt(fileid, varid), &
         last_wrt(fileid, varid), last_wrt_chk(fileid, varid), do_write)

    ! 5.0 histwrite called

    IF (do_oper.OR.do_write) THEN

       !- 5.1 Get the sizes of the data we will handle

       IF (datasz_in(fileid, varid, 1) <= 0) THEN
          !--- There is the risk here that the user has over-sized the array.
          !--- But how can we catch this ?
          !--- In the worst case we will do impossible operations
          !--- on part of the data !
          datasz_in(fileid, varid, 1) = SIZE(pdata, DIM=1)
          datasz_in(fileid, varid, 2) = SIZE(pdata, DIM=2)
          datasz_in(fileid, varid, 3) = -1
       ENDIF

       !- 5.2 The maximum size of the data will give the size of the buffer

       IF (datasz_max(fileid, varid) <= 0) THEN
          largebuf = .FALSE.
          DO io=1, nbopp(fileid, varid)
             IF (INDEX(fuchnbout, sopps(fileid, varid, io)) > 0) THEN
                largebuf = .TRUE.
             ENDIF
          ENDDO
          IF (largebuf) THEN
             datasz_max(fileid, varid) = &
                  scsize(fileid, varid, 1) &
                  *scsize(fileid, varid, 2) &
                  *scsize(fileid, varid, 3)
          ELSE
             datasz_max(fileid, varid) = &
                  datasz_in(fileid, varid, 1) &
                  *datasz_in(fileid, varid, 2)
          ENDIF
       ENDIF

       IF (.NOT.ALLOCATED(buff_tmp)) THEN
          ALLOCATE (buff_tmp(datasz_max(fileid, varid)))
          buff_tmp_sz = datasz_max(fileid, varid)
       ELSE IF (datasz_max(fileid, varid) > buff_tmp_sz) THEN
          DEALLOCATE (buff_tmp)
          ALLOCATE (buff_tmp(datasz_max(fileid, varid)))
          buff_tmp_sz = datasz_max(fileid, varid)
       ENDIF

       !- We have to do the first operation anyway.
       !- Thus we do it here and change the ranke
       !- of the data at the same time. This should speed up things.

       nbpt_in(1:2) = datasz_in(fileid, varid, 1:2)
       nbpt_out = datasz_max(fileid, varid)
       CALL mathop (sopps(fileid, varid, 1), nbpt_in, pdata, &
            missing_val, nbindex, nindex, &
            scal(fileid, varid, 1), nbpt_out, buff_tmp)
       CALL histwrite_real (datasz_max, fileid, varid, itau, nbpt_out, &
            buff_tmp, nbindex, nindex, do_oper, do_write)
    ENDIF

    ! 6.0 Manage time steps

    IF ((TRIM(tmp_opp) /= "once").AND.(TRIM(tmp_opp) /= "never")) THEN
       last_opp_chk(fileid, varid) = itau
       last_wrt_chk(fileid, varid) = itau
    ELSE
       last_opp_chk(fileid, varid) = -99
       last_wrt_chk(fileid, varid) = -99
    ENDIF

  END SUBROUTINE histwrite_r2d

  !************************************************************************

  SUBROUTINE histwrite_r3d (fileid, varname, itau, pdata)

    INTEGER, INTENT(IN):: fileid, itau
    REAL, DIMENSION(:, :, :), INTENT(IN):: pdata
    CHARACTER(LEN=*), INTENT(IN):: varname

    integer nbindex, nindex(size(pdata))
    LOGICAL:: do_oper, do_write, largebuf
    INTEGER:: varid, io, nbpt_in(1:3), nbpt_out
    REAL, ALLOCATABLE, SAVE:: buff_tmp(:)
    INTEGER, SAVE:: buff_tmp_sz
    CHARACTER(LEN=7):: tmp_opp

    !--------------------------------------------------------------------

    nbindex = size(nindex)
    nindex = 0

    ! 1.0 Try to catch errors like specifying the wrong file ID.
    !     Thanks Marine for showing us what errors users can make !

    IF ( (fileid < 1).OR.(fileid > nb_files) ) THEN
       CALL histerr (3, "histwrite", &
            'Illegal file ID in the histwrite of variable', varname, ' ')
    ENDIF

    ! 1.1 Find the id of the variable to be written and the real time

    CALL histvar_seq (fileid, varname, varid)

    ! 2.0 do nothing for never operation

    tmp_opp = topp(fileid, varid)

    IF (TRIM(tmp_opp) == "never") THEN
       last_opp_chk(fileid, varid) = -99
       last_wrt_chk(fileid, varid) = -99
    ENDIF

    ! 3.0 We check if we need to do an operation

    IF (last_opp_chk(fileid, varid) == itau) THEN
       CALL histerr (3, "histwrite", &
            'This variable as already been analysed at the present', &
            'time step', ' ')
    ENDIF

    CALL isittime &
         (itau, date0(fileid), deltat(fileid), freq_opp(fileid, varid), &
         last_opp(fileid, varid), last_opp_chk(fileid, varid), do_oper)

    ! 4.0 We check if we need to write the data

    IF (last_wrt_chk(fileid, varid) == itau) THEN
       CALL histerr (3, "histwrite", &
            'This variable as already been written for the present', &
            'time step', ' ')
    ENDIF

    CALL isittime &
         (itau, date0(fileid), deltat(fileid), freq_wrt(fileid, varid), &
         last_wrt(fileid, varid), last_wrt_chk(fileid, varid), do_write)

    ! 5.0 histwrite called

    IF (do_oper.OR.do_write) THEN

       !- 5.1 Get the sizes of the data we will handle

       IF (datasz_in(fileid, varid, 1) <= 0) THEN
          !--- There is the risk here that the user has over-sized the array.
          !--- But how can we catch this ?
          !--- In the worst case we will do impossible operations
          !--- on part of the data !
          datasz_in(fileid, varid, 1) = SIZE(pdata, DIM=1)
          datasz_in(fileid, varid, 2) = SIZE(pdata, DIM=2)
          datasz_in(fileid, varid, 3) = SIZE(pdata, DIM=3)
       ENDIF

       !- 5.2 The maximum size of the data will give the size of the buffer

       IF (datasz_max(fileid, varid) <= 0) THEN
          largebuf = .FALSE.
          DO io =1, nbopp(fileid, varid)
             IF (INDEX(fuchnbout, sopps(fileid, varid, io)) > 0) THEN
                largebuf = .TRUE.
             ENDIF
          ENDDO
          IF (largebuf) THEN
             datasz_max(fileid, varid) = &
                  scsize(fileid, varid, 1) &
                  *scsize(fileid, varid, 2) &
                  *scsize(fileid, varid, 3)
          ELSE
             datasz_max(fileid, varid) = &
                  datasz_in(fileid, varid, 1) &
                  *datasz_in(fileid, varid, 2) &
                  *datasz_in(fileid, varid, 3)
          ENDIF
       ENDIF

       IF (.NOT.ALLOCATED(buff_tmp)) THEN
          ALLOCATE (buff_tmp(datasz_max(fileid, varid)))
          buff_tmp_sz = datasz_max(fileid, varid)
       ELSE IF (datasz_max(fileid, varid) > buff_tmp_sz) THEN
          DEALLOCATE (buff_tmp)
          ALLOCATE (buff_tmp(datasz_max(fileid, varid)))
          buff_tmp_sz = datasz_max(fileid, varid)
       ENDIF

       !- We have to do the first operation anyway.
       !- Thus we do it here and change the ranke
       !- of the data at the same time. This should speed up things.

       nbpt_in(1:3) = datasz_in(fileid, varid, 1:3)
       nbpt_out = datasz_max(fileid, varid)
       CALL mathop (sopps(fileid, varid, 1), nbpt_in, pdata, &
            missing_val, nbindex, nindex, &
            scal(fileid, varid, 1), nbpt_out, buff_tmp)
       CALL histwrite_real(datasz_max, fileid, varid, itau, nbpt_out, &
            buff_tmp, nbindex, nindex, do_oper, do_write)
    ENDIF

    ! 6.0 Manage time steps

    IF ((TRIM(tmp_opp) /= "once").AND.(TRIM(tmp_opp) /= "never")) THEN
       last_opp_chk(fileid, varid) = itau
       last_wrt_chk(fileid, varid) = itau
    ELSE
       last_opp_chk(fileid, varid) = -99
       last_wrt_chk(fileid, varid) = -99
    ENDIF

  END SUBROUTINE histwrite_r3d

END MODULE histwrite_m
