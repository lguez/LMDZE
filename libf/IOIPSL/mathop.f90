MODULE mathop_m

  ! From mathelp.f90, version 2.0 2004/04/05 14:47:50

  implicit none

  PRIVATE
  PUBLIC :: mathop

  INTERFACE mathop
     MODULE PROCEDURE mathop_r11, mathop_r21, mathop_r31
  END INTERFACE

  !- Variables used to detect and identify the operations
  CHARACTER(LEN=120), SAVE :: indexfu = &
       'gather, scatter, fill, coll, undef, only'

CONTAINS

  SUBROUTINE mathop_r11 &
       &  (fun, nb, work_in, miss_val, nb_index, nindex, scal, nb_max, work_out)
    !- This subroutines gives an interface to the various operation
    !- which are allowed. The interface is general enough to allow its use
    !- for other cases.

    !- INPUT

    !- fun      : function to be applied to the vector of data
    !- nb       : Length of input vector
    !- work_in  : Input vector of data (REAL)
    !- miss_val : The value of the missing data flag (it has to be a
    !-            maximum value, in f90 : huge( a real ))
    !- nb_index : Length of index vector
    !- nindex   : Vector of indices
    !- scal     : A scalar value for vector/scalar operations
    !- nb_max   : maximum length of output vector

    !- OUTPUT

    !- nb_max   : Actual length of output variable
    !- work_out : Output vector after the operation was applied
    USE errioipsl, ONLY : histerr

    CHARACTER(LEN=7) :: fun
    INTEGER :: nb, nb_max, nb_index
    INTEGER :: nindex(nb_index)
    REAL :: work_in(nb), scal, miss_val
    REAL :: work_out(nb_max)

    INTEGER :: ierr

    INTRINSIC SIN, COS, TAN, ASIN, ACOS, ATAN, EXP, ALOG, SQRT, ABS
    !---------------------------------------------------------------------
    ierr = 0

    IF (scal >= miss_val-1.) THEN
       IF (INDEX(indexfu, fun(1:LEN_TRIM(fun))) == 0) THEN
          SELECT CASE (fun)
          CASE('sin')
             ierr = ma_sin_r11(nb, work_in, nb_max, work_out)
          CASE('cos')
             ierr = ma_cos_r11(nb, work_in, nb_max, work_out)
          CASE('tan')
             ierr = ma_tan_r11(nb, work_in, nb_max, work_out)
          CASE('asin')
             ierr = ma_asin_r11(nb, work_in, nb_max, work_out)
          CASE('acos')
             ierr = ma_acos_r11(nb, work_in, nb_max, work_out)
          CASE('atan')
             ierr = ma_atan_r11(nb, work_in, nb_max, work_out)
          CASE('exp')
             ierr = ma_exp_r11(nb, work_in, nb_max, work_out)
          CASE('log')
             ierr = ma_alog_r11(nb, work_in, nb_max, work_out)
          CASE('sqrt')
             ierr = ma_sqrt_r11(nb, work_in, nb_max, work_out)
          CASE('chs')
             ierr = ma_chs_r11(nb, work_in, nb_max, work_out)
          CASE('abs')
             ierr = ma_abs_r11(nb, work_in, nb_max, work_out)
          CASE('cels')
             ierr = ma_cels_r11(nb, work_in, nb_max, work_out)
          CASE('kelv')
             ierr = ma_kelv_r11(nb, work_in, nb_max, work_out)
          CASE('deg')
             ierr = ma_deg_r11(nb, work_in, nb_max, work_out)
          CASE('rad')
             ierr = ma_rad_r11(nb, work_in, nb_max, work_out)
          CASE('ident')
             ierr = ma_ident_r11(nb, work_in, nb_max, work_out)
          CASE DEFAULT
             CALL histerr(3, "mathop", &
                  &        'scalar variable undefined and no indexing', &
                  &        'but still unknown function', fun)
          END SELECT
          IF (ierr > 0) THEN
             CALL histerr(3, "mathop", &
                  &        'Error while executing a simple function', fun, ' ')
          ENDIF
       ELSE
          SELECT CASE (fun)
          CASE('gather')
             ierr = ma_fugath_r11(nb, work_in, nb_index, nindex, &
                  &                            miss_val, nb_max, work_out)
          CASE('scatter')
             IF (nb_index > nb) THEN
                work_out(1:nb_max) = miss_val
                ierr=1
             ELSE
                ierr = ma_fuscat_r11(nb, work_in, nb_index, nindex, &
                     &                              miss_val, nb_max, work_out)
             ENDIF
          CASE('coll')
             ierr = ma_fucoll_r11(nb, work_in, nb_index, nindex, &
                  &                            miss_val, nb_max, work_out)
          CASE('fill')
             ierr = ma_fufill_r11(nb, work_in, nb_index, nindex, &
                  &                            miss_val, nb_max, work_out)
          CASE('undef')
             ierr = ma_fuundef_r11(nb, work_in, nb_index, nindex, &
                  &                            miss_val, nb_max, work_out)
          CASE('only')
             ierr = ma_fuonly_r11(nb, work_in, nb_index, nindex, &
                  &                            miss_val, nb_max, work_out)
          CASE DEFAULT
             CALL histerr(3, "mathop", &
                  &        'scalar variable undefined and indexing', &
                  &        'was requested but with unknown function', fun)
          END SELECT
          IF (ierr > 0) THEN
             CALL histerr(3, "mathop_r11", &
                  &        'Error while executing an indexing function', fun, ' ')
          ENDIF
       ENDIF
    ELSE
       SELECT CASE (fun)
       CASE('fumin')
          ierr = ma_fumin_r11(nb, work_in, scal, nb_max, work_out)
       CASE('fumax')
          ierr = ma_fumax_r11(nb, work_in, scal, nb_max, work_out)
       CASE('add')
          ierr = ma_add_r11(nb, work_in, scal, nb_max, work_out)
       CASE('subi')
          ierr = ma_subi_r11(nb, work_in, scal, nb_max, work_out)
       CASE('sub')
          ierr = ma_sub_r11(nb, work_in, scal, nb_max, work_out)
       CASE('mult')
          ierr = ma_mult_r11(nb, work_in, scal, nb_max, work_out)
       CASE('div')
          ierr = ma_div_r11(nb, work_in, scal, nb_max, work_out)
       CASE('divi')
          ierr = ma_divi_r11(nb, work_in, scal, nb_max, work_out)
       CASE('power')
          ierr = ma_power_r11(nb, work_in, scal, nb_max, work_out)
       CASE DEFAULT
          CALL histerr(3, "mathop", &
               &      'Unknown operation with a scalar', fun, ' ')
       END SELECT
       IF (ierr > 0) THEN
          CALL histerr(3, "mathop", &
               &      'Error while executing a scalar function', fun, ' ')
       ENDIF
    ENDIF
    !------------------------
  END SUBROUTINE mathop_r11
  !-
  !=== FUNCTIONS (only one argument)
  !-
  INTEGER FUNCTION ma_sin_r11(nb, x, nbo, y)
    INTEGER :: nb, nbo, i
    REAL :: x(nb), y(nbo)
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = SIN(x(i))
    ENDDO

    nbo = nb
    ma_sin_r11 = 0
    !----------------------
  END FUNCTION ma_sin_r11

  !************************************************

  INTEGER FUNCTION ma_cos_r11(nb, x, nbo, y)
    INTEGER :: nb, nbo, i
    REAL :: x(nb), y(nbo)
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = COS(x(i))
    ENDDO

    nbo = nb
    ma_cos_r11 = 0
    !----------------------
  END FUNCTION ma_cos_r11

  !************************************************

  INTEGER FUNCTION ma_tan_r11(nb, x, nbo, y)
    INTEGER :: nb, nbo, i
    REAL :: x(nb), y(nbo)
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = TAN(x(i))
    ENDDO

    nbo = nb
    ma_tan_r11 = 0
    !----------------------
  END FUNCTION ma_tan_r11

  !************************************************

  INTEGER FUNCTION ma_asin_r11(nb, x, nbo, y)
    INTEGER :: nb, nbo, i
    REAL :: x(nb), y(nbo)
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = ASIN(x(i))
    ENDDO

    nbo = nb
    ma_asin_r11 = 0
    !-----------------------
  END FUNCTION ma_asin_r11

  !************************************************

  INTEGER FUNCTION ma_acos_r11(nb, x, nbo, y)
    INTEGER :: nb, nbo, i
    REAL :: x(nb), y(nbo)
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = ACOS(x(i))
    ENDDO

    nbo = nb
    ma_acos_r11 = 0
    !-----------------------
  END FUNCTION ma_acos_r11

  !************************************************

  INTEGER FUNCTION ma_atan_r11(nb, x, nbo, y)
    INTEGER :: nb, nbo, i
    REAL :: x(nb), y(nbo)
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = ATAN(x(i))
    ENDDO

    nbo = nb
    ma_atan_r11 = 0
    !-----------------------
  END FUNCTION ma_atan_r11

  !************************************************

  INTEGER FUNCTION ma_exp_r11(nb, x, nbo, y)
    INTEGER :: nb, nbo, i
    REAL :: x(nb), y(nbo)
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = EXP(x(i))
    ENDDO

    nbo = nb
    ma_exp_r11 = 0
    !----------------------
  END FUNCTION ma_exp_r11

  !************************************************

  INTEGER FUNCTION ma_alog_r11(nb, x, nbo, y)
    INTEGER :: nb, nbo, i
    REAL :: x(nb), y(nbo)
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = alog(x(i))
    ENDDO

    nbo = nb
    ma_alog_r11 = 0
    !-----------------------
  END FUNCTION ma_alog_r11

  !************************************************

  INTEGER FUNCTION ma_sqrt_r11(nb, x, nbo, y)
    INTEGER :: nb, nbo, i
    REAL :: x(nb), y(nbo)
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = SQRT(x(i))
    ENDDO

    nbo = nb
    ma_sqrt_r11 = 0
    !-----------------------
  END FUNCTION ma_sqrt_r11

  !************************************************

  INTEGER FUNCTION ma_abs_r11(nb, x, nbo, y)
    INTEGER :: nb, nbo, i
    REAL :: x(nb), y(nbo)
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = ABS(x(i))
    ENDDO

    nbo = nb
    ma_abs_r11 = 0
    !----------------------
  END FUNCTION ma_abs_r11

  !************************************************

  INTEGER FUNCTION ma_chs_r11(nb, x, nbo, y)
    INTEGER :: nb, nbo, i
    REAL :: x(nb), y(nbo)
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = x(i)*(-1.)
    ENDDO

    nbo = nb
    ma_chs_r11 = 0
    !----------------------
  END FUNCTION ma_chs_r11

  !************************************************

  INTEGER FUNCTION ma_cels_r11(nb, x, nbo, y)
    INTEGER :: nb, nbo, i
    REAL :: x(nb), y(nbo)
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = x(i)-273.15
    ENDDO

    nbo = nb
    ma_cels_r11 = 0
    !-----------------------
  END FUNCTION ma_cels_r11

  !************************************************

  INTEGER FUNCTION ma_kelv_r11(nb, x, nbo, y)
    INTEGER :: nb, nbo, i
    REAL :: x(nb), y(nbo)
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = x(i)+273.15
    ENDDO

    nbo = nb
    ma_kelv_r11 = 0
    !-----------------------
  END FUNCTION ma_kelv_r11

  !************************************************

  INTEGER FUNCTION ma_deg_r11(nb, x, nbo, y)
    INTEGER :: nb, nbo, i
    REAL :: x(nb), y(nbo)
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = x(i)*57.29577951
    ENDDO

    nbo = nb
    ma_deg_r11 = 0
    !-----------------------
  END FUNCTION ma_deg_r11

  !************************************************

  INTEGER FUNCTION ma_rad_r11(nb, x, nbo, y)
    INTEGER :: nb, nbo, i
    REAL :: x(nb), y(nbo)
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = x(i)*0.01745329252
    ENDDO

    nbo = nb
    ma_rad_r11 = 0
    !----------------------
  END FUNCTION ma_rad_r11

  !************************************************

  INTEGER FUNCTION ma_ident_r11(nb, x, nbo, y)
    INTEGER :: nb, nbo, i
    REAL :: x(nb), y(nbo)
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = x(i)
    ENDDO

    nbo = nb
    ma_ident_r11 = 0
    !------------------------
  END FUNCTION ma_ident_r11
  !-
  !=== OPERATIONS (two argument)
  !-
  INTEGER FUNCTION ma_add_r11(nb, x, s, nbo, y)
    INTEGER :: nb, nbo
    REAL :: x(nb), s, y(nbo)

    INTEGER :: i
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = x(i)+s
    ENDDO

    nbo = nb
    ma_add_r11 = 0
    !-----------------------
  END FUNCTION ma_add_r11

  !************************************************

  INTEGER FUNCTION ma_sub_r11(nb, x, s, nbo, y)
    INTEGER :: nb, nbo
    REAL :: x(nb), s, y(nbo)

    INTEGER :: i
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = x(i)-s
    ENDDO

    nbo = nb
    ma_sub_r11 = 0
    !----------------------
  END FUNCTION ma_sub_r11

  !************************************************

  INTEGER FUNCTION ma_subi_r11(nb, x, s, nbo, y)
    INTEGER :: nb, nbo
    REAL :: x(nb), s, y(nbo)

    INTEGER :: i
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) =  s-x(i)
    ENDDO

    nbo = nb
    ma_subi_r11 = 0
    !-----------------------
  END FUNCTION ma_subi_r11

  !************************************************

  INTEGER FUNCTION ma_mult_r11(nb, x, s, nbo, y)
    INTEGER :: nb, nbo
    REAL :: x(nb), s, y(nbo)

    INTEGER :: i
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = x(i)*s
    ENDDO

    nbo = nb
    ma_mult_r11 = 0
    !-----------------------
  END FUNCTION ma_mult_r11

  !************************************************

  INTEGER FUNCTION ma_div_r11(nb, x, s, nbo, y)
    INTEGER :: nb, nbo
    REAL :: x(nb), s, y(nbo)

    INTEGER :: i
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = x(i)/s
    ENDDO

    nbo = nb
    ma_div_r11 = 0
    !-----------------------
  END FUNCTION ma_div_r11

  !************************************************

  INTEGER FUNCTION ma_divi_r11(nb, x, s, nbo, y)
    INTEGER :: nb, nbo
    REAL :: x(nb), s, y(nbo)

    INTEGER :: i
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = s/x(i)
    ENDDO

    nbo = nb
    ma_divi_r11 = 0
    !-----------------------
  END FUNCTION ma_divi_r11

  !************************************************

  INTEGER FUNCTION ma_power_r11(nb, x, s, nbo, y)
    INTEGER :: nb, nbo
    REAL :: x(nb), s, y(nbo)

    INTEGER :: i
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = x(i)**s
    ENDDO

    nbo = nb
    ma_power_r11 = 0
    !-----------------------
  END FUNCTION ma_power_r11

  !************************************************

  INTEGER FUNCTION ma_fumin_r11(nb, x, s, nbo, y)
    INTEGER :: nb, nbo
    REAL :: x(nb), s, y(nbo)

    INTEGER :: i
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = MIN(x(i), s)
    ENDDO

    nbo = nb
    ma_fumin_r11 = 0
    !------------------------
  END FUNCTION ma_fumin_r11

  !************************************************

  INTEGER FUNCTION ma_fumax_r11(nb, x, s, nbo, y)
    INTEGER :: nb, nbo
    REAL :: x(nb), s, y(nbo)

    INTEGER :: i
    !---------------------------------------------------------------------
    DO i=1, nb
       y(i) = MAX(x(i), s)
    ENDDO

    nbo = nb
    ma_fumax_r11 = 0
    !------------------------
  END FUNCTION ma_fumax_r11

  !************************************************

  INTEGER FUNCTION ma_fuscat_r11(nb, x, nbi, ind, miss_val, nbo, y)
    INTEGER :: nb, nbo, nbi
    INTEGER :: ind(nbi)
    REAL :: x(nb), miss_val, y(nbo)

    INTEGER :: i, ii, ipos
    !---------------------------------------------------------------------
    ma_fuscat_r11 = 0

    y(1:nbo) = miss_val

    IF (nbi <= nb) THEN
       ipos = 0
       DO i=1, nbi
          IF (ind(i) <= nbo .AND. ind(i) > 0) THEN
             ipos = ipos+1
             y(ind(i)) = x(ipos)
          ELSE
             IF (ind(i) > nbo) ma_fuscat_r11  = ma_fuscat_r11+1
          ENDIF
       ENDDO
       !-- Repeat the data if needed
       IF (MINVAL(ind) < 0) THEN
          DO i=1, nbi
             IF (ind(i) <= 0) THEN
                DO ii=1, ABS(ind(i))-1
                   IF (ind(i+1)+ii <= nbo) THEN
                      y(ind(i+1)+ii) = y(ind(i+1))
                   ELSE
                      ma_fuscat_r11  = ma_fuscat_r11+1
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDIF
    ELSE
       ma_fuscat_r11  = 1
    ENDIF
    !-------------------------
  END FUNCTION ma_fuscat_r11

  !************************************************

  INTEGER FUNCTION ma_fugath_r11(nb, x, nbi, ind, miss_val, nbo, y)
    INTEGER :: nb, nbo, nbi
    INTEGER :: ind(nbi)
    REAL :: x(nb), miss_val, y(nbo)

    INTEGER :: i, ipos
    !---------------------------------------------------------------------
    IF (nbi <= nbo) THEN
       ma_fugath_r11 = 0
       y(1:nbo) = miss_val
       ipos = 0
       DO i=1, nbi
          IF (ipos+1 <= nbo) THEN
             IF (ind(i) > 0) THEN
                ipos = ipos+1
                y(ipos) = x(ind(i))
             ENDIF
          ELSE
             IF (ipos+1 > nbo) ma_fugath_r11  = ma_fugath_r11+1
          ENDIF
       ENDDO
    ELSE
       ma_fugath_r11 = 1
    ENDIF

    nbo = ipos
    !-------------------------
  END FUNCTION ma_fugath_r11

  !************************************************

  INTEGER FUNCTION ma_fufill_r11(nb, x, nbi, ind, miss_val, nbo, y)
    INTEGER :: nb, nbo, nbi
    INTEGER :: ind(nbi)
    REAL :: x(nb), miss_val, y(nbo)

    INTEGER :: i, ii, ipos
    !---------------------------------------------------------------------
    ma_fufill_r11 = 0

    IF (nbi <= nb) THEN
       ipos = 0
       DO i=1, nbi
          IF (ind(i) <= nbo .AND. ind(i) > 0) THEN
             ipos = ipos+1
             y(ind(i)) = x(ipos)
          ELSE
             IF (ind(i) > nbo) ma_fufill_r11  = ma_fufill_r11+1
          ENDIF
       ENDDO
       !-- Repeat the data if needed
       IF (MINVAL(ind) < 0) THEN
          DO i=1, nbi
             IF (ind(i) <= 0) THEN
                DO ii=1, ABS(ind(i))-1
                   IF (ind(i+1)+ii <= nbo) THEN
                      y(ind(i+1)+ii) = y(ind(i+1))
                   ELSE
                      ma_fufill_r11  = ma_fufill_r11+1
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDIF
    ELSE
       ma_fufill_r11  = 1
    ENDIF
    !-------------------------
  END FUNCTION ma_fufill_r11

  !************************************************

  INTEGER FUNCTION ma_fucoll_r11(nb, x, nbi, ind, miss_val, nbo, y)
    INTEGER :: nb, nbo, nbi
    INTEGER :: ind(nbi)
    REAL :: x(nb), miss_val, y(nbo)

    INTEGER :: i, ipos
    !---------------------------------------------------------------------
    IF (nbi <= nbo) THEN
       ma_fucoll_r11 = 0
       ipos = 0
       DO i=1, nbi
          IF (ipos+1 <= nbo) THEN
             IF (ind(i) > 0) THEN
                ipos = ipos+1
                y(ipos) = x(ind(i))
             ENDIF
          ELSE
             IF (ipos+1 > nbo) ma_fucoll_r11  = ma_fucoll_r11+1
          ENDIF
       ENDDO
    ELSE
       ma_fucoll_r11 = 1
    ENDIF

    nbo = ipos
    !-------------------------
  END FUNCTION ma_fucoll_r11

  !************************************************

  INTEGER FUNCTION ma_fuundef_r11(nb, x, nbi, ind, miss_val, nbo, y)
    INTEGER :: nb, nbo, nbi
    INTEGER :: ind(nbi)
    REAL :: x(nb), miss_val, y(nbo)

    INTEGER :: i
    !---------------------------------------------------------------------
    IF (nbi <= nbo .AND. nbo == nb) THEN
       ma_fuundef_r11 = 0
       DO i=1, nbo
          y(i) = x(i)
       ENDDO
       DO i=1, nbi
          IF (ind(i) <= nbo .AND. ind(i) > 0) THEN
             y(ind(i)) =  miss_val
          ELSE
             IF (ind(i) > nbo) ma_fuundef_r11  = ma_fuundef_r11+1
          ENDIF
       ENDDO
    ELSE
       ma_fuundef_r11 = 1
    ENDIF
    !--------------------------
  END FUNCTION ma_fuundef_r11

  !************************************************

  INTEGER FUNCTION ma_fuonly_r11(nb, x, nbi, ind, miss_val, nbo, y)
    INTEGER :: nb, nbo, nbi
    INTEGER :: ind(nbi)
    REAL :: x(nb), miss_val, y(nbo)

    INTEGER :: i
    !---------------------------------------------------------------------
    IF (     (nbi <= nbo).AND.(nbo == nb) &
         &    .AND.ALL(ind(1:nbi) <= nbo) ) THEN
       ma_fuonly_r11 = 0
       y(1:nbo) = miss_val
       DO i=1, nbi
          IF (ind(i) > 0) THEN
             y(ind(i)) =  x(ind(i))
          ENDIF
       ENDDO
    ELSE
       ma_fuonly_r11 = 1
    ENDIF
    !-------------------------
  END FUNCTION ma_fuonly_r11

  !************************************************

  !************************************************

  SUBROUTINE mathop_r21 &
       &  (fun, nb, work_in, miss_val, nb_index, nindex, scal, nb_max, work_out)
    !- This subroutines gives an interface to the various operations
    !- which are allowed. The interface is general enough to allow its use
    !- for other cases.

    !- INPUT

    !- fun      : function to be applied to the vector of data
    !- nb       : Length of input vector
    !- work_in  : Input vector of data (REAL)
    !- miss_val : The value of the missing data flag (it has to be a
    !-            maximum value, in f90 : huge( a real ))
    !- nb_index : Length of index vector
    !- nindex   : Vector of indices
    !- scal     : A scalar value for vector/scalar operations
    !- nb_max   : maximum length of output vector

    !- OUTPUT

    !- nb_max   : Actual length of output variable
    !- work_out : Output vector after the operation was applied
    USE errioipsl, ONLY : histerr

    CHARACTER(LEN=7) :: fun
    INTEGER :: nb(2), nb_max, nb_index
    INTEGER :: nindex(nb_index)
    REAL :: work_in(nb(1), nb(2)), scal, miss_val
    REAL :: work_out(nb_max)

    INTEGER :: ierr

    INTRINSIC SIN, COS, TAN, ASIN, ACOS, ATAN, EXP, ALOG, SQRT, ABS
    !---------------------------------------------------------------------
    ierr = 0

    IF (scal >= miss_val-1.) THEN
       IF (INDEX(indexfu, fun(1:LEN_TRIM(fun))) == 0) THEN
          SELECT CASE (fun)
          CASE('sin')
             ierr = ma_sin_r21(nb, work_in, nb_max, work_out)
          CASE('cos')
             ierr = ma_cos_r21(nb, work_in, nb_max, work_out)
          CASE('tan')
             ierr = ma_tan_r21(nb, work_in, nb_max, work_out)
          CASE('asin')
             ierr = ma_asin_r21(nb, work_in, nb_max, work_out)
          CASE('acos')
             ierr = ma_acos_r21(nb, work_in, nb_max, work_out)
          CASE('atan')
             ierr = ma_atan_r21(nb, work_in, nb_max, work_out)
          CASE('exp')
             ierr = ma_exp_r21(nb, work_in, nb_max, work_out)
          CASE('log')
             ierr = ma_alog_r21(nb, work_in, nb_max, work_out)
          CASE('sqrt')
             ierr = ma_sqrt_r21(nb, work_in, nb_max, work_out)
          CASE('chs')
             ierr = ma_chs_r21(nb, work_in, nb_max, work_out)
          CASE('abs')
             ierr = ma_abs_r21(nb, work_in, nb_max, work_out)
          CASE('cels')
             ierr = ma_cels_r21(nb, work_in, nb_max, work_out)
          CASE('kelv')
             ierr = ma_kelv_r21(nb, work_in, nb_max, work_out)
          CASE('deg')
             ierr = ma_deg_r21(nb, work_in, nb_max, work_out)
          CASE('rad')
             ierr = ma_rad_r21(nb, work_in, nb_max, work_out)
          CASE('ident')
             ierr = ma_ident_r21(nb, work_in, nb_max, work_out)
          CASE DEFAULT
             CALL histerr(3, "mathop", &
                  &        'scalar variable undefined and no indexing', &
                  &        'but still unknown function', fun)
          END SELECT
          IF (ierr > 0) THEN
             CALL histerr(3, "mathop", &
                  &        'Error while executing a simple function', fun, ' ')
          ENDIF
       ELSE
          SELECT CASE (fun)
          CASE('gather')
             ierr = ma_fugath_r21(nb, work_in, nb_index, nindex, &
                  &                            miss_val, nb_max, work_out)
          CASE('scatter')
             IF (nb_index > (nb(1)*nb(2)) ) THEN
                work_out(1:nb_max) = miss_val
                ierr=1
             ELSE
                ierr = ma_fuscat_r21(nb, work_in, nb_index, nindex, &
                     &                             miss_val, nb_max, work_out)
             ENDIF
          CASE('coll')
             ierr = ma_fucoll_r21(nb, work_in, nb_index, nindex, &
                  &                           miss_val, nb_max, work_out)
          CASE('fill')
             ierr = ma_fufill_r21(nb, work_in, nb_index, nindex, &
                  &                           miss_val, nb_max, work_out)
          CASE('undef')
             ierr = ma_fuundef_r21(nb, work_in, nb_index, nindex, &
                  &                           miss_val, nb_max, work_out)
          CASE('only')
             ierr = ma_fuonly_r21(nb, work_in, nb_index, nindex, &
                  &                           miss_val, nb_max, work_out)
          CASE DEFAULT
             CALL histerr(3, "mathop", &
                  &        'scalar variable undefined and indexing', &
                  &        'was requested but with unknown function', fun)
          END SELECT
          IF (ierr > 0) THEN
             CALL histerr(3, "mathop_r21", &
                  &        'Error while executing an indexing function', fun, ' ')
          ENDIF
       ENDIF
    ELSE
       SELECT CASE (fun)
       CASE('fumin')
          ierr = ma_fumin_r21(nb, work_in, scal, nb_max, work_out)
       CASE('fumax')
          ierr = ma_fumax_r21(nb, work_in, scal, nb_max, work_out)
       CASE('add')
          ierr = ma_add_r21(nb, work_in, scal, nb_max, work_out)
       CASE('subi')
          ierr = ma_subi_r21(nb, work_in, scal, nb_max, work_out)
       CASE('sub')
          ierr = ma_sub_r21(nb, work_in, scal, nb_max, work_out)
       CASE('mult')
          ierr = ma_mult_r21(nb, work_in, scal, nb_max, work_out)
       CASE('div')
          ierr = ma_div_r21(nb, work_in, scal, nb_max, work_out)
       CASE('divi')
          ierr = ma_divi_r21(nb, work_in, scal, nb_max, work_out)
       CASE('power')
          ierr = ma_power_r21(nb, work_in, scal, nb_max, work_out)
       CASE DEFAULT
          CALL histerr(3, "mathop", &
               &      'Unknown operation with a scalar', fun, ' ')
       END SELECT
       IF (ierr > 0) THEN
          CALL histerr(3, "mathop", &
               &      'Error while executing a scalar function', fun, ' ')
       ENDIF
    ENDIF
    !------------------------
  END SUBROUTINE mathop_r21
  !-
  !=== FUNCTIONS (only one argument)
  !-
  INTEGER FUNCTION ma_sin_r21(nb, x, nbo, y)
    INTEGER :: nb(2), nbo, i, j, ij
    REAL :: x(nb(1), nb(2)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = SIN(x(i, j))
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_sin_r21 = 0
    !----------------------
  END FUNCTION ma_sin_r21

  !************************************************

  INTEGER FUNCTION ma_cos_r21(nb, x, nbo, y)
    INTEGER :: nb(2), nbo, i, j, ij
    REAL :: x(nb(1), nb(2)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = COS(x(i, j))
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_cos_r21 = 0
    !----------------------
  END FUNCTION ma_cos_r21

  !************************************************

  INTEGER FUNCTION ma_tan_r21(nb, x, nbo, y)
    INTEGER :: nb(2), nbo, i, j, ij
    REAL :: x(nb(1), nb(2)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = TAN(x(i, j))
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_tan_r21 = 0
    !----------------------
  END FUNCTION ma_tan_r21

  !************************************************

  INTEGER FUNCTION ma_asin_r21(nb, x, nbo, y)
    INTEGER :: nb(2), nbo, i, j, ij
    REAL :: x(nb(1), nb(2)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = ASIN(x(i, j))
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_asin_r21 = 0
    !-----------------------
  END FUNCTION ma_asin_r21

  !************************************************

  INTEGER FUNCTION ma_acos_r21(nb, x, nbo, y)
    INTEGER :: nb(2), nbo, i, j, ij
    REAL :: x(nb(1), nb(2)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = ACOS(x(i, j))
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_acos_r21 = 0
    !-----------------------
  END FUNCTION ma_acos_r21

  !************************************************

  INTEGER FUNCTION ma_atan_r21(nb, x, nbo, y)
    INTEGER :: nb(2), nbo, i, j, ij
    REAL :: x(nb(1), nb(2)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = ATAN(x(i, j))
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_atan_r21 = 0
    !-----------------------
  END FUNCTION ma_atan_r21

  !************************************************

  INTEGER FUNCTION ma_exp_r21(nb, x, nbo, y)
    INTEGER :: nb(2), nbo, i, j, ij
    REAL :: x(nb(1), nb(2)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = EXP(x(i, j))
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_exp_r21 = 0
    !----------------------
  END FUNCTION ma_exp_r21

  !************************************************

  INTEGER FUNCTION ma_alog_r21(nb, x, nbo, y)
    INTEGER :: nb(2), nbo, i, j, ij
    REAL :: x(nb(1), nb(2)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = ALOG(x(i, j))
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_alog_r21 = 0
    !-----------------------
  END FUNCTION ma_alog_r21

  !************************************************

  INTEGER FUNCTION ma_sqrt_r21(nb, x, nbo, y)
    INTEGER :: nb(2), nbo, i, j, ij
    REAL :: x(nb(1), nb(2)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = SQRT(x(i, j))
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_sqrt_r21 = 0
    !-----------------------
  END FUNCTION ma_sqrt_r21

  !************************************************

  INTEGER FUNCTION ma_abs_r21(nb, x, nbo, y)
    INTEGER :: nb(2), nbo, i, j, ij
    REAL :: x(nb(1), nb(2)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = ABS(x(i, j))
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_abs_r21 = 0
    !----------------------
  END FUNCTION ma_abs_r21

  !************************************************

  INTEGER FUNCTION ma_chs_r21(nb, x, nbo, y)
    INTEGER :: nb(2), nbo, i, j, ij
    REAL :: x(nb(1), nb(2)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = x(i, j)*(-1.)
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_chs_r21 = 0
    !----------------------
  END FUNCTION ma_chs_r21

  !************************************************

  INTEGER FUNCTION ma_cels_r21(nb, x, nbo, y)
    INTEGER :: nb(2), nbo, i, j, ij
    REAL :: x(nb(1), nb(2)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = x(i, j)-273.15
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_cels_r21 = 0
    !-----------------------
  END FUNCTION ma_cels_r21

  !************************************************

  INTEGER FUNCTION ma_kelv_r21(nb, x, nbo, y)
    INTEGER :: nb(2), nbo, i, j, ij
    REAL :: x(nb(1), nb(2)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = x(i, j)+273.15
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_kelv_r21 = 0
    !-----------------------
  END FUNCTION ma_kelv_r21

  !************************************************

  INTEGER FUNCTION ma_deg_r21(nb, x, nbo, y)
    INTEGER :: nb(2), nbo, i, j, ij
    REAL :: x(nb(1), nb(2)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = x(i, j)*57.29577951
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_deg_r21 = 0
    !----------------------
  END FUNCTION ma_deg_r21

  !************************************************

  INTEGER FUNCTION ma_rad_r21(nb, x, nbo, y)
    INTEGER :: nb(2), nbo, i, j, ij
    REAL :: x(nb(1), nb(2)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = x(i, j)*0.01745329252
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_rad_r21 = 0
    !----------------------
  END FUNCTION ma_rad_r21

  !************************************************

  INTEGER FUNCTION ma_ident_r21(nb, x, nbo, y)
    INTEGER :: nb(2), nbo, i, j, ij
    REAL :: x(nb(1), nb(2)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = x(i, j)
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_ident_r21 = 0
    !------------------------
  END FUNCTION ma_ident_r21
  !-
  !=== OPERATIONS (two argument)
  !-
  INTEGER FUNCTION ma_add_r21(nb, x, s, nbo, y)
    INTEGER :: nb(2), nbo
    REAL :: x(nb(1), nb(2)), s, y(nbo)

    INTEGER :: i, j, ij
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = x(i, j)+s
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_add_r21 = 0
    !----------------------
  END FUNCTION ma_add_r21

  !************************************************

  INTEGER FUNCTION ma_sub_r21(nb, x, s, nbo, y)
    INTEGER :: nb(2), nbo
    REAL :: x(nb(1), nb(2)), s, y(nbo)

    INTEGER :: i, j, ij
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = x(i, j)-s
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_sub_r21 = 0
    !----------------------
  END FUNCTION ma_sub_r21

  !************************************************

  INTEGER FUNCTION ma_subi_r21(nb, x, s, nbo, y)
    INTEGER :: nb(2), nbo
    REAL :: x(nb(1), nb(2)), s, y(nbo)

    INTEGER :: i, j, ij
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) =  s-x(i, j)
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_subi_r21 = 0
    !-----------------------
  END FUNCTION ma_subi_r21

  !************************************************

  INTEGER FUNCTION ma_mult_r21(nb, x, s, nbo, y)
    INTEGER :: nb(2), nbo
    REAL :: x(nb(1), nb(2)), s, y(nbo)

    INTEGER :: i, j, ij
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = x(i, j)*s
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_mult_r21 = 0
    !-----------------------
  END FUNCTION ma_mult_r21

  !************************************************

  INTEGER FUNCTION ma_div_r21(nb, x, s, nbo, y)
    INTEGER :: nb(2), nbo
    REAL :: x(nb(1), nb(2)), s, y(nbo)

    INTEGER :: i, j, ij
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = x(i, j)/s
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_div_r21 = 0
    !----------------------
  END FUNCTION ma_div_r21

  !************************************************

  INTEGER FUNCTION ma_divi_r21(nb, x, s, nbo, y)
    INTEGER :: nb(2), nbo
    REAL :: x(nb(1), nb(2)), s, y(nbo)

    INTEGER :: i, j, ij
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = s/x(i, j)
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_divi_r21 = 0
    !-----------------------
  END FUNCTION ma_divi_r21

  !************************************************

  INTEGER FUNCTION ma_power_r21(nb, x, s, nbo, y)
    INTEGER :: nb(2), nbo
    REAL :: x(nb(1), nb(2)), s, y(nbo)

    INTEGER :: i, j, ij
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = x(i, j) ** s
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_power_r21 = 0
    !------------------------
  END FUNCTION ma_power_r21

  !************************************************

  INTEGER FUNCTION ma_fumin_r21(nb, x, s, nbo, y)
    INTEGER :: nb(2), nbo
    REAL :: x(nb(1), nb(2)), s, y(nbo)

    INTEGER :: i, j, ij
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = MIN(x(i, j), s)
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_fumin_r21 = 0
    !------------------------
  END FUNCTION ma_fumin_r21

  !************************************************

  INTEGER FUNCTION ma_fumax_r21(nb, x, s, nbo, y)
    INTEGER :: nb(2), nbo
    REAL :: x(nb(1), nb(2)), s, y(nbo)

    INTEGER :: i, j, ij
    !---------------------------------------------------------------------
    ij = 0
    DO j=1, nb(2)
       DO i=1, nb(1)
          ij = ij+1
          y(ij) = MAX(x(i, j), s)
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)
    ma_fumax_r21 = 0
    !------------------------
  END FUNCTION ma_fumax_r21

  !************************************************

  INTEGER FUNCTION ma_fuscat_r21(nb, x, nbi, ind, miss_val, nbo, y)
    INTEGER :: nb(2), nbo, nbi
    INTEGER :: ind(nbi)
    REAL :: x(nb(1), nb(2)), miss_val, y(nbo)

    INTEGER :: i, j, ij, ii, ipos
    !---------------------------------------------------------------------
    ma_fuscat_r21 = 0

    y(1:nbo) = miss_val

    IF (nbi <= nb(1)*nb(2)) THEN
       ipos = 0
       DO ij=1, nbi
          IF (ind(ij) <= nbo .AND. ind(ij) > 0) THEN
             ipos = ipos+1
             j = ((ipos-1)/nb(1))+1
             i = (ipos-(j-1)*nb(1))
             y(ind(ij)) = x(i, j)
          ELSE
             IF (ind(ij) > nbo) ma_fuscat_r21  = ma_fuscat_r21+1
          ENDIF
       ENDDO
       !-- Repeat the data if needed
       IF (MINVAL(ind) < 0) THEN
          DO i=1, nbi
             IF (ind(i) <= 0) THEN
                DO ii=1, ABS(ind(i))-1
                   IF (ind(i+1)+ii <= nbo) THEN
                      y(ind(i+1)+ii) = y(ind(i+1))
                   ELSE
                      ma_fuscat_r21  = ma_fuscat_r21+1
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDIF
    ELSE
       ma_fuscat_r21  = 1
    ENDIF
    !-------------------------
  END FUNCTION ma_fuscat_r21

  !************************************************

  INTEGER FUNCTION ma_fugath_r21(nb, x, nbi, ind, miss_val, nbo, y)
    INTEGER :: nb(2), nbo, nbi
    INTEGER :: ind(nbi)
    REAL :: x(nb(1), nb(2)), miss_val, y(nbo)

    INTEGER :: i, j, ij, ipos
    !---------------------------------------------------------------------
    IF (nbi <= nbo) THEN
       ma_fugath_r21 = 0
       y(1:nbo) = miss_val
       ipos = 0
       DO ij=1, nbi
          IF (ipos+1 <= nbo) THEN
             IF (ind(ij) > 0) THEN
                j = ((ind(ij)-1)/nb(1))+1
                i = (ind(ij)-(j-1)*nb(1))
                ipos = ipos+1
                y(ipos) = x(i, j)
             ENDIF
          ELSE
             IF (ipos+1 > nbo) ma_fugath_r21  = ma_fugath_r21+1
          ENDIF
       ENDDO
    ELSE
       ma_fugath_r21 = 1
    ENDIF
    nbo = ipos
    !-------------------------
  END FUNCTION ma_fugath_r21

  !************************************************

  INTEGER FUNCTION ma_fufill_r21(nb, x, nbi, ind, miss_val, nbo, y)
    INTEGER :: nb(2), nbo, nbi
    INTEGER :: ind(nbi)
    REAL :: x(nb(1), nb(2)), miss_val, y(nbo)

    INTEGER :: i, j, ij, ii, ipos
    !---------------------------------------------------------------------
    ma_fufill_r21 = 0

    IF (nbi <= nb(1)*nb(2)) THEN
       ipos = 0
       DO ij=1, nbi
          IF (ind(ij) <= nbo .AND. ind(ij) > 0) THEN
             ipos = ipos+1
             j = ((ipos-1)/nb(1))+1
             i = (ipos-(j-1)*nb(1))
             y(ind(ij)) = x(i, j)
          ELSE
             IF (ind(ij) > nbo) ma_fufill_r21  = ma_fufill_r21+1
          ENDIF
       ENDDO
       !-- Repeat the data if needed
       IF (MINVAL(ind) < 0) THEN
          DO i=1, nbi
             IF (ind(i) <= 0) THEN
                DO ii=1, ABS(ind(i))-1
                   IF (ind(i+1)+ii <= nbo) THEN
                      y(ind(i+1)+ii) = y(ind(i+1))
                   ELSE
                      ma_fufill_r21  = ma_fufill_r21+1
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDIF
    ELSE
       ma_fufill_r21  = 1
    ENDIF
    !-------------------------
  END FUNCTION ma_fufill_r21

  !************************************************

  INTEGER FUNCTION ma_fucoll_r21(nb, x, nbi, ind, miss_val, nbo, y)
    INTEGER :: nb(2), nbo, nbi
    INTEGER :: ind(nbi)
    REAL :: x(nb(1), nb(2)), miss_val, y(nbo)

    INTEGER :: i, j, ij, ipos
    !---------------------------------------------------------------------
    IF (nbi <= nbo) THEN
       ma_fucoll_r21 = 0
       ipos = 0
       DO ij=1, nbi
          IF (ipos+1 <= nbo) THEN
             IF (ind(ij) > 0) THEN
                j = ((ind(ij)-1)/nb(1))+1
                i = (ind(ij)-(j-1)*nb(1))
                ipos = ipos+1
                y(ipos) = x(i, j)
             ENDIF
          ELSE
             IF (ipos+1 > nbo) ma_fucoll_r21  = ma_fucoll_r21+1
          ENDIF
       ENDDO
    ELSE
       ma_fucoll_r21 = 1
    ENDIF
    nbo = ipos
    !-------------------------
  END FUNCTION ma_fucoll_r21

  !************************************************

  INTEGER FUNCTION ma_fuundef_r21(nb, x, nbi, ind, miss_val, nbo, y)
    INTEGER :: nb(2), nbo, nbi
    INTEGER :: ind(nbi)
    REAL :: x(nb(1), nb(2)), miss_val, y(nbo)

    INTEGER :: i, j, ij
    !---------------------------------------------------------------------
    IF (nbi <= nbo .AND. nbo == nb(1)*nb(2)) THEN
       ma_fuundef_r21 = 0
       DO ij=1, nbo
          j = ((ij-1)/nb(1))+1
          i = (ij-(j-1)*nb(1))
          y(ij) = x(i, j)
       ENDDO
       DO i=1, nbi
          IF (ind(i) <= nbo .AND. ind(i) > 0) THEN
             y(ind(i)) =  miss_val
          ELSE
             IF (ind(i) > nbo) ma_fuundef_r21 = ma_fuundef_r21+1
          ENDIF
       ENDDO
    ELSE
       ma_fuundef_r21 = 1
    ENDIF
    !--------------------------
  END FUNCTION ma_fuundef_r21

  !************************************************

  INTEGER FUNCTION ma_fuonly_r21(nb, x, nbi, ind, miss_val, nbo, y)
    INTEGER :: nb(2), nbo, nbi
    INTEGER :: ind(nbi)
    REAL :: x(nb(1), nb(2)), miss_val, y(nbo)

    INTEGER :: i, j, ij
    !---------------------------------------------------------------------
    IF (     (nbi <= nbo).AND.(nbo == nb(1)*nb(2)) &
         &    .AND.ALL(ind(1:nbi) <= nbo) ) THEN
       ma_fuonly_r21 = 0
       y(1:nbo) = miss_val
       DO ij=1, nbi
          IF (ind(ij) > 0) THEN
             j = ((ind(ij)-1)/nb(1))+1
             i = (ind(ij)-(j-1)*nb(1))
             y(ind(ij)) =  x(i, j)
          ENDIF
       ENDDO
    ELSE
       ma_fuonly_r21 = 1
    ENDIF
    !-------------------------
  END FUNCTION ma_fuonly_r21

  !************************************************

  !************************************************

  SUBROUTINE mathop_r31 &
       &  (fun, nb, work_in, miss_val, nb_index, nindex, scal, nb_max, work_out)
    !- This subroutines gives an interface to the various operations
    !- which are allowed. The interface is general enough to allow its use
    !- for other cases.

    !- INPUT

    !- fun      : function to be applied to the vector of data
    !- nb       : Length of input vector
    !- work_in  : Input vector of data (REAL)
    !- miss_val : The value of the missing data flag (it has to be a
    !-            maximum value, in f90 : huge( a real ))
    !- nb_index : Length of index vector
    !- nindex   : Vector of indices
    !- scal     : A scalar value for vector/scalar operations
    !- nb_max   : maximum length of output vector

    !- OUTPUT

    !- nb_max   : Actual length of output variable
    !- work_out : Output vector after the operation was applied
    USE errioipsl, ONLY : histerr

    CHARACTER(LEN=7) :: fun
    INTEGER :: nb(3), nb_max, nb_index
    INTEGER :: nindex(nb_index)
    REAL :: work_in(nb(1), nb(2), nb(3)), scal, miss_val
    REAL :: work_out(nb_max)

    INTEGER :: ierr

    INTRINSIC SIN, COS, TAN, ASIN, ACOS, ATAN, EXP, ALOG, SQRT, ABS
    !---------------------------------------------------------------------
    ierr = 0

    IF (scal >= miss_val-1.) THEN
       IF (INDEX(indexfu, fun(1:LEN_TRIM(fun))) == 0) THEN
          SELECT CASE (fun)
          CASE('sin')
             ierr = ma_sin_r31(nb, work_in, nb_max, work_out)
          CASE('cos')
             ierr = ma_cos_r31(nb, work_in, nb_max, work_out)
          CASE('tan')
             ierr = ma_tan_r31(nb, work_in, nb_max, work_out)
          CASE('asin')
             ierr = ma_asin_r31(nb, work_in, nb_max, work_out)
          CASE('acos')
             ierr = ma_acos_r31(nb, work_in, nb_max, work_out)
          CASE('atan')
             ierr = ma_atan_r31(nb, work_in, nb_max, work_out)
          CASE('exp')
             ierr = ma_exp_r31(nb, work_in, nb_max, work_out)
          CASE('log')
             ierr = ma_alog_r31(nb, work_in, nb_max, work_out)
          CASE('sqrt')
             ierr = ma_sqrt_r31(nb, work_in, nb_max, work_out)
          CASE('chs')
             ierr = ma_chs_r31(nb, work_in, nb_max, work_out)
          CASE('abs')
             ierr = ma_abs_r31(nb, work_in, nb_max, work_out)
          CASE('cels')
             ierr = ma_cels_r31(nb, work_in, nb_max, work_out)
          CASE('kelv')
             ierr = ma_kelv_r31(nb, work_in, nb_max, work_out)
          CASE('deg')
             ierr = ma_deg_r31(nb, work_in, nb_max, work_out)
          CASE('rad')
             ierr = ma_rad_r31(nb, work_in, nb_max, work_out)
          CASE('ident')
             ierr = ma_ident_r31(nb, work_in, nb_max, work_out)
          CASE DEFAULT
             CALL histerr(3, "mathop", &
                  &        'scalar variable undefined and no indexing', &
                  &        'but still unknown function', fun)
          END SELECT
          IF (ierr > 0) THEN
             CALL histerr(3, "mathop", &
                  &        'Error while executing a simple function', fun, ' ')
          ENDIF
       ELSE
          SELECT CASE (fun)
          CASE('gather')
             ierr = ma_fugath_r31(nb, work_in, nb_index, nindex, &
                  &                           miss_val, nb_max, work_out)
          CASE('scatter')
             IF (nb_index > (nb(1)*nb(2)*nb(3))) THEN
                work_out(1:nb_max) = miss_val
                ierr=1
             ELSE
                ierr = ma_fuscat_r31(nb, work_in, nb_index, nindex, &
                     &                             miss_val, nb_max, work_out)
             ENDIF
          CASE('coll')
             ierr = ma_fucoll_r31(nb, work_in, nb_index, nindex, &
                  &                           miss_val, nb_max, work_out)
          CASE('fill')
             ierr = ma_fufill_r31(nb, work_in, nb_index, nindex, &
                  &                           miss_val, nb_max, work_out)
          CASE('undef')
             ierr = ma_fuundef_r31(nb, work_in, nb_index, nindex, &
                  &                           miss_val, nb_max, work_out)
          CASE('only')
             ierr = ma_fuonly_r31(nb, work_in, nb_index, nindex, &
                  &                           miss_val, nb_max, work_out)
          CASE DEFAULT
             CALL histerr(3, "mathop", &
                  &        'scalar variable undefined and indexing', &
                  &        'was requested but with unknown function', fun)
          END SELECT
          IF (ierr > 0) THEN
             CALL histerr(3, "mathop_r31", &
                  &        'Error while executing an indexing function', fun, ' ')
          ENDIF
       ENDIF
    ELSE
       SELECT CASE (fun)
       CASE('fumin')
          ierr = ma_fumin_r31(nb, work_in, scal, nb_max, work_out)
       CASE('fumax')
          ierr = ma_fumax_r31(nb, work_in, scal, nb_max, work_out)
       CASE('add')
          ierr = ma_add_r31(nb, work_in, scal, nb_max, work_out)
       CASE('subi')
          ierr = ma_subi_r31(nb, work_in, scal, nb_max, work_out)
       CASE('sub')
          ierr = ma_sub_r31(nb, work_in, scal, nb_max, work_out)
       CASE('mult')
          ierr = ma_mult_r31(nb, work_in, scal, nb_max, work_out)
       CASE('div')
          ierr = ma_div_r31(nb, work_in, scal, nb_max, work_out)
       CASE('divi')
          ierr = ma_divi_r31(nb, work_in, scal, nb_max, work_out)
       CASE('power')
          ierr = ma_power_r31(nb, work_in, scal, nb_max, work_out)
       CASE DEFAULT
          CALL histerr(3, "mathop", &
               &      'Unknown operation with a scalar', fun, ' ')
       END SELECT
       IF (ierr > 0) THEN
          CALL histerr(3, "mathop", &
               &      'Error while executing a scalar function', fun, ' ')
       ENDIF
    ENDIF
    !------------------------
  END SUBROUTINE mathop_r31
  !-
  !=== FUNCTIONS (only one argument)
  !-
  INTEGER FUNCTION ma_sin_r31(nb, x, nbo, y)
    INTEGER :: nb(3), nbo, i, j, k, ij
    REAL :: x(nb(1), nb(2), nb(3)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = SIN(x(i, j, k))
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_sin_r31 = 0
    !----------------------
  END FUNCTION ma_sin_r31

  !************************************************

  INTEGER FUNCTION ma_cos_r31(nb, x, nbo, y)
    INTEGER :: nb(3), nbo, i, j, k, ij
    REAL :: x(nb(1), nb(2), nb(3)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = COS(x(i, j, k))
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_cos_r31 = 0
    !----------------------
  END FUNCTION ma_cos_r31

  !************************************************

  INTEGER FUNCTION ma_tan_r31(nb, x, nbo, y)
    INTEGER :: nb(3), nbo, i, j, k, ij
    REAL :: x(nb(1), nb(2), nb(3)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = TAN(x(i, j, k))
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_tan_r31 = 0
    !----------------------
  END FUNCTION ma_tan_r31

  !************************************************

  INTEGER FUNCTION ma_asin_r31(nb, x, nbo, y)
    INTEGER :: nb(3), nbo, i, j, k, ij
    REAL :: x(nb(1), nb(2), nb(3)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = ASIN(x(i, j, k))
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_asin_r31 = 0
    !-----------------------
  END FUNCTION ma_asin_r31

  !************************************************

  INTEGER FUNCTION ma_acos_r31(nb, x, nbo, y)
    INTEGER :: nb(3), nbo, i, j, k, ij
    REAL :: x(nb(1), nb(2), nb(3)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = ACOS(x(i, j, k))
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_acos_r31 = 0
    !-----------------------
  END FUNCTION ma_acos_r31

  !************************************************

  INTEGER FUNCTION ma_atan_r31(nb, x, nbo, y)
    INTEGER :: nb(3), nbo, i, j, k, ij
    REAL :: x(nb(1), nb(2), nb(3)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = ATAN(x(i, j, k))
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_atan_r31 = 0
    !-----------------------
  END FUNCTION ma_atan_r31

  !************************************************

  INTEGER FUNCTION ma_exp_r31(nb, x, nbo, y)
    INTEGER :: nb(3), nbo, i, j, k, ij
    REAL :: x(nb(1), nb(2), nb(3)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = EXP(x(i, j, k))
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_exp_r31 = 0
    !----------------------
  END FUNCTION ma_exp_r31

  !************************************************

  INTEGER FUNCTION ma_alog_r31(nb, x, nbo, y)
    INTEGER :: nb(3), nbo, i, j, k, ij
    REAL :: x(nb(1), nb(2), nb(3)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = ALOG(x(i, j, k))
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_alog_r31 = 0
    !-----------------------
  END FUNCTION ma_alog_r31

  !************************************************

  INTEGER FUNCTION ma_sqrt_r31(nb, x, nbo, y)
    INTEGER :: nb(3), nbo, i, j, k, ij
    REAL :: x(nb(1), nb(2), nb(3)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = SQRT(x(i, j, k))
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_sqrt_r31 = 0
    !-----------------------
  END FUNCTION ma_sqrt_r31

  !************************************************

  INTEGER FUNCTION ma_abs_r31(nb, x, nbo, y)
    INTEGER :: nb(3), nbo, i, j, k, ij
    REAL :: x(nb(1), nb(2), nb(3)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = ABS(x(i, j, k))
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_abs_r31 = 0
    !----------------------
  END FUNCTION ma_abs_r31

  !************************************************

  INTEGER FUNCTION ma_chs_r31(nb, x, nbo, y)
    INTEGER :: nb(3), nbo, i, j, k, ij
    REAL :: x(nb(1), nb(2), nb(3)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = x(i, j, k)*(-1.)
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_chs_r31 = 0
    !----------------------
  END FUNCTION ma_chs_r31

  !************************************************

  INTEGER FUNCTION ma_cels_r31(nb, x, nbo, y)
    INTEGER :: nb(3), nbo, i, j, k, ij
    REAL :: x(nb(1), nb(2), nb(3)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = x(i, j, k)-273.15
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_cels_r31 = 0
    !-----------------------
  END FUNCTION ma_cels_r31

  !************************************************

  INTEGER FUNCTION ma_kelv_r31(nb, x, nbo, y)
    INTEGER :: nb(3), nbo, i, j, k, ij
    REAL :: x(nb(1), nb(2), nb(3)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = x(i, j, k)+273.15
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_kelv_r31 = 0
    !-----------------------
  END FUNCTION ma_kelv_r31

  !************************************************

  INTEGER FUNCTION ma_deg_r31(nb, x, nbo, y)
    INTEGER :: nb(3), nbo, i, j, k, ij
    REAL :: x(nb(1), nb(2), nb(3)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = x(i, j, k)*57.29577951
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_deg_r31 = 0
    !----------------------
  END FUNCTION ma_deg_r31

  !************************************************

  INTEGER FUNCTION ma_rad_r31(nb, x, nbo, y)
    INTEGER :: nb(3), nbo, i, j, k, ij
    REAL :: x(nb(1), nb(2), nb(3)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = x(i, j, k)*0.01745329252
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_rad_r31 = 0
    !----------------------
  END FUNCTION ma_rad_r31

  !************************************************

  INTEGER FUNCTION ma_ident_r31(nb, x, nbo, y)
    INTEGER :: nb(3), nbo, i, j, k, ij
    REAL :: x(nb(1), nb(2), nb(3)), y(nbo)
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = x(i, j, k)
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_ident_r31 = 0
    !------------------------
  END FUNCTION ma_ident_r31
  !-
  !=== OPERATIONS (two argument)
  !-
  INTEGER FUNCTION ma_add_r31(nb, x, s, nbo, y)
    INTEGER :: nb(3), nbo
    REAL :: x(nb(1), nb(2), nb(3)), s, y(nbo)

    INTEGER :: i, j, k, ij
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = x(i, j, k)+s
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_add_r31 = 0
    !----------------------
  END FUNCTION ma_add_r31

  !************************************************

  INTEGER FUNCTION ma_sub_r31(nb, x, s, nbo, y)
    INTEGER :: nb(3), nbo
    REAL :: x(nb(1), nb(2), nb(3)), s, y(nbo)

    INTEGER :: i, j, k, ij
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = x(i, j, k)-s
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_sub_r31 = 0
    !----------------------
  END FUNCTION ma_sub_r31

  !************************************************

  INTEGER FUNCTION ma_subi_r31(nb, x, s, nbo, y)
    INTEGER :: nb(3), nbo
    REAL :: x(nb(1), nb(2), nb(3)), s, y(nbo)

    INTEGER :: i, j, k, ij
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) =  s-x(i, j, k)
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_subi_r31 = 0
    !-----------------------
  END FUNCTION ma_subi_r31

  !************************************************

  INTEGER FUNCTION ma_mult_r31(nb, x, s, nbo, y)
    INTEGER :: nb(3), nbo
    REAL :: x(nb(1), nb(2), nb(3)), s, y(nbo)

    INTEGER :: i, j, k, ij
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = x(i, j, k)*s
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_mult_r31 = 0
    !-----------------------
  END FUNCTION ma_mult_r31

  !************************************************

  INTEGER FUNCTION ma_div_r31(nb, x, s, nbo, y)
    INTEGER :: nb(3), nbo
    REAL :: x(nb(1), nb(2), nb(3)), s, y(nbo)

    INTEGER :: i, j, k, ij
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = x(i, j, k)/s
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_div_r31 = 0
    !----------------------
  END FUNCTION ma_div_r31

  !************************************************

  INTEGER FUNCTION ma_divi_r31(nb, x, s, nbo, y)
    INTEGER :: nb(3), nbo
    REAL :: x(nb(1), nb(2), nb(3)), s, y(nbo)

    INTEGER :: i, j, k, ij
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = s/x(i, j, k)
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_divi_r31 = 0
    !-----------------------
  END FUNCTION ma_divi_r31

  !************************************************

  INTEGER FUNCTION ma_power_r31(nb, x, s, nbo, y)
    INTEGER :: nb(3), nbo
    REAL :: x(nb(1), nb(2), nb(3)), s, y(nbo)

    INTEGER :: i, j, k, ij
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = x(i, j, k)**s
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_power_r31 = 0
    !------------------------
  END FUNCTION ma_power_r31

  !************************************************

  INTEGER FUNCTION ma_fumin_r31(nb, x, s, nbo, y)
    INTEGER :: nb(3), nbo
    REAL :: x(nb(1), nb(2), nb(3)), s, y(nbo)

    INTEGER :: i, j, k, ij
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = MIN(x(i, j, k), s)
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_fumin_r31 = 0
    !------------------------
  END FUNCTION ma_fumin_r31

  !************************************************

  INTEGER FUNCTION ma_fumax_r31(nb, x, s, nbo, y)
    INTEGER :: nb(3), nbo
    REAL :: x(nb(1), nb(2), nb(3)), s, y(nbo)

    INTEGER :: i, j, k, ij
    !---------------------------------------------------------------------
    ij = 0
    DO k=1, nb(3)
       DO j=1, nb(2)
          DO i=1, nb(1)
             ij = ij+1
             y(ij) = MAX(x(i, j, k), s)
          ENDDO
       ENDDO
    ENDDO

    nbo = nb(1)*nb(2)*nb(3)
    ma_fumax_r31 = 0
    !------------------------
  END FUNCTION ma_fumax_r31

  !************************************************

  INTEGER FUNCTION ma_fuscat_r31(nb, x, nbi, ind, miss_val, nbo, y)
    INTEGER :: nb(3), nbo, nbi
    INTEGER :: ind(nbi)
    REAL :: x(nb(1), nb(2), nb(3)), miss_val, y(nbo)

    INTEGER :: i, j, k, ij, ii, ipos, ipp, isb
    !---------------------------------------------------------------------
    ma_fuscat_r31 = 0

    y(1:nbo) = miss_val

    IF (nbi <= nb(1)*nb(2)*nb(3)) THEN
       ipos = 0
       isb = nb(1)*nb(2)
       DO ij=1, nbi
          IF (ind(ij) <= nbo .AND. ind(ij) > 0) THEN
             ipos = ipos+1
             k = ((ipos-1)/isb)+1
             ipp = ipos-(k-1)*isb
             j = ((ipp-1)/nb(1))+1
             i = (ipp-(j-1)*nb(1))
             y(ind(ij)) = x(i, j, k)
          ELSE
             IF (ind(ij) > nbo) ma_fuscat_r31  = ma_fuscat_r31+1
          ENDIF
       ENDDO
       !-- Repeat the data if needed
       IF (MINVAL(ind) < 0) THEN
          DO i=1, nbi
             IF (ind(i) <= 0) THEN
                DO ii=1, ABS(ind(i))-1
                   IF (ind(i+1)+ii <= nbo) THEN
                      y(ind(i+1)+ii) = y(ind(i+1))
                   ELSE
                      ma_fuscat_r31  = ma_fuscat_r31+1
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDIF
    ELSE
       ma_fuscat_r31  = 1
    ENDIF
    !-------------------------
  END FUNCTION ma_fuscat_r31

  !************************************************

  INTEGER FUNCTION ma_fugath_r31(nb, x, nbi, ind, miss_val, nbo, y)
    INTEGER :: nb(3), nbo, nbi
    INTEGER :: ind(nbi)
    REAL :: x(nb(1), nb(2), nb(3)), miss_val, y(nbo)

    INTEGER :: i, j, k, ij, ipos, ipp, isb
    !---------------------------------------------------------------------
    IF (nbi <= nbo) THEN
       ma_fugath_r31 = 0
       y(1:nbo) = miss_val
       ipos = 0
       isb = nb(1)*nb(2)
       DO ij=1, nbi
          IF (ipos+1 <= nbo) THEN
             IF (ind(ij) > 0) THEN
                k = ((ind(ij)-1)/isb)+1
                ipp = ind(ij)-(k-1)*isb
                j = ((ipp-1)/nb(1))+1
                i = (ipp-(j-1)*nb(1))
                ipos = ipos+1
                y(ipos) = x(i, j, k)
             ENDIF
          ELSE
             IF (ipos+1 > nbo) ma_fugath_r31  = ma_fugath_r31+1
          ENDIF
       ENDDO
    ELSE
       ma_fugath_r31 = 1
    ENDIF
    nbo = ipos
    !-------------------------
  END FUNCTION ma_fugath_r31

  !************************************************

  INTEGER FUNCTION ma_fufill_r31(nb, x, nbi, ind, miss_val, nbo, y)
    INTEGER :: nb(3), nbo, nbi
    INTEGER :: ind(nbi)
    REAL :: x(nb(1), nb(2), nb(3)), miss_val, y(nbo)

    INTEGER :: i, j, k, ij, ii, ipos, ipp, isb
    !---------------------------------------------------------------------
    ma_fufill_r31 = 0
    IF (nbi <= nb(1)*nb(2)*nb(3)) THEN
       ipos = 0
       isb = nb(1)*nb(2)
       DO ij=1, nbi
          IF (ind(ij) <= nbo .AND. ind(ij) > 0) THEN
             ipos = ipos+1
             k = ((ipos-1)/isb)+1
             ipp = ipos-(k-1)*isb
             j = ((ipp-1)/nb(1))+1
             i = (ipp-(j-1)*nb(1))
             y(ind(ij)) = x(i, j, k)
          ELSE
             IF (ind(ij) > nbo) ma_fufill_r31  = ma_fufill_r31+1
          ENDIF
       ENDDO
       !-- Repeat the data if needed
       IF (MINVAL(ind) < 0) THEN
          DO i=1, nbi
             IF (ind(i) <= 0) THEN
                DO ii=1, ABS(ind(i))-1
                   IF (ind(i+1)+ii <= nbo) THEN
                      y(ind(i+1)+ii) = y(ind(i+1))
                   ELSE
                      ma_fufill_r31  = ma_fufill_r31+1
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDIF
    ELSE
       ma_fufill_r31  = 1
    ENDIF
    !-------------------------
  END FUNCTION ma_fufill_r31

  !************************************************

  INTEGER FUNCTION ma_fucoll_r31(nb, x, nbi, ind, miss_val, nbo, y)
    INTEGER :: nb(3), nbo, nbi
    INTEGER :: ind(nbi)
    REAL :: x(nb(1), nb(2), nb(3)), miss_val, y(nbo)

    INTEGER :: i, j, k, ij, ipos, ipp, isb
    !---------------------------------------------------------------------
    IF (nbi <= nbo) THEN
       ma_fucoll_r31 = 0
       ipos = 0
       isb = nb(1)*nb(2)
       DO ij=1, nbi
          IF (ipos+1 <= nbo) THEN
             IF (ind(ij) > 0) THEN
                k = ((ind(ij)-1)/isb)+1
                ipp = ind(ij)-(k-1)*isb
                j = ((ipp-1)/nb(1))+1
                i = (ipp-(j-1)*nb(1))
                ipos = ipos+1
                y(ipos) = x(i, j, k)
             ENDIF
          ELSE
             IF (ipos+1 > nbo) ma_fucoll_r31  = ma_fucoll_r31+1
          ENDIF
       ENDDO
    ELSE
       ma_fucoll_r31 = 1
    ENDIF
    nbo = ipos
    !-------------------------
  END FUNCTION ma_fucoll_r31

  !************************************************

  INTEGER FUNCTION ma_fuundef_r31(nb, x, nbi, ind, miss_val, nbo, y)
    INTEGER :: nb(3), nbo, nbi
    INTEGER :: ind(nbi)
    REAL :: x(nb(1), nb(2), nb(3)), miss_val, y(nbo)

    INTEGER :: i, j, k, ij, ipp, isb
    !---------------------------------------------------------------------
    IF (nbi <= nbo .AND. nbo == nb(1)*nb(2)*nb(3)) THEN
       ma_fuundef_r31 = 0
       isb = nb(1)*nb(2)
       DO ij=1, nbo
          k = ((ij-1)/isb)+1
          ipp = ij-(k-1)*isb
          j = ((ipp-1)/nb(1))+1
          i = (ipp-(j-1)*nb(1))
          y(ij) = x(i, j, k)
       ENDDO
       DO i=1, nbi
          IF (ind(i) <= nbo .AND. ind(i) > 0) THEN
             y(ind(i)) =  miss_val
          ELSE
             IF (ind(i) > nbo) ma_fuundef_r31  = ma_fuundef_r31+1
          ENDIF
       ENDDO
    ELSE
       ma_fuundef_r31 = 1
    ENDIF
    !--------------------------
  END FUNCTION ma_fuundef_r31

  !************************************************

  INTEGER FUNCTION ma_fuonly_r31(nb, x, nbi, ind, miss_val, nbo, y)
    INTEGER :: nb(3), nbo, nbi
    INTEGER :: ind(nbi)
    REAL :: x(nb(1), nb(2), nb(3)), miss_val, y(nbo)

    INTEGER :: i, j, k, ij, ipp, isb
    !---------------------------------------------------------------------
    IF (     (nbi <= nbo).AND.(nbo == nb(1)*nb(2)*nb(3)) &
         &    .AND.ALL(ind(1:nbi) <= nbo) ) THEN
       ma_fuonly_r31 = 0
       y(1:nbo) = miss_val
       isb = nb(1)*nb(2)
       DO ij=1, nbi
          IF (ind(ij) > 0) THEN
             k = ((ind(ij)-1)/isb)+1
             ipp = ind(ij)-(k-1)*isb
             j = ((ipp-1)/nb(1))+1
             i = (ipp-(j-1)*nb(1))
             y(ind(ij)) = x(i, j, k)
          ENDIF
       ENDDO
    ELSE
       ma_fuonly_r31 = 1
    ENDIF
    !-------------------------
  END FUNCTION ma_fuonly_r31

END MODULE mathop_m
