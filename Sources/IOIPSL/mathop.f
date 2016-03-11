MODULE mathop_m

  ! From mathelp.f90, version 2.0 2004/04/05 14:47:50

  implicit none

  PRIVATE
  PUBLIC mathop

  INTERFACE mathop
     MODULE PROCEDURE mathop_r11, mathop_r21, mathop_r31
  END INTERFACE mathop

  ! Variables used to detect and identify the operations
  CHARACTER(LEN=120):: indexfu = 'gather, scatter, fill, coll, undef, only'

CONTAINS

  SUBROUTINE mathop_r11(fun, nb, work_in, miss_val, nb_index, nindex, scal, &
       nb_max, work_out)

    ! This subroutine gives an interface to the various operations
    ! which are allowed. The interface is general enough to allow its
    ! use for other cases.

    ! INPUT

    ! fun      : function to be applied to the vector of data
    ! nb       : Length of input vector
    ! work_in  : Input vector of data (REAL)
    ! miss_val : The value of the missing data flag (it has to be a
    !            maximum value, in f90 : huge( a real ))
    ! nb_index : Length of index vector
    ! nindex   : Vector of indices
    ! scal     : A scalar value for vector/scalar operations
    ! nb_max   : maximum length of output vector

    ! OUTPUT

    ! nb_max   : Actual length of output variable
    ! work_out : Output vector after the operation was applied
    USE errioipsl, ONLY : histerr
    USE mathop2, ONLY : ma_abs_r11, ma_acos_r11, ma_add_r11, ma_alog_r11, &
         ma_asin_r11, ma_atan_r11, ma_cels_r11, ma_chs_r11, ma_cos_r11, &
         ma_deg_r11, ma_divi_r11, ma_div_r11, ma_exp_r11, ma_fucoll_r11, &
         ma_fufill_r11, ma_fugath_r11, ma_fumax_r11, ma_fumin_r11, &
         ma_fuonly_r11, ma_fuscat_r11, ma_fuundef_r11, ma_ident_r11, &
         ma_kelv_r11, ma_mult_r11, ma_power_r11, ma_rad_r11, ma_sin_r11, &
         ma_sqrt_r11, ma_subi_r11, ma_sub_r11, ma_tan_r11

    CHARACTER(LEN=7) :: fun
    INTEGER :: nb, nb_max, nb_index
    INTEGER :: nindex(nb_index)
    REAL :: work_in(nb), scal, miss_val
    REAL :: work_out(nb_max)

    INTEGER :: ierr

    !--------------------------------------------------------------------
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
        'scalar variable undefined and no indexing', &
        'but still unknown function', fun)
          END SELECT
          IF (ierr > 0) THEN
             CALL histerr(3, "mathop", &
        'Error while executing a simple function', fun, ' ')
          ENDIF
       ELSE
          SELECT CASE (fun)
          CASE('gather')
             ierr = ma_fugath_r11(nb, work_in, nb_index, nindex, &
                            miss_val, nb_max, work_out)
          CASE('scatter')
             IF (nb_index > nb) THEN
                work_out(1:nb_max) = miss_val
                ierr=1
             ELSE
                ierr = ma_fuscat_r11(nb, work_in, nb_index, nindex, &
                              miss_val, nb_max, work_out)
             ENDIF
          CASE('coll')
             ierr = ma_fucoll_r11(nb, work_in, nb_index, nindex, nb_max, &
                  work_out)
          CASE('fill')
             ierr = ma_fufill_r11(nb, work_in, nb_index, nindex, nb_max, &
                  work_out)
          CASE('undef')
             ierr = ma_fuundef_r11(nb, work_in, nb_index, nindex, &
                            miss_val, nb_max, work_out)
          CASE('only')
             ierr = ma_fuonly_r11(nb, work_in, nb_index, nindex, &
                            miss_val, nb_max, work_out)
          CASE DEFAULT
             CALL histerr(3, "mathop", &
        'scalar variable undefined and indexing', &
        'was requested but with unknown function', fun)
          END SELECT
          IF (ierr > 0) THEN
             CALL histerr(3, "mathop_r11", &
        'Error while executing an indexing function', fun, ' ')
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
      'Unknown operation with a scalar', fun, ' ')
       END SELECT
       IF (ierr > 0) THEN
          CALL histerr(3, "mathop", &
      'Error while executing a scalar function', fun, ' ')
       ENDIF
    ENDIF
    !-----------------------
  END SUBROUTINE mathop_r11
  !
  SUBROUTINE mathop_r21(fun, nb, work_in, miss_val, nb_index, nindex, scal, &
       nb_max, work_out)

    ! This subroutines gives an interface to the various operations
    ! which are allowed. The interface is general enough to allow its use
    ! for other cases.

    ! INPUT

    ! fun      : function to be applied to the vector of data
    ! nb       : Length of input vector
    ! work_in  : Input vector of data (REAL)
    ! miss_val : The value of the missing data flag (it has to be a
    !            maximum value, in f90 : huge( a real ))
    ! nb_index : Length of index vector
    ! nindex   : Vector of indices
    ! scal     : A scalar value for vector/scalar operations
    ! nb_max   : maximum length of output vector

    ! OUTPUT

    ! nb_max   : Actual length of output variable
    ! work_out : Output vector after the operation was applied
    USE errioipsl, ONLY : histerr
    USE mathop2, ONLY : ma_abs_r21, ma_acos_r21, ma_add_r21, ma_alog_r21, &
         ma_asin_r21, ma_atan_r21, ma_cels_r21, ma_chs_r21, ma_cos_r21, &
         ma_deg_r21, ma_divi_r21, ma_div_r21, ma_exp_r21, ma_fucoll_r21, &
         ma_fufill_r21, ma_fugath_r21, ma_fumax_r21, ma_fumin_r21, &
         ma_fuonly_r21, ma_fuscat_r21, ma_fuundef_r21, ma_ident_r21, &
         ma_kelv_r21, ma_mult_r21, ma_power_r21, ma_rad_r21, ma_sin_r21, &
         ma_sqrt_r21, ma_subi_r21, ma_sub_r21, ma_tan_r21

    CHARACTER(LEN=7) :: fun
    INTEGER :: nb(2), nb_max, nb_index
    INTEGER :: nindex(nb_index)
    REAL :: work_in(nb(1), nb(2)), scal, miss_val
    REAL :: work_out(nb_max)

    INTEGER :: ierr

    !--------------------------------------------------------------------
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
        'scalar variable undefined and no indexing', &
        'but still unknown function', fun)
          END SELECT
          IF (ierr > 0) THEN
             CALL histerr(3, "mathop", &
        'Error while executing a simple function', fun, ' ')
          ENDIF
       ELSE
          SELECT CASE (fun)
          CASE('gather')
             ierr = ma_fugath_r21(nb, work_in, nb_index, nindex, &
                            miss_val, nb_max, work_out)
          CASE('scatter')
             IF (nb_index > (nb(1)*nb(2)) ) THEN
                work_out(1:nb_max) = miss_val
                ierr=1
             ELSE
                ierr = ma_fuscat_r21(nb, work_in, nb_index, nindex, &
                             miss_val, nb_max, work_out)
             ENDIF
          CASE('coll')
             ierr = ma_fucoll_r21(nb, work_in, nb_index, nindex, nb_max, &
                  work_out)
          CASE('fill')
             ierr = ma_fufill_r21(nb, work_in, nb_index, nindex, nb_max, &
                  work_out)
          CASE('undef')
             ierr = ma_fuundef_r21(nb, work_in, nb_index, nindex, &
                           miss_val, nb_max, work_out)
          CASE('only')
             ierr = ma_fuonly_r21(nb, work_in, nb_index, nindex, &
                           miss_val, nb_max, work_out)
          CASE DEFAULT
             CALL histerr(3, "mathop", &
        'scalar variable undefined and indexing', &
        'was requested but with unknown function', fun)
          END SELECT
          IF (ierr > 0) THEN
             CALL histerr(3, "mathop_r21", &
        'Error while executing an indexing function', fun, ' ')
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
      'Unknown operation with a scalar', fun, ' ')
       END SELECT
       IF (ierr > 0) THEN
          CALL histerr(3, "mathop", &
      'Error while executing a scalar function', fun, ' ')
       ENDIF
    ENDIF
    !-----------------------
  END SUBROUTINE mathop_r21
  !
  SUBROUTINE mathop_r31(fun, nb, work_in, miss_val, nb_index, nindex, scal, &
       nb_max, work_out)

    ! This subroutines gives an interface to the various operations
    ! which are allowed. The interface is general enough to allow its use
    ! for other cases.

    ! INPUT

    ! fun      : function to be applied to the vector of data
    ! nb       : Length of input vector
    ! work_in  : Input vector of data (REAL)
    ! miss_val : The value of the missing data flag (it has to be a
    !            maximum value, in f90 : huge( a real ))
    ! nb_index : Length of index vector
    ! nindex   : Vector of indices
    ! scal     : A scalar value for vector/scalar operations
    ! nb_max   : maximum length of output vector

    ! OUTPUT

    ! nb_max   : Actual length of output variable
    ! work_out : Output vector after the operation was applied
    USE errioipsl, ONLY : histerr
    USE mathop2, ONLY : ma_abs_r31, ma_acos_r31, ma_add_r31, ma_alog_r31, &
         ma_asin_r31, ma_atan_r31, ma_cels_r31, ma_chs_r31, ma_cos_r31, &
         ma_deg_r31, ma_divi_r31, ma_div_r31, ma_exp_r31, ma_fucoll_r31, &
         ma_fufill_r31, ma_fugath_r31, ma_fumax_r31, ma_fumin_r31, &
         ma_fuonly_r31, ma_fuscat_r31, ma_fuundef_r31, ma_ident_r31, &
         ma_kelv_r31, ma_mult_r31, ma_power_r31, ma_rad_r31, ma_sin_r31, &
         ma_sqrt_r31, ma_subi_r31, ma_sub_r31, ma_tan_r31

    CHARACTER(LEN=7) :: fun
    INTEGER :: nb(3), nb_max, nb_index
    INTEGER :: nindex(nb_index)
    REAL :: work_in(nb(1), nb(2), nb(3)), scal, miss_val
    REAL :: work_out(nb_max)

    INTEGER :: ierr

    !--------------------------------------------------------------------
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
                  'scalar variable undefined and no indexing', &
                  'but still unknown function', fun)
          END SELECT
          IF (ierr > 0) THEN
             CALL histerr(3, "mathop", &
                  'Error while executing a simple function', fun, ' ')
          ENDIF
       ELSE
          SELECT CASE (fun)
          CASE('gather')
             ierr = ma_fugath_r31(nb, work_in, nb_index, nindex, &
                  miss_val, nb_max, work_out)
          CASE('scatter')
             IF (nb_index > (nb(1)*nb(2)*nb(3))) THEN
                work_out(1:nb_max) = miss_val
                ierr=1
             ELSE
                ierr = ma_fuscat_r31(nb, work_in, nb_index, nindex, &
                     miss_val, nb_max, work_out)
             ENDIF
          CASE('coll')
             ierr = ma_fucoll_r31(nb, work_in, nb_index, nindex, nb_max, &
                  work_out)
          CASE('fill')
             ierr = ma_fufill_r31(nb, work_in, nb_index, nindex, nb_max, &
                  work_out)
          CASE('undef')
             ierr = ma_fuundef_r31(nb, work_in, nb_index, nindex, &
                  miss_val, nb_max, work_out)
          CASE('only')
             ierr = ma_fuonly_r31(nb, work_in, nb_index, nindex, &
                  miss_val, nb_max, work_out)
          CASE DEFAULT
             CALL histerr(3, "mathop", &
                  'scalar variable undefined and indexing', &
                  'was requested but with unknown function', fun)
          END SELECT
          IF (ierr > 0) THEN
             CALL histerr(3, "mathop_r31", &
                  'Error while executing an indexing function', fun, ' ')
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
               'Unknown operation with a scalar', fun, ' ')
       END SELECT
       IF (ierr > 0) THEN
          CALL histerr(3, "mathop", &
               'Error while executing a scalar function', fun, ' ')
       ENDIF
    ENDIF

  END SUBROUTINE mathop_r31

END MODULE mathop_m
