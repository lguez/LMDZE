MODULE mathelp

  ! From mathelp.f90, version 2.0 2004/04/05 14:47:50

  USE errioipsl, ONLY : histerr
  USE stringop

  implicit none

  PRIVATE
  PUBLIC :: mathop, moycum, trans_buff, buildop

  INTERFACE mathop
     MODULE PROCEDURE mathop_r11, mathop_r21, mathop_r31
  END INTERFACE

  !- Variables used to detect and identify the operations
  CHARACTER(LEN=80), SAVE :: seps='( ) , + - / * ^', ops = '+ - * / ^', &
       mima = 'min max'
  CHARACTER(LEN=250), SAVE :: funcs = &
       'sin cos tan asin acos atan exp log sqrt chs abs ' &
       //'cels kelv deg rad gather scatter fill coll undef only ident'
  CHARACTER(LEN=120), SAVE :: indexfu = &
       'gather, scatter, fill, coll, undef, only'

CONTAINS

  SUBROUTINE buildop (str, ex_topps, topp, nbops_max, &
       &                  missing_val, opps, scal, nbops)
    !- This subroutine decomposes the input string in the elementary
    !- functions which need to be applied to the vector of data.
    !- This vector is represented by X in the string.
    !- This subroutine is the driver of the decomposition and gets
    !- the time operation but then call decoop for the other operations
    !- INPUT

    !- str      : String containing the operations
    !- ex_toops : The time operations that can be expected
    !-            within the string

    !- OUTPUT

    CHARACTER(LEN=80) :: str
    CHARACTER(LEN=*) :: ex_topps
    CHARACTER(LEN=7) :: topp
    INTEGER :: nbops_max, nbops
    CHARACTER(LEN=7) :: opps(nbops_max)
    REAL :: scal(nbops_max), missing_val

    CHARACTER(LEN=80) :: new_str
    INTEGER :: leng, ind_opb, ind_clb

    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    IF (check) WRITE(*, *) 'buildop : Some preliminary cleaning'

    leng = LEN_TRIM(str)
    IF ( str(1:1) == '(' .AND. str(leng:leng) == ')' ) THEN
       str = str(2:leng-1)
       leng = leng-2
    ENDIF

    IF (check) &
         &  WRITE(*, *) 'buildop : Starting to test the various options'

    IF (leng <= 5 .AND. INDEX(ex_topps, str(1:leng)) > 0) THEN
       IF (check) WRITE(*, *) 'buildop : Time operation only'
       nbops = 0
       topp = str(1:leng)
    ELSE
       IF (check) THEN
          WRITE(*, *) 'buildop : Time operation and something else'
       ENDIF
       !--
       ind_opb = INDEX(str(1:leng), '(')
       IF (ind_opb > 0) THEN
          IF (INDEX(ex_topps, str(1:ind_opb-1)) > 0) THEN
             IF (check) THEN
                WRITE(*, '(2a)') &
                     &          ' buildop : Extract time operation from : ', str
             ENDIF
             topp = str(1:ind_opb-1)
             ind_clb = INDEX(str(1:leng), ')', BACK=.TRUE.)
             new_str = str(ind_opb+1:ind_clb-1)
             IF (check) THEN
                WRITE(*, '(2a, 2I3)') &
                     &          ' buildop : Call decoop ', new_str, ind_opb, ind_clb
             ENDIF
             CALL decoop (new_str, nbops_max, missing_val, opps, scal, nbops)
          ELSE
             CALL histerr(3, 'buildop', &
                  &        'time opperation does not exist', str(1:ind_opb-1), ' ')
          ENDIF
       ELSE
          CALL histerr(3, 'buildop', &
               &      'some long opperation exists but wihout parenthesis', &
               &      str(1:leng), ' ')
       ENDIF
    ENDIF

    IF (check) THEN
       DO leng=1, nbops
          WRITE(*, *) &
               &      'buildop : i -- opps, scal : ', leng, opps(leng), scal(leng)
       ENDDO
    ENDIF
    !---------------------
  END SUBROUTINE buildop

  !************************************************

  SUBROUTINE decoop (pstr, nbops_max, missing_val, opps, scal, nbops)
    CHARACTER(LEN=80) :: pstr
    INTEGER :: nbops_max, nbops
    CHARACTER(LEN=7) :: opps(nbops_max)
    REAL :: scal(nbops_max), missing_val

    CHARACTER(LEN=1) :: f_char(2), s_char(2)
    INTEGER :: nbsep, f_pos(2), s_pos(2)
    CHARACTER(LEN=20) :: opp_str, scal_str
    CHARACTER(LEN=80) :: str
    INTEGER :: xpos, leng, ppos, epos, int_tmp
    CHARACTER(LEN=3) :: tl, dl
    CHARACTER(LEN=10) :: fmt

    LOGICAL :: check = .FALSE., prio
    !---------------------------------------------------------------------
    IF (check) WRITE(*, '(2a)') ' decoop : Incoming string : ', pstr

    nbops = 0
    str = pstr

    CALL findsep (str, nbsep, f_char, f_pos, s_char, s_pos)
    IF (check) WRITE(*, *) 'decoop : Out of findsep', nbsep
    DO WHILE (nbsep > 0)
       xpos = INDEX(str, 'X')
       leng = LEN_TRIM(str)
       nbops = nbops+1
       !--
       IF (check) THEN
          WRITE(*, *) 'decoop : str -->', str(1:leng)
          WRITE(*, *) s_char(1), '-', f_char(1), '|', f_char(2), '-', s_char(2)
          WRITE(*, *) s_pos(1), '-', f_pos(1), '|', f_pos(2), '-', s_pos(2)
       ENDIF
       !--
       IF (nbops > nbops_max-1) THEN
          CALL histerr(3, 'decoop', 'Expression too complex', str, ' ')
       ENDIF
       !--
       IF (check) WRITE(*, *) 'decoop : --', nbops, ' ', str(1:leng)
       !---
       !-- Start the analysis of the syntax. 3 types of constructs
       !-- are recognized.  They are scanned sequentialy
       !---
       IF (nbsep == 1) THEN
          IF (check) WRITE(*, *) 'decoop : Only one operation'
          IF (INDEX(ops, f_char(1)) > 0) THEN
             !------ Type : scal+X
             IF (f_char(1) == '-' .OR. f_char(1) == '/') THEN
                opp_str = f_char(1)//'I'
             ELSE
                opp_str = f_char(1)
             ENDIF
             scal_str = str(s_pos(1)+1:f_pos(1)-1)
             str = 'X'
          ELSE IF (INDEX(ops, f_char(2)) > 0) THEN
             !------ Type : X+scal
             opp_str = f_char(2)
             scal_str = str(f_pos(2)+1:s_pos(2)-1)
             str = 'X'
          ELSE
             CALL histerr(3, 'decoop', &
                  &        'Unknown operations of type X+scal', f_char(1), pstr)
          ENDIF
       ELSE
          IF (check) WRITE(*, *) 'decoop : More complex operation'
          IF ( f_char(1) == '(' .AND. f_char(2) == ')' ) THEN
             !------ Type : sin(X)
             opp_str = str(s_pos(1)+1:f_pos(1)-1)
             scal_str = '?'
             str = str(1:s_pos(1))//'X'//str(f_pos(2)+1:leng)
          ELSE IF (    (f_char(1) == '(' .AND. f_char(2) == ', ')&
               &             .OR.(f_char(1) == ', ' .AND. f_char(2) == ')')) THEN
             !------ Type : max(X, scal) or max(scal, X)
             IF (f_char(1) == '(' .AND. s_char(2) == ')') THEN
                !-------- Type : max(X, scal)
                opp_str = str(f_pos(1)-3:f_pos(1)-1)
                scal_str = str(f_pos(2)+1:s_pos(2)-1)
                str = str(1:f_pos(1)-4)//'X'//str(s_pos(2)+1:leng)
             ELSE IF (f_char(1) == ', ' .AND. s_char(1) == '(') THEN
                !-------- Type : max(scal, X)
                opp_str = str(s_pos(1)-3:s_pos(1)-1)
                scal_str = str(s_pos(1)+1:f_pos(1)-1)
                str = str(1:s_pos(1)-4)//'X'//str(f_pos(2)+1:leng)
             ELSE
                CALL histerr(3, 'decoop', 'Syntax error 1', str, ' ')
             ENDIF
          ELSE
             prio = (f_char(2) == '*').OR.(f_char(2) == '^')
             IF (     (INDEX(ops, f_char(1)) > 0) &
                  &          .AND.(xpos-f_pos(1) == 1).AND.(.NOT.prio) ) THEN
                !-------- Type : ... scal+X ...
                IF (f_char(1) == '-' .OR. f_char(1) == '/') THEN
                   opp_str = f_char(1)//'I'
                ELSE
                   opp_str = f_char(1)
                ENDIF
                scal_str = str(s_pos(1)+1:f_pos(1)-1)
                str = str(1:s_pos(1))//'X'//str(f_pos(1)+2:leng)
             ELSE IF (     (INDEX(ops, f_char(2)) > 0) &
                  &               .AND.(f_pos(2)-xpos == 1) ) THEN
                !-------- Type : ... X+scal ...
                opp_str = f_char(2)
                scal_str = str(f_pos(2)+1:s_pos(2)-1)
                str = str(1:f_pos(2)-2)//'X'//str(s_pos(2):leng)
             ELSE
                CALL histerr(3, 'decoop', 'Syntax error 2', str, ' ')
             ENDIF
          ENDIF
       ENDIF
       !---
       IF (check) WRITE(*, *) 'decoop : Finished syntax, str = ', TRIM(str)
       !---
       !-- Now that the different components of the operation are identified
       !-- we transform them into what is going to be used in the program
       !---
       IF (INDEX(scal_str, '?') > 0) THEN
          IF (INDEX(funcs, opp_str(1:LEN_TRIM(opp_str))) > 0) THEN
             opps(nbops) = opp_str(1:LEN_TRIM(opp_str))
             scal(nbops) =  missing_val
          ELSE
             CALL histerr(3, 'decoop', &
                  &        'Unknown function', opp_str(1:LEN_TRIM(opp_str)), ' ')
          ENDIF
       ELSE
          leng = LEN_TRIM(opp_str)
          IF (INDEX(mima, opp_str(1:leng)) > 0) THEN
             opps(nbops) = 'fu'//opp_str(1:leng)
          ELSE
             IF (INDEX(opp_str(1:leng), '+') > 0) THEN
                opps(nbops) = 'add'
             ELSE IF (INDEX(opp_str(1:leng), '-I') > 0) THEN
                opps(nbops) = 'subi'
             ELSE IF (INDEX(opp_str(1:leng), '-') > 0) THEN
                opps(nbops) = 'sub'
             ELSE IF (INDEX(opp_str(1:leng), '*') > 0) THEN
                opps(nbops) = 'mult'
             ELSE IF (INDEX(opp_str(1:leng), '/') > 0) THEN
                opps(nbops) = 'div'
             ELSE IF (INDEX(opp_str(1:leng), '/I') > 0) THEN
                opps(nbops) = 'divi'
             ELSE IF (INDEX(opp_str(1:leng), '^') > 0) THEN
                opps(nbops) = 'power'
             ELSE
                CALL histerr(3, 'decoop', &
                     &          'Unknown operation', opp_str(1:leng), ' ')
             ENDIF
          ENDIF
          !-----
          leng = LEN_TRIM(scal_str)
          ppos = INDEX(scal_str, '.')
          epos = INDEX(scal_str, 'e')
          IF (epos == 0) epos = INDEX(scal_str, 'E')
          !-----
          !---- Try to catch a few errors
          !-----
          IF (INDEX(ops, scal_str) > 0) THEN
             CALL histerr(3, 'decoop', &
                  &        'Strange scalar you have here ', scal_str, pstr)
          ENDIF
          IF (epos > 0) THEN
             WRITE(tl, '(I3.3)') leng
             WRITE(dl, '(I3.3)') epos-ppos-1
             fmt='(e'//tl//'.'//dl//')'
             READ(scal_str, fmt) scal(nbops)
          ELSE IF (ppos > 0) THEN
             WRITE(tl, '(I3.3)') leng
             WRITE(dl, '(I3.3)') leng-ppos
             fmt='(f'//tl//'.'//dl//')'
             READ(scal_str, fmt) scal(nbops)
          ELSE
             WRITE(tl, '(I3.3)') leng
             fmt = '(I'//tl//')'
             READ(scal_str, fmt) int_tmp
             scal(nbops) = REAL(int_tmp)
          ENDIF
       ENDIF
       IF (check) WRITE(*, *) 'decoop : Finished interpretation'
       CALL findsep(str, nbsep, f_char, f_pos, s_char, s_pos)
    ENDDO
    !--------------------
  END SUBROUTINE decoop

  !************************************************

  SUBROUTINE findsep (str, nbsep, f_char, f_pos, s_char, s_pos)
    !- Subroutine finds all separators in a given string
    !- It returns the following information about str :
    !-   f_char : The first separation character
    !-            (1 for before and 2 for after)
    !-   f_pos  : The position of the first separator
    !-   s_char : The second separation character
    !-            (1 for before and 2 for after)
    !-   s_pos  : The position of the second separator
    CHARACTER(LEN=80) :: str
    INTEGER :: nbsep
    CHARACTER(LEN=1), DIMENSION(2) :: f_char, s_char
    INTEGER, DIMENSION(2) :: f_pos, s_pos

    CHARACTER(LEN=70) :: str_tmp
    LOGICAL :: f_found, s_found
    INTEGER :: ind, xpos, leng, i

    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    IF (check) WRITE(*, *) 'findsep : call cleanstr: ', TRIM(str)

    CALL cleanstr(str)

    IF (check) WRITE(*, *) 'findsep : out of cleanstr: ', TRIM(str)

    xpos = INDEX(str, 'X')
    leng = LEN_TRIM(str)

    f_pos(1:2) = (/ 0, leng+1 /)
    f_char(1:2) = (/ '?', '?' /)
    s_pos(1:2) = (/ 0, leng+1 /)
    s_char(1:2) = (/ '?', '?' /)

    nbsep = 0

    f_found = .FALSE.
    s_found = .FALSE.
    IF (xpos > 1) THEN
       DO i=xpos-1, 1, -1
          ind = INDEX(seps, str(i:i))
          IF (ind > 0) THEN
             IF (.NOT.f_found) THEN
                f_char(1) = str(i:i)
                f_pos(1) = i
                nbsep = nbsep+1
                f_found = .TRUE.
             ELSE IF (.NOT.s_found) THEN
                s_char(1) = str(i:i)
                s_pos(1) = i
                nbsep = nbsep+1
                s_found = .TRUE.
             ENDIF
          ENDIF
       ENDDO
    ENDIF

    f_found = .FALSE.
    s_found = .FALSE.
    IF (xpos < leng) THEN
       DO i=xpos+1, leng
          ind = INDEX(seps, str(i:i))
          IF (ind > 0) THEN
             IF (.NOT.f_found) THEN
                f_char(2) = str(i:i)
                f_pos(2) = i
                nbsep = nbsep+1
                f_found = .TRUE.
             ELSE IF (.NOT.s_found) THEN
                s_char(2) = str(i:i)
                s_pos(2) = i
                nbsep = nbsep+1
                s_found = .TRUE.
             ENDIF
          ENDIF
       ENDDO
    ENDIF

    IF (nbsep > 4) THEN
       WRITE(str_tmp, '("number :", I3)') nbsep
       CALL histerr(3, 'findsep', &
            &    'How can I find that many separators', str_tmp, str)
    ENDIF

    IF (check) WRITE(*, *) 'Finished findsep : ', nbsep, leng
    !---------------------
  END SUBROUTINE findsep

  !************************************************

  SUBROUTINE cleanstr(str)
    !- We clean up the string by taking out the extra () and puting
    !- everything in lower case except for the X describing the variable
    CHARACTER(LEN=80) :: str

    INTEGER :: ind, leng, ic, it
    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    leng = LEN_TRIM(str)
    CALL strlowercase(str)

    ind = INDEX(str, 'x')
    IF (check) THEN
       WRITE (*, *) 'cleanstr 1.0 : ind = ', ind, &
            &               ' str = ', str(1:leng), '---'
    ENDIF

    ! If the character before the x is not a letter then we can assume
    ! that it is the variable and promote it to a capital letter

    DO WHILE (ind > 0)
       ic = 0
       IF (ind > 1) ic = IACHAR(str(ind-1:ind-1))
       IF (ic < 97 .OR. ic > 122) THEN
          str(ind:ind) = 'X'
       ENDIF
       it = INDEX(str(ind+1:leng), 'x')
       IF (it > 0) THEN
          ind = ind+it
       ELSE
          ind = it
       ENDIF
    ENDDO

    IF (check) WRITE (*, *) 'cleanstr 2.0 : str = ', str(1:leng), '---'

    IF ( str(1:1) == '(' .AND. str(leng:leng) == ')' ) THEN
       str = str(2:leng-1)
    ENDIF

    IF (check) WRITE (*, *) 'cleanstr 3.0 : str = ', str(1:leng), '---'

    leng = LEN_TRIM(str)
    ind = INDEX(str, '((X))')
    IF (ind > 0) THEN
       str=str(1:ind-1)//'(X)'//str(ind+5:leng)//'  '
    ENDIF

    IF (check) WRITE (*, *) 'cleanstr 4.0 : str = ', str(1:leng), '---'

    leng = LEN_TRIM(str)
    ind = INDEX(str, '(X)')
    IF (ind > 0 .AND. ind+3 < leng) THEN
       IF (      (INDEX(seps, str(ind-1:ind-1)) > 0) &
            &      .AND. (INDEX(seps, str(ind+3:ind+3)) > 0) ) THEN
          str=str(1:ind-1)//'X'//str(ind+3:leng)//'  '
       ENDIF
    ENDIF

    IF (check) WRITE (*, *) 'cleanstr 5.0 : str = ', str(1:leng), '---'

    leng = LEN_TRIM(str)
    ind = INDEX(str(1:leng), ' ')
    DO WHILE (ind > 0)
       str=str(1:ind-1)//str(ind+1:leng)//' '
       leng = LEN_TRIM(str)
       ind = INDEX(str(1:leng), ' ')
    ENDDO

    IF (check) WRITE (*, *) 'cleanstr 6.0 : str = ', str(1:leng), '---'
    !----------------------
  END SUBROUTINE cleanstr

  !************************************************

  !************************************************

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

  !************************************************

  SUBROUTINE moycum (opp, np, px, py, pwx)
    !- Does time operations
    CHARACTER(LEN=7) :: opp
    INTEGER :: np
    REAL, DIMENSION(:) :: px, py
    INTEGER :: pwx
    !---------------------------------------------------------------------
    IF (pwx /= 0) THEN
       IF      (opp == 'ave') THEN
          px(1:np)=(px(1:np)*pwx+py(1:np))/REAL(pwx+1)
       ELSE IF (opp == 't_sum') THEN
          px(1:np)=px(1:np)+py(1:np)
       ELSE IF ( (opp == 'l_min').OR.(opp == 't_min') ) THEN
          px(1:np)=MIN(px(1:np), py(1:np))
       ELSE IF ( (opp == 'l_max').OR.(opp == 't_max') ) THEN
          px(1:np)=MAX(px(1:np), py(1:np))
       ELSE
          CALL histerr(3, "moycum", 'Unknown time operation', opp, ' ')
       ENDIF
    ELSE
       IF      (opp == 'l_min') THEN
          px(1:np)=MIN(px(1:np), py(1:np))
       ELSE IF (opp == 'l_max') THEN
          px(1:np)=MAX(px(1:np), py(1:np))
       ELSE
          px(1:np)=py(1:np)
       ENDIF
    ENDIF
    !--------------------
  END SUBROUTINE moycum

  !************************************************

  SUBROUTINE trans_buff (ox, sx, oy, sy, oz, sz, xsz, ysz, zsz, v3d, sl, v1d)
    !- This subroutine extracts from the full 3D variable the slab of
    !- data that will be used later. Perhaps there are hardware routines
    !- for this task on some computers. This routine will be obsolete in
    !- a F90 environnement

    !- INPUT
    !- ox  : Origin of slab of data in X
    !- sx  : Size of slab in X
    !- oy  : Origin of slab of data in Y
    !- sy  : Size of slab in Y
    !- oz  : Origin of slab of data in Z
    !- sz  : Size of slab in Z
    !- xsz, ysz, zsz : 3 sizes of full variable v3d
    !- v3d : The full 3D variable
    !- sl  : size of variable v1d
    !- v1d : The 1D variable containing the slab

    INTEGER :: ox, sx, oy, sy, oz, sz
    INTEGER :: xsz, ysz, zsz
    INTEGER :: sl
    REAL :: v3d(xsz, ysz, zsz)
    REAL :: v1d(sl)

    !---------------------------------------------------------------------

    V1d(:sx*sy*sz) = reshape(V3d(ox:ox-1+sx, oy:oy-1+sy, oz:oz-1+sz), &
         SHAPE=(/sx*sy*sz/))

  END SUBROUTINE trans_buff

END MODULE mathelp
