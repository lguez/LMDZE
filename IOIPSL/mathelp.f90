MODULE mathelp

  ! From mathelp.f90, version 2.0 2004/04/05 14:47:50

  implicit none

  PRIVATE
  PUBLIC :: moycum, trans_buff, buildop

  !- Variables used to detect and identify the operations
  CHARACTER(LEN=80), SAVE :: seps='( ) , + - / * ^'

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

    USE errioipsl, ONLY : histerr

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
    USE errioipsl, ONLY : histerr

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
    CHARACTER(LEN=80), SAVE :: ops = '+ - * / ^'
    CHARACTER(LEN=80), SAVE :: mima = 'min max'
    CHARACTER(LEN=250), SAVE :: funcs = &
         'sin cos tan asin acos atan exp log sqrt chs abs ' &
         //'cels kelv deg rad gather scatter fill coll undef only ident'
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
    USE errioipsl, ONLY : histerr

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
    use strlowercase_m, only: strlowercase

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

  SUBROUTINE moycum (opp, np, px, py, pwx)
    !- Does time operations
    USE errioipsl, ONLY : histerr

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
