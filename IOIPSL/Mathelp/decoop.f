module decoop_m

  implicit none

contains

  SUBROUTINE decoop (pstr, nbops_max, missing_val, opps, scal, nbops)

    USE errioipsl, ONLY : histerr
    use findsep_m, only: findsep

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
                  'Unknown operations of type X+scal', f_char(1), pstr)
          ENDIF
       ELSE
          IF (check) WRITE(*, *) 'decoop : More complex operation'
          IF ( f_char(1) == '(' .AND. f_char(2) == ')' ) THEN
             !------ Type : sin(X)
             opp_str = str(s_pos(1)+1:f_pos(1)-1)
             scal_str = '?'
             str = str(1:s_pos(1))//'X'//str(f_pos(2)+1:leng)
          ELSE IF (    (f_char(1) == '(' .AND. f_char(2) == ', ')&
               .OR.(f_char(1) == ', ' .AND. f_char(2) == ')')) THEN
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
                  .AND.(xpos-f_pos(1) == 1).AND.(.NOT.prio) ) THEN
                !-------- Type : ... scal+X ...
                IF (f_char(1) == '-' .OR. f_char(1) == '/') THEN
                   opp_str = f_char(1)//'I'
                ELSE
                   opp_str = f_char(1)
                ENDIF
                scal_str = str(s_pos(1)+1:f_pos(1)-1)
                str = str(1:s_pos(1))//'X'//str(f_pos(1)+2:leng)
             ELSE IF (     (INDEX(ops, f_char(2)) > 0) &
                  .AND.(f_pos(2)-xpos == 1) ) THEN
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
                  'Unknown function', opp_str(1:LEN_TRIM(opp_str)), ' ')
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
                     'Unknown operation', opp_str(1:leng), ' ')
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
                  'Strange scalar you have here ', scal_str, pstr)
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

  END SUBROUTINE decoop

end module decoop_m
