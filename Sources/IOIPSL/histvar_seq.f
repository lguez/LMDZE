module histvar_seq_m

  implicit none

contains

  SUBROUTINE histvar_seq (pfid, pvarname, pvid)

    ! This subroutine optimized the search for the variable in the table.
    ! In a first phase it will learn the succession of the variables
    ! called and then it will use the table to guess what comes next.
    ! It is the best solution to avoid lengthy searches through array
    ! vectors.

    ! ARGUMENTS :

    ! pfid  : id of the file on which we work
    ! pvarname : The name of the variable we are looking for
    ! pvid     : The var id we found

    USE find_str_m, ONLY: find_str
    USE errioipsl, ONLY : histerr
    use histcom_var

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

    !--------------------------------------------------------------------
    nb = nb_var(pfid)

    IF (learning(pfid)) THEN

       !- 1.0 We compute the length over which we are going
       !-     to check the overlap

       IF (overlap(pfid) <= 0) THEN
          IF (nb_var(pfid) > 6) THEN
             overlap(pfid) = nb_var(pfid)/3*2
          ELSE
             overlap(pfid) = nb_var(pfid)
          ENDIF
       ENDIF

       !- 1.1 Find the position of this string

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

       !- 1.2 If we have not given up we store the position
       !-     in the sequence of calls

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

          !--- 1.3 Check if we have found the right overlap

          IF (varseq_len(pfid) .GE. overlap(pfid)*2) THEN

             !----- We skip a few variables if needed as they could come
             !----- from the initialisation of the model.

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

       !- 2.0 Now we know how the calls to histwrite are sequenced
       !-     and we can get a guess at the var ID

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

          !--- We only keep the new position if we have found the variable
          !--- this way. This way an out of sequence call to histwrite does
          !--- not defeat the process.

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

end module histvar_seq_m
