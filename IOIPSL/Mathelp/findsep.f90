module findsep_m

  implicit none

contains

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
    use mathelp, only: seps
    use cleanstr_m, only: cleanstr

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

end module findsep_m
