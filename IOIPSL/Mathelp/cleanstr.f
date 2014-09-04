module cleanstr_m

  implicit none

contains

  SUBROUTINE cleanstr(str)
    !- We clean up the string by taking out the extra () and puting
    !- everything in lower case except for the X describing the variable
    use strlowercase_m, only: strlowercase
    use mathelp, only: seps

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

end module cleanstr_m
