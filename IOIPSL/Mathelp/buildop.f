module buildop_m

  implicit none

contains

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
    use decoop_m, only: decoop

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

end module buildop_m
