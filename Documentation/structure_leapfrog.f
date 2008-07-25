  itau = 0
  ...
1 CONTINUE
  ... bloc 1
  forward = .TRUE.
  leapf   = .FALSE.
  ... bloc 2
2 CONTINUE
  ... bloc 3
  IF (.NOT. purmats) THEN
     IF (forward .OR. leapf) THEN
        itau= itau + 1
        ... bloc 4
     end IF
     ... bloc 5
     IF (MOD(itau, iperiod) == 0) THEN
        GO TO 1
     ELSE IF (MOD(itau-1, iperiod) == 0) THEN
        IF (forward) THEN
           forward = .FALSE.
           leapf = .FALSE.
           GO TO 2
        ELSE
           leapf =  .TRUE.
           ... bloc 6
           GO TO 2
        END IF
     ELSE
        leapf = .TRUE.
        ... bloc 7
        GO TO 2
     END IF
  ELSE
     IF (forward)  THEN
        itau =  itau + 1
        ... bloc 8
        forward =  .FALSE.
        ... bloc 9
        GO TO 2
     ELSE
        ... bloc 10
        forward = .TRUE.
        GO TO 1
     ENDIF
  END IF
