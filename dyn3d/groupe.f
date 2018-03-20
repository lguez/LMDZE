module groupe_m

  IMPLICIT NONE

contains

  SUBROUTINE groupe(pbaru, pbarv, pbarum, pbarvm, wm)

    ! From dyn3d/groupe.F, v 1.1.1.1 2004/05/19 12:53:06

    ! sous-programme servant a fitlrer les champs de flux de masse aux    
    ! poles en "regroupant" les mailles 2 par 2 puis 4 par 4 etc. au fur  
    ! et a mesure qu'on se rapproche du pole.                             

    ! en entree: pbaru et pbarv                                     
    ! en sortie:  pbarum, pbarvm et wm.                                    

    ! remarque, le wm est recalcule a partir des pbaru pbarv et on n'a    
    ! donc pas besoin de w en entree.

    USE dimensions
    USE paramet_m
    USE comconst
    USE disvert_m
    USE comgeom
    use vitvert_m, only: vitvert

    INTEGER, PARAMETER:: ngroup = 3

    REAL pbaru(iip1, jjp1, llm), pbarv(iip1, jjm, llm)
    REAL, intent(out):: pbarum(iip1, jjp1, llm)
    real pbarvm(iip1, jjm, llm)
    REAL wm(iip1, jjp1, llm)

    REAL zconvm(iip1, jjp1, llm), zconvmm(iip1, jjp1, llm)
    REAL uu
    INTEGER i, j, l
    LOGICAL:: firstcall = .TRUE.

    !------------------------------------------------------

    IF (firstcall) THEN
       IF (mod(iim, 2**ngroup) /= 0) then
          print *, 'groupe: bad iim'
          STOP 1
       end IF
       firstcall = .FALSE.
    END IF

    ! Champs 1D                                                           

    CALL convflu(pbaru, pbarv, llm, zconvm)

    CALL scopy(ijp1llm, zconvm, 1, zconvmm, 1)
    CALL scopy(ijmllm, pbarv, 1, pbarvm, 1)

    CALL groupeun(jjp1, llm, zconvmm)
    CALL groupeun(jjm, llm, pbarvm)

    ! Champs 3D                                                           

    DO l = 1, llm
       DO j = 2, jjm
          uu = pbaru(iim, j, l)
          DO i = 1, iim
             uu = uu + pbarvm(i, j, l) - pbarvm(i, j-1, l) - zconvmm(i, j, l)
             pbarum(i, j, l) = uu
          END DO
          pbarum(iip1, j, l) = pbarum(1, j, l)
       END DO
    END DO
    pbarum(:, 1, :) = 0
    pbarum(:, jjm + 1, :) = 0

    ! integration de la convergence de masse de haut  en bas
    DO l = 1, llm
       DO j = 1, jjp1
          DO i = 1, iip1
             zconvmm(i, j, l) = zconvmm(i, j, l)
          END DO
       END DO
    END DO
    DO l = llm - 1, 1, -1
       DO j = 1, jjp1
          DO i = 1, iip1
             zconvmm(i, j, l) = zconvmm(i, j, l) + zconvmm(i, j, l+1)
          END DO
       END DO
    END DO

    wm = vitvert(zconvmm)

  END SUBROUTINE groupe

end module groupe_m
