module groupe_m

  IMPLICIT NONE

contains

  SUBROUTINE groupe(pbaru, pbarv, pbarum, pbarvm, wm)

    ! From dyn3d/groupe.F, v 1.1.1.1 2004/05/19 12:53:06

    ! Sous-programme servant \`a filtrer les champs de flux de masse
    ! aux p\^oles en "regroupant" les mailles 2 par 2 puis 4 par 4
    ! etc. au fur et \`a mesure qu'on se rapproche du p\^ole.

    ! Remarque : wm est recalcul\'e \`a partir de pbaru et pbarv, et
    ! on n'a donc pas besoin de w en entr\'ee.

    use conf_gcm_m, only: ngroup
    use convflu_m, only: convflu
    use groupeun_m, only: groupeun
    USE dimensions, only: iim, jjm, llm
    USE paramet_m, only: iip1, jjp1
    use vitvert_m, only: vitvert

    REAL, intent(in):: pbaru(:, :, :) ! (iip1, jjp1, llm)
    REAL, intent(in):: pbarv(:, :, :) ! (iip1, jjm, llm)
    REAL, intent(out):: pbarum(iip1, jjp1, llm), pbarvm(iip1, jjm, llm)
    REAL, intent(out):: wm(iip1, jjp1, llm)

    ! Local:
    REAL zconvm(iip1, jjp1, llm)
    REAL uu
    INTEGER i, j, l
    
    !------------------------------------------------------

    ! Champs 1D                                                           
    CALL convflu(pbaru, pbarv, llm, zconvm)
    pbarvm = pbarv
    CALL groupeun(zconvm)
    CALL groupeun(pbarvm)

    ! Champs 3D                                                           

    DO l = 1, llm
       DO j = 2, jjm
          uu = pbaru(iim, j, l)
          
          DO i = 1, iim
             uu = uu + pbarvm(i, j, l) - pbarvm(i, j-1, l) - zconvm(i, j, l)
             pbarum(i, j, l) = uu
          END DO
          
          pbarum(iip1, j, l) = pbarum(1, j, l)
       END DO
    END DO
    
    pbarum(:, 1, :) = 0
    pbarum(:, jjm + 1, :) = 0

    ! Int\'egration de la convergence de masse de haut en bas :
    
    DO l = 1, llm
       DO j = 1, jjp1
          DO i = 1, iip1
             zconvm(i, j, l) = zconvm(i, j, l)
          END DO
       END DO
    END DO
    
    DO l = llm - 1, 1, -1
       DO j = 1, jjp1
          DO i = 1, iip1
             zconvm(i, j, l) = zconvm(i, j, l) + zconvm(i, j, l+1)
          END DO
       END DO
    END DO

    wm = vitvert(zconvm)

  END SUBROUTINE groupe

end module groupe_m
