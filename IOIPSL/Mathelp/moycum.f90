module moycum_m

  implicit none

contains

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

end module moycum_m
