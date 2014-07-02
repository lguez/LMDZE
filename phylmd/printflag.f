SUBROUTINE printflag(radpas, ok_journe, ok_instan, ok_region)

  ! From phylmd/printflag.F, v 1.1.1.1 2004/05/19 12:53:09
  ! Auteur : P. Le Van

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: radpas
  LOGICAL, INTENT (IN) :: ok_journe, ok_instan, ok_region

  !--------------------------------------------------

  PRINT *, 'Choix des principales clés de la physique'

  PRINT 8, radpas
  PRINT 4, ok_journe, ok_instan, ok_region

4 FORMAT ('ok_journe= ', L3, ', ok_instan = ', L3, ', ok_region = ', L3)
8 FORMAT ('radpas = ', I4)

END SUBROUTINE printflag
