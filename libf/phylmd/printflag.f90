SUBROUTINE printflag(radpas, ok_ocean, ok_oasis, ok_journe, ok_instan, &
     ok_region)

  ! From phylmd/printflag.F, v 1.1.1.1 2004/05/19 12:53:09
  ! Auteur : P. Le Van

  USE clesphys2, ONLY: cycle_diurne, iflag_con, nbapp_rad, new_oliq, &
       ok_limitvrai, ok_orodr, ok_orolf, soil_model

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: radpas
  LOGICAL, INTENT (IN) :: ok_ocean, ok_oasis, ok_journe, ok_instan, ok_region

  !--------------------------------------------------

  PRINT *, 'Choix des principales cles de la physique'
  PRINT 10, cycle_diurne, soil_model

  select case (iflag_con)
  case (1)
     PRINT *, 'Shema convection LMD'
  case (2)
     PRINT *, 'Shema convection Tiedtke'
  case (3)
     PRINT *, 'Shema convection CCM'
  END select

  PRINT 11, new_oliq, ok_orodr, ok_orolf
  PRINT 7, ok_limitvrai
  PRINT 12, nbapp_rad
  PRINT 8, radpas
  PRINT 5, ok_ocean, ok_oasis
  PRINT 4, ok_journe, ok_instan, ok_region

4 FORMAT ('ok_journe= ', L3, ', ok_instan = ', L3, ', ok_region = ', L3)
5 FORMAT ('ok_ocean = ', L3, ', ok_oasis = ', L3)
7 FORMAT ('ok_limitvrai   = ', L3)
8 FORMAT ('radpas = ', I4)
10 FORMAT ('Cycle_diurne = ', L3, ', Soil_model = ', L3)
11 FORMAT ('new_oliq = ', L3, ', Ok_orodr = ', L3, ', Ok_orolf = ', L3)
12 FORMAT ('Nb d appels /jour des routines de rayonn. = ', I4)

END SUBROUTINE printflag
