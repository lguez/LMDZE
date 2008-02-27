SUBROUTINE caladvtrac(q, pbaru, pbarv, p, masse, dq, teta, pk)

  ! From dyn3d/caladvtrac.F, version 1.3 2005/04/13 08:58:34

  ! Authors : F. Hourdin, P. Le Van, F. Forget, F. Codron
  ! F. Codron (10/99) : ajout humidité spécifique pour eau vapeur
  ! Schéma de Van Leer

  use dimens_m, only: iim, jjm, llm, nqmx
  use paramet_m, only: ip1jmp1
  use comconst, only: dtvr
  use conf_gcm_m, only: iapp_tracvl

  IMPLICIT NONE

  REAL pbaru(ip1jmp1, llm), pbarv((iim + 1) * jjm, llm), masse(ip1jmp1, llm)
  REAL p(ip1jmp1, llm + 1), q(ip1jmp1, llm, nqmx)

  real, intent(out):: dq(ip1jmp1, llm, 2)
  ! (n'est utilisé et dimensionné que pour l'eau vapeur et liquide)

  REAL teta(ip1jmp1, llm), pk(ip1jmp1, llm)

  ! Local:

  EXTERNAL  advtrac, qminimum
  INTEGER l, iq, iapptrac
  REAL finmasse(ip1jmp1, llm), dtvrtrac

  !------------------------------------------------

  dq(:, :, :) = q(:, :, :2) ! initialisation

  ! Advection:
  CALL advtrac(pbaru, pbarv, p, masse, q, iapptrac, teta, pk)

  IF (iapptrac == iapp_tracvl) THEN
     ! Calcul  de deltap  qu'on stocke dans finmasse
     forall (l = 1:llm) finmasse(:, l) =  p(:, l) - p(:, l+1) 

     ! On appelle "qminimum" uniquement  pour l'eau vapeur et liquide
     CALL qminimum(q, 2, finmasse)

     finmasse(:, :) = masse(:, :)
     CALL filtreg(finmasse, jjm + 1, llm, -2, 2, .TRUE., 1)

     ! Calcul de "dq" pour l'eau, pour le passer à la physique
     dtvrtrac = iapp_tracvl * dtvr
     DO iq = 1, 2
        dq(:, :, iq) = (q(:, :, iq) - dq(:, :, iq)) * finmasse(:, :) &
             /  dtvrtrac
     ENDDO
  ELSE
     dq(:, :, :)  = 0.
  ENDIF

END SUBROUTINE caladvtrac
