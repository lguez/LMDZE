module caladvtrac_m

  IMPLICIT NONE

contains

  SUBROUTINE caladvtrac(q, pbaru, pbarv, p, masse, dq, teta, pk)

    ! From dyn3d/caladvtrac.F, version 1.3 2005/04/13 08:58:34

    ! Authors : F. Hourdin, P. Le Van, F. Forget, F. Codron
    ! F. Codron (10/99) : ajout humidité spécifique pour eau vapeur
    ! Schéma de Van Leer

    use advtrac_m, only: advtrac
    use comconst, only: dtvr
    use conf_gcm_m, only: iapp_tracvl
    use dimens_m, only: iim, jjm, llm, nqmx
    use filtreg_m, only: filtreg
    use paramet_m, only: ip1jmp1

    REAL pbaru(ip1jmp1, llm), pbarv((iim + 1) * jjm, llm)
    real masse(iim + 1, jjm + 1, llm)
    REAL, intent(in):: p(iim + 1, jjm + 1, llm + 1)
    real, intent(inout):: q(:, :, :, :) ! (iim + 1, jjm + 1, llm, nqmx)

    real, intent(out):: dq(iim + 1, jjm + 1, llm, 2)
    ! (n'est utilisé et dimensionné que pour l'eau vapeur et liquide)

    REAL, intent(in):: teta(ip1jmp1, llm)
    real pk(ip1jmp1, llm)

    ! Local:

    EXTERNAL qminimum
    INTEGER l, iq, iapptrac
    REAL finmasse(iim + 1, jjm + 1, llm), dtvrtrac

    !------------------------------------------------

    dq = q(:, :, :, :2) ! initialisation

    ! Advection:
    CALL advtrac(pbaru, pbarv, p, masse, q, iapptrac, teta, pk)

    IF (iapptrac == iapp_tracvl) THEN
       ! Calcul  de deltap  qu'on stocke dans finmasse
       forall (l = 1:llm) finmasse(:, :, l) =  p(:, :, l) - p(:, :, l+1) 

       ! On appelle "qminimum" uniquement  pour l'eau vapeur et liquide
       CALL qminimum(q, 2, finmasse)

       finmasse = masse
       CALL filtreg(finmasse, jjm + 1, llm, -2, 2, .TRUE.)

       ! Calcul de "dq" pour l'eau, pour le passer à la physique
       dtvrtrac = iapp_tracvl * dtvr
       DO iq = 1, 2
          dq(:, :, :, iq) = (q(:, :, :, iq) - dq(:, :, :, iq)) * finmasse &
               /  dtvrtrac
       ENDDO
    ELSE
       dq = 0.
    ENDIF

  END SUBROUTINE caladvtrac

end module caladvtrac_m
