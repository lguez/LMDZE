module caladvtrac_m

  IMPLICIT NONE

contains

  SUBROUTINE caladvtrac(q, pbaru, pbarv, p3d, masse, teta, pk)

    ! From dyn3d/caladvtrac.F, version 1.3 2005/04/13 08:58:34
    ! Authors: F. Hourdin, P. Le Van, F. Forget, F. Codron
    ! F. Codron (10/99) : ajout humidit\'e sp\'ecifique pour eau vapeur
    ! Sch\'ema de Van Leer

    ! Calcul des tendances advection des traceurs (dont l'humidit\'e)

    use advtrac_m, only: advtrac
    use conf_gcm_m, only: iapp_tracvl
    use dimensions, only: iim, jjm, llm
    use qminimum_m, only: qminimum

    real, intent(inout):: q(:, :, :, :) ! (iim + 1, jjm + 1, llm, nqmx)
    REAL, intent(in):: pbaru(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(in):: pbarv(:, :, :) ! (iim + 1, jjm, llm)
    REAL, intent(in):: p3d(:, :, :) ! (iim + 1, jjm + 1, llm + 1)
    real, intent(in), dimension(:, :, :):: masse, teta, pk
    ! (iim + 1, jjm + 1, llm)

    ! Local:
    INTEGER l, iapptrac
    REAL finmasse(iim + 1, jjm + 1, llm)

    !------------------------------------------------

    ! Advection:
    CALL advtrac(pbaru, pbarv, p3d, masse, q, iapptrac, teta, pk)

    IF (iapptrac == iapp_tracvl) THEN
       forall (l = 1:llm) finmasse(:, :, l) =  p3d(:, :, l) - p3d(:, :, l+1) 

       ! Uniquement pour l'eau vapeur et liquide:
       CALL qminimum(q, 2, finmasse)
    ENDIF

  END SUBROUTINE caladvtrac

end module caladvtrac_m
