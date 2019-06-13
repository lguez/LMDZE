module qminimum_m

  IMPLICIT none

contains

  SUBROUTINE qminimum(q, nq, deltap)

    ! From LMDZ4/libf/dyn3d/qminimum.F, version 1.1.1.1 2004/05/19 12:53:05

    ! Objet : traiter les valeurs trop petites (même négatives) pour
    ! l'eau vapeur et l'eau liquide

    USE dimensions, ONLY: llm
    USE paramet_m, ONLY: ip1jmp1

    INTEGER, intent(in):: nq
    REAL, intent(inout):: q(ip1jmp1, llm, nq)
    real, intent(in):: deltap(ip1jmp1, llm)

    ! Local:

    INTEGER, PARAMETER:: iq_vap = 1 ! indice pour l'eau vapeur
    INTEGER, PARAMETER:: iq_liq = 2 ! indice pour l'eau liquide

    ! Il est souhaitable mais non obligatoire que les valeurs des
    ! paramètres seuil_vap, seuil_liq soient identiques à celles de la
    ! routine ADDFI.
    REAL, PARAMETER:: seuil_vap = 1e-10 ! seuil pour l'eau vapeur
    REAL, PARAMETER:: seuil_liq = 1e-11 ! seuil pour l'eau liquide

    INTEGER i, k
    REAL zx_defau, zx_abc, zx_pump(ip1jmp1), pompe
    INTEGER:: imprim = 0

    !-------------------------------------------------------------------

    ! Quand l'eau liquide est trop petite (ou négative), on prend
    ! l'eau vapeur de la même couche et on la convertit en eau liquide
    ! (sans changer la température !).
    DO k = 1, llm
       DO i = 1, ip1jmp1
          zx_defau = MAX(seuil_liq - q(i, k, iq_liq), 0.0)
          q(i, k, iq_vap) = q(i, k, iq_vap) - zx_defau
          q(i, k, iq_liq) = q(i, k, iq_liq) + zx_defau
       end DO
    end DO

    ! Quand l'eau vapeur est trop faible (ou négative), on complète
    ! le défaut en prenant de l'eau vapeur de la couche au-dessous.
    DO k = llm, 2, -1
       DO i = 1, ip1jmp1
          zx_abc = deltap(i, k)/deltap(i, k-1)
          zx_defau = MAX(seuil_vap - q(i, k, iq_vap), 0.0)
          q(i, k-1, iq_vap) = q(i, k-1, iq_vap) - zx_defau * zx_abc
          q(i, k, iq_vap) = q(i, k, iq_vap) + zx_defau 
       ENDDO
    ENDDO

    ! Quand il s'agit de la première couche au dessus du sol, on doit
    ! imprimer un message d'avertissement (saturation possible).

    DO i = 1, ip1jmp1
       zx_pump(i) = MAX(0., seuil_vap - q(i, 1, iq_vap))
       q(i, 1, iq_vap) = MAX(q(i, 1, iq_vap), seuil_vap)
    ENDDO
    pompe = SUM(zx_pump)
    IF (imprim <= 500 .AND. pompe > 0.) THEN
       print *, "Attention : on pompe de l'eau au sol, pompe = ", pompe
       DO i = 1, ip1jmp1
          IF (zx_pump(i) > 0.) THEN
             imprim = imprim + 1
          ENDIF
       ENDDO
    ENDIF

  END SUBROUTINE qminimum

end module qminimum_m
