module fluxstokenc_m

  IMPLICIT NONE

contains

  SUBROUTINE fluxstokenc(pbaru, pbarv, masse, teta, phi, phis, time_step, itau)

    ! Author: F. Hourdin

    USE histwrite_m, ONLY: histwrite
    USE dimens_m, ONLY: jjm, llm, nqmx
    USE paramet_m, ONLY: iip1, ijmllm, ijp1llm, ip1jm, ip1jmp1, jjp1
    USE comgeom, ONLY: aire
    USE tracstoke, ONLY: istdyn, istphy

    REAL pbaru(ip1jmp1, llm), pbarv(ip1jm, llm)
    REAL masse(ip1jmp1, llm)
    real, intent(in):: phi(ip1jmp1, llm)
    real, intent(in):: teta(ip1jmp1, llm)
    REAL, intent(in):: phis(ip1jmp1)
    REAL, intent(in):: time_step
    INTEGER, INTENT (IN):: itau

    ! Variables local to the procedure:

    REAL, SAVE:: pbaruc(ip1jmp1, llm), pbarvc(ip1jm, llm)
    REAL, SAVE:: massem(ip1jmp1, llm)
    real, SAVE:: tetac(ip1jmp1, llm), phic(ip1jmp1, llm)

    REAL pbarug(ip1jmp1, llm), pbarvg(iip1, jjm, llm), wg(ip1jmp1, llm)
    REAL tst(1), ist(1), istp(1)
    INTEGER ij, l
    INTEGER, save:: fluxid, fluxvid
    integer fluxdid

    !-------------------------------------------------------------

    IF (itau == 0) THEN
       CALL initfluxsto(time_step, istdyn*time_step, istdyn*time_step, nqmx, &
            fluxid, fluxvid, fluxdid)
       CALL histwrite(fluxid, 'phis', 1, phis)
       CALL histwrite(fluxid, 'aire', 1, aire)
       tst(1) = time_step
       CALL histwrite(fluxdid, 'dtvr', 1, tst)
       ist(1) = istdyn
       CALL histwrite(fluxdid, 'istdyn', 1, ist)
       istp(1) = istphy
       CALL histwrite(fluxdid, 'istphy', 1, istp)

       CALL initial0(ijp1llm, phic)
       CALL initial0(ijp1llm, tetac)
       CALL initial0(ijp1llm, pbaruc)
       CALL initial0(ijmllm, pbarvc)
    END IF

    !   accumulation des flux de masse horizontaux
    DO l = 1, llm
       DO ij = 1, ip1jmp1
          pbaruc(ij, l) = pbaruc(ij, l) + pbaru(ij, l)
          tetac(ij, l) = tetac(ij, l) + teta(ij, l)
          phic(ij, l) = phic(ij, l) + phi(ij, l)
       END DO
       DO ij = 1, ip1jm
          pbarvc(ij, l) = pbarvc(ij, l) + pbarv(ij, l)
       END DO
    END DO

    !   selection de la masse instantannee des mailles avant le transport.
    IF (itau == 0) THEN
       CALL scopy(ip1jmp1*llm, masse, 1, massem, 1)
    END IF

    IF (mod(itau + 1, istdyn) == 0) THEN
       ! on advecte a ce pas de temps
       !    normalisation
       DO l = 1, llm
          DO ij = 1, ip1jmp1
             pbaruc(ij, l) = pbaruc(ij, l)/float(istdyn)
             tetac(ij, l) = tetac(ij, l)/float(istdyn)
             phic(ij, l) = phic(ij, l)/float(istdyn)
          END DO
          DO ij = 1, ip1jm
             pbarvc(ij, l) = pbarvc(ij, l)/float(istdyn)
          END DO
       END DO

       !   traitement des flux de masse avant advection.
       !     1. calcul de w
       !     2. groupement des mailles pres du pole.

       CALL groupe(massem, pbaruc, pbarvc, pbarug, pbarvg, wg)

       CALL histwrite(fluxid, 'masse', itau, massem)
       CALL histwrite(fluxid, 'pbaru', itau, pbarug)
       CALL histwrite(fluxvid, 'pbarv', itau, pbarvg)
       CALL histwrite(fluxid, 'w', itau, wg)
       CALL histwrite(fluxid, 'teta', itau, tetac)
       CALL histwrite(fluxid, 'phi', itau, phic)
    END IF

  END SUBROUTINE fluxstokenc

end module fluxstokenc_m
