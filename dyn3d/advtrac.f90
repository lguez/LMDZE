module advtrac_m

  IMPLICIT NONE

contains

  SUBROUTINE advtrac(pbaru, pbarv, p3d, masse, q, iapptrac, teta, pk)

    ! From dyn3d/advtrac.F, version 1.4 2005/04/13 08:58:34
    ! Author: F. Hourdin

    USE comconst, ONLY : dtvr
    USE conf_gcm_m, ONLY : iapp_tracvl
    USE dimensions, ONLY : jjm, llm, nqmx
    use groupe_m, only: groupe
    USE infotrac_init_m, ONLY : iadv
    use massbar_m, only: massbar
    USE paramet_m, ONLY : iip1, iip2, ijmllm, ijp1llm, ip1jm, ip1jmp1, llmp1
    use vlsplt_m, only: vlsplt
    use vlspltqs_m, only: vlspltqs

    REAL, intent(in):: pbaru(ip1jmp1, llm), pbarv(ip1jm, llm)
    REAL, intent(in):: p3d(:, :, :) ! (iim + 1, jjm + 1, llm + 1)
    real, intent(in):: masse(ip1jmp1, llm)
    REAL, intent(inout):: q(ip1jmp1, llm, nqmx)
    INTEGER, intent(out):: iapptrac
    real, intent(in):: teta(ip1jmp1, llm)
    REAL, intent(in):: pk(ip1jmp1, llm)

    ! Variables locales

    REAL massebx(ip1jmp1, llm), masseby(ip1jm, llm)
    REAL, save, allocatable:: pbaruc(:, :) ! (ip1jmp1, llm)
    REAL, save, allocatable:: pbarvc(:, :) ! (ip1jm, llm)
    REAL, save, allocatable:: massem(:, :) ! (ip1jmp1, llm)
    real zdp(ip1jmp1)
    REAL pbarug(ip1jmp1, llm), pbarvg(ip1jm, llm), wg(ip1jmp1, llm)

    INTEGER:: iadvtr = 0
    INTEGER ij, l, iq
    REAL zdpmin, zdpmax
    EXTERNAL minmax

    INTEGER indice, n
    ! Pas de temps adaptatif pour que CFL < 1 
    REAL dtbon
    logical:: first_call = .true.

    !-----------------------------------------------------------

    if (first_call) then
       allocate(pbaruc(ip1jmp1, llm), pbarvc(ip1jm, llm))
       allocate(massem(ip1jmp1, llm))
       first_call = .false.
    end if

    IF (iadvtr==0) THEN
       CALL initial0(ijp1llm, pbaruc)
       CALL initial0(ijmllm, pbarvc)
    END IF

    ! accumulation des flux de masse horizontaux
    DO l = 1, llm
       DO ij = 1, ip1jmp1
          pbaruc(ij, l) = pbaruc(ij, l) + pbaru(ij, l)
       END DO
       DO ij = 1, ip1jm
          pbarvc(ij, l) = pbarvc(ij, l) + pbarv(ij, l)
       END DO
    END DO

    ! selection de la masse instantannee des mailles avant le transport.
    IF (iadvtr==0) massem = masse

    iadvtr = iadvtr + 1
    iapptrac = iadvtr

    ! Test pour savoir si on advecte a ce pas de temps
    IF (iadvtr == iapp_tracvl) THEN
       ! traitement des flux de masse avant advection.
       ! 1. calcul de w
       ! 2. groupement des mailles pres du pole.

       CALL groupe(pbaruc, pbarvc, pbarug, pbarvg, wg)

       ! test sur l'eventuelle creation de valeurs negatives de la masse
       DO l = 1, llm - 1
          DO ij = iip2 + 1, ip1jm
             zdp(ij) = pbarug(ij-1, l) - pbarug(ij, l) - pbarvg(ij-iip1, l) + &
                  pbarvg(ij, l) + wg(ij, l+1) - wg(ij, l)
          END DO
          CALL scopy(jjm-1, zdp(iip1+iip1), iip1, zdp(iip2), iip1)
          DO ij = iip2, ip1jm
             zdp(ij) = zdp(ij)*dtvr/massem(ij, l)
          END DO

          CALL minmax(ip1jm-iip1, zdp(iip2), zdpmin, zdpmax)

          IF (max(abs(zdpmin), abs(zdpmax))>0.5) THEN
             PRINT *, 'WARNING DP/P l=', l, ' MIN:', zdpmin, ' MAX:', zdpmax
          END IF
       END DO

       ! Advection proprement dite

       ! Calcul des moyennes bas\'ees sur la masse
       CALL massbar(massem, massebx, masseby)

       ! Appel des sous programmes d'advection

       DO iq = 1, nqmx
          select case (iadv(iq))
          case (10)
             ! Schema de Van Leer I MUSCL
             CALL vlsplt(q(:, :, iq), 2., massem, wg, pbarug, pbarvg, dtvr)
          case (12)
             ! Schema de Frederic Hourdin
             ! Pas de temps adaptatif
             CALL adaptdt(dtbon, n, pbarug, massem)
             IF (n>1) THEN
                WRITE (*, *) 'WARNING horizontal dt=', dtbon, 'dtvr=', dtvr, &
                     'n=', n
             END IF
             DO indice = 1, n
                CALL advn(q(1, 1, iq), massem, wg, pbarug, pbarvg, dtbon, 1)
             END DO
          case (13)
             ! Pas de temps adaptatif
             CALL adaptdt(dtbon, n, pbarug, massem)
             IF (n>1) THEN
                WRITE (*, *) 'WARNING horizontal dt=', dtbon, 'dtvr=', dtvr, &
                     'n=', n
             END IF
             DO indice = 1, n
                CALL advn(q(1, 1, iq), massem, wg, pbarug, pbarvg, dtbon, 2)
             END DO
          case (14)
             ! Schema "pseudo amont" + test sur humidite specifique
             ! pour la vapeur d'eau. F. Codron
             CALL vlspltqs(q(1, 1, 1), 2., massem, wg, pbarug, pbarvg, dtvr, &
                  p3d, pk, teta)
          END select
       END DO

       ! on reinitialise a zero les flux de masse cumules
       iadvtr = 0
    END IF

  END SUBROUTINE advtrac

end module advtrac_m
