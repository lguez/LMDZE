module vlspltqs_m

  IMPLICIT NONE

contains

  SUBROUTINE vlspltqs(q, pente_max, masse, w, pbaru, pbarv, pdt, p, pk, teta)

    ! From LMDZ4/libf/dyn3d/vlspltqs.F, version 1.2 2005/02/24 12:16:57 fairhead

    ! Authors: P. Le Van, F. Hourdin, F. Forget, F. Codron

    ! Schéma d'advection "pseudo amont"
    ! + test sur humidité spécifique : Q advecté < Qsat aval
    ! (F. Codron, 10/99)

    ! q, pbaru, pbarv, w sont des arguments d'entree pour le sous-programme

    ! pente_max facteur de limitation des pentes: 2 en général
    ! 0 pour un schéma amont
    ! pbaru, pbarv, w flux de masse en u , v , w
    ! pdt pas de temps

    ! teta température potentielle, p pression aux interfaces,
    ! pk exner au milieu des couches nécessaire pour calculer Qsat

    USE dimensions, ONLY : iim, llm
    use FCTTRE, only: foeew
    USE paramet_m, ONLY : iip1, iip2, ijp1llm, ip1jm, ip1jmp1, llmp1
    use SUPHEC_M, only: rtt, rcpd
    use vlyqs_m, only: vlyqs

    ! Arguments:

    REAL masse(ip1jmp1, llm), pente_max
    REAL, intent(in):: pbaru(ip1jmp1, llm), pbarv(ip1jm, llm)
    REAL q(ip1jmp1, llm)
    REAL w(ip1jmp1, llm)
    real, intent(in):: pdt
    REAL, intent(in):: p(ip1jmp1, llmp1)
    real, intent(in):: teta(ip1jmp1, llm)
    real, intent(in):: pk(ip1jmp1, llm)

    ! Local

    INTEGER ij, l

    REAL qsat(ip1jmp1, llm)
    REAL zm(ip1jmp1, llm)
    REAL mu(ip1jmp1, llm)
    REAL mv(ip1jm, llm)
    REAL mw(ip1jmp1, llm+1)
    REAL zq(ip1jmp1, llm)
    REAL zzpbar, zzw

    !--pour rapport de melange saturant--

    REAL retv, r2es, play
    logical zdelta
    REAL tempe(ip1jmp1)

    !------------------------------------------------------------------

    r2es = 380.11733
    retv = 0.6077667

    !-- Calcul de Qsat en chaque point
    !-- approximation: au milieu des couches play(l)=(p(l)+p(l+1))/2
    ! pour eviter une exponentielle.
    DO l = 1, llm
       DO ij = 1, ip1jmp1
          tempe(ij) = teta(ij, l) * pk(ij, l) /rcpd
       ENDDO
       DO ij = 1, ip1jmp1
          zdelta = rtt > tempe(ij)
          play = 0.5*(p(ij, l)+p(ij, l+1))
          qsat(ij, l) = MIN(0.5, r2es* FOEEW(tempe(ij), zdelta) / play)
          qsat(ij, l) = qsat(ij, l) / (1. - retv * qsat(ij, l))
       ENDDO
    ENDDO

    zzpbar = 0.5 * pdt
    zzw = pdt
    DO l=1, llm
       DO ij = iip2, ip1jm
          mu(ij, l)=pbaru(ij, l) * zzpbar
       ENDDO
       DO ij=1, ip1jm
          mv(ij, l)=pbarv(ij, l) * zzpbar
       ENDDO
       DO ij=1, ip1jmp1
          mw(ij, l)=w(ij, l) * zzw
       ENDDO
    ENDDO

    DO ij=1, ip1jmp1
       mw(ij, llm+1)=0.
    ENDDO

    CALL SCOPY(ijp1llm, q, 1, zq, 1)
    CALL SCOPY(ijp1llm, masse, 1, zm, 1)

    ! call minmaxq(zq, qmin, qmax, 'avant vlxqs ')
    call vlxqs(zq, pente_max, zm, mu, qsat)

    ! call minmaxq(zq, qmin, qmax, 'avant vlyqs ')

    call vlyqs(zq, pente_max, zm, mv, qsat)

    ! call minmaxq(zq, qmin, qmax, 'avant vlz ')

    call vlz(zq, pente_max, zm, mw)

    ! call minmaxq(zq, qmin, qmax, 'avant vlyqs ')
    ! call minmaxq(zm, qmin, qmax, 'M avant vlyqs ')

    call vlyqs(zq, pente_max, zm, mv, qsat)

    ! call minmaxq(zq, qmin, qmax, 'avant vlxqs ')
    ! call minmaxq(zm, qmin, qmax, 'M avant vlxqs ')

    call vlxqs(zq, pente_max, zm, mu, qsat)

    ! call minmaxq(zq, qmin, qmax, 'apres vlxqs ')
    ! call minmaxq(zm, qmin, qmax, 'M apres vlxqs ')

    DO l=1, llm
       DO ij=1, ip1jmp1
          q(ij, l)=zq(ij, l)
       ENDDO
       DO ij=1, ip1jm+1, iip1
          q(ij+iim, l)=q(ij, l)
       ENDDO
    ENDDO

  END SUBROUTINE vlspltqs

end module vlspltqs_m
