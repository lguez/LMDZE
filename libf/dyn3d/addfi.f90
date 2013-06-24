module addfi_m

  IMPLICIT NONE

contains

  SUBROUTINE addfi(nq, ucov, vcov, teta, pq, pps, dufi, dvfi, pdhfi, &
       pdqfi, pdpfi)

    ! From dyn3d/addfi.F, v 1.1.1.1 2004/05/19 12:53:06

    ! Addition of the physical tendencies

    USE comconst, ONLY: dtphys
    USE comgeom, ONLY: aire, apoln, apols
    USE dimens_m, ONLY: iim, llm
    USE paramet_m, ONLY: iip1, iip2, ip1jm, ip1jmp1

    INTEGER, intent(in):: nq

    REAL, intent(inout):: ucov(ip1jmp1, llm), vcov(ip1jm, llm)
    ! first and second components of the covariant velocity

    REAL, intent(inout):: teta(ip1jmp1, llm) ! potential temperature
    real, intent(inout):: pq(ip1jmp1, llm, nq), pps(ip1jmp1)
    REAL, intent(in):: dufi(ip1jmp1, llm), dvfi(ip1jm, llm) ! tendencies
    real, intent(in):: pdhfi(ip1jmp1, llm) ! tendency
    REAL, intent(in):: pdqfi(ip1jmp1, llm, nq), pdpfi(ip1jmp1)

    ! Local variables :
    REAL xpn(iim), xps(iim), tpn, tps
    INTEGER j, k, iq, ij
    REAL, PARAMETER:: qtestw = 1e-15, qtestt = 1e-40

    !-----------------------------------------------------------------------

    DO k = 1, llm
       DO j = 1, ip1jmp1
          teta(j, k)= teta(j, k) + pdhfi(j, k) * dtphys
       ENDDO
    ENDDO

    DO k = 1, llm
       DO ij = 1, iim
          xpn(ij) = aire(ij) * teta(ij , k)
          xps(ij) = aire(ij+ip1jm) * teta(ij+ip1jm, k)
       ENDDO
       tpn = SUM(xpn)/ apoln
       tps = SUM(xps)/ apols

       DO ij = 1, iip1
          teta(ij , k) = tpn
          teta(ij+ip1jm, k) = tps
       ENDDO
    ENDDO

    DO k = 1, llm
       DO j = iip2, ip1jm
          ucov(j, k)= ucov(j, k) + dufi(j, k) * dtphys
       ENDDO
    ENDDO

    DO k = 1, llm
       DO j = 1, ip1jm
          vcov(j, k)= vcov(j, k) + dvfi(j, k) * dtphys
       ENDDO
    ENDDO

    DO j = 1, ip1jmp1
       pps(j) = pps(j) + pdpfi(j) * dtphys
    ENDDO

    DO iq = 1, 2
       DO k = 1, llm
          DO j = 1, ip1jmp1
             pq(j, k, iq)= pq(j, k, iq) + pdqfi(j, k, iq) * dtphys
             pq(j, k, iq)= MAX(pq(j, k, iq), qtestw)
          ENDDO
       ENDDO
    ENDDO

    DO iq = 3, nq
       DO k = 1, llm
          DO j = 1, ip1jmp1
             pq(j, k, iq)= pq(j, k, iq) + pdqfi(j, k, iq) * dtphys
             pq(j, k, iq)= MAX(pq(j, k, iq), qtestt)
          ENDDO
       ENDDO
    ENDDO

    DO ij = 1, iim
       xpn(ij) = aire(ij) * pps(ij)
       xps(ij) = aire(ij+ip1jm) * pps(ij+ip1jm)
    ENDDO
    tpn = SUM(xpn)/apoln
    tps = SUM(xps)/apols

    DO ij = 1, iip1
       pps (ij) = tpn
       pps (ij+ip1jm) = tps
    ENDDO

    DO iq = 1, nq
       DO k = 1, llm
          DO ij = 1, iim
             xpn(ij) = aire(ij) * pq(ij , k, iq)
             xps(ij) = aire(ij+ip1jm) * pq(ij+ip1jm, k, iq)
          ENDDO
          tpn = SUM(xpn)/apoln
          tps = SUM(xps)/apols

          DO ij = 1, iip1
             pq (ij , k, iq) = tpn
             pq (ij+ip1jm, k, iq) = tps
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE addfi

end module addfi_m
