module addfi_m

  IMPLICIT NONE

contains

  SUBROUTINE addfi(ucov, vcov, teta, q, ps, dufi, dvfi, dtetafi, dqfi, dpfi)

    ! From dyn3d/addfi.F, v 1.1.1.1 2004/05/19 12:53:06

    ! Addition of the physical tendencies

    USE comconst, ONLY: dtphys
    USE comgeom, ONLY: aire, apoln, apols
    USE dimens_m, ONLY: iim, jjm, llm, nqmx

    ! First and second components of the covariant velocity:
    REAL, intent(inout):: ucov((iim + 1) * (jjm + 1), llm)
    REAL, intent(inout):: vcov((iim + 1) * jjm, llm)

    REAL, intent(inout):: teta((iim + 1) * (jjm + 1), llm)
    ! potential temperature

    real, intent(inout):: q((iim + 1) * (jjm + 1), llm, nqmx)
    real, intent(inout):: ps((iim + 1) * (jjm + 1))

    ! Tendencies:
    REAL, intent(in):: dufi((iim + 1) * (jjm + 1), llm)
    REAL, intent(in):: dvfi((iim + 1) * jjm, llm)
    real, intent(in):: dtetafi((iim + 1) * (jjm + 1), llm)
    REAL, intent(in):: dqfi((iim + 1) * (jjm + 1), llm, nqmx)
    REAL, intent(in):: dpfi((iim + 1) * (jjm + 1))

    ! Local variables :
    REAL xpn(iim), xps(iim), tpn, tps
    INTEGER j, k, iq, ij
    REAL, PARAMETER:: qtestw = 1e-15, qtestt = 1e-40

    !-----------------------------------------------------------------------

    teta = teta + dtetafi * dtphys

    DO k = 1, llm
       DO ij = 1, iim
          xpn(ij) = aire(ij) * teta(ij , k)
          xps(ij) = aire(ij+(iim + 1) * jjm) * teta(ij+(iim + 1) * jjm, k)
       ENDDO
       tpn = SUM(xpn)/ apoln
       tps = SUM(xps)/ apols

       DO ij = 1, iim + 1
          teta(ij , k) = tpn
          teta(ij+(iim + 1) * jjm, k) = tps
       ENDDO
    ENDDO

    DO k = 1, llm
       DO j = iim + 2, (iim + 1) * jjm
          ucov(j, k)= ucov(j, k) + dufi(j, k) * dtphys
       ENDDO
    ENDDO

    vcov = vcov + dvfi * dtphys
    ps = ps + dpfi * dtphys

    DO iq = 1, 2
       DO k = 1, llm
          DO j = 1, (iim + 1) * (jjm + 1)
             q(j, k, iq)= q(j, k, iq) + dqfi(j, k, iq) * dtphys
             q(j, k, iq)= MAX(q(j, k, iq), qtestw)
          ENDDO
       ENDDO
    ENDDO

    DO iq = 3, nqmx
       DO k = 1, llm
          DO j = 1, (iim + 1) * (jjm + 1)
             q(j, k, iq)= q(j, k, iq) + dqfi(j, k, iq) * dtphys
             q(j, k, iq)= MAX(q(j, k, iq), qtestt)
          ENDDO
       ENDDO
    ENDDO

    DO ij = 1, iim
       xpn(ij) = aire(ij) * ps(ij)
       xps(ij) = aire(ij+(iim + 1) * jjm) * ps(ij+(iim + 1) * jjm)
    ENDDO
    tpn = SUM(xpn)/apoln
    tps = SUM(xps)/apols

    DO ij = 1, iim + 1
       ps(ij) = tpn
       ps(ij+(iim + 1) * jjm) = tps
    ENDDO

    DO iq = 1, nqmx
       DO k = 1, llm
          DO ij = 1, iim
             xpn(ij) = aire(ij) * q(ij , k, iq)
             xps(ij) = aire(ij+(iim + 1) * jjm) * q(ij+(iim + 1) * jjm, k, iq)
          ENDDO
          tpn = SUM(xpn)/apoln
          tps = SUM(xps)/apols

          DO ij = 1, iim + 1
             q(ij , k, iq) = tpn
             q(ij+(iim + 1) * jjm, k, iq) = tps
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE addfi

end module addfi_m
