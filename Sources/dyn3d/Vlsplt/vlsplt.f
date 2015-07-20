module vlsplt_m

  IMPLICIT NONE

contains

  SUBROUTINE vlsplt(q, pente_max, masse, w, pbaru, pbarv, pdt)

    ! From LMDZ4/libf/dyn3d/vlsplt.F, version 1.2 2005/02/24 12:16:57 fairhead

    ! Authors: P. Le Van, F. Hourdin, F. Forget

    ! Sch\'ema d'advection "pseudo-amont".

    USE dimens_m, ONLY: iim, llm
    USE paramet_m, ONLY: iip1, iip2, ijp1llm, ip1jm, ip1jmp1
    use vlx_m, only: vlx

    REAL, intent(inout):: q(ip1jmp1, llm)

    REAL, intent(in):: pente_max
    ! facteur de limitation des pentes, 2 en general

    real, intent(in):: masse(ip1jmp1, llm)
    REAL, intent(in):: w(ip1jmp1, llm) ! flux de masse

    REAL, intent(in):: pbaru( ip1jmp1, llm ), pbarv( ip1jm, llm)
    ! flux de masse en u, v

    real, intent(in):: pdt ! pas de temps

    ! Local:

    INTEGER ij, l
    REAL zm(ip1jmp1, llm)
    REAL mu(ip1jmp1, llm)
    REAL mv(ip1jm, llm)
    REAL mw(ip1jmp1, llm+1)
    REAL zzpbar, zzw

    !---------------------------------------------------------------

    zzpbar = 0.5 * pdt
    zzw = pdt
    DO l = 1, llm
       DO ij = iip2, ip1jm
          mu(ij, l) = pbaru(ij, l) * zzpbar
       ENDDO
       DO ij = 1, ip1jm
          mv(ij, l) = pbarv(ij, l) * zzpbar
       ENDDO
       DO ij = 1, ip1jmp1
          mw(ij, l) = w(ij, l) * zzw
       ENDDO
    ENDDO

    DO ij = 1, ip1jmp1
       mw(ij, llm+1) = 0.
    ENDDO

    zm = masse

    call vlx(q, pente_max, zm, mu)
    call vly(q, pente_max, zm, mv)
    call vlz(q, pente_max, zm, mw)
    call vly(q, pente_max, zm, mv)
    call vlx(q, pente_max, zm, mu)

    DO ij = 1, ip1jm + 1, iip1
       q(ij + iim, :) = q(ij, :)
    ENDDO

  END SUBROUTINE vlsplt

end module vlsplt_m
