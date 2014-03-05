SUBROUTINE vlsplt(q, pente_max, masse, w, pbaru, pbarv, pdt)

  ! From LMDZ4/libf/dyn3d/vlsplt.F, version 1.2 2005/02/24 12:16:57 fairhead
  ! Auteurs: P. Le Van, F. Hourdin, F. Forget

  ! Shema d'advection "pseudo amont".
  ! q, pbaru, pbarv, w sont des arguments d'entree pour le sous-programme.
  ! pente_max facteur de limitation des pentes : 2 en general
  ! 0 pour un schema amont
  ! pbaru, pbarv, w flux de masse en u, v, w
  ! pdt pas de temps

  USE dimens_m, ONLY: iim, llm
  USE paramet_m, ONLY: iip1, iip2, ijp1llm, ip1jm, ip1jmp1

  IMPLICIT NONE

  ! Arguments:

  REAL masse(ip1jmp1, llm), pente_max
  REAL, intent(in):: pbaru( ip1jmp1, llm ), pbarv( ip1jm, llm)
  REAL, intent(inout):: q(ip1jmp1, llm)
  REAL w(ip1jmp1, llm)
  real, intent(in):: pdt

  ! Local

  INTEGER i, ij, l, j, ii
  INTEGER ijlqmin, iqmin, jqmin, lqmin

  REAL zm(ip1jmp1, llm), newmasse
  REAL mu(ip1jmp1, llm)
  REAL mv(ip1jm, llm)
  REAL mw(ip1jmp1, llm+1)
  REAL zq(ip1jmp1, llm), zz
  REAL dqx(ip1jmp1, llm), dqy(ip1jmp1, llm), dqz(ip1jmp1, llm)
  REAL second, temps0, temps1, temps2, temps3
  REAL ztemps1, ztemps2, ztemps3
  REAL zzpbar, zzw
  LOGICAL testcpu
  SAVE testcpu
  SAVE temps1, temps2, temps3
  INTEGER iminn, imaxx

  REAL qmin, qmax
  DATA qmin, qmax/0., 1.e33/
  DATA testcpu/.false./
  DATA temps1, temps2, temps3/0., 0., 0./

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

  zq = q
  zm = masse

  call vlx(zq, pente_max, zm, mu)

  call vly(zq, pente_max, zm, mv)
  call vlz(zq, pente_max, zm, mw)

  call vly(zq, pente_max, zm, mv)

  call vlx(zq, pente_max, zm, mu)

  q = zq
  DO ij = 1, ip1jm+1, iip1
     q(ij+iim, :) = q(ij, :)
  ENDDO

END SUBROUTINE vlsplt
