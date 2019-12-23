SUBROUTINE advn(q, masse, w, pbaru, pbarv, pdt, mode)

  ! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/advn.F,v 1.1.1.1 2004/05/19
  ! 12:53:06 lmdzadmin Exp $

  ! Auteur : F. Hourdin

  ! ********************************************************************
  ! Shema  d'advection " pseudo amont " .
  ! ********************************************************************
  ! q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....

  ! pbaru,pbarv,w flux de masse en u ,v ,w
  ! pdt pas de temps

  ! --------------------------------------------------------------------
  USE dimensions
  USE paramet_m
  USE comconst
  USE disvert_m
  USE conf_gcm_m
  USE comgeom
  IMPLICIT NONE



  ! Arguments:
  ! ----------
  INTEGER mode
  REAL masse(ip1jmp1, llm)
  REAL, INTENT (IN) :: pbaru(ip1jmp1, llm), pbarv(ip1jm, llm)
  REAL q(ip1jmp1, llm)
  REAL w(ip1jmp1, llm), pdt

  ! Local
  ! ---------

  INTEGER ij, l
  REAL zm(ip1jmp1, llm)
  REAL mu(ip1jmp1, llm)
  REAL mv(ip1jm, llm)
  REAL mw(ip1jmp1, llm+1)
  REAL zq(ip1jmp1, llm), qpn, qps
  REAL zqg(ip1jmp1, llm), zqd(ip1jmp1, llm)
  REAL zqs(ip1jmp1, llm), zqn(ip1jmp1, llm)
  REAL zqh(ip1jmp1, llm), zqb(ip1jmp1, llm)
  REAL ssum
  REAL zzpbar, zzw

  zzpbar = 0.5*pdt
  zzw = pdt

  DO l = 1, llm
    DO ij = iip2, ip1jm
      mu(ij, l) = pbaru(ij, l)*zzpbar
    END DO
    DO ij = 1, ip1jm
      mv(ij, l) = pbarv(ij, l)*zzpbar
    END DO
    DO ij = 1, ip1jmp1
      mw(ij, l) = w(ij, l)*zzw
    END DO
  END DO

  DO ij = 1, ip1jmp1
    mw(ij, llm+1) = 0.
  END DO

  DO l = 1, llm
    qpn = 0.
    qps = 0.
    DO ij = 1, iim
      qpn = qpn + q(ij, l)*masse(ij, l)
      qps = qps + q(ip1jm+ij, l)*masse(ip1jm+ij, l)
    END DO
    qpn = qpn/ssum(iim, masse(1,l), 1)
    qps = qps/ssum(iim, masse(ip1jm+1,l), 1)
    DO ij = 1, iip1
      q(ij, l) = qpn
      q(ip1jm+ij, l) = qps
    END DO
  END DO

  DO ij = 1, ip1jmp1
    mw(ij, llm+1) = 0.
  END DO
  DO l = 1, llm
    DO ij = 1, ip1jmp1
      zq(ij, l) = q(ij, l)
      zm(ij, l) = masse(ij, l)
    END DO
  END DO

  ! call minmaxq(zq,qmin,qmax,'avant vlx     ')
  CALL advnqx(zq, zqg, zqd)
  CALL advnx(zq, zqg, zqd, zm, mu, mode)
  CALL advnqy(zq, zqs, zqn)
  CALL advny(zq, zqs, zqn, zm, mv)
  CALL advnqz(zq, zqh, zqb)
  CALL advnz(zq, zqh, zqb, zm, mw)
  ! call vlz(zq,0.,zm,mw)
  CALL advnqy(zq, zqs, zqn)
  CALL advny(zq, zqs, zqn, zm, mv)
  CALL advnqx(zq, zqg, zqd)
  CALL advnx(zq, zqg, zqd, zm, mu, mode)
  ! call minmaxq(zq,qmin,qmax,'apres vlx     ')

  DO l = 1, llm
    DO ij = 1, ip1jmp1
      q(ij, l) = zq(ij, l)
    END DO
    DO ij = 1, ip1jm + 1, iip1
      q(ij+iim, l) = q(ij, l)
    END DO
  END DO

  RETURN
END SUBROUTINE advn

