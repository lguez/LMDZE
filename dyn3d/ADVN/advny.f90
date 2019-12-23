SUBROUTINE advny(q, qs, qn, masse, v_m)

  ! Auteur : F. Hourdin

  ! ********************************************************************
  ! Shema  d'advection " pseudo amont " .
  ! ********************************************************************
  ! nq,iq,q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....


  ! --------------------------------------------------------------------
  USE dimensions
  USE paramet_m
  USE comgeom
  USE conf_gcm_m
  IMPLICIT NONE



  ! Arguments:
  ! ----------
  REAL masse(ip1jmp1, llm)
  REAL v_m(ip1jm, llm)
  REAL q(ip1jmp1, llm), qn(ip1jmp1, llm), qs(ip1jmp1, llm)

  ! Local
  ! ---------

  INTEGER ij, l

  REAL new_m, zdq, zz
  REAL zsigs(ip1jmp1), zsign(ip1jmp1), zsig
  REAL v_mq(ip1jm, llm)
  REAL convpn, convps, convmpn, convmps, massen, masses
  REAL zm, zq, zsigm, zsigp, zqm, zqp
  REAL ssum
  REAL prec
  SAVE prec

  DATA prec/1.E-15/

  DO l = 1, llm
    DO ij = 1, ip1jmp1
      zdq = qn(ij, l) - qs(ij, l)
      IF (abs(zdq)>prec) THEN
        zsign(ij) = (q(ij,l)-qs(ij,l))/zdq
        zsigs(ij) = 1. - zsign(ij)
      ELSE
        zsign(ij) = 0.5
        zsigs(ij) = 0.5
      END IF
    END DO

    ! calcul de la pente maximum dans la maille en valeur absolue

    DO ij = 1, ip1jm
      IF (v_m(ij,l)>=0.) THEN
        zsigp = zsign(ij+iip1)
        zsigm = zsigs(ij+iip1)
        zqp = qn(ij+iip1, l)
        zqm = qs(ij+iip1, l)
        zm = masse(ij+iip1, l)
        zq = q(ij+iip1, l)
      ELSE
        zsigm = zsign(ij)
        zsigp = zsigs(ij)
        zqm = qn(ij, l)
        zqp = qs(ij, l)
        zm = masse(ij, l)
        zq = q(ij, l)
      END IF
      zsig = abs(v_m(ij,l))/zm
      IF (zsig==0.) zsigp = 0.1
      IF (zsig<=zsigp) THEN
        v_mq(ij, l) = v_m(ij, l)*(zqp-0.5*zsig/zsigp*(zqp-zq))
      ELSE
        zz = 0.5*(zsig-zsigp)/zsigm
        v_mq(ij, l) = sign(zm, v_m(ij,l))*(0.5*(zq+zqp)*zsigp+(zsig-zsigp)*( &
          zq+zz*(zqm-zq)))
      END IF
    END DO
  END DO

  DO l = 1, llm
    DO ij = iip2, ip1jm
      new_m = masse(ij, l) + v_m(ij, l) - v_m(ij-iip1, l)
      q(ij, l) = (q(ij,l)*masse(ij,l)+v_mq(ij,l)-v_mq(ij-iip1,l))/new_m
      masse(ij, l) = new_m
    END DO
    ! .-. ancienne version
    convpn = ssum(iim, v_mq(1,l), 1)
    convmpn = ssum(iim, v_m(1,l), 1)
    massen = ssum(iim, masse(1,l), 1)
    new_m = massen + convmpn
    q(1, l) = (q(1,l)*massen+convpn)/new_m
    DO ij = 1, iip1
      q(ij, l) = q(1, l)
      masse(ij, l) = new_m*aire(ij)/apoln
    END DO

    convps = -ssum(iim, v_mq(ip1jm-iim,l), 1)
    convmps = -ssum(iim, v_m(ip1jm-iim,l), 1)
    masses = ssum(iim, masse(ip1jm+1,l), 1)
    new_m = masses + convmps
    q(ip1jm+1, l) = (q(ip1jm+1,l)*masses+convps)/new_m
    DO ij = ip1jm + 1, ip1jmp1
      q(ij, l) = q(ip1jm+1, l)
      masse(ij, l) = new_m*aire(ij)/apols
    END DO
  END DO

  RETURN
END SUBROUTINE advny
