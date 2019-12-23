SUBROUTINE advnz(q, qh, qb, masse, w_m)

  ! Auteurs:   F.Hourdin

  ! ********************************************************************
  ! Shema  d'advection " pseudo amont " .
  ! b designe le bas et h le haut
  ! il y a une correspondance entre le b en z et le d en x
  ! ********************************************************************


  ! --------------------------------------------------------------------
  USE dimensions
  USE paramet_m
  USE comgeom
  USE conf_gcm_m
  IMPLICIT NONE



  ! Arguments:
  ! ----------
  REAL masse(ip1jmp1, llm)
  REAL w_m(ip1jmp1, llm+1)
  REAL q(ip1jmp1, llm), qb(ip1jmp1, llm), qh(ip1jmp1, llm)


  ! Local
  ! ---------

  INTEGER ij, l

  REAL new_m, zdq, zz
  REAL zsigh(ip1jmp1, llm), zsigb(ip1jmp1, llm), zsig
  REAL w_mq(ip1jmp1, llm+1)
  REAL zm, zq, zsigm, zsigp, zqm, zqp
  REAL prec
  SAVE prec

  DATA prec/1.E-13/

  DO l = 1, llm
    DO ij = 1, ip1jmp1
      zdq = qb(ij, l) - qh(ij, l)
      IF (abs(zdq)>prec) THEN
        zsigb(ij, l) = (q(ij,l)-qh(ij,l))/zdq
        zsigh(ij, l) = 1. - zsigb(ij, l)
        zsigb(ij, l) = min(max(zsigb(ij,l),0.), 1.)
      ELSE
        zsigb(ij, l) = 0.5
        zsigh(ij, l) = 0.5
      END IF
    END DO
  END DO

  ! calcul de la pente maximum dans la maille en valeur absolue
  DO l = 2, llm
    DO ij = 1, ip1jmp1
      IF (w_m(ij,l)>=0.) THEN
        zsigp = zsigb(ij, l)
        zsigm = zsigh(ij, l)
        zqp = qb(ij, l)
        zqm = qh(ij, l)
        zm = masse(ij, l)
        zq = q(ij, l)
      ELSE
        zsigm = zsigb(ij, l-1)
        zsigp = zsigh(ij, l-1)
        zqm = qb(ij, l-1)
        zqp = qh(ij, l-1)
        zm = masse(ij, l-1)
        zq = q(ij, l-1)
      END IF
      zsig = abs(w_m(ij,l))/zm
      IF (zsig==0.) zsigp = 0.1
      IF (zsig<=zsigp) THEN
        w_mq(ij, l) = w_m(ij, l)*(zqp-0.5*zsig/zsigp*(zqp-zq))
      ELSE
        zz = 0.5*(zsig-zsigp)/zsigm
        w_mq(ij, l) = sign(zm, w_m(ij,l))*(0.5*(zq+zqp)*zsigp+(zsig-zsigp)*( &
          zq+zz*(zqm-zq)))
      END IF
    END DO
  END DO

  DO ij = 1, ip1jmp1
    w_mq(ij, llm+1) = 0.
    w_mq(ij, 1) = 0.
  END DO

  DO l = 1, llm
    DO ij = 1, ip1jmp1
      new_m = masse(ij, l) + w_m(ij, l+1) - w_m(ij, l)
      q(ij, l) = (q(ij,l)*masse(ij,l)+w_mq(ij,l+1)-w_mq(ij,l))/new_m
      masse(ij, l) = new_m
    END DO
  END DO

END SUBROUTINE advnz
