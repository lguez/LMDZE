
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/limy.F,v 1.1.1.1 2004/05/19
! 12:53:07 lmdzadmin Exp $

SUBROUTINE limy(s0, sy, sm, pente_max)

  ! Auteurs:   P.Le Van, F.Hourdin, F.Forget

  ! ********************************************************************
  ! Shema  d'advection " pseudo amont " .
  ! ********************************************************************
  ! q,w sont des arguments d'entree  pour le s-pg ....
  ! dq 	       sont des arguments de sortie pour le s-pg ....


  ! --------------------------------------------------------------------
  USE comconst
  use comgeom, only: aire
  USE conf_gcm_m
  USE dimens_m
  USE disvert_m
  USE dynetat0_m, only: rlonv, rlonu
  USE nr_util, ONLY: pi
  USE paramet_m

  IMPLICIT NONE



  ! Arguments:
  ! ----------
  REAL pente_max
  REAL s0(ip1jmp1, llm), sy(ip1jmp1, llm), sm(ip1jmp1, llm)

  ! Local
  ! ---------

  INTEGER i, ij, l

  REAL q(ip1jmp1, llm)
  REAL airej2, airejjm, airescb(iim), airesch(iim)
  REAL sigv, dyq(ip1jmp1), dyqv(ip1jm)
  REAL adyqv(ip1jm), dyqmax(ip1jmp1)
  REAL qbyv(ip1jm, llm)

  REAL qpns, qpsn, apn, aps, dyn1, dys1, dyn2, dys2
  LOGICAL extremum, first
  SAVE first

  REAL convpn, convps, convmpn, convmps
  REAL sinlon(iip1), sinlondlon(iip1)
  REAL coslon(iip1), coslondlon(iip1)
  SAVE sinlon, coslon, sinlondlon, coslondlon


  REAL ssum
  INTEGER ismax, ismin
  EXTERNAL ssum, convflu, ismin, ismax

  DATA first/.TRUE./

  IF (first) THEN
    PRINT *, 'SCHEMA AMONT NOUVEAU'
    first = .FALSE.
    DO i = 2, iip1
      coslon(i) = cos(rlonv(i))
      sinlon(i) = sin(rlonv(i))
      coslondlon(i) = coslon(i)*(rlonu(i)-rlonu(i-1))/pi
      sinlondlon(i) = sinlon(i)*(rlonu(i)-rlonu(i-1))/pi
    END DO
    coslon(1) = coslon(iip1)
    coslondlon(1) = coslondlon(iip1)
    sinlon(1) = sinlon(iip1)
    sinlondlon(1) = sinlondlon(iip1)
  END IF



  DO l = 1, llm

    DO ij = 1, ip1jmp1
      q(ij, l) = s0(ij, l)/sm(ij, l)
      dyq(ij) = sy(ij, l)/sm(ij, l)
    END DO

    ! --------------------------------
    ! CALCUL EN LATITUDE
    ! --------------------------------

    ! On commence par calculer la valeur du traceur moyenne sur le premier
    ! cercle
    ! de latitude autour du pole (qpns pour le pole nord et qpsn pour
    ! le pole nord) qui sera utilisee pour evaluer les pentes au pole.

    airej2 = ssum(iim, aire(iip2), 1)
    airejjm = ssum(iim, aire(ip1jm-iim), 1)
    DO i = 1, iim
      airescb(i) = aire(i+iip1)*q(i+iip1, l)
      airesch(i) = aire(i+ip1jm-iip1)*q(i+ip1jm-iip1, l)
    END DO
    qpns = ssum(iim, airescb, 1)/airej2
    qpsn = ssum(iim, airesch, 1)/airejjm

    ! calcul des pentes aux points v

    DO ij = 1, ip1jm
      dyqv(ij) = q(ij, l) - q(ij+iip1, l)
      adyqv(ij) = abs(dyqv(ij))
    END DO

    ! calcul des pentes aux points scalaires

    DO ij = iip2, ip1jm
      dyqmax(ij) = min(adyqv(ij-iip1), adyqv(ij))
      dyqmax(ij) = pente_max*dyqmax(ij)
    END DO

    ! calcul des pentes aux poles

    ! calcul des pentes limites aux poles

    ! cas ou on a un extremum au pole

    ! if(dyqv(ismin(iim,dyqv,1))*dyqv(ismax(iim,dyqv,1)).le.0.)
    ! &   apn=0.
    ! if(dyqv(ismax(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1)*
    ! &   dyqv(ismin(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1).le.0.)
    ! &   aps=0.

    ! limitation des pentes aux poles
    ! do ij=1,iip1
    ! dyq(ij)=apn*dyq(ij)
    ! dyq(ip1jm+ij)=aps*dyq(ip1jm+ij)
    ! enddo

    ! test
    ! do ij=1,iip1
    ! dyq(iip1+ij)=0.
    ! dyq(ip1jm+ij-iip1)=0.
    ! enddo
    ! do ij=1,ip1jmp1
    ! dyq(ij)=dyq(ij)*cos(rlatu((ij-1)/iip1+1))
    ! enddo

    IF (dyqv(ismin(iim,dyqv,1))*dyqv(ismax(iim,dyqv,1))<=0.) THEN
      DO ij = 1, iip1
        dyqmax(ij) = 0.
      END DO
    ELSE
      DO ij = 1, iip1
        dyqmax(ij) = pente_max*abs(dyqv(ij))
      END DO
    END IF

    IF (dyqv(ismax(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1)*dyqv(ismin(iim, &
        dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1)<=0.) THEN
      DO ij = ip1jm + 1, ip1jmp1
        dyqmax(ij) = 0.
      END DO
    ELSE
      DO ij = ip1jm + 1, ip1jmp1
        dyqmax(ij) = pente_max*abs(dyqv(ij-iip1))
      END DO
    END IF

    ! calcul des pentes limitees

    DO ij = 1, ip1jmp1
      IF (dyqv(ij)*dyqv(ij-iip1)>0.) THEN
        dyq(ij) = sign(min(abs(dyq(ij)),dyqmax(ij)), dyq(ij))
      ELSE
        dyq(ij) = 0.
      END IF
    END DO

    DO ij = 1, ip1jmp1
      sy(ij, l) = dyq(ij)*sm(ij, l)
    END DO

  END DO ! fin de la boucle sur les couches verticales

  RETURN
END SUBROUTINE limy
