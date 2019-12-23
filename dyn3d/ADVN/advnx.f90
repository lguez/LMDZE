SUBROUTINE advnx(q, qg, qd, masse, u_m, mode)

  ! Auteur : F. Hourdin

  ! ********************************************************************
  ! Shema  d'advection " pseudo amont " .
  ! ********************************************************************
  ! nq,iq,q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....


  ! --------------------------------------------------------------------
  USE dimensions
  USE paramet_m
  USE comconst
  USE disvert_m
  USE conf_gcm_m
  IMPLICIT NONE



  ! Arguments:
  ! ----------
  INTEGER mode
  REAL masse(ip1jmp1, llm)
  REAL u_m(ip1jmp1, llm)
  REAL q(ip1jmp1, llm), qd(ip1jmp1, llm), qg(ip1jmp1, llm)

  ! Local
  ! ---------

  INTEGER i, j, ij, l, indu(ip1jmp1), niju, iju, ijq
  INTEGER n0, nl(llm)

  REAL new_m, zu_m, zdq, zz
  REAL zsigg(ip1jmp1, llm), zsigd(ip1jmp1, llm), zsig
  REAL u_mq(ip1jmp1, llm)

  REAL zm, zq, zsigm, zsigp, zqm, zqp, zu

  LOGICAL ladvplus(ip1jmp1, llm)

  REAL prec
  SAVE prec

  DATA prec/1.E-15/

  DO l = 1, llm
    DO ij = iip2, ip1jm
      zdq = qd(ij, l) - qg(ij, l)
      IF (abs(zdq)>prec) THEN
        zsigd(ij, l) = (q(ij,l)-qg(ij,l))/zdq
        zsigg(ij, l) = 1. - zsigd(ij, l)
      ELSE
        zsigd(ij, l) = 0.5
        zsigg(ij, l) = 0.5
        qd(ij, l) = q(ij, l)
        qg(ij, l) = q(ij, l)
      END IF
    END DO
  END DO

  ! calcul de la pente maximum dans la maille en valeur absolue

  DO l = 1, llm
    DO ij = iip2, ip1jm - 1
      IF (u_m(ij,l)>=0.) THEN
        zsigp = zsigd(ij, l)
        zsigm = zsigg(ij, l)
        zqp = qd(ij, l)
        zqm = qg(ij, l)
        zm = masse(ij, l)
        zq = q(ij, l)
      ELSE
        zsigm = zsigd(ij+1, l)
        zsigp = zsigg(ij+1, l)
        zqm = qd(ij+1, l)
        zqp = qg(ij+1, l)
        zm = masse(ij+1, l)
        zq = q(ij+1, l)
      END IF
      zu = abs(u_m(ij,l))
      ladvplus(ij, l) = zu > zm
      zsig = zu/zm
      IF (zsig==0.) zsigp = 0.1
      IF (mode==1) THEN
        IF (zsig<=zsigp) THEN
          u_mq(ij, l) = u_m(ij, l)*zqp
        ELSE IF (mode==1) THEN
          u_mq(ij, l) = sign(zm, u_m(ij,l))*(zsigp*zqp+(zsig-zsigp)*zqm)
        END IF
      ELSE
        IF (zsig<=zsigp) THEN
          u_mq(ij, l) = u_m(ij, l)*(zqp-0.5*zsig/zsigp*(zqp-zq))
        ELSE
          zz = 0.5*(zsig-zsigp)/zsigm
          u_mq(ij, l) = sign(zm, u_m(ij,l))*(0.5*(zq+zqp)*zsigp+(zsig-zsigp)* &
            (zq+zz*(zqm-zq)))
        END IF
      END IF
    END DO
  END DO

  DO l = 1, llm
    DO ij = iip1 + iip1, ip1jm, iip1
      u_mq(ij, l) = u_mq(ij-iim, l)
      ladvplus(ij, l) = ladvplus(ij-iim, l)
    END DO
  END DO

  ! =================================================================
  ! SCHEMA SEMI-LAGRAGIEN EN X DANS LES REGIONS POLAIRES
  ! =================================================================
  ! tris des regions a traiter
  n0 = 0
  DO l = 1, llm
    nl(l) = 0
    DO ij = iip2, ip1jm
      IF (ladvplus(ij,l)) THEN
        nl(l) = nl(l) + 1
        u_mq(ij, l) = 0.
      END IF
    END DO
    n0 = n0 + nl(l)
  END DO

  IF (n0>1) THEN
    IF (prt_level>9) PRINT *, &
      'Nombre de points pour lesquels on advect plus que le', &
      'contenu de la maille : ', n0

    DO l = 1, llm
      IF (nl(l)>0) THEN
        iju = 0
        ! indicage des mailles concernees par le traitement special
        DO ij = iip2, ip1jm
          IF (ladvplus(ij,l) .AND. mod(ij,iip1)/=0) THEN
            iju = iju + 1
            indu(iju) = ij
          END IF
        END DO
        niju = iju

        ! traitement des mailles
        DO iju = 1, niju
          ij = indu(iju)
          j = (ij-1)/iip1 + 1
          zu_m = u_m(ij, l)
          u_mq(ij, l) = 0.
          IF (zu_m>0.) THEN
            ijq = ij
            i = ijq - (j-1)*iip1
            ! accumulation pour les mailles completements advectees
            DO WHILE (zu_m>masse(ijq,l))
              u_mq(ij, l) = u_mq(ij, l) + q(ijq, l)*masse(ijq, l)
              zu_m = zu_m - masse(ijq, l)
              i = mod(i-2+iim, iim) + 1
              ijq = (j-1)*iip1 + i
            END DO
            ! MODIFS SPECIFIQUES DU SCHEMA
            ! ajout de la maille non completement advectee
            zsig = zu_m/masse(ijq, l)
            IF (zsig<=zsigd(ijq,l)) THEN
              u_mq(ij, l) = u_mq(ij, l) + zu_m*(qd(ijq,l)-0.5*zsig/zsigd(ijq, &
                l)*(qd(ijq,l)-q(ijq,l)))
            ELSE
              ! u_mq(ij,l)=u_mq(ij,l)+zu_m*q(ijq,l)
              ! goto 8888
              zz = 0.5*(zsig-zsigd(ijq,l))/zsigg(ijq, l)
              IF (.NOT. (zz>0. .AND. zz<=0.5)) THEN
                PRINT *, 'probleme2 au point ij=', ij, '  l=', l
                PRINT *, 'zz=', zz
                STOP
              END IF
              u_mq(ij, l) = u_mq(ij, l) + masse(ijq, l)*(0.5*(q(ijq, &
                l)+qd(ijq,l))*zsigd(ijq,l)+(zsig-zsigd(ijq,l))*(q(ijq, &
                l)+zz*(qg(ijq,l)-q(ijq,l))))
            END IF
          ELSE
            ijq = ij + 1
            i = ijq - (j-1)*iip1
            ! accumulation pour les mailles completements advectees
            DO WHILE (-zu_m>masse(ijq,l))
              u_mq(ij, l) = u_mq(ij, l) - q(ijq, l)*masse(ijq, l)
              zu_m = zu_m + masse(ijq, l)
              i = mod(i, iim) + 1
              ijq = (j-1)*iip1 + i
            END DO
            ! ajout de la maille non completement advectee
            ! 2eme MODIF SPECIFIQUE
            zsig = -zu_m/masse(ij+1, l)
            IF (zsig<=zsigg(ijq,l)) THEN
              u_mq(ij, l) = u_mq(ij, l) + zu_m*(qg(ijq,l)-0.5*zsig/zsigg(ijq, &
                l)*(qg(ijq,l)-q(ijq,l)))
            ELSE
              ! u_mq(ij,l)=u_mq(ij,l)+zu_m*q(ijq,l)
              ! goto 9999
              zz = 0.5*(zsig-zsigg(ijq,l))/zsigd(ijq, l)
              IF (.NOT. (zz>0. .AND. zz<=0.5)) THEN
                PRINT *, 'probleme22 au point ij=', ij, '  l=', l
                PRINT *, 'zz=', zz
                STOP
              END IF
              u_mq(ij, l) = u_mq(ij, l) - masse(ijq, l)*(0.5*(q(ijq, &
                l)+qg(ijq,l))*zsigg(ijq,l)+(zsig-zsigg(ijq,l))*(q(ijq, &
                l)+zz*(qd(ijq,l)-q(ijq,l))))
            END IF
            ! fin de la modif
          END IF
        END DO
      END IF
    END DO
  END IF ! n0.gt.0

  ! bouclage en latitude
  DO l = 1, llm
    DO ij = iip1 + iip1, ip1jm, iip1
      u_mq(ij, l) = u_mq(ij-iim, l)
    END DO
  END DO

  ! =================================================================
  ! CALCUL DE LA CONVERGENCE DES FLUX
  ! =================================================================

  DO l = 1, llm
    DO ij = iip2 + 1, ip1jm
      new_m = masse(ij, l) + u_m(ij-1, l) - u_m(ij, l)
      q(ij, l) = (q(ij,l)*masse(ij,l)+u_mq(ij-1,l)-u_mq(ij,l))/new_m
      masse(ij, l) = new_m
    END DO
    ! Modif Fred 22 03 96 correction d'un bug (les scopy ci-dessous)
    DO ij = iip1 + iip1, ip1jm, iip1
      q(ij-iim, l) = q(ij, l)
      masse(ij-iim, l) = masse(ij, l)
    END DO
  END DO

  RETURN
END SUBROUTINE advnx
