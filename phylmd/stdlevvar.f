module stdlevvar_m

 IMPLICIT NONE

contains

  SUBROUTINE stdlevvar(klon, knon, nsrf, zxli, u1, v1, t1, q1, z1, ts1, &
       qsurf, rugos, psol, pat1, t_2m, q_2m, t_10m, q_10m, u_10m, ustar)

    ! From LMDZ4/libf/phylmd/stdlevvar.F90, version 1.3 2005/05/25 13:10:09

    use coefcdrag_m, only: coefcdrag
    USE suphec_m, ONLY: rg, rkappa

    ! Objet : calcul de la température et de l'humidité relative à 2 m
    ! et du module du vent à 10 m à partir des relations de
    ! Dyer-Businger et des équations de Louis.

    ! Reference: Hess, Colman and McAvaney (1995) 

    ! Author: I. Musat, 01.07.2002

    INTEGER, intent(in):: klon
    ! dimension de la grille physique (= nb_pts_latitude X nb_pts_longitude)

    INTEGER, intent(in):: knon
    ! knon----input-I- nombre de points pour un type de surface
    INTEGER, intent(in):: nsrf
    ! nsrf----input-I- indice pour le type de surface; voir indicesol.inc
    LOGICAL, intent(in):: zxli
    ! zxli----input-L- TRUE si calcul des cdrags selon Laurent Li
    REAL, dimension(klon), intent(in):: u1
    ! u1------input-R- vent zonal au 1er niveau du modele
    REAL, dimension(klon), intent(in):: v1
    ! v1------input-R- vent meridien au 1er niveau du modele
    REAL, dimension(klon), intent(in):: t1
    ! t1------input-R- temperature de l'air au 1er niveau du modele
    REAL, dimension(klon), intent(in):: q1
    ! q1------input-R- humidite relative au 1er niveau du modele
    REAL, dimension(klon), intent(in):: z1
    ! z1------input-R- geopotentiel au 1er niveau du modele
    REAL, dimension(klon), intent(in):: ts1
    ! ts1-----input-R- temperature de l'air a la surface
    REAL, dimension(klon), intent(in):: qsurf
     ! qsurf---input-R- humidite relative a la surface
    REAL, dimension(klon), intent(in):: rugos
    ! rugos---input-R- rugosite
   REAL, dimension(klon), intent(in):: psol
    ! psol----input-R- pression au sol
   REAL, dimension(klon), intent(in):: pat1
    ! pat1----input-R- pression au 1er niveau du modele

    REAL, dimension(klon), intent(out):: t_2m
    ! t_2m---output-R- temperature de l'air a 2m
    REAL, dimension(klon), intent(out):: q_2m
    ! q_2m---output-R- humidite relative a 2m
    REAL, dimension(klon), intent(out):: t_10m
    ! t_10m--output-R- temperature de l'air a 10m
    REAL, dimension(klon), intent(out):: q_10m
    ! q_10m--output-R- humidite specifique a 10m
    REAL, dimension(klon), intent(out):: u_10m
    ! u_10m--output-R- vitesse du vent a 10m
    REAL, intent(out):: ustar(klon) ! u*

    ! Local:

    ! RKAR : constante de von Karman
    REAL, PARAMETER:: RKAR=0.40
    ! niter : nombre iterations calcul "corrector"
    INTEGER, parameter:: niter=2, ncon=niter-1

    ! Variables locales
    INTEGER i, n
    REAL zref
    REAL, dimension(klon):: speed
    ! tpot : temperature potentielle
    REAL, dimension(klon):: tpot
    REAL, dimension(klon):: zri1, cdran
    REAL cdram(klon), cdrah(klon)
    ! ri1 : nb. de Richardson entre la surface --> la 1ere couche
    REAL, dimension(klon):: ri1 
    REAL, dimension(klon):: testar, qstar
    REAL, dimension(klon):: zdte, zdq 
    ! lmon : longueur de Monin-Obukhov selon Hess, Colman and McAvaney 
    DOUBLE PRECISION, dimension(klon):: lmon
    DOUBLE PRECISION, parameter:: eps=1.0D-20
    REAL, dimension(klon):: delu, delte, delq
    REAL, dimension(klon):: u_zref, te_zref, q_zref 
    REAL, dimension(klon):: temp, pref
    LOGICAL okri
    REAL, dimension(klon):: u_zref_p, temp_p, q_zref_p
    !convertgence
    REAL, dimension(klon):: te_zref_con, q_zref_con
    REAL, dimension(klon):: u_zref_c, temp_c, q_zref_c
    REAL, dimension(klon):: ok_pred, ok_corr

    !------------------------------------------------------------------------- 

    DO i=1, knon
       speed(i)=SQRT(u1(i)**2+v1(i)**2)
       ri1(i) = 0.0
    ENDDO

    okri=.FALSE.
    CALL coefcdrag(klon, knon, nsrf, zxli, speed, t1, q1, z1, psol, ts1, &
         qsurf, rugos, okri, ri1, cdram, cdrah, cdran, zri1, pref) 

    ! Star variables 

    DO i = 1, knon
       ri1(i) = zri1(i)
       tpot(i) = t1(i)* (psol(i)/pat1(i))**RKAPPA
       ustar(i) = sqrt(cdram(i) * speed(i) * speed(i))
       zdte(i) = tpot(i) - ts1(i)
       zdq(i) = max(q1(i), 0.0) - max(qsurf(i), 0.0)

       zdte(i) = sign(max(abs(zdte(i)), 1.e-10), zdte(i))

       testar(i) = (cdrah(i) * zdte(i) * speed(i))/ustar(i)
       qstar(i) = (cdrah(i) * zdq(i) * speed(i))/ustar(i)
       lmon(i) = (ustar(i) * ustar(i) * tpot(i))/ &
            (RKAR * RG * testar(i))
    ENDDO

    ! First aproximation of variables at zref  
    zref = 2.0
    CALL screenp(klon, knon, nsrf, speed, tpot, q1, &
         ts1, qsurf, rugos, lmon, &
         ustar, testar, qstar, zref, &
         delu, delte, delq)

    DO i = 1, knon
       u_zref(i) = delu(i)
       q_zref(i) = max(qsurf(i), 0.0) + delq(i)
       te_zref(i) = ts1(i) + delte(i)
       temp(i) = te_zref(i) * (psol(i)/pat1(i))**(-RKAPPA)
       q_zref_p(i) = q_zref(i)
       temp_p(i) = temp(i)
    ENDDO

    ! Iteration of the variables at the reference level zref :
    ! corrector calculation ; see Hess & McAvaney, 1995

    DO n = 1, niter
       okri=.TRUE.
       CALL screenc(klon, knon, nsrf, zxli, &
            u_zref, temp, q_zref, zref, &
            ts1, qsurf, rugos, psol, & 
            ustar, testar, qstar, okri, ri1, &
            pref, delu, delte, delq) 

       DO i = 1, knon
          u_zref(i) = delu(i)
          q_zref(i) = delq(i) + max(qsurf(i), 0.0)
          te_zref(i) = delte(i) + ts1(i) 

          ! return to normal temperature

          temp(i) = te_zref(i) * (psol(i)/pref(i))**(-RKAPPA)

          IF(n == ncon) THEN
             te_zref_con(i) = te_zref(i)
             q_zref_con(i) = q_zref(i)
          ENDIF
       ENDDO
    ENDDO

    ! verifier le critere de convergence : 0.25% pour te_zref et 5% pour qe_zref

    DO i = 1, knon
       q_zref_c(i) = q_zref(i)
       temp_c(i) = temp(i)

       ok_pred(i)=0.
       ok_corr(i)=1.

       t_2m(i) = temp_p(i) * ok_pred(i) + temp_c(i) * ok_corr(i)
       q_2m(i) = q_zref_p(i) * ok_pred(i) + q_zref_c(i) * ok_corr(i)
    ENDDO

    ! First aproximation of variables at zref  

    zref = 10.0
    CALL screenp(klon, knon, nsrf, speed, tpot, q1, &
         ts1, qsurf, rugos, lmon, &
         ustar, testar, qstar, zref, &
         delu, delte, delq)

    DO i = 1, knon
       u_zref(i) = delu(i)
       q_zref(i) = max(qsurf(i), 0.0) + delq(i)
       te_zref(i) = ts1(i) + delte(i)
       temp(i) = te_zref(i) * (psol(i)/pat1(i))**(-RKAPPA)
       u_zref_p(i) = u_zref(i)
    ENDDO

    ! Iteration of the variables at the reference level zref:
    ! corrector ; see Hess & McAvaney, 1995

    DO n = 1, niter
       okri=.TRUE.
       CALL screenc(klon, knon, nsrf, zxli, &
            u_zref, temp, q_zref, zref, &
            ts1, qsurf, rugos, psol, &
            ustar, testar, qstar, okri, ri1, &
            pref, delu, delte, delq)

       DO i = 1, knon
          u_zref(i) = delu(i)
          q_zref(i) = delq(i) + max(qsurf(i), 0.0)
          te_zref(i) = delte(i) + ts1(i)
          temp(i) = te_zref(i) * (psol(i)/pref(i))**(-RKAPPA)
       ENDDO
    ENDDO

    DO i = 1, knon
       u_zref_c(i) = u_zref(i)

       u_10m(i) = u_zref_p(i) * ok_pred(i) + u_zref_c(i) * ok_corr(i)

       q_zref_c(i) = q_zref(i)
       temp_c(i) = temp(i)
       t_10m(i) = temp_p(i) * ok_pred(i) + temp_c(i) * ok_corr(i)
       q_10m(i) = q_zref_p(i) * ok_pred(i) + q_zref_c(i) * ok_corr(i)
    ENDDO

  END subroutine stdlevvar

end module stdlevvar_m
