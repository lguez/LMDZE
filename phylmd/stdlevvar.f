module stdlevvar_m

 IMPLICIT NONE

contains

  SUBROUTINE stdlevvar(nsrf, u1, v1, t1, q1, z1, ts1, qsurf, rugos, psol, &
       pat1, t_2m, q_2m, t_10m, q_10m, wind10m, ustar)

    ! From LMDZ4/libf/phylmd/stdlevvar.F90, version 1.3, 2005/05/25 13:10:09

    ! Objet : calcul de la température et de l'humidité relative à 2 m
    ! et du module du vent à 10 m à partir des relations de
    ! Dyer-Businger et des équations de Louis.

    ! Reference: Hess, Colman and McAvaney (1995) 

    ! Author: I. Musat, July 1st, 2002

    use nr_util, only: assert_eq

    use coefcdrag_m, only: coefcdrag
    USE dimphy, ONLY: klon
    USE suphec_m, ONLY: rg, rkappa
    use screenc_m, only: screenc
    use screenp_m, only: screenp

    INTEGER, intent(in):: nsrf ! indice pour le type de surface
    REAL, intent(in):: u1(:) ! (knon) vent zonal au 1er niveau du modele
    REAL, intent(in):: v1(:) ! (knon) vent meridien au 1er niveau du modele
    REAL, intent(in):: t1(:) ! (knon) temperature de l'air au 1er
                             ! niveau du modele
    REAL, intent(in):: q1(klon) ! humidite relative au 1er niveau du modele
    REAL, intent(in):: z1 (klon) ! geopotentiel au 1er niveau du modele
    REAL, intent(in):: ts1(klon) ! temperature de l'air a la surface
    REAL, intent(in):: qsurf(klon) ! humidite relative a la surface
    REAL, intent(in):: rugos(klon) ! rugosite
    REAL, intent(in):: psol(klon) ! pression au sol
    REAL, intent(in):: pat1(klon) ! pression au 1er niveau du modele
    REAL, intent(out):: t_2m(klon) ! temperature de l'air a 2m
    REAL, intent(out):: q_2m(klon) ! humidite relative a 2m
    REAL, intent(out):: t_10m(klon) ! temperature de l'air a 10m
    REAL, intent(out):: q_10m(klon) ! humidite specifique a 10m
    REAL, intent(out):: wind10m(:) ! (knon) norme du vent \`a 10m
    REAL, intent(out):: ustar(:) ! (knon) u*

    ! Local:
    INTEGER knon ! nombre de points pour un type de surface
    REAL, PARAMETER:: RKAR = 0.4 ! constante de von Karman
    INTEGER, parameter:: niter = 2 ! nombre iterations calcul "corrector"
    INTEGER i, n
    REAL zref
    REAL, dimension(klon):: speed
    ! tpot : temperature potentielle
    REAL, dimension(klon):: tpot
    REAL cdram(klon), cdrah(klon)
    REAL, dimension(klon):: testar, qstar
    REAL, dimension(klon):: zdte, zdq 
    ! lmon : longueur de Monin-Obukhov selon Hess, Colman and McAvaney 
    DOUBLE PRECISION, dimension(klon):: lmon
    REAL, dimension(klon):: delu, delte, delq
    REAL, dimension(klon):: u_zref, te_zref, q_zref 
    REAL, dimension(klon):: temp
    real pref(size(u1)) ! (knon)

    !------------------------------------------------------------------------- 

    knon = assert_eq([size(u1), size(v1), size(t1), size(wind10m), &
         size(ustar)], "stdlevvar knon")
    
    DO i=1, knon
       speed(i)=SQRT(u1(i)**2+v1(i)**2)
    ENDDO

    CALL coefcdrag(nsrf, speed(:knon), t1(:knon), q1(:knon), z1(:knon), &
         psol(:knon), ts1, qsurf, rugos, cdram, cdrah) 

    ! Star variables 

    DO i = 1, knon
       tpot(i) = t1(i)* (psol(i)/pat1(i))**RKAPPA
       ustar(i) = sqrt(cdram(i) * speed(i) * speed(i))
       zdte(i) = tpot(i) - ts1(i)
       zdq(i) = max(q1(i), 0.0) - max(qsurf(i), 0.0)

       zdte(i) = sign(max(abs(zdte(i)), 1.e-10), zdte(i))

       testar(i) = (cdrah(i) * zdte(i) * speed(i))/ustar(i)
       qstar(i) = (cdrah(i) * zdq(i) * speed(i))/ustar(i)
       lmon(i) = (ustar(i) * ustar(i) * tpot(i)) / (RKAR * RG * testar(i))
    ENDDO

    ! First aproximation of variables at zref  
    zref = 2.0
    CALL screenp(knon, speed, tpot, q1, ts1, qsurf, rugos, lmon, ustar, &
         testar, qstar, zref, delu, delte, delq)

    DO i = 1, knon
       u_zref(i) = delu(i)
       q_zref(i) = max(qsurf(i), 0.0) + delq(i)
       te_zref(i) = ts1(i) + delte(i)
       temp(i) = te_zref(i) * (psol(i)/pat1(i))**(-RKAPPA)
    ENDDO

    ! Iteration of the variables at the reference level zref :
    ! corrector calculation ; see Hess & McAvaney, 1995

    DO n = 1, niter
       CALL screenc(klon, knon, nsrf, u_zref, temp, q_zref, zref, ts1, &
            qsurf, rugos, psol, ustar, testar, qstar, pref, delu, delte, delq) 

       DO i = 1, knon
          u_zref(i) = delu(i)
          q_zref(i) = delq(i) + max(qsurf(i), 0.0)
          te_zref(i) = delte(i) + ts1(i) 

          ! return to normal temperature
          temp(i) = te_zref(i) * (psol(i)/pref(i))**(-RKAPPA)
       ENDDO
    ENDDO

    ! verifier le critere de convergence : 0.25% pour te_zref et 5% pour qe_zref

    DO i = 1, knon
       t_2m(i) = temp(i)
       q_2m(i) = q_zref(i)
    ENDDO

    ! First aproximation of variables at zref  

    zref = 10.
    CALL screenp(knon, speed, tpot, q1, ts1, qsurf, rugos, lmon, ustar, &
         testar, qstar, zref, delu, delte, delq)

    DO i = 1, knon
       u_zref(i) = delu(i)
       q_zref(i) = max(qsurf(i), 0.0) + delq(i)
       te_zref(i) = ts1(i) + delte(i)
       temp(i) = te_zref(i) * (psol(i)/pat1(i))**(-RKAPPA)
    ENDDO

    ! Iteration of the variables at the reference level zref:
    ! corrector ; see Hess & McAvaney, 1995

    DO n = 1, niter
       CALL screenc(klon, knon, nsrf, u_zref, temp, q_zref, zref, ts1, &
            qsurf, rugos, psol, ustar, testar, qstar, pref, delu, delte, delq)

       DO i = 1, knon
          u_zref(i) = delu(i)
          q_zref(i) = delq(i) + max(qsurf(i), 0.0)
          te_zref(i) = delte(i) + ts1(i)
          temp(i) = te_zref(i) * (psol(i)/pref(i))**(-RKAPPA)
       ENDDO
    ENDDO

    DO i = 1, knon
       wind10m(i) = u_zref(i)
       t_10m(i) = temp(i)
       q_10m(i) = q_zref(i)
    ENDDO

  END subroutine stdlevvar

end module stdlevvar_m
