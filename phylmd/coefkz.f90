module coefkz_m

  IMPLICIT none

contains

  SUBROUTINE coefkz(nsrf, knon, paprs, pplay, ksta, ksta_ter, ts, rugos, u, v, &
       t, q, qsurf, coefm, coefh)

    ! Authors: F. Hourdin, M. Forichon, Z. X. Li (LMD/CNRS)
    ! date: 1993/09/22
    ! Objet : calculer le coefficient de frottement du sol ("Cdrag") et les
    ! coefficients d'échange turbulent dans l'atmosphère.

    USE indicesol, ONLY: is_oce
    USE dimphy, ONLY: klev, klon
    USE conf_gcm_m, ONLY: prt_level
    USE suphec_m, ONLY: rcpd, rd, retv, rg, rkappa, rlstt, rlvtt, rtt
    USE yoethf_m, ONLY: r2es, r5ies, r5les, rvtmp2
    USE fcttre, ONLY: dqsatl, dqsats, foede, foeew, qsatl, qsats, thermcep
    USE conf_phys_m, ONLY: iflag_pbl
    use clcdrag_m, only: clcdrag

    ! Arguments:

    integer, intent(in):: nsrf ! indicateur de la nature du sol
    INTEGER, intent(in):: knon ! nombre de points a traiter

    REAL, intent(in):: paprs(klon, klev+1)
    ! pression a chaque intercouche (en Pa)

    real, intent(in):: pplay(klon, klev)
    ! pression au milieu de chaque couche (en Pa)

    REAL, intent(in):: ksta, ksta_ter
    REAL, intent(in):: ts(klon) ! temperature du sol (en Kelvin)
    REAL, intent(in):: rugos(klon) ! longeur de rugosite (en m)
    REAL, intent(in):: u(klon, klev), v(klon, klev) ! wind
    REAL, intent(in):: t(klon, klev) ! temperature (K)
    real, intent(in):: q(klon, klev) ! vapeur d'eau (kg/kg)
    real, intent(in):: qsurf(klon) 
    REAL, intent(out):: coefm(:, :) ! (knon, klev) coefficient, vitesse

    real, intent(out):: coefh(:, :) ! (knon, klev) 
    ! coefficient, chaleur et humidité

    ! Local:

    INTEGER itop(knon) ! numero de couche du sommet de la couche limite

    ! Quelques constantes et options:

    REAL, PARAMETER:: cepdu2 =0.1**2
    REAL, PARAMETER:: CKAP = 0.4
    REAL, PARAMETER:: cb = 5.
    REAL, PARAMETER:: cc = 5.
    REAL, PARAMETER:: cd = 5.
    REAL, PARAMETER:: clam = 160.
    REAL, PARAMETER:: ratqs = 0.05 ! largeur de distribution de vapeur d'eau
    
    LOGICAL, PARAMETER:: richum = .TRUE.
    ! utilise le nombre de Richardson humide

    REAL, PARAMETER:: ric = 0.4 ! nombre de Richardson critique
    REAL, PARAMETER:: prandtl = 0.4

    REAL kstable ! diffusion minimale (situation stable)
    REAL, PARAMETER:: mixlen = 35. ! constante contrôlant longueur de mélange
    INTEGER, PARAMETER:: isommet = klev ! sommet de la couche limite

    LOGICAL, PARAMETER:: tvirtu = .TRUE.
    ! calculer Ri d'une maniere plus performante

    LOGICAL, PARAMETER:: opt_ec = .FALSE.
    ! formule du Centre Europeen dans l'atmosphere

    INTEGER i, k, kk
    REAL zgeop(klon, klev)
    REAL zmgeom(klon)
    REAL ri(klon)
    REAL l2(klon)

    REAL u1(klon), v1(klon), t1(klon), q1(klon), z1(klon)

    REAL zdphi, zdu2, ztvd, ztvu, cdn
    REAL scf
    REAL zt, zq, zdelta, zcvm5, zcor, zqs, zfr, zdqs
    REAL z2geomf, zalh2, alm2, zscfh, scfm
    REAL, PARAMETER:: t_coup = 273.15
    REAL gamt(2:klev) ! contre-gradient pour la chaleur sensible: Kelvin/metre

    !--------------------------------------------------------------------

    ! Prescrire la valeur de contre-gradient
    if (iflag_pbl.eq.1) then
       DO k = 3, klev
          gamt(k) = -1.0E-03
       ENDDO
       gamt(2) = -2.5E-03
    else
       DO k = 2, klev
          gamt(k) = 0.0
       ENDDO
    ENDIF

    IF ( nsrf .NE. is_oce ) THEN
       kstable = ksta_ter
    ELSE
       kstable = ksta
    ENDIF

    ! Calculer les géopotentiels de chaque couche
    DO i = 1, knon
       zgeop(i, 1) = RD * t(i, 1) / (0.5 * (paprs(i, 1) + pplay(i, 1))) &
            * (paprs(i, 1) - pplay(i, 1))
    ENDDO
    DO k = 2, klev
       DO i = 1, knon
          zgeop(i, k) = zgeop(i, k-1) &
               + RD * 0.5*(t(i, k-1)+t(i, k)) / paprs(i, k) &
               * (pplay(i, k-1)-pplay(i, k))
       ENDDO
    ENDDO

    ! Calculer le frottement au sol (Cdrag)

    DO i = 1, knon
       u1(i) = u(i, 1)
       v1(i) = v(i, 1)
       t1(i) = t(i, 1)
       q1(i) = q(i, 1)
       z1(i) = zgeop(i, 1)
    ENDDO

    CALL clcdrag(klon, knon, nsrf, .false., u1, v1, t1, q1, z1, ts, qsurf, &
         rugos, coefm(:, 1), coefh(:, 1)) 

    ! Calculer les coefficients turbulents dans l'atmosphere

    itop = isommet

    loop_vertical: DO k = 2, isommet
       loop_horiz: DO i = 1, knon
          zdu2 = MAX(cepdu2, (u(i, k)-u(i, k-1))**2 &
               +(v(i, k)-v(i, k-1))**2)
          zmgeom(i) = zgeop(i, k)-zgeop(i, k-1)
          zdphi =zmgeom(i) / 2.0
          zt = (t(i, k)+t(i, k-1)) * 0.5
          zq = (q(i, k)+q(i, k-1)) * 0.5

          ! calculer Qs et dQs/dT:

          IF (thermcep) THEN
             zdelta = MAX(0., SIGN(1., RTT-zt))
             zcvm5 = R5LES*RLVTT/RCPD/(1.0+RVTMP2*zq)*(1.-zdelta)  &
                  + R5IES*RLSTT/RCPD/(1.0+RVTMP2*zq)*zdelta
             zqs = R2ES * FOEEW(zt, zdelta) / pplay(i, k)
             zqs = MIN(0.5, zqs)
             zcor = 1./(1.-RETV*zqs)
             zqs = zqs*zcor
             zdqs = FOEDE(zt, zdelta, zcvm5, zqs, zcor)
          ELSE
             IF (zt  <  t_coup) THEN
                zqs = qsats(zt) / pplay(i, k)
                zdqs = dqsats(zt, zqs)
             ELSE
                zqs = qsatl(zt) / pplay(i, k)
                zdqs = dqsatl(zt, zqs)
             ENDIF
          ENDIF

          ! calculer la fraction nuageuse (processus humide):

          zfr = (zq+ratqs*zq-zqs) / (2.0*ratqs*zq)
          zfr = MAX(0.0, MIN(1.0, zfr))
          IF (.NOT.richum) zfr = 0.0

          !  calculer le nombre de Richardson:

          IF (tvirtu) THEN
             ztvd =( t(i, k) &
                  + zdphi/RCPD/(1.+RVTMP2*zq) &
                  *( (1.-zfr) + zfr*(1.+RLVTT*zqs/RD/zt)/(1.+zdqs) ) &
                  )*(1.+RETV*q(i, k))
             ztvu =( t(i, k-1) &
                  - zdphi/RCPD/(1.+RVTMP2*zq) &
                  *( (1.-zfr) + zfr*(1.+RLVTT*zqs/RD/zt)/(1.+zdqs) ) &
                  )*(1.+RETV*q(i, k-1))
             ri(i) =zmgeom(i)*(ztvd-ztvu)/(zdu2*0.5*(ztvd+ztvu))
             ri(i) = ri(i) &
                  + zmgeom(i)*zmgeom(i)/RG*gamt(k) &
                  *(paprs(i, k)/101325.0)**RKAPPA &
                  /(zdu2*0.5*(ztvd+ztvu))
          ELSE
             ! calcul de Ridchardson compatible LMD5
             ri(i) =(RCPD*(t(i, k)-t(i, k-1)) &
                  -RD*0.5*(t(i, k)+t(i, k-1))/paprs(i, k) &
                  *(pplay(i, k)-pplay(i, k-1)) &
                  )*zmgeom(i)/(zdu2*0.5*RCPD*(t(i, k-1)+t(i, k)))
             ri(i) = ri(i) + &
                  zmgeom(i)*zmgeom(i)*gamt(k)/RG &
                  *(paprs(i, k)/101325.0)**RKAPPA &
                  /(zdu2*0.5*(t(i, k-1)+t(i, k)))
          ENDIF

          ! finalement, les coefficients d'echange sont obtenus:

          cdn = SQRT(zdu2) / zmgeom(i) * RG

          IF (opt_ec) THEN
             z2geomf = zgeop(i, k-1)+zgeop(i, k)
             alm2 = (0.5*ckap/RG*z2geomf &
                  /(1.+0.5*ckap/rg/clam*z2geomf))**2
             zalh2 = (0.5*ckap/rg*z2geomf &
                  /(1.+0.5*ckap/RG/(clam*SQRT(1.5*cd))*z2geomf))**2
             IF (ri(i) < 0.) THEN
                ! situation instable
                scf = ((zgeop(i, k)/zgeop(i, k-1))**(1./3.)-1.)**3 &
                     / (zmgeom(i)/RG)**3 / (zgeop(i, k-1)/RG)
                scf = SQRT(-ri(i)*scf)
                scfm = 1.0 / (1.0+3.0*cb*cc*alm2*scf)
                zscfh = 1.0 / (1.0+3.0*cb*cc*zalh2*scf)
                coefm(i, k) = cdn * alm2 * (1. - 2. * cb * ri(i) * scfm)
                coefh(i, k) = cdn*zalh2*(1.-3.0*cb*ri(i)*zscfh)
             ELSE
                ! situation stable
                scf = SQRT(1.+cd*ri(i))
                coefm(i, k) = cdn * alm2 / (1. + 2. * cb * ri(i) / scf)
                coefh(i, k) = cdn*zalh2/(1.+3.0*cb*ri(i)*scf)
             ENDIF
          ELSE
             l2(i) = (mixlen*MAX(0.0, (paprs(i, k)-paprs(i, itop(i)+1)) &
                  /(paprs(i, 2)-paprs(i, itop(i)+1)) ))**2
             coefm(i, k) = sqrt(max(cdn**2 * (ric - ri(i)) / ric, kstable))
             coefm(i, k)= l2(i) * coefm(i, k)
             coefh(i, k) = coefm(i, k) / prandtl ! h et m different
          ENDIF
       ENDDO loop_horiz
    ENDDO loop_vertical

    ! Au-delà du sommet, pas de diffusion turbulente :
    forall (i = 1: knon)
       coefh(i, itop(i) + 1:) = 0.
       coefm(i, itop(i) + 1:) = 0.
    END forall

  END SUBROUTINE coefkz

end module coefkz_m
