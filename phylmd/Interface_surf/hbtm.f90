module HBTM_m

  IMPLICIT none

contains

  SUBROUTINE HBTM(paprs, pplay, t2m, q2m, ustar, flux_t, flux_q, u, v, t, q, &
       pblh, cape, EauLiq, ctei, pblT, therm, plcl)

    ! D'apr\'es Holstag et Boville et Troen et Mahrt
    ! JAS 47 BLM

    ! Algorithme th\'ese Anne Mathieu. Crit\'ere d'entra\^inement
    ! Peter Duynkerke (JAS 50). Written by: Anne MATHIEU and Alain
    ! LAHELLEC, 22nd November 1999.

    ! Modifications : d\'ecembre 99 passage th \`a niveau plus bas. Voir fixer
    ! la prise du th \`a z / Lambda = -.2 (max Ray)
    ! Autre algorithme : entra\^inement ~ Theta + v =constante
    ! mais comment ? The ?
    ! On peut fixer q \`a 0.7 qsat (cf. non adiabatique) d'où T2 et The2.
    ! Voir aussi KE pblh = niveau The_e ou l = env.

    ! Adaptation \`a LMDZ version coupl\'ee. Pour le moment on fait
    ! passer en argument les grandeurs de surface : flux, t, q2m. On
    ! va utiliser syst\'ematiquement les grandeurs \`a 2 m mais on
    ! garde la possibilit\'e de changer si besoin (jusqu'\`a pr\'esent
    ! la forme de HB avec le premier niveau mod\`ele \'etait
    ! conserv\'ee).

    USE dimphy, ONLY: klev, klon
    USE fcttre, ONLY: foeew
    USE suphec_m, ONLY: rcpd, rd, retv, rg, rkappa, rtt
    USE yoethf_m, ONLY: r2es, rvtmp2

    ! Arguments:

    REAL, intent(in):: paprs(klon, klev+1) ! pression a inter-couche (Pa)
    REAL, intent(in):: pplay(klon, klev) ! pression au milieu de couche (Pa)
    REAL, intent(in):: t2m(klon) ! temperature a 2 m
    REAL, intent(in):: q2m(klon) ! q a 2 et 10m
    REAL, intent(in):: ustar(:) ! (knon)
    REAL, intent(in):: flux_t(:), flux_q(:) ! (knon) flux à la surface
    REAL, intent(in):: u(:, :) ! (knon, klev) vitesse U (m/s)
    REAL, intent(in):: v(:, :) ! (knon, klev) vitesse V (m/s)
    REAL, intent(in):: t(:, :) ! (knon, klev) temperature (K)
    REAL, intent(in):: q(:, :) ! (knon, klev) vapeur d'eau (kg/kg)
    REAL, intent(out):: pblh(:) ! (knon)
    REAL, INTENT(OUT):: Cape(:) ! (knon) Cape du thermique
    REAL, INTENT(OUT):: EauLiq(:) ! (knon) Eau liqu integr du thermique

    REAL, INTENT(OUT):: ctei(:) ! (knon)
    ! Critere d'instab d'entrainmt des nuages de

    REAL, INTENT(OUT):: pblT(:) ! (knon)
    REAL, INTENT(OUT):: therm(:) ! (knon) ! thermal virtual temperature excess
    REAL, INTENT(OUT):: plcl(:) ! (knon)

    ! Local:
    
    INTEGER knon ! nombre de points a calculer
    REAL vk
    ! Von Karman => passer a .41 ! cf U.Olgstrom
    PARAMETER (vk=0.35)
    REAL ricr
    PARAMETER (ricr=0.4)
    ! a
    REAL onet
    PARAMETER (onet=1.0/3.0)
    REAL zkmin
    PARAMETER (zkmin=0.01)
    REAL betam
    ! pour Phim / h dans la S.L stable
    PARAMETER (betam=15.0)
    ! z/OBL<>1
    REAL sffrac
    ! S.L. = z/h < .1
    PARAMETER (sffrac=0.1)
    REAL binm
    PARAMETER (binm=betam*sffrac)

    REAL q_star, t_star
    ! Lambert correlations T' q' avec T* q*
    REAL b1, b2, b212, b2sr
    PARAMETER (b1=70., b2=20.)

    REAL z(klon, klev)

    REAL zref
    ! Niveau de ref a 2m peut eventuellement
    PARAMETER (zref=2.)
    ! etre choisi a 10m

    INTEGER i, k
    REAL zxt
    ! surface kinematic heat flux [mK/s]
    REAL khfs(klon)
    ! sfc kinematic constituent flux [m/s]
    REAL kqfs(klon)
    ! surface virtual heat flux
    REAL heatv(klon)
    ! bulk Richardon no. mais en Theta_v
    REAL rhino(klon, klev)
    ! pts w/unstbl pbl (positive virtual ht flx)
    LOGICAL unstbl(klon)
    LOGICAL check(klon) ! Richardson number > critical
    ! flag de prolongerment cape pour pt Omega
    LOGICAL omegafl(klon)

    ! Monin-Obukhov lengh
    REAL obklen(klon)

    REAL zdu2
    ! Algorithme thermique
    REAL s(klon, klev) ! [P/Po]^Kappa milieux couches
    ! total water of thermal
    REAL qT_th(klon)
    ! T thermique niveau precedent
    REAL qsatbef(klon)
    ! le thermique est sature
    LOGICAL Zsat(klon)
    REAL zthvd, zthvu, qqsat
    REAL t2

    ! inverse phi function for momentum
    REAL phiminv(klon)
    ! turbulent velocity scale for momentum
    REAL wm(klon)
    ! current level height + one level up
    REAL zp(klon)
    REAL zcor

    REAL pblmin

    !-----------------------------------------------------------------

    knon = size(pblh)

    ! initialisations
    q_star = 0
    t_star = 0
    therm = 0.

    b212=sqrt(b1*b2)
    b2sr=sqrt(b2)

    ! Calculer les hauteurs de chaque couche
    ! (geopotentielle Int_dp/ro = Int_[Rd.T.dp/p] z = geop/g)
    ! pourquoi ne pas utiliser Phi/RG ?
    DO i = 1, knon
       z(i, 1) = RD * t(i, 1) / (0.5*(paprs(i, 1)+pplay(i, 1))) &
            * (paprs(i, 1)-pplay(i, 1)) / RG
       s(i, 1) = (pplay(i, 1)/paprs(i, 1))**RKappa
    ENDDO
    ! s(k) = [pplay(k)/ps]^kappa
    ! + + + + + + + + + pplay <-> s(k) t dp=pplay(k-1)-pplay(k)
    ! ----------------- paprs <-> sig(k)
    ! + + + + + + + + + pplay <-> s(k-1)
    ! + + + + + + + + + pplay <-> s(1) t dp=paprs-pplay z(1)
    ! ----------------- paprs <-> sig(1)

    DO k = 2, klev
       DO i = 1, knon
          z(i, k) = z(i, k-1) &
               + RD * 0.5*(t(i, k-1)+t(i, k)) / paprs(i, k) &
               * (pplay(i, k-1)-pplay(i, k)) / RG
          s(i, k) = (pplay(i, k) / paprs(i, 1))**RKappa
       ENDDO
    ENDDO

    ! Determination des grandeurs de surface
    DO i = 1, knon
       ! Niveau de ref choisi a 2m
       zxt = t2m(i)

       ! convention >0 vers le bas ds lmdz
       khfs(i) = - flux_t(i)*zxt*Rd / (RCPD*paprs(i, 1))
       kqfs(i) = - flux_q(i)*zxt*Rd / paprs(i, 1)
       ! verifier que khfs et kqfs sont bien de la forme w'l'
       heatv(i) = khfs(i) + 0.608*zxt*kqfs(i)
       ! a comparer aussi aux sorties de clqh : flux_T/RoCp et flux_q/RoLv
       ! Theta et qT du thermique sans exces (interpolin vers surf)
       ! chgt de niveau du thermique (jeudi 30/12/1999)
       ! (interpolation lineaire avant integration phi_h)
       qT_th(i) = q2m(i)
    ENDDO

    DO i = 1, knon
       ! Global Richardson
       rhino(i, 1) = 0.0
       check(i) = .TRUE.
       ! on initialise pblh a l'altitude du 1er niv
       pblh(i) = z(i, 1)
       plcl(i) = 6000.
       ! Lambda = -u*^3 / (alpha.g.kvon.<w'Theta'v>
       obklen(i) = -t(i, 1)*ustar(i)**3/(RG*vk*heatv(i))
    ENDDO

    ! PBL height calculation: Search for level of pbl. Scan upward
    ! until the Richardson number between the first level and the
    ! current level exceeds the "critical" value.  (bonne idee Nu de
    ! separer le Ric et l'exces de temp du thermique)
    DO k = 2, klev
       DO i = 1, knon
          IF (check(i)) THEN
             ! pourquoi / niveau 1 (au lieu du sol) et le terme en u*^2 ?
             zdu2 = u(i, k)**2+v(i, k)**2
             zdu2 = max(zdu2, 1.0e-20)
             ! Theta_v environnement
             zthvd=t(i, k)/s(i, k)*(1.+RETV*q(i, k))

             ! therm Theta_v sans exces (avec hypothese fausse de H&B, sinon,
             ! passer par Theta_e et virpot)
             zthvu = T2m(i)*(1.+RETV*qT_th(i))
             ! Le Ri par Theta_v
             ! On a nveau de ref a 2m ???
             rhino(i, k) = (z(i, k)-zref)*RG*(zthvd-zthvu) &
                  /(zdu2*0.5*(zthvd+zthvu))

             IF (rhino(i, k).GE.ricr) THEN
                pblh(i) = z(i, k-1) + (z(i, k-1)-z(i, k)) * &
                     (ricr-rhino(i, k-1))/(rhino(i, k-1)-rhino(i, k))
                ! test04
                pblh(i) = pblh(i) + 100.
                pblT(i) = t(i, k-1) + (t(i, k)-t(i, k-1)) * &
                     (pblh(i)-z(i, k-1))/(z(i, k)-z(i, k-1))
                check(i) = .FALSE.
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    ! Set pbl height to maximum value where computation exceeds number of
    ! layers allowed
    DO i = 1, knon
       if (check(i)) pblh(i) = z(i, klev)
    ENDDO

    ! Improve estimate of pbl height for the unstable points.
    ! Find unstable points (sensible heat flux is upward):
    DO i = 1, knon
       IF (heatv(i) > 0.) THEN
          unstbl(i) = .TRUE.
          check(i) = .TRUE.
       ELSE
          unstbl(i) = .FALSE.
          check(i) = .FALSE.
       ENDIF
    ENDDO

    ! For the unstable case, compute velocity scale and the
    ! convective temperature excess:
    DO i = 1, knon
       IF (check(i)) THEN
          phiminv(i) = (1.-binm*pblh(i)/obklen(i))**onet

          ! CALCUL DE wm
          ! Ici on considerera que l'on est dans la couche de surf jusqu'a 100
          ! On prend svt couche de surface=0.1*h mais on ne connait pas h
          ! Dans la couche de surface
          wm(i)= ustar(i)*phiminv(i)

          ! forme Mathieu :
          q_star = kqfs(i)/wm(i)
          t_star = khfs(i)/wm(i)

          therm(i) = sqrt( b1*(1.+2.*RETV*qT_th(i))*t_star**2 &
               + (RETV*T2m(i))**2*b2*q_star**2 &
               + max(0., 2.*RETV*T2m(i)*b212*q_star*t_star))

          ! Theta et qT du thermique (forme H&B) avec exces
          ! (attention, on ajoute therm(i) qui est virtuelle ...)
          ! pourquoi pas sqrt(b1)*t_star ?
          qT_th(i) = qT_th(i) + b2sr*q_star
          ! new on diff\`ere le calcul de Theta_e
          rhino(i, 1) = 0.
       ENDIF
    ENDDO

    ! Improve pblh estimate for unstable conditions using the
    ! convective temperature excess :
    DO k = 2, klev
       DO i = 1, knon
          IF (check(i)) THEN
             zdu2 = u(i, k)**2 + v(i, k)**2
             zdu2 = max(zdu2, 1e-20)
             ! Theta_v environnement
             zthvd=t(i, k)/s(i, k)*(1.+RETV*q(i, k))

             ! et therm Theta_v (avec hypothese de constance de H&B,
             zthvu = T2m(i)*(1.+RETV*qT_th(i)) + therm(i)

             ! Le Ri par Theta_v
             ! Niveau de ref 2m
             rhino(i, k) = (z(i, k)-zref)*RG*(zthvd-zthvu) &
                  /(zdu2*0.5*(zthvd+zthvu))

             IF (rhino(i, k).GE.ricr) THEN
                pblh(i) = z(i, k-1) + (z(i, k-1)-z(i, k)) * &
                     (ricr-rhino(i, k-1))/(rhino(i, k-1)-rhino(i, k))
                ! test04
                pblh(i) = pblh(i) + 100.
                pblT(i) = t(i, k-1) + (t(i, k)-t(i, k-1)) * &
                     (pblh(i)-z(i, k-1))/(z(i, k)-z(i, k-1))
                check(i) = .FALSE.
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    ! Set pbl height to maximum value where computation exceeds number of
    ! layers allowed
    DO i = 1, knon
       if (check(i)) pblh(i) = z(i, klev)
    ENDDO

    ! PBL height must be greater than some minimum mechanical mixing depth
    ! Several investigators have proposed minimum mechanical mixing depth
    ! relationships as a function of the local friction velocity, u*. We
    ! make use of a linear relationship of the form h = c u* where c=700.
    ! The scaling arguments that give rise to this relationship most often
    ! represent the coefficient c as some constant over the local coriolis
    ! parameter. Here we make use of the experimental results of Koracin
    ! and Berkowicz (1988) [BLM, Vol 43] for wich they recommend 0.07/f
    ! where f was evaluated at 39.5 N and 52 N. Thus we use a typical mid
    ! latitude value for f so that c = 0.07/f = 700.
    DO i = 1, knon
       pblmin = 700. * ustar(i)
       pblh(i) = MAX(pblh(i), pblmin)
       ! par exemple :
       pblT(i) = t(i, 2) + (t(i, 3)-t(i, 2)) * &
            (pblh(i)-z(i, 2))/(z(i, 3)-z(i, 2))
    ENDDO

    ! pblh is now available; do preparation for diffusivity calculation:
    DO i = 1, knon
       check(i) = .TRUE.
       Zsat(i) = .FALSE.
       ! omegafl utilise pour prolongement CAPE
       omegafl(i) = .FALSE.
       Cape(i) = 0.
       EauLiq(i) = 0.
       CTEI(i) = 0.

       ! Do additional preparation for unstable cases only, set temperature
       ! and moisture perturbations depending on stability.
       ! Remarque : les formule sont prises dans leur forme CS
       IF (unstbl(i)) THEN
          ! Niveau de ref du thermique
          zxt=(T2m(i)-zref*0.5*RG/RCPD/(1.+RVTMP2*qT_th(i))) &
               *(1.+RETV*qT_th(i))
          phiminv(i) = (1. - binm*pblh(i)/obklen(i))**onet
          wm(i) = ustar(i)*phiminv(i)
       ENDIF
    ENDDO

    ! Main level loop to compute the diffusivities and
    ! counter-gradient terms:
    loop_level: DO k = 2, klev
       ! Find levels within boundary layer:
       DO i = 1, knon
          zp(i) = z(i, k)
          IF (zkmin == 0. .AND. zp(i) > pblh(i)) zp(i) = pblh(i)
       ENDDO

       ! For all layers, compute integral info and CTEI
       DO i = 1, knon
          if (check(i) .or. omegafl(i)) then
             if (.not. Zsat(i)) then
                T2 = T2m(i) * s(i, k)
                ! thermodyn functions
                qqsat= r2es * FOEEW(T2, RTT >= T2) / pplay(i, k)
                qqsat=MIN(0.5, qqsat)
                zcor=1./(1.-retv*qqsat)
                qqsat=qqsat*zcor

                if (qqsat < qT_th(i)) then
                   ! on calcule lcl
                   if (k == 2) then
                      plcl(i) = z(i, k)
                   else
                      plcl(i) = z(i, k-1) + (z(i, k-1)-z(i, k)) &
                           * (qT_th(i)-qsatbef(i)) / (qsatbef(i)-qqsat)
                   endif
                   Zsat(i) = .true.
                endif
             endif
             qsatbef(i) = qqsat
             ! cette ligne a deja ete faite normalement ?
          endif
       ENDDO
    end DO loop_level

  END SUBROUTINE HBTM

end module HBTM_m
