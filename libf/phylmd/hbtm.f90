module HBTM_m

  IMPLICIT none

contains

  SUBROUTINE HBTM(knon, paprs, pplay, t2m, t10m, q2m, q10m, ustar, flux_t, &
       flux_q, u, v, t, q, pblh, cape, EauLiq, ctei, pblT, therm, trmb1, &
       trmb2, trmb3, plcl)

    use dimens_m
    use dimphy
    use SUPHEC_M
    use yoethf_m
    use fcttre

    ! D'apres Holstag & Boville et Troen & Mahrt
    ! JAS 47 BLM
    ! Algorithme th�se Anne Mathieu
    ! Crit�re d'entra�nement Peter Duynkerke (JAS 50)
    ! written by: Anne MATHIEU and Alain LAHELLEC, 22nd November 1999
    ! features : implem. exces Mathieu

    ! modifications : decembre 99 passage th a niveau plus bas. voir fixer
    ! la prise du th a z/Lambda = -.2 (max Ray)
    ! Autre algo : entrainement ~ Theta+v =cste mais comment=>The?
    ! on peut fixer q a .7 qsat (cf. non adiabatique) => T2 et The2
    ! voir aussi //KE pblh = niveau The_e ou l = env.

    ! fin therm a la HBTM passage a forme Mathieu 12/09/2001

    ! Adaptation a LMDZ version couplee
    ! Pour le moment on fait passer en argument les grandeurs de surface :
    ! flux, t, q2m, t, q10m, on va utiliser systematiquement les grandeurs a 2m
    ! mais on garde la possibilit� de changer si besoin est (jusqu'� pr�sent
    ! la forme de HB avec le 1er niveau modele etait conservee)

    REAL RLvCp, REPS
    ! Arguments:

    ! nombre de points a calculer
    INTEGER, intent(in):: knon

    REAL, intent(in):: t2m(klon) ! temperature a 2 m
    real t10m(klon) ! temperature a 10 m
    ! q a 2 et 10m
    REAL q2m(klon), q10m(klon)
    REAL ustar(klon)
    ! pression a inter-couche (Pa)
    REAL paprs(klon, klev+1)
    ! pression au milieu de couche (Pa)
    REAL pplay(klon, klev)
    ! Flux
    REAL flux_t(klon, klev), flux_q(klon, klev)
    ! vitesse U (m/s)
    REAL u(klon, klev)
    ! vitesse V (m/s)
    REAL v(klon, klev)
    ! temperature (K)
    REAL t(klon, klev)
    ! vapeur d'eau (kg/kg)
    REAL q(klon, klev)

    INTEGER isommet
    ! limite max sommet pbl
    PARAMETER (isommet=klev)
    REAL vk
    ! Von Karman => passer a .41 ! cf U.Olgstrom
    PARAMETER (vk=0.35)
    REAL ricr
    PARAMETER (ricr=0.4)
    REAL fak
    ! b calcul du Prandtl et de dTetas
    PARAMETER (fak=8.5)
    REAL fakn
    ! a
    PARAMETER (fakn=7.2)
    REAL onet
    PARAMETER (onet=1.0/3.0)
    REAL t_coup
    PARAMETER(t_coup=273.15)
    REAL zkmin
    PARAMETER (zkmin=0.01)
    REAL betam
    ! pour Phim / h dans la S.L stable
    PARAMETER (betam=15.0)
    REAL betah
    PARAMETER (betah=15.0)
    REAL betas
    ! Phit dans la S.L. stable (mais 2 formes /
    PARAMETER (betas=5.0)
    ! z/OBL<>1
    REAL sffrac
    ! S.L. = z/h < .1
    PARAMETER (sffrac=0.1)
    REAL binm
    PARAMETER (binm=betam*sffrac)
    REAL binh
    PARAMETER (binh=betah*sffrac)
    REAL ccon
    PARAMETER (ccon=fak*sffrac*vk)

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
    ! stable pbl with levels within pbl
    LOGICAL stblev(klon)
    ! unstbl pbl with levels within pbl
    LOGICAL unslev(klon)
    ! unstb pbl w/lvls within srf pbl lyr
    LOGICAL unssrf(klon)
    ! unstb pbl w/lvls in outer pbl lyr
    LOGICAL unsout(klon)
    LOGICAL check(klon) ! Richardson number > critical
    ! flag de prolongerment cape pour pt Omega
    LOGICAL omegafl(klon)
    REAL pblh(klon)
    REAL pblT(klon)
    REAL plcl(klon)

    ! Monin-Obukhov lengh
    REAL obklen(klon)

    REAL zdu2
    ! thermal virtual temperature excess
    REAL therm(klon)
    REAL trmb1(klon), trmb2(klon), trmb3(klon)
    ! Algorithme thermique
    REAL s(klon, klev) ! [P/Po]^Kappa milieux couches
    ! equivalent potential temperature of therma
    REAL The_th(klon)
    ! total water of thermal
    REAL qT_th(klon)
    ! T thermique niveau precedent
    REAL Tbef(klon)
    REAL qsatbef(klon)
    ! le thermique est sature
    LOGICAL Zsat(klon)
    ! Cape du thermique
    REAL Cape(klon)
    ! Cape locale
    REAL Kape(klon)
    ! Eau liqu integr du thermique
    REAL EauLiq(klon)
    ! Critere d'instab d'entrainmt des nuages de
    REAL ctei(klon)
    REAL the1, the2, aa, zthvd, zthvu, xintpos, qqsat
    REAL a1, a2, a3
    REAL xhis, rnum, th1, thv1, thv2, ql2
    REAL qsat2, qT1, q2, t1, t2, xnull
    REAL quadsat, spblh, reduc

    ! inverse phi function for momentum
    REAL phiminv(klon)
    ! inverse phi function for heat
    REAL phihinv(klon)
    ! turbulent velocity scale for momentum
    REAL wm(klon)
    ! k*ustar*pblh
    REAL fak1(klon)
    ! k*wm*pblh
    REAL fak2(klon)
    ! fakn*wstr/wm
    REAL fak3(klon)
    ! level eddy diffusivity for momentum
    REAL pblk(klon)
    ! Prandtl number for eddy diffusivities
    REAL pr(klon)
    ! zmzp / Obukhov length
    REAL zl(klon)
    ! zmzp / pblh
    REAL zh(klon)
    ! (1-(zmzp/pblh))**2
    REAL zzh(klon)
    ! w*, convective velocity scale
    REAL wstr(klon)
    ! current level height
    REAL zm(klon)
    ! current level height + one level up
    REAL zp(klon)
    REAL zcor, zdelta, zcvm5

    REAL fac, pblmin, zmzp, term

    !-----------------------------------------------------------------

    ! initialisations
    q_star = 0
    t_star = 0

    b212=sqrt(b1*b2)
    b2sr=sqrt(b2)

    ! Initialisation
    RLvCp = RLVTT/RCPD
    REPS = RD/RV

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
       khfs(i) = - flux_t(i, 1)*zxt*Rd / (RCPD*paprs(i, 1))
       kqfs(i) = - flux_q(i, 1)*zxt*Rd / (paprs(i, 1))
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
       trmb1(i) = 0.
       trmb2(i) = 0.
       trmb3(i) = 0.
    ENDDO

    ! PBL height calculation: Search for level of pbl. Scan upward
    ! until the Richardson number between the first level and the
    ! current level exceeds the "critical" value.  (bonne idee Nu de
    ! separer le Ric et l'exces de temp du thermique)
    fac = 100.
    DO k = 2, isommet
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
       if (check(i)) pblh(i) = z(i, isommet)
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

          a1=b1*(1.+2.*RETV*qT_th(i))*t_star**2
          a2=(RETV*T2m(i))**2*b2*q_star**2
          a3=2.*RETV*T2m(i)*b212*q_star*t_star
          aa=a1+a2+a3

          therm(i) = sqrt( b1*(1.+2.*RETV*qT_th(i))*t_star**2 &
               + (RETV*T2m(i))**2*b2*q_star**2 &
               + max(0., 2.*RETV*T2m(i)*b212*q_star*t_star))

          ! Theta et qT du thermique (forme H&B) avec exces
          ! (attention, on ajoute therm(i) qui est virtuelle ...)
          ! pourquoi pas sqrt(b1)*t_star ?
          qT_th(i) = qT_th(i) + b2sr*q_star
          ! new on differre le calcul de Theta_e
          rhino(i, 1) = 0.
       ENDIF
    ENDDO

    ! Improve pblh estimate for unstable conditions using the
    ! convective temperature excess :
    DO k = 2, isommet
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
       if (check(i)) pblh(i) = z(i, isommet)
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
       Kape(i) = 0.
       EauLiq(i) = 0.
       CTEI(i) = 0.
       pblk(i) = 0.0
       fak1(i) = ustar(i)*pblh(i)*vk

       ! Do additional preparation for unstable cases only, set temperature
       ! and moisture perturbations depending on stability.
       ! Remarque : les formule sont prises dans leur forme CS
       IF (unstbl(i)) THEN
          ! Niveau de ref du thermique
          zxt=(T2m(i)-zref*0.5*RG/RCPD/(1.+RVTMP2*qT_th(i))) &
               *(1.+RETV*qT_th(i))
          phiminv(i) = (1. - binm*pblh(i)/obklen(i))**onet
          phihinv(i) = sqrt(1. - binh*pblh(i)/obklen(i))
          wm(i) = ustar(i)*phiminv(i)
          fak2(i) = wm(i)*pblh(i)*vk
          wstr(i) = (heatv(i)*RG*pblh(i)/zxt)**onet
          fak3(i) = fakn*wstr(i)/wm(i)
       ENDIF
       ! Computes Theta_e for thermal (all cases : to be modified)
       ! attention ajout therm(i) = virtuelle
       The_th(i) = T2m(i) + therm(i) + RLvCp*qT_th(i)
    ENDDO

    ! Main level loop to compute the diffusivities and
    ! counter-gradient terms:
    DO k = 2, isommet
       ! Find levels within boundary layer:
       DO i = 1, knon
          unslev(i) = .FALSE.
          stblev(i) = .FALSE.
          zm(i) = z(i, k-1)
          zp(i) = z(i, k)
          IF (zkmin == 0. .AND. zp(i) > pblh(i)) zp(i) = pblh(i)
          IF (zm(i) < pblh(i)) THEN
             zmzp = 0.5*(zm(i) + zp(i))
             zh(i) = zmzp/pblh(i)
             zl(i) = zmzp/obklen(i)
             zzh(i) = 0.
             IF (zh(i) <= 1.) zzh(i) = (1. - zh(i))**2

             ! stblev for points zm < plbh and stable and neutral
             ! unslev for points zm < plbh and unstable
             IF (unstbl(i)) THEN
                unslev(i) = .TRUE.
             ELSE
                stblev(i) = .TRUE.
             ENDIF
          ENDIF
       ENDDO

       ! Stable and neutral points; set diffusivities; counter-gradient
       ! terms zero for stable case:
       DO i = 1, knon
          IF (stblev(i)) THEN
             IF (zl(i) <= 1.) THEN
                pblk(i) = fak1(i)*zh(i)*zzh(i)/(1. + betas*zl(i))
             ELSE
                pblk(i) = fak1(i)*zh(i)*zzh(i)/(betas + zl(i))
             ENDIF
          ENDIF
       ENDDO

       ! unssrf, unstable within surface layer of pbl
       ! unsout, unstable within outer layer of pbl
       DO i = 1, knon
          unssrf(i) = .FALSE.
          unsout(i) = .FALSE.
          IF (unslev(i)) THEN
             IF (zh(i) < sffrac) THEN
                unssrf(i) = .TRUE.
             ELSE
                unsout(i) = .TRUE.
             ENDIF
          ENDIF
       ENDDO

       ! Unstable for surface layer; counter-gradient terms zero
       DO i = 1, knon
          IF (unssrf(i)) THEN
             term = (1. - betam*zl(i))**onet
             pblk(i) = fak1(i)*zh(i)*zzh(i)*term
             pr(i) = term/sqrt(1. - betah*zl(i))
          ENDIF
       ENDDO

       ! Unstable for outer layer; counter-gradient terms non-zero:
       DO i = 1, knon
          IF (unsout(i)) THEN
             pblk(i) = fak2(i)*zh(i)*zzh(i)
             pr(i) = phiminv(i)/phihinv(i) + ccon*fak3(i)/fak
          ENDIF
       ENDDO

       ! For all layers, compute integral info and CTEI
       DO i = 1, knon
          if (check(i).or.omegafl(i)) then
             if (.not.Zsat(i)) then
                T2 = T2m(i) * s(i, k)
                ! thermodyn functions
                zdelta=MAX(0., SIGN(1., RTT - T2))
                qqsat= r2es * FOEEW(T2, zdelta) / pplay(i, k)
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
                   Tbef(i) = T2
                endif
             endif
             qsatbef(i) = qqsat
             ! cette ligne a deja ete faite normalement ?
          endif
       ENDDO
    end DO

  END SUBROUTINE HBTM

end module HBTM_m
