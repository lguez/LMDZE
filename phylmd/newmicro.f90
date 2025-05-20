module newmicro_m

  IMPLICIT none

contains

  SUBROUTINE newmicro (paprs, play, t, cldliq, clc, cldtau, clemi, cldh, cldl, &
       cldm, cldt, ctlwp, flwp, fiwp, flwc, fiwc)

    ! From LMDZ4/libf/phylmd/newmicro.F, version 1.2 2004/06/03 09:22:43

    ! Authors: Z. X. Li (LMD/CNRS), Johannes Quaas
    ! Date: 1993/09/10
    ! Objet : calcul de l'\'epaisseur optique et de l'\'emissivit\'e des nuages.

    USE conf_phys_m, ONLY: rad_chau1, rad_chau2
    USE dimphy, ONLY: klev, klon
    USE histwrite_phy_m, ONLY: histwrite_phy
    USE suphec_m, ONLY: rg

    REAL, intent(in):: paprs(:, :) ! (klon, klev+1)
    ! pression pour chaque inter-couche, en Pa

    real, intent(in):: play(:, :) ! (klon, klev)
    REAL, intent(in):: t(:, :) ! (klon, klev) temperature

    REAL, intent(in):: cldliq(:, :) ! (klon, klev)
    ! mass fraction of liquid water in atmosphere

    REAL, intent(inout):: clc(:, :) ! (klon, klev)
    ! couverture nuageuse pour le rayonnement (0 à 1)

    REAL, intent(out):: cldtau(:, :) ! (klon, klev)
    ! \'epaisseur optique des nuages
    
    REAL, intent(out):: clemi(:, :) ! (klon, klev)
    ! \'emissivit\'e des nuages (0 à 1)

    REAL, intent(out):: cldh(:), cldl(:), cldm(:), cldt(:) ! (klon)
    REAL, intent(out):: ctlwp(:) ! (klon)
    REAL, intent(out):: flwp(:), fiwp(:) ! (klon)
    REAL, intent(out):: flwc(:, :), fiwc(:, :) ! (klon, klev)

    ! Local:

    REAL re(klon, klev)
    ! cloud droplet effective radius multiplied by fl (micro m)

    REAL fl(klon, klev) 
    ! Denominator to re, introduced to avoid problems in the averaging
    ! of the output. fl is the fraction of liquid water clouds within
    ! a grid cell.

    REAL, PARAMETER:: cetahb = 0.45, cetamb = 0.8
    INTEGER i, k
    REAL zflwp ! liquid water path, in micrometers
    real fice ! fraction of ice in cloud
    REAL rad_chaud ! effective radius of liquid cloud droplets, in micrometers
    REAL, PARAMETER:: coef_chau = 0.13
    REAL, PARAMETER:: seuil_neb = 0.001, t_glace = 258.
    real tc, rei, zfiwp
    real k_ice
    real, parameter:: k_ice0 = 0.005 ! units=m2 / g
    real, parameter:: DF = 1.66 ! diffusivity factor

    !-----------------------------------------------------------------

    ! Calculer l'\'epaisseur optique et l'\'emissivit\'e des nuages

    loop_horizontal: DO i = 1, klon
       flwp(i) = 0.
       fiwp(i) = 0.

       loop_vertical: DO k = 1, klev
          clc(i, k) = MAX(clc(i, k), seuil_neb)

          ! liquid/ice cloud water paths:

          ! Linear transition:
          fice = 1. - (t(i, k) - t_glace) / (273.13 - t_glace)
          fice = MIN(MAX(fice, 0.), 1.)

          zflwp = 1000. * (1. - fice) * cldliq(i, k) / clc(i, k) &
               * (paprs(i, k) - paprs(i, k + 1)) / RG
          zfiwp = 1000. * fice * cldliq(i, k) / clc(i, k) &
               * (paprs(i, k) - paprs(i, k + 1)) / RG

          flwp(i) = flwp(i) + (1. - fice) * cldliq(i, k) &
               * (paprs(i, k) - paprs(i, k + 1)) / RG
          fiwp(i) = fiwp(i) &
               + fice * cldliq(i, k) * (paprs(i, k) - paprs(i, k + 1)) / RG

          ! Total Liquid/Ice water content
          flwc(i, k) = (1.-fice) * cldliq(i, k)
          fiwc(i, k) = fice * cldliq(i, k)
          ! In-Cloud Liquid/Ice water content

          ! effective cloud droplet radius (microns):

          ! for liquid water clouds: 
          rad_chaud = merge(rad_chau2, rad_chau1, k <= 3)
          
          ! For output diagnostics

          ! Cloud droplet effective radius (micro m)

          ! we multiply here with f * xl (fraction of liquid water
          ! clouds in the grid cell) to avoid problems in the
          ! averaging of the output.
          ! In the output of IOIPSL, derive the real cloud droplet 
          ! effective radius as re/fl

          fl(i, k) = clc(i, k) * (1.-fice) 
          re(i, k) = rad_chaud * fl(i, k)

          ! for ice clouds: as a function of the ambiant temperature
          ! (formula used by Iacobellis and Somerville (2000), with an 
          ! asymptotical value of 3.5 microns at T<-81.4 C added to be 
          ! consistent with observations of Heymsfield et al. 1986):
          tc = t(i, k)-273.15
          rei = merge(3.5, 0.71 * tc + 61.29, tc <= -81.4)

          ! Cloud optical thickness. For liquid clouds, traditional
          ! formula (e. g. Liou 2002 k0795 § 8.4.5.2). For ice clouds,
          ! Ebert and Curry (1992).
          if (zfiwp == 0. .or. rei <= 0.) rei = 1. 
          cldtau(i, k) = 3. / 2. * zflwp / rad_chaud &
               + zfiwp * (3.448e-03 + 2.431 / rei)

          ! cloud infrared emissivity:

          ! (the broadband infrared absorption coefficient is parameterized
          ! as a function of the effective cld droplet radius)

          ! Ebert and Curry (1992) formula as used by Kiehl & Zender (1995):
          k_ice = k_ice0 + 1. / rei

          clemi(i, k) = 1. - EXP(- coef_chau * zflwp - DF * k_ice * zfiwp)

          if (clc(i, k) <= seuil_neb) then
             clc(i, k) = 0.
             cldtau(i, k) = 0.
             clemi(i, k) = 0.
          end if
       ENDDO loop_vertical
    ENDDO loop_horizontal

    ! COMPUTE CLOUD LIQUID PATH AND TOTAL CLOUDINESS

    DO i = 1, klon
       cldt(i)=1.
       cldh(i)=1.
       cldm(i) = 1.
       cldl(i) = 1.
       ctlwp(i) = 0.
    ENDDO

    DO k = klev, 1, -1
       DO i = 1, klon
          ctlwp(i) = ctlwp(i) &
               + cldliq(i, k) * (paprs(i, k) - paprs(i, k + 1)) / RG
          cldt(i) = cldt(i) * (1.-clc(i, k))
          if (play(i, k) <= cetahb * paprs(i, 1)) &
               cldh(i) = cldh(i) * (1. - clc(i, k))
          if (play(i, k) > cetahb * paprs(i, 1) .AND. &
               play(i, k) <= cetamb * paprs(i, 1)) &
               cldm(i) = cldm(i) * (1.-clc(i, k))
          if (play(i, k) > cetamb * paprs(i, 1)) &
               cldl(i) = cldl(i) * (1. - clc(i, k))
       ENDDO
    ENDDO

    DO i = 1, klon
       cldt(i)=1.-cldt(i)
       cldh(i)=1.-cldh(i)
       cldm(i)=1.-cldm(i)
       cldl(i)=1.-cldl(i)
    ENDDO

    CALL histwrite_phy("re", re)
    CALL histwrite_phy("fl", fl)

  END SUBROUTINE newmicro

end module newmicro_m
