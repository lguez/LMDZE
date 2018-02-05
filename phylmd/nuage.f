module nuage_m

  IMPLICIT none

contains

  SUBROUTINE nuage (paprs, pplay, t, pqlwp, pclc, pcltau, pclemi, pch, pcl, &
       pcm, pct, pctlwp)

    ! From LMDZ4/libf/phylmd/nuage.F, version 1.1.1.1, 2004/05/19 12:53:07

    use dimphy, only: klon, klev
    use SUPHEC_M, only: rg

    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 1993/09/10
    ! Objet: Calculer \'epaisseur optique et \'emissivit\'e des nuages

    ! Arguments:

    REAL, intent(in):: paprs(klon, klev+1)
    real, intent(in):: pplay(klon, klev)
    REAL, intent(in):: t(klon, klev) ! temperature

    REAL, intent(in):: pqlwp(klon, klev)
    ! eau liquide nuageuse dans l'atmosphere (kg/kg)

    REAL, intent(inout):: pclc(klon, klev)
    ! couverture nuageuse pour le rayonnement (0 \`a 1)

    REAL, intent(out):: pcltau(klon, klev) ! \'epaisseur optique des nuages
    real, intent(out):: pclemi(klon, klev) ! \'emissivit\'e des nuages (0 \`a 1)
    REAL pch(klon), pcl(klon), pcm(klon), pct(klon), pctlwp(klon)

    ! Local:

    LOGICAL lo

    REAL cetahb, cetamb
    PARAMETER (cetahb = 0.45, cetamb = 0.80)

    INTEGER i, k
    REAL zflwp, zfice

    REAL radius, rad_froid, rad_chaud, rad_chau1, rad_chau2
    PARAMETER (rad_chau1=13.0, rad_chau2=9.0, rad_froid=35.0)
    REAL coef, coef_froi, coef_chau
    PARAMETER (coef_chau=0.13, coef_froi=0.09)
    REAL seuil_neb, t_glace
    PARAMETER (seuil_neb=0.001, t_glace=273.0-15.0)
    INTEGER nexpo ! exponentiel pour glace/eau
    PARAMETER (nexpo=6)

    !--------------------------------------------------------------------
    
    ! Calculer l'epaisseur optique et l'emmissivite des nuages

    DO k = 1, klev
       DO i = 1, klon
          rad_chaud = rad_chau1
          IF (k <= 3) rad_chaud = rad_chau2

          pclc(i, k) = MAX(pclc(i, k), seuil_neb)
          zflwp = 1000.*pqlwp(i, k)/RG/pclc(i, k) &
               *(paprs(i, k)-paprs(i, k+1))
          zfice = 1.0 - (t(i, k)-t_glace) / (273.13-t_glace)
          zfice = MIN(MAX(zfice, 0.0), 1.0)
          zfice = zfice**nexpo

          radius = rad_chaud * (1.-zfice) + rad_froid * zfice
          coef = coef_chau * (1.-zfice) + coef_froi * zfice
          pcltau(i, k) = 3.0/2.0 * zflwp / radius
          pclemi(i, k) = 1.0 - EXP(- coef * zflwp)
          lo = (pclc(i, k) <= seuil_neb)
          IF (lo) pclc(i, k) = 0.0
          IF (lo) pcltau(i, k) = 0.0
          IF (lo) pclemi(i, k) = 0.0
       END DO
    END DO

    ! COMPUTE CLOUD LIQUID PATH AND TOTAL CLOUDINESS

    DO i = 1, klon
       pct(i)=1.0
       pch(i)=1.0
       pcm(i) = 1.0
       pcl(i) = 1.0
       pctlwp(i) = 0.0
    END DO

    DO k = klev, 1, -1
       DO i = 1, klon
          pctlwp(i) = pctlwp(i) &
               + pqlwp(i, k)*(paprs(i, k)-paprs(i, k+1))/RG
          pct(i) = pct(i)*(1.0-pclc(i, k))
          if (pplay(i, k) <= cetahb*paprs(i, 1)) &
               pch(i) = pch(i)*(1.0-pclc(i, k))
          if (pplay(i, k) > cetahb*paprs(i, 1) .AND. &
               pplay(i, k) <= cetamb*paprs(i, 1)) &
               pcm(i) = pcm(i)*(1.0-pclc(i, k))
          if (pplay(i, k) > cetamb*paprs(i, 1)) &
               pcl(i) = pcl(i)*(1.0-pclc(i, k))
       END DO
    END DO

    DO i = 1, klon
       pct(i)=1.-pct(i)
       pch(i)=1.-pch(i)
       pcm(i)=1.-pcm(i)
       pcl(i)=1.-pcl(i)
    END DO

  END SUBROUTINE nuage

end module nuage_m
