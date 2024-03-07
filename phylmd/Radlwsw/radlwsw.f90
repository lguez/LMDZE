module radlwsw_m

  IMPLICIT none

contains

  SUBROUTINE radlwsw(dist, mu0, fract, paprs, play, tsol, albedo, t_seri, &
       q_seri, wo, cldfra, cldemi, cldtau, heat, heat0, cool, cool0, radsol, &
       topsw, toplw, solsw, sollw, sollwdown, topsw0, toplw0, solsw0, sollw0, &
       lwdn0, lwdn, lwup0, lwup, swdn0, swdn, swup0, swup)

    ! From LMDZ4/libf/phylmd/radlwsw.F, version 1.4, 2005/06/06 13:16:33
    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 1996/07/19

    ! Objet : interface entre le mod\`ele et les rayonnements solaire
    ! et infrarouge

    USE clesphys, ONLY: solaire
    USE dimphy, ONLY: klev, klon
    use lw_m, only: lw
    USE suphec_m, ONLY: rg
    use sw_m, only: sw
    USE yoethf_m, ONLY: rvtmp2

    real, intent(in):: dist ! distance Terre-Soleil, en ua
    real, intent(in):: mu0(:) ! (klon) cosinus de l'angle zenithal
    real, intent(in):: fract(:) ! (klon)  duree d'ensoleillement normalisee

    real, intent(in):: paprs(:, :) ! (klon, klev + 1)
    ! pression a inter-couche (Pa)

    real, intent(in):: play(:, :) ! (klon, klev)
    ! pression au milieu de couche (Pa)

    real, intent(in):: tsol(:) ! (klon)  temperature du sol (en K)
    real, intent(in):: albedo(:) ! (klon)  albedo du sol (entre 0 et 1)
    real, intent(in):: t_seri(:, :) ! (klon, klev) temperature (K)
    real, intent(in):: q_seri(:, :) ! (klon, klev) vapeur d'eau (en kg/kg)

    real, intent(in):: wo(:, :) ! (klon, klev)
    ! column-density of ozone in a layer, in kilo-Dobsons

    real, intent(in):: cldfra(:, :) ! (klon, klev)
    ! fraction nuageuse (entre 0 et 1)

    real, intent(in):: cldemi(:, :) ! (klon, klev)
    ! emissivite des nuages dans l'IR (entre 0 et 1)

    real, intent(in):: cldtau(:, :) ! (klon, klev)
    ! \'epaisseur optique des nuages dans le visible (present-day value)

    real, intent(out):: heat(:, :) ! (klon, klev)
    ! échauffement atmosphérique (visible) (K/jour)

    real, intent(out):: heat0(:, :) ! (klon, klev) chauffage solaire ciel clair

    real, intent(out):: cool(:, :) ! (klon, klev)
    ! refroidissement dans l'IR (K/jour)

    real, intent(out):: cool0(:, :) ! (klon, klev)
    ! refroidissement infrarouge ciel clair

    real, intent(out):: radsol(:) ! (klon)
    ! bilan radiatif net au sol (W m-2), positif vers le bas

    real, intent(out):: topsw(:) ! (klon)
    ! flux solaire net au sommet de l'atmosph`ere

    real, intent(out):: toplw(:) ! (klon)
    ! rayonnement infrarouge montant au sommet de l'atmosphère

    real, intent(out):: solsw(:) ! (klon) flux solaire net à la surface

    real, intent(out):: sollw(:) ! (klon)
    ! rayonnement infrarouge net à la surface

    real, intent(out):: sollwdown(:) ! (klon)
    real, intent(out):: topsw0(:) ! (klon)
    real, intent(out):: toplw0(:) ! (klon)
    real, intent(out):: solsw0(:), sollw0(:) ! (klon)
    REAL, intent(out):: lwdn0(:, :), lwdn(:, :) ! (klon, klev + 1)
    REAL, intent(out):: lwup0(:, :), lwup(:, :) ! (klon, klev + 1)
    REAL, intent(out):: swdn0(:, :), swdn(:, :) ! (klon, klev + 1)
    REAL, intent(out):: swup0(:, :), swup(:, :) ! (klon, klev + 1)

    ! Local:

    DOUBLE PRECISION ZFSUP(KLON, KLEV + 1)
    DOUBLE PRECISION ZFSDN(KLON, KLEV + 1)
    DOUBLE PRECISION ZFSUP0(KLON, KLEV + 1)
    DOUBLE PRECISION ZFSDN0(KLON, KLEV + 1)

    DOUBLE PRECISION ZFLUP(KLON, KLEV + 1)
    DOUBLE PRECISION ZFLDN(KLON, KLEV + 1)
    DOUBLE PRECISION ZFLUP0(KLON, KLEV + 1)
    DOUBLE PRECISION ZFLDN0(KLON, KLEV + 1)

    DOUBLE PRECISION zx_alpha1, zx_alpha2
    INTEGER k, kk, i
    DOUBLE PRECISION PSCT

    DOUBLE PRECISION PALBD(klon, 2), PALBP(klon, 2)
    DOUBLE PRECISION PEMIS(klon), PDT0(klon), PVIEW(klon)
    DOUBLE PRECISION PPSOL(klon), PDP(klon, klev)
    DOUBLE PRECISION PTL(klon, klev + 1), PPMB(klon, klev + 1)
    DOUBLE PRECISION PTAVE(klon, klev)
    DOUBLE PRECISION PWV(klon, klev), PQS(klon, klev)
    DOUBLE PRECISION POZON(klon, klev) ! mass fraction of ozone
    DOUBLE PRECISION PAER(klon, klev, 5) ! AEROSOLS' OPTICAL THICKNESS
    DOUBLE PRECISION PCLDLD(klon, klev)
    DOUBLE PRECISION PCLDLU(klon, klev)
    DOUBLE PRECISION PCLDSW(klon, klev)
    DOUBLE PRECISION PTAU(klon, 2, klev)
    DOUBLE PRECISION POMEGA(klon, 2, klev)
    DOUBLE PRECISION PCG(klon, 2, klev)

    DOUBLE PRECISION zfract(klon), zrmu0(klon)

    DOUBLE PRECISION zheat(klon, klev), zcool(klon, klev)
    DOUBLE PRECISION zheat0(klon, klev), zcool0(klon, klev)
    DOUBLE PRECISION ztopsw(klon), ztoplw(klon)
    DOUBLE PRECISION zsolsw(klon), zsollw(klon)
    DOUBLE PRECISION zsollwdown(klon)

    DOUBLE PRECISION ztopsw0(klon), ztoplw0(klon)
    DOUBLE PRECISION zsolsw0(klon), zsollw0(klon)
    DOUBLE PRECISION zznormcp

    ! The following quantities are needed for the aerosol radiative forcings:
    DOUBLE PRECISION ztopswad(klon), zsolswad(klon)
    ! Aerosol direct forcing at TOA and surface

    real, parameter:: dobson_u = 2.1415e-05 ! Dobson unit, in kg m-2

    !----------------------------------------------------------------------

    heat = 0.
    cool = 0.
    heat0 = 0.
    cool0 = 0.
    PSCT = solaire / dist**2

    DO i = 1, klon
       zfract(i) = fract(i)
       zrmu0(i) = mu0(i)
       PALBD(i, 1) = albedo(i)
       PALBD(i, 2) = albedo(i)
       PALBP(i, 1) = albedo(i)
       PALBP(i, 2) = albedo(i)
       ! cf. JLD pour etre en accord avec ORCHIDEE il faut mettre
       ! PEMIS(i) = 0.96
       PEMIS(i) = 1.
       PVIEW(i) = 1.66
       PPSOL(i) = paprs(i, 1)
       zx_alpha1 = (paprs(i, 1)-play(i, 2)) &
            / (play(i, 1)-play(i, 2))
       zx_alpha2 = 1. - zx_alpha1
       PTL(i, 1) = t_seri(i, 1) * zx_alpha1 + t_seri(i, 2) * zx_alpha2
       PTL(i, klev + 1) = t_seri(i, klev)
       PDT0(i) = tsol(i) - PTL(i, 1)
    ENDDO
    DO k = 2, klev
       DO i = 1, klon
          PTL(i, k) = (t_seri(i, k) + t_seri(i, k-1))*0.5
       ENDDO
    ENDDO
    DO k = 1, klev
       DO i = 1, klon
          PDP(i, k) = paprs(i, k)-paprs(i, k + 1)
          PTAVE(i, k) = t_seri(i, k)
          PWV(i, k) = MAX(q_seri(i, k), 1e-12)
          PQS(i, k) = PWV(i, k)
          POZON(i, k) = wo(i, k) * RG * dobson_u * 1e3 &
               / (paprs(i, k) - paprs(i, k + 1))
          PCLDLD(i, k) = cldfra(i, k)*cldemi(i, k)
          PCLDLU(i, k) = cldfra(i, k)*cldemi(i, k)
          PCLDSW(i, k) = cldfra(i, k)
          PTAU(i, 1, k) = MAX(cldtau(i, k), 1e-05)
          ! (1e-12 serait instable)
          PTAU(i, 2, k) = MAX(cldtau(i, k), 1e-05)
          ! (pour 32-bit machines)
          POMEGA(i, 1, k) = 0.9999 - 5e-04 * EXP(-0.5 * PTAU(i, 1, k))
          POMEGA(i, 2, k) = 0.9988 - 2.5e-03 * EXP(-0.05 * PTAU(i, 2, k))
          PCG(i, 1, k) = 0.865
          PCG(i, 2, k) = 0.910
       ENDDO
    ENDDO

    DO k = 1, klev + 1
       DO i = 1, klon
          PPMB(i, k) = paprs(i, k)/100.
       ENDDO
    ENDDO

    DO kk = 1, 5
       DO k = 1, klev
          DO i = 1, klon
             PAER(i, k, kk) = 1E-15
          ENDDO
       ENDDO
    ENDDO

    CALL LW(PPMB, PDP, PDT0, PEMIS, PTL, PTAVE, PWV, POZON, PAER, PCLDLD, &
         PCLDLU, PVIEW, zcool, zcool0, ztoplw, zsollw, ztoplw0, zsollw0, &
         zsollwdown, ZFLUP, ZFLDN, ZFLUP0, ZFLDN0)
    CALL SW(PSCT, zrmu0, zfract, PPMB, PDP, PPSOL, PALBD, PALBP, PTAVE, &
         PWV, PQS, POZON, PCLDSW, PTAU, POMEGA, PCG, zheat, zheat0, ztopsw, &
         zsolsw, ztopsw0, zsolsw0, ZFSUP, ZFSDN, ZFSUP0, ZFSDN0, ztopswad, &
         zsolswad)

    DO i = 1, klon
       radsol(i) = zsolsw(i) + zsollw(i)
       topsw(i) = ztopsw(i)
       toplw(i) = ztoplw(i)
       solsw(i) = zsolsw(i)
       sollw(i) = zsollw(i)
       sollwdown(i) = zsollwdown(i)

       DO k = 1, klev + 1
          lwdn0(i, k) = ZFLDN0(i, k)
          lwdn(i, k) = ZFLDN(i, k)
          lwup0(i, k) = ZFLUP0(i, k)
          lwup(i, k) = ZFLUP(i, k)
       ENDDO

       topsw0(i) = ztopsw0(i)
       toplw0(i) = ztoplw0(i)
       solsw0(i) = zsolsw0(i)
       sollw0(i) = zsollw0(i)

       DO k = 1, klev + 1
          swdn0(i, k) = ZFSDN0(i, k)
          swdn(i, k) = ZFSDN(i, k)
          swup0(i, k) = ZFSUP0(i, k)
          swup(i, k) = ZFSUP(i, k)
       ENDDO
    ENDDO

    DO k = 1, klev
       DO i = 1, klon
          ! scale factor to take into account the difference
          ! between dry air and water vapour specific heat capacity
          zznormcp = 1. + RVTMP2 * PWV(i, k)
          heat(i, k) = zheat(i, k) / zznormcp
          cool(i, k) = zcool(i, k)/zznormcp
          heat0(i, k) = zheat0(i, k)/zznormcp
          cool0(i, k) = zcool0(i, k)/zznormcp
       ENDDO
    ENDDO

  END SUBROUTINE radlwsw

end module radlwsw_m
