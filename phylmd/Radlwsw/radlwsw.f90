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
    real, intent(in):: mu0(:) ! (klon) cosinus de l'angle z\'enithal
    real, intent(in):: fract(:) ! (klon) dur\'ee d'ensoleillement normalis\'ee

    real, intent(in):: paprs(:, :) ! (klon, klev + 1)
    ! pression \`a inter-couche (Pa)

    real, intent(in):: play(:, :) ! (klon, klev)
    ! pression au milieu de couche (Pa)

    real, intent(in):: tsol(:) ! (klon)  temp\'erature du sol (en K)
    real, intent(in):: albedo(:) ! (klon)  alb\'edo du sol (entre 0 et 1)
    real, intent(in):: t_seri(:, :) ! (klon, klev) temp\'erature (K)
    real, intent(in):: q_seri(:, :) ! (klon, klev) vapeur d'eau (en kg/kg)

    real, intent(in):: wo(:, :) ! (klon, klev)
    ! column-density of ozone in a layer, in kilo-Dobsons

    real, intent(in):: cldfra(:, :) ! (klon, klev)
    ! fraction nuageuse (entre 0 et 1)

    real, intent(in):: cldemi(:, :) ! (klon, klev)
    ! \'emissivit\'e des nuages dans l'infrarouge (entre 0 et 1)

    real, intent(in):: cldtau(:, :) ! (klon, klev)
    ! \'epaisseur optique des nuages dans le visible (present-day value)

    real, intent(out):: heat(:, :) ! (klon, klev)
    ! \'echauffement atmosph\'erique (visible) (K/jour)

    real, intent(out):: heat0(:, :) ! (klon, klev) chauffage solaire ciel clair

    real, intent(out):: cool(:, :) ! (klon, klev)
    ! refroidissement dans l'infrarouge (K/jour)

    real, intent(out):: cool0(:, :) ! (klon, klev)
    ! refroidissement infrarouge ciel clair

    real, intent(out):: radsol(:) ! (klon)
    ! bilan radiatif net au sol (W m-2), positif vers le bas

    real, intent(out):: topsw(:) ! (klon)
    ! flux solaire net au sommet de l'atmosph`ere

    real, intent(out):: toplw(:) ! (klon)
    ! rayonnement infrarouge montant au sommet de l'atmosph\`ere

    real, intent(out):: solsw(:) ! (klon) flux solaire net \`a la surface

    real, intent(out):: sollw(:) ! (klon)
    ! rayonnement infrarouge net \`a la surface

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
    DOUBLE PRECISION alpha1
    INTEGER k, i
    DOUBLE PRECISION SCT
    DOUBLE PRECISION ALBD(klon, 2), ALBP(klon, 2)
    DOUBLE PRECISION EMIS(klon), DT0(klon), VIEW(klon)
    DOUBLE PRECISION PSOL(klon), PDP(klon, klev)
    DOUBLE PRECISION TL(klon, klev + 1), PPMB(klon, klev + 1)
    DOUBLE PRECISION PTAVE(klon, klev)
    DOUBLE PRECISION PWV(klon, klev), PQS(klon, klev)
    DOUBLE PRECISION POZON(klon, klev) ! mass fraction of ozone
    DOUBLE PRECISION PCLDLD(klon, klev)
    DOUBLE PRECISION PCLDLU(klon, klev)
    DOUBLE PRECISION PCLDSW(klon, klev)
    DOUBLE PRECISION PTAU(klon, 2, klev)
    DOUBLE PRECISION POMEGA(klon, 2, klev)
    DOUBLE PRECISION PCG(klon, 2, klev)
    DOUBLE PRECISION rmu0(klon)
    DOUBLE PRECISION zheat(klon, klev), zcool(klon, klev)
    DOUBLE PRECISION zheat0(klon, klev), zcool0(klon, klev)
    DOUBLE PRECISION ztopsw(klon), ztoplw(klon)
    DOUBLE PRECISION zsolsw(klon), zsollw(klon)
    DOUBLE PRECISION zsollwdown(klon)
    DOUBLE PRECISION ztopsw0(klon), ztoplw0(klon)
    DOUBLE PRECISION zsolsw0(klon), zsollw0(klon)
    DOUBLE PRECISION zznormcp
    real, parameter:: dobson_u = 2.1415e-05 ! Dobson unit, in kg m-2

    !----------------------------------------------------------------------

    heat = 0.
    cool = 0.
    heat0 = 0.
    cool0 = 0.
    SCT = solaire / dist**2
    rmu0 = mu0
    ALBD(:, 1) = albedo
    ALBD(:, 2) = albedo
    ALBP(:, 1) = albedo
    ALBP(:, 2) = albedo
    ! cf. JLD pour etre en accord avec ORCHIDEE il faut mettre
    ! EMIS = 0.96
    EMIS = 1.
    VIEW = 1.66
    PSOL = paprs(:, 1)

    DO i = 1, klon
       alpha1 = (paprs(i, 1) - play(i, 2)) / (play(i, 1) - play(i, 2))
       TL(i, 1) = t_seri(i, 1) * alpha1 + t_seri(i, 2) * (1. - alpha1)
    ENDDO

    TL(:, klev + 1) = t_seri(:, klev)
    DT0 = tsol - TL(:, 1)
    forall (k = 2:klev) TL(:, k) = (t_seri(:, k) + t_seri(:, k - 1)) * 0.5

    DO k = 1, klev
       DO i = 1, klon
          PDP(i, k) = paprs(i, k) - paprs(i, k + 1)
          PTAVE(i, k) = t_seri(i, k)
          PWV(i, k) = MAX(q_seri(i, k), 1e-12)
          PQS(i, k) = PWV(i, k)
          POZON(i, k) = wo(i, k) * RG * dobson_u * 1e3 &
               / (paprs(i, k) - paprs(i, k + 1))
          PCLDLD(i, k) = cldfra(i, k) * cldemi(i, k)
          PCLDLU(i, k) = cldfra(i, k) * cldemi(i, k)
          PCLDSW(i, k) = cldfra(i, k)
          PTAU(i, 1, k) = MAX(cldtau(i, k), 1e-05)
          ! (1e-12 serait instable)
          PTAU(i, 2, k) = MAX(cldtau(i, k), 1e-05)
          ! (pour 32-bit machines)
          POMEGA(i, 1, k) = 0.9999 - 5e-04 * EXP(- 0.5 * PTAU(i, 1, k))
          POMEGA(i, 2, k) = 0.9988 - 2.5e-03 * EXP(- 0.05 * PTAU(i, 2, k))
          PCG(i, 1, k) = 0.865
          PCG(i, 2, k) = 0.910
       ENDDO
    ENDDO

    DO k = 1, klev + 1
       DO i = 1, klon
          PPMB(i, k) = paprs(i, k)/100.
       ENDDO
    ENDDO

    CALL LW(PPMB, PDP, DT0, EMIS, TL, PTAVE, PWV, POZON, PCLDLD, &
         PCLDLU, VIEW, zcool, zcool0, ztoplw, zsollw, ztoplw0, zsollw0, &
         zsollwdown, ZFLUP, ZFLDN, ZFLUP0, ZFLDN0)
    CALL SW(SCT, rmu0, dble(fract), PPMB, PDP, PSOL, ALBD, ALBP, PTAVE, &
         PWV, PQS, POZON, PCLDSW, PTAU, POMEGA, PCG, zheat, zheat0, ztopsw, &
         zsolsw, ztopsw0, zsolsw0, ZFSUP, ZFSDN, ZFSUP0, ZFSDN0)

    radsol = zsolsw + zsollw
    topsw = ztopsw
    toplw = ztoplw
    solsw = zsolsw
    sollw = zsollw
    sollwdown = zsollwdown
    lwdn0 = ZFLDN0
    lwdn = ZFLDN
    lwup0 = ZFLUP0
    lwup = ZFLUP
    topsw0 = ztopsw0
    toplw0 = ztoplw0
    solsw0 = zsolsw0
    sollw0 = zsollw0
    swdn0 = ZFSDN0
    swdn = ZFSDN
    swup0 = ZFSUP0
    swup = ZFSUP

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
