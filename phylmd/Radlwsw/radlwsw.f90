module radlwsw_m

  IMPLICIT none

contains

  SUBROUTINE radlwsw(dist, mu0, fract, paprs, play, tsol, albsol, t_seri, &
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
    real, intent(in):: albsol(:) ! (klon)  alb\'edo du sol (entre 0 et 1)
    real, intent(in):: t_seri(:, :) ! (klon, klev) temp\'erature (K)
    real, intent(in):: q_seri(:, :) ! (klon, klev) vapeur d'eau (en kg / kg)

    real, intent(in):: wo(:, :) ! (klon, klev)
    ! column-density of ozone in a layer, in kilo-Dobsons

    real, intent(in):: cldfra(:, :) ! (klon, klev)
    ! fraction nuageuse (entre 0 et 1)

    real, intent(in):: cldemi(:, :) ! (klon, klev)
    ! \'emissivit\'e des nuages dans l'infrarouge (entre 0 et 1)

    real, intent(in):: cldtau(:, :) ! (klon, klev)
    ! \'epaisseur optique des nuages dans le visible (present-day value)

    real, intent(out):: heat(:, :) ! (klon, klev)
    ! \'echauffement atmosph\'erique (visible) (K / jour)

    real, intent(out):: heat0(:, :) ! (klon, klev) chauffage solaire ciel clair

    real, intent(out):: cool(:, :) ! (klon, klev)
    ! refroidissement dans l'infrarouge (K / jour)

    real, intent(out):: cool0(:, :) ! (klon, klev)
    ! refroidissement infrarouge ciel clair

    real, intent(out):: radsol(:) ! (klon)
    ! bilan radiatif net au sol (W m-2), positif vers le bas

    real, intent(out):: topsw(:) ! (klon)
    ! net incoming shortwave flux at top of atmosphere

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
    DOUBLE PRECISION FSUP(KLON, KLEV + 1)
    DOUBLE PRECISION FSDN(KLON, KLEV + 1)
    DOUBLE PRECISION FSUP0(KLON, KLEV + 1)
    DOUBLE PRECISION FSDN0(KLON, KLEV + 1)
    DOUBLE PRECISION FLUP(KLON, KLEV + 1)
    DOUBLE PRECISION FLDN(KLON, KLEV + 1)
    DOUBLE PRECISION FLUP0(KLON, KLEV + 1)
    DOUBLE PRECISION FLDN0(KLON, KLEV + 1)
    DOUBLE PRECISION alpha1
    INTEGER k, i
    DOUBLE PRECISION ALBD(klon, 2) ! second dimension is spectral band
    DOUBLE PRECISION EMIS(klon), VIEW(klon)
    DOUBLE PRECISION DP(klon, klev)
    DOUBLE PRECISION TL(klon, klev + 1), PMB(klon, klev + 1)
    DOUBLE PRECISION TAVE(klon, klev)
    DOUBLE PRECISION WV(klon, klev)
    DOUBLE PRECISION OZON(klon, klev) ! mass fraction of ozone
    DOUBLE PRECISION CLDLD(klon, klev)
    DOUBLE PRECISION TAU(klon, 2, klev)
    DOUBLE PRECISION OMEGA(klon, 2, klev)
    DOUBLE PRECISION CG(klon, 2, klev)
    DOUBLE PRECISION zheat(klon, klev), zcool(klon, klev)
    DOUBLE PRECISION zheat0(klon, klev), zcool0(klon, klev)
    DOUBLE PRECISION ztopsw(klon), ztoplw(klon)
    DOUBLE PRECISION zsolsw(klon), zsollw(klon)
    DOUBLE PRECISION zsollwdown(klon)
    DOUBLE PRECISION ztopsw0(klon), ztoplw0(klon)
    DOUBLE PRECISION zsolsw0(klon), zsollw0(klon)
    DOUBLE PRECISION zznormcp
    real, parameter:: dobson_u = 2.1415e-5 ! Dobson unit, in kg m-2

    !----------------------------------------------------------------------

    ALBD = spread(albsol, 2, 2)
    EMIS = 1.
    VIEW = 1.66

    DO i = 1, klon
       alpha1 = (paprs(i, 1) - play(i, 2)) / (play(i, 1) - play(i, 2))
       TL(i, 1) = t_seri(i, 1) * alpha1 + t_seri(i, 2) * (1. - alpha1)
    ENDDO

    TL(:, klev + 1) = t_seri(:, klev)
    forall (k = 2:klev) TL(:, k) = (t_seri(:, k) + t_seri(:, k - 1)) * 0.5
    forall (k = 1:klev) DP(:, k) = paprs(:, k) - paprs(:, k + 1)
    OZON = wo * RG * dobson_u * 1e3 / real(DP)
    TAVE = t_seri
    WV = MAX(q_seri, 1e-12)
    CLDLD = cldfra * cldemi

    TAU(:, 1, :) = MAX(cldtau, 1e-5)
    ! (1e-12 serait instable pour machines 32 bits)

    TAU(:, 2, :) = TAU(:, 1, :)
    OMEGA(:, 1, :) = 0.9999 - 5e-4 * EXP(- 0.5 * TAU(:, 1, :))
    OMEGA(:, 2, :) = 0.9988 - 2.5e-3 * EXP(- 0.05 * TAU(:, 2, :))
    CG(:, 1, :) = 0.865
    CG(:, 2, :) = 0.910
    PMB = paprs / 100.
    CALL LW(PMB, DP, tsol - TL(:, 1), EMIS, TL, TAVE, WV, OZON, CLDLD, CLDLd, &
         VIEW, zcool, zcool0, ztoplw, zsollw, ztoplw0, zsollw0, zsollwdown, &
         FLUP, FLDN, FLUP0, FLDN0)
    CALL SW(dble(solaire / dist**2), dble(mu0), dble(fract), PMB, DP, &
         dble(paprs(:, 1)), ALBD, ALBD, TAVE, WV, wv, OZON, dble(cldfra), TAU, &
         OMEGA, CG, zheat, zheat0, ztopsw, zsolsw, ztopsw0, zsolsw0, FSUP, &
         FSDN, FSUP0, FSDN0)
    radsol = zsolsw + zsollw
    topsw = ztopsw
    toplw = ztoplw
    solsw = zsolsw
    sollw = zsollw
    sollwdown = zsollwdown
    lwdn0 = FLDN0
    lwdn = FLDN
    lwup0 = FLUP0
    lwup = FLUP
    topsw0 = ztopsw0
    toplw0 = ztoplw0
    solsw0 = zsolsw0
    sollw0 = zsollw0
    swdn0 = FSDN0
    swdn = FSDN
    swup0 = FSUP0
    swup = FSUP

    DO k = 1, klev
       DO i = 1, klon
          ! scale factor to take into account the difference
          ! between dry air and water vapour specific heat capacity
          zznormcp = 1. + RVTMP2 * WV(i, k)
          heat(i, k) = zheat(i, k) / zznormcp
          cool(i, k) = zcool(i, k) / zznormcp
          heat0(i, k) = zheat0(i, k) / zznormcp
          cool0(i, k) = zcool0(i, k) / zznormcp
       ENDDO
    ENDDO

  END SUBROUTINE radlwsw

end module radlwsw_m
