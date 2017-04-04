module radlwsw_m

  IMPLICIT none

contains

  SUBROUTINE radlwsw(dist, mu0, fract, paprs, play, tsol, albedo, t, q, wo, &
       cldfra, cldemi, cldtaupd, heat, heat0, cool, cool0, radsol, albpla, &
       topsw, toplw, solsw, sollw, sollwdown, topsw0, toplw0, solsw0, sollw0, &
       lwdn0, lwdn, lwup0, lwup, swdn0, swdn, swup0, swup, ok_ade, topswad, &
       solswad)

    ! From LMDZ4/libf/phylmd/radlwsw.F, version 1.4, 2005/06/06 13:16:33
    ! Author: Z. X. Li (LMD/CNRS) 
    ! Date: 1996/07/19

    ! Objet : interface entre le modèle et les rayonnements solaire et
    ! infrarouge

    ! ATTENTION: swad has to be interpreted in the following manner:
    ! not ok_ade zero
    ! ok_ade aerosol direct forcing is F_{AD} = topsw - topswad

    USE clesphys, ONLY: solaire
    USE dimphy, ONLY: klev, klon
    use lw_m, only: lw
    USE raddim, ONLY: kdlon
    USE suphec_m, ONLY: rg
    use sw_m, only: sw
    USE yoethf_m, ONLY: rvtmp2

    real, intent(in):: dist ! distance astronomique terre-soleil
    real, intent(in):: mu0(klon) ! cosinus de l'angle zenithal
    real, intent(in):: fract(klon) ! duree d'ensoleillement normalisee
    real, intent(in):: paprs(klon, klev+1) ! pression a inter-couche (Pa)
    real, intent(in):: play(klon, klev) ! pression au milieu de couche (Pa)
    real, intent(in):: tsol(klon) ! temperature du sol (en K)
    real, intent(in):: albedo(klon) ! albedo du sol (entre 0 et 1)
    real, intent(in):: t(klon, klev) ! temperature (K)
    real, intent(in):: q(klon, klev) ! vapeur d'eau (en kg/kg)

    real, intent(in):: wo(klon, klev)
    ! column-density of ozone in a layer, in kilo-Dobsons

    real, intent(in):: cldfra(klon, klev) ! fraction nuageuse (entre 0 et 1)

    real, intent(in):: cldemi(klon, klev)
    ! emissivite des nuages dans l'IR (entre 0 et 1)

    real, intent(in):: cldtaupd(klon, klev)
    ! epaisseur optique des nuages dans le visible (present-day value)

    real, intent(out):: heat(klon, klev)
    ! échauffement atmosphérique (visible) (K/jour)

    real, intent(out):: heat0(klon, klev) ! chauffage solaire ciel clair
    real, intent(out):: cool(klon, klev) ! refroidissement dans l'IR (K/jour)

    real, intent(out):: cool0(klon, klev)
    ! refroidissement infrarouge ciel clair

    real, intent(out):: radsol(klon)
    ! bilan radiatif net au sol (W/m**2) (+ vers le bas)

    real, intent(out):: albpla(klon) ! albedo planetaire (entre 0 et 1)
    real, intent(out):: topsw(klon) ! flux solaire net au sommet de l'atm.

    real, intent(out):: toplw(klon)
    ! rayonnement infrarouge montant au sommet de l'atmosphère

    real, intent(out):: solsw(klon) ! flux solaire net à la surface

    real, intent(out):: sollw(klon)
    ! rayonnement infrarouge montant à la surface

    real, intent(out):: sollwdown(klon)
    real, intent(out):: topsw0(klon)
    real, intent(out):: toplw0(klon)
    real, intent(out):: solsw0(klon), sollw0(klon)
    REAL, intent(out):: lwdn0(klon, klev+1), lwdn(klon, klev+1)
    REAL, intent(out):: lwup0(klon, klev+1), lwup(klon, klev+1)
    REAL, intent(out):: swdn0(klon, klev+1), swdn(klon, klev+1)
    REAL, intent(out):: swup0(klon, klev+1), swup(klon, klev+1)

    logical, intent(in):: ok_ade ! apply the Aerosol Direct Effect

    real, intent(out):: topswad(klon), solswad(klon)
    ! aerosol direct forcing at TOA and surface
    ! rayonnement solaire net absorb\'e

    ! Local:

    DOUBLE PRECISION ZFSUP(KDLON, KLEV+1)
    DOUBLE PRECISION ZFSDN(KDLON, KLEV+1)
    DOUBLE PRECISION ZFSUP0(KDLON, KLEV+1)
    DOUBLE PRECISION ZFSDN0(KDLON, KLEV+1)

    DOUBLE PRECISION ZFLUP(KDLON, KLEV+1)
    DOUBLE PRECISION ZFLDN(KDLON, KLEV+1)
    DOUBLE PRECISION ZFLUP0(KDLON, KLEV+1)
    DOUBLE PRECISION ZFLDN0(KDLON, KLEV+1)

    DOUBLE PRECISION zx_alpha1, zx_alpha2
    INTEGER k, kk, i, iof, nb_gr
    DOUBLE PRECISION PSCT

    DOUBLE PRECISION PALBD(kdlon, 2), PALBP(kdlon, 2)
    DOUBLE PRECISION PEMIS(kdlon), PDT0(kdlon), PVIEW(kdlon)
    DOUBLE PRECISION PPSOL(kdlon), PDP(kdlon, klev)
    DOUBLE PRECISION PTL(kdlon, klev+1), PPMB(kdlon, klev+1)
    DOUBLE PRECISION PTAVE(kdlon, klev)
    DOUBLE PRECISION PWV(kdlon, klev), PQS(kdlon, klev)
    DOUBLE PRECISION POZON(kdlon, klev) ! mass fraction of ozone
    DOUBLE PRECISION PAER(kdlon, klev, 5) ! AEROSOLS' OPTICAL THICKNESS
    DOUBLE PRECISION PCLDLD(kdlon, klev)
    DOUBLE PRECISION PCLDLU(kdlon, klev)
    DOUBLE PRECISION PCLDSW(kdlon, klev)
    DOUBLE PRECISION PTAU(kdlon, 2, klev)
    DOUBLE PRECISION POMEGA(kdlon, 2, klev)
    DOUBLE PRECISION PCG(kdlon, 2, klev)

    DOUBLE PRECISION zfract(kdlon), zrmu0(kdlon)

    DOUBLE PRECISION zheat(kdlon, klev), zcool(kdlon, klev)
    DOUBLE PRECISION zheat0(kdlon, klev), zcool0(kdlon, klev)
    DOUBLE PRECISION ztopsw(kdlon), ztoplw(kdlon)
    DOUBLE PRECISION zsolsw(kdlon), zsollw(kdlon), zalbpla(kdlon)
    DOUBLE PRECISION zsollwdown(kdlon)

    DOUBLE PRECISION ztopsw0(kdlon), ztoplw0(kdlon)
    DOUBLE PRECISION zsolsw0(kdlon), zsollw0(kdlon)
    DOUBLE PRECISION zznormcp

    !jq the following quantities are needed for the aerosol radiative forcings
    DOUBLE PRECISION ztopswad(kdlon), zsolswad(kdlon) 
    ! Aerosol direct forcing at TOA and surface

    real, parameter:: dobson_u = 2.1415e-05 ! Dobson unit, in kg m-2

    !----------------------------------------------------------------------

    nb_gr = klon / kdlon
    IF (nb_gr * kdlon /= klon) THEN
       PRINT *, "kdlon mauvais :", klon, kdlon, nb_gr
       stop 1
    ENDIF

    heat = 0.
    cool = 0.
    heat0 = 0.
    cool0 = 0.
    PSCT = solaire / dist**2

    loop_iof: DO iof = 0, klon - kdlon, kdlon
       DO i = 1, kdlon
          zfract(i) = fract(iof+i)
          zrmu0(i) = mu0(iof+i)
          PALBD(i, 1) = albedo(iof+i)
          PALBD(i, 2) = albedo(iof+i)
          PALBP(i, 1) = albedo(iof+i)
          PALBP(i, 2) = albedo(iof+i)
          ! cf. JLD pour etre en accord avec ORCHIDEE il faut mettre
          ! PEMIS(i) = 0.96
          PEMIS(i) = 1. 
          PVIEW(i) = 1.66
          PPSOL(i) = paprs(iof+i, 1)
          zx_alpha1 = (paprs(iof+i, 1)-play(iof+i, 2))  &
               / (play(iof+i, 1)-play(iof+i, 2))
          zx_alpha2 = 1. - zx_alpha1
          PTL(i, 1) = t(iof+i, 1) * zx_alpha1 + t(iof+i, 2) * zx_alpha2
          PTL(i, klev+1) = t(iof+i, klev)
          PDT0(i) = tsol(iof+i) - PTL(i, 1)
       ENDDO
       DO k = 2, klev
          DO i = 1, kdlon
             PTL(i, k) = (t(iof+i, k)+t(iof+i, k-1))*0.5
          ENDDO
       ENDDO
       DO k = 1, klev
          DO i = 1, kdlon
             PDP(i, k) = paprs(iof+i, k)-paprs(iof+i, k+1)
             PTAVE(i, k) = t(iof+i, k)
             PWV(i, k) = MAX (q(iof+i, k), 1e-12)
             PQS(i, k) = PWV(i, k)
             POZON(i, k) = wo(iof+i, k) * RG * dobson_u * 1e3 &
                  / (paprs(iof+i, k) - paprs(iof+i, k+1))
             PCLDLD(i, k) = cldfra(iof+i, k)*cldemi(iof+i, k)
             PCLDLU(i, k) = cldfra(iof+i, k)*cldemi(iof+i, k)
             PCLDSW(i, k) = cldfra(iof+i, k)
             PTAU(i, 1, k) = MAX(cldtaupd(iof+i, k), 1e-05)
             ! (1e-12 serait instable)
             PTAU(i, 2, k) = MAX(cldtaupd(iof+i, k), 1e-05)
             ! (pour 32-bit machines)
             POMEGA(i, 1, k) = 0.9999 - 5e-04 * EXP(-0.5 * PTAU(i, 1, k))
             POMEGA(i, 2, k) = 0.9988 - 2.5e-03 * EXP(-0.05 * PTAU(i, 2, k))
             PCG(i, 1, k) = 0.865
             PCG(i, 2, k) = 0.910
          ENDDO
       ENDDO

       DO k = 1, klev+1
          DO i = 1, kdlon
             PPMB(i, k) = paprs(iof+i, k)/100.
          ENDDO
       ENDDO

       DO kk = 1, 5
          DO k = 1, klev
             DO i = 1, kdlon
                PAER(i, k, kk) = 1E-15
             ENDDO
          ENDDO
       ENDDO

       CALL LW(PPMB, PDP, PDT0, PEMIS, PTL, PTAVE, PWV, POZON, PAER, PCLDLD, &
            PCLDLU, PVIEW, zcool, zcool0, ztoplw, zsollw, ztoplw0, zsollw0, &
            zsollwdown, ZFLUP, ZFLDN, ZFLUP0, ZFLDN0)
       CALL SW(PSCT, zrmu0, zfract, PPMB, PDP, PPSOL, PALBD, PALBP, PTAVE, &
            PWV, PQS, POZON, PCLDSW, PTAU, POMEGA, PCG, zheat, zheat0, &
            zalbpla, ztopsw, zsolsw, ztopsw0, zsolsw0, ZFSUP, ZFSDN, ZFSUP0, &
            ZFSDN0, ztopswad, zsolswad, ok_ade)

       DO i = 1, kdlon
          radsol(iof+i) = zsolsw(i) + zsollw(i)
          topsw(iof+i) = ztopsw(i)
          toplw(iof+i) = ztoplw(i)
          solsw(iof+i) = zsolsw(i)
          sollw(iof+i) = zsollw(i)
          sollwdown(iof+i) = zsollwdown(i)

          DO k = 1, klev+1
             lwdn0 ( iof+i, k)   = ZFLDN0 ( i, k)
             lwdn  ( iof+i, k)   = ZFLDN  ( i, k)
             lwup0 ( iof+i, k)   = ZFLUP0 ( i, k)
             lwup  ( iof+i, k)   = ZFLUP  ( i, k)
          ENDDO

          topsw0(iof+i) = ztopsw0(i)
          toplw0(iof+i) = ztoplw0(i)
          solsw0(iof+i) = zsolsw0(i)
          sollw0(iof+i) = zsollw0(i)
          albpla(iof+i) = zalbpla(i)

          DO k = 1, klev+1
             swdn0 ( iof+i, k)   = ZFSDN0 ( i, k)
             swdn  ( iof+i, k)   = ZFSDN  ( i, k)
             swup0 ( iof+i, k)   = ZFSUP0 ( i, k)
             swup  ( iof+i, k)   = ZFSUP  ( i, k)
          ENDDO
       ENDDO
       ! transform the aerosol forcings, if they have to be calculated
       IF (ok_ade) THEN
          DO i = 1, kdlon
             topswad(iof+i) = ztopswad(i)
             solswad(iof+i) = zsolswad(i)
          ENDDO
       ELSE
          DO i = 1, kdlon
             topswad(iof+i) = 0.
             solswad(iof+i) = 0.
          ENDDO
       ENDIF

       DO k = 1, klev
          DO i = 1, kdlon
             ! scale factor to take into account the difference
             ! between dry air and water vapour specific heat capacity
             zznormcp = 1. + RVTMP2 * PWV(i, k)
             heat(iof+i, k) = zheat(i, k) / zznormcp
             cool(iof+i, k) = zcool(i, k)/zznormcp
             heat0(iof+i, k) = zheat0(i, k)/zznormcp
             cool0(iof+i, k) = zcool0(i, k)/zznormcp
          ENDDO
       ENDDO
    end DO loop_iof

  END SUBROUTINE radlwsw

end module radlwsw_m
