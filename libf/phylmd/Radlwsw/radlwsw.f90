module radlwsw_m

  IMPLICIT none

contains

  SUBROUTINE radlwsw(dist, rmu0, fract, paprs, pplay, tsol, albedo, alblw, &
       t, q, wo, cldfra, cldemi, cldtaupd, heat, heat0, cool, cool0, radsol, &
       albpla, topsw, toplw, solsw, sollw, sollwdown, topsw0, toplw0, solsw0, &
       sollw0, lwdn0, lwdn, lwup0, lwup, swdn0, swdn, swup0, swup, ok_ade, &
       ok_aie, tau_ae, piz_ae, cg_ae, topswad, solswad, cldtaupi, topswai, &
       solswai)

    ! From LMDZ4/libf/phylmd/radlwsw.F, version 1.4 2005/06/06 13:16:33
    ! Author: Z. X. Li (LMD/CNRS) 
    ! Date: 1996/07/19

    ! Objet : interface entre le modèle et les rayonnements solaire et
    ! infrarouge

    ! ATTENTION: swai and swad have to be interpreted in the following manner:

    ! not ok_ade and not ok_aie
    ! both are zero

    ! ok_ade and not ok_aie
    ! aerosol direct forcing is F_{AD} = topsw - topswad
    ! indirect is zero

    ! not ok_ade and ok_aie
    ! aerosol indirect forcing is F_{AI} = topsw - topswai
    ! direct is zero

    ! ok_ade and ok_aie
    ! aerosol indirect forcing is F_{AI} = topsw - topswai
    ! aerosol direct forcing is F_{AD} = topswai - topswad

    USE clesphys, ONLY: bug_ozone, solaire
    USE dimphy, ONLY: klev, klon
    use lw_m, only: lw
    USE raddim, ONLY: kdlon
    USE suphec_m, ONLY: rg
    use sw_m, only: sw
    USE yoethf_m, ONLY: rvtmp2
        
    ! Arguments:

    real rmu0(klon), fract(klon), dist
    ! dist-----input-R- distance astronomique terre-soleil
    ! rmu0-----input-R- cosinus de l'angle zenithal
    ! fract----input-R- duree d'ensoleillement normalisee

    real, intent(in):: paprs(klon, klev+1)
    ! paprs----input-R- pression a inter-couche (Pa)
    real, intent(in):: pplay(klon, klev)
    ! pplay----input-R- pression au milieu de couche (Pa)
    real albedo(klon), alblw(klon), tsol(klon)
    ! albedo---input-R- albedo du sol (entre 0 et 1)
    ! tsol-----input-R- temperature du sol (en K)
    real, intent(in):: t(klon, klev)
    ! t--------input-R- temperature (K)
    real q(klon, klev)
    ! q--------input-R- vapeur d'eau (en kg/kg)
    real, intent(in):: wo(klon, klev)
    ! wo-------input-R- contenu en ozone (en kg/kg) correction MPL 100505
    real cldfra(klon, klev), cldemi(klon, klev)
    ! cldfra---input-R- fraction nuageuse (entre 0 et 1)
    ! cldemi---input-R- emissivite des nuages dans l'IR (entre 0 et 1)

    real cldtaupd(klon, klev)
    ! input-R- epaisseur optique des nuages dans le visible (present-day value)

    real, intent(out):: heat(klon, klev)
    ! échauffement atmosphérique (visible) (K/jour)

    real cool(klon, klev)
    ! cool-----output-R- refroidissement dans l'IR (K/jour)
    real heat0(klon, klev), cool0(klon, klev)
    real radsol(klon), topsw(klon)
    ! radsol---output-R- bilan radiatif net au sol (W/m**2) (+ vers le bas)
    ! topsw----output-R- flux solaire net au sommet de l'atm.

    real, intent(out):: toplw(klon)
    ! rayonnement infrarouge montant au sommet de l'atmosphère

    real solsw(klon), sollw(klon), albpla(klon)
    ! solsw----output-R- flux solaire net a la surface
    ! sollw----output-R- ray. IR montant a la surface
    ! albpla---output-R- albedo planetaire (entre 0 et 1)
    real topsw0(klon), solsw0(klon), sollw0(klon)
    real, intent(out):: toplw0(klon)
    real sollwdown(klon)
    !IM output 3D 
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
    DOUBLE PRECISION PWV(kdlon, klev), PQS(kdlon, klev), POZON(kdlon, klev)
    DOUBLE PRECISION PAER(kdlon, klev, 5)
    DOUBLE PRECISION PCLDLD(kdlon, klev)
    DOUBLE PRECISION PCLDLU(kdlon, klev)
    DOUBLE PRECISION PCLDSW(kdlon, klev)
    DOUBLE PRECISION PTAU(kdlon, 2, klev)
    DOUBLE PRECISION POMEGA(kdlon, 2, klev)
    DOUBLE PRECISION PCG(kdlon, 2, klev)

    DOUBLE PRECISION zfract(kdlon), zrmu0(kdlon), zdist

    DOUBLE PRECISION zheat(kdlon, klev), zcool(kdlon, klev)
    DOUBLE PRECISION zheat0(kdlon, klev), zcool0(kdlon, klev)
    DOUBLE PRECISION ztopsw(kdlon), ztoplw(kdlon)
    DOUBLE PRECISION zsolsw(kdlon), zsollw(kdlon), zalbpla(kdlon)
    DOUBLE PRECISION zsollwdown(kdlon)

    DOUBLE PRECISION ztopsw0(kdlon), ztoplw0(kdlon)
    DOUBLE PRECISION zsolsw0(kdlon), zsollw0(kdlon)
    DOUBLE PRECISION zznormcp
    !IM output 3D: SWup, SWdn, LWup, LWdn
    REAL swdn(klon, klev+1), swdn0(klon, klev+1)
    REAL swup(klon, klev+1), swup0(klon, klev+1)
    REAL lwdn(klon, klev+1), lwdn0(klon, klev+1)
    REAL lwup(klon, klev+1), lwup0(klon, klev+1)

    !jq the following quantities are needed for the aerosol radiative forcings

    real topswad(klon), solswad(klon)
    ! output: aerosol direct forcing at TOA and surface
    ! topswad---output-R- ray. solaire absorbe au sommet de l'atm. (aerosol dir)
    ! solswad---output-R- ray. solaire net absorbe a la surface (aerosol dir)

    real topswai(klon), solswai(klon)
    ! output: aerosol indirect forcing atTOA and surface
    ! topswai---output-R- ray. solaire absorbe au sommet de l'atm. (aerosol ind)
    ! solswai---output-R- ray. solaire net absorbe a la surface (aerosol ind)

    real tau_ae(klon, klev, 2), piz_ae(klon, klev, 2), cg_ae(klon, klev, 2)
    ! input-R- aerosol optical properties (calculated in aeropt.F)

    real cldtaupi(klon, klev)
    ! cloud optical thickness for pre-industrial aerosol concentrations
    ! (i.e. with a smaller droplet concentration and thus larger droplet radii)
    ! -input-R- epaisseur optique des nuages dans le visible
    ! calculated for pre-industrial (pi) aerosol concentrations,
    ! i.e. with smaller droplet concentration, thus larger droplets,
    ! thus generally cdltaupi cldtaupd it is needed for the
    ! diagnostics of the aerosol indirect radiative forcing

    logical ok_ade, ok_aie 
    ! switches whether to use aerosol direct (indirect) effects or not
    ! ok_ade---input-L- apply the Aerosol Direct Effect or not?
    ! ok_aie---input-L- apply the Aerosol Indirect Effect or not?

    double precision tauae(kdlon, klev, 2) ! aer opt properties
    double precision pizae(kdlon, klev, 2)
    double precision cgae(kdlon, klev, 2)

    DOUBLE PRECISION PTAUA(kdlon, 2, klev)
    ! present-day value of cloud opt thickness (PTAU is pre-industrial
    ! value), local use

    DOUBLE PRECISION POMEGAA(kdlon, 2, klev) ! dito for single scatt albedo

    DOUBLE PRECISION ztopswad(kdlon), zsolswad(kdlon) 
    ! Aerosol direct forcing at TOAand surface

    DOUBLE PRECISION ztopswai(kdlon), zsolswai(kdlon) ! dito, indirect

    !----------------------------------------------------------------------

    tauae = 0.
    pizae = 0.
    cgae = 0.

    nb_gr = klon / kdlon
    IF (nb_gr * kdlon /= klon) THEN
       PRINT *, "kdlon mauvais :", klon, kdlon, nb_gr
       stop 1
    ENDIF
    
    heat = 0.
    cool = 0.
    heat0 = 0.
    cool0 = 0.
    zdist = dist
    PSCT = solaire / zdist / zdist

    loop_iof: DO iof = 0, klon - kdlon, kdlon
       DO i = 1, kdlon
          zfract(i) = fract(iof+i)
          zrmu0(i) = rmu0(iof+i)
          PALBD(i, 1) = albedo(iof+i)
          PALBD(i, 2) = alblw(iof+i)
          PALBP(i, 1) = albedo(iof+i)
          PALBP(i, 2) = alblw(iof+i)
          ! cf. JLD pour etre en accord avec ORCHIDEE il faut mettre
          ! PEMIS(i) = 0.96
          PEMIS(i) = 1.0 
          PVIEW(i) = 1.66
          PPSOL(i) = paprs(iof+i, 1)
          zx_alpha1 = (paprs(iof+i, 1)-pplay(iof+i, 2))  &
               / (pplay(iof+i, 1)-pplay(iof+i, 2))
          zx_alpha2 = 1.0 - zx_alpha1
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
             PWV(i, k) = MAX (q(iof+i, k), 1.0e-12)
             PQS(i, k) = PWV(i, k)
             ! wo:    cm.atm (epaisseur en cm dans la situation standard)
             ! POZON: kg/kg
             IF (bug_ozone) then
                POZON(i, k) = MAX(wo(iof+i, k), 1.0e-12)*RG/46.6968 &
                     /(paprs(iof+i, k)-paprs(iof+i, k+1)) &
                     *(paprs(iof+i, 1)/101325.0)
             ELSE
                ! le calcul qui suit est maintenant fait dans ozonecm (MPL)
                POZON(i, k) = wo(i, k)
             ENDIF
             PCLDLD(i, k) = cldfra(iof+i, k)*cldemi(iof+i, k)
             PCLDLU(i, k) = cldfra(iof+i, k)*cldemi(iof+i, k)
             PCLDSW(i, k) = cldfra(iof+i, k)
             PTAU(i, 1, k) = MAX(cldtaupi(iof+i, k), 1.0e-05)
             ! (1e-12 serait instable)
             PTAU(i, 2, k) = MAX(cldtaupi(iof+i, k), 1.0e-05)
             ! (pour 32-bit machines)
             POMEGA(i, 1, k) = 0.9999 - 5.0e-04 * EXP(-0.5 * PTAU(i, 1, k))
             POMEGA(i, 2, k) = 0.9988 - 2.5e-03 * EXP(-0.05 * PTAU(i, 2, k))
             PCG(i, 1, k) = 0.865
             PCG(i, 2, k) = 0.910

             ! Introduced for aerosol indirect forcings.  The
             ! following values use the cloud optical thickness
             ! calculated from present-day aerosol concentrations
             ! whereas the quantities without the "A" at the end are
             ! for pre-industial (natural-only) aerosol concentrations
             PTAUA(i, 1, k) = MAX(cldtaupd(iof+i, k), 1.0e-05)
             ! (1e-12 serait instable)
             PTAUA(i, 2, k) = MAX(cldtaupd(iof+i, k), 1.0e-05)
             ! (pour 32-bit machines)
             POMEGAA(i, 1, k) = 0.9999 - 5.0e-04 * EXP(-0.5 * PTAUA(i, 1, k))
             POMEGAA(i, 2, k) = 0.9988 - 2.5e-03 * EXP(-0.05 * PTAUA(i, 2, k))
             !jq-end
          ENDDO
       ENDDO

       DO k = 1, klev+1
          DO i = 1, kdlon
             PPMB(i, k) = paprs(iof+i, k)/100.0
          ENDDO
       ENDDO

       DO kk = 1, 5
          DO k = 1, klev
             DO i = 1, kdlon
                PAER(i, k, kk) = 1.0E-15
             ENDDO
          ENDDO
       ENDDO

       DO k = 1, klev
          DO i = 1, kdlon
             tauae(i, k, 1) = tau_ae(iof+i, k, 1)
             pizae(i, k, 1) = piz_ae(iof+i, k, 1)
             cgae(i, k, 1) =cg_ae(iof+i, k, 1)
             tauae(i, k, 2) = tau_ae(iof+i, k, 2)
             pizae(i, k, 2) = piz_ae(iof+i, k, 2)
             cgae(i, k, 2) =cg_ae(iof+i, k, 2)
          ENDDO
       ENDDO

       CALL LW(PPMB, PDP, PPSOL, PDT0, PEMIS, PTL, PTAVE, PWV, POZON, PAER, &
            PCLDLD, PCLDLU, PVIEW, zcool, zcool0, ztoplw, zsollw, ztoplw0, &
            zsollw0, zsollwdown, ZFLUP, ZFLDN, ZFLUP0, ZFLDN0)
       CALL SW(PSCT, zrmu0, zfract, PPMB, PDP, PPSOL, PALBD, PALBP, PTAVE, &
            PWV, PQS, POZON, PAER, PCLDSW, PTAU, POMEGA, PCG, zheat, zheat0, &
            zalbpla, ztopsw, zsolsw, ztopsw0, zsolsw0, ZFSUP, ZFSDN, ZFSUP0, &
            ZFSDN0, tauae, pizae, cgae, PTAUA, POMEGAA, ztopswad, zsolswad, &
            ztopswai, zsolswai, ok_ade, ok_aie)

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
             topswad(iof+i) = 0.0
             solswad(iof+i) = 0.0
          ENDDO
       ENDIF
       IF (ok_aie) THEN
          DO i = 1, kdlon
             topswai(iof+i) = ztopswai(i)
             solswai(iof+i) = zsolswai(i)
          ENDDO
       ELSE
          DO i = 1, kdlon
             topswai(iof+i) = 0.0
             solswai(iof+i) = 0.0
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
