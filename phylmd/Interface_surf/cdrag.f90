module cdrag_m

  IMPLICIT NONE

contains

  SUBROUTINE cdrag(nsrf, speed, t, q, zgeop, psol, ts, qsurf, rugos, cdragm, &
       cdragh, pref)

    ! From LMDZ4/libf/phylmd/clcdrag.F90 and
    ! LMDZ4/libf/phylmd/coefcdrag.F90, version 1.1.1.1, 2004 May 19th

    ! Objet : calcul des drag coefficients au sol pour le moment et
    ! les flux de chaleurs sensible et latente et calcul de la
    ! pression au niveau de r\'ef\'erence.

    ! Ionela MUSAT, July, 1st, 2002

    ! Louis, J. F., Tiedtke, M. and Geleyn, J. F., 1982. A short
    ! history of the operational PBL parametrization at
    ! ECMWF. Workshop on boundary layer parametrization, November
    ! 1981, ECMWF, Reading, England. Page: 19. Equations in Table 1.

    ! Miller, M. J., A. C. M. Beljaars, T. N. Palmer, 1992. The
    ! sensitivity of the ECMWF model to the parameterization of
    ! evaporation from the tropical oceans. J. Climate, 5:418-434.

    ! Library:
    use nr_util, only: assert_eq

    use clesphys, only: f_cdrag_oce, f_cdrag_ter
    use indicesol, only: is_oce
    use SUPHEC_M, only: rcpd, rd, retv, rg
    USE yoethf_m, ONLY: rvtmp2

    INTEGER, intent(in):: nsrf ! indice pour le type de surface

    REAL, intent(in):: speed(:) ! (knon)
    ! norm of the wind in the first model layer, in m s-1

    REAL, intent(in):: t(:) ! (knon)
    ! temperature de l'air au 1er niveau du modele

    REAL, intent(in):: q(:) ! (knon) ! humidite de l'air au 1er niveau du modele

    REAL, intent(in):: zgeop(:) ! (knon)
    ! g\'eopotentiel au 1er niveau du mod\`ele

    REAL, intent(in) :: psol(:) ! (knon) pression au sol
    REAL, intent(in):: ts(:) ! (knon) temperature de l'air a la surface
    REAL, intent(in):: qsurf(:) ! (knon) humidite de l'air a la surface
    REAL, intent(in):: rugos(:) ! (knon) rugosit\'e
    REAL, intent(out):: cdragm(:) ! (knon) drag coefficient for momentum

    REAL, intent(out):: cdragh(:) ! (knon)
    ! drag coefficient for latent and sensible heat fluxes

    REAL, intent(out), optional:: pref(:) ! (knon) pression au niveau zgeop / RG

    ! Local:

    REAL, PARAMETER:: ckap = 0.4, cb = 5., cc = 5., cd = 5.
    REAL, PARAMETER:: cepdu2 = 0.1**2 ! in m2 s-2
    real, parameter:: f_ri_cd_min = 0.1
    INTEGER i, knon
    REAL du2 ! in m2 s-2
    real tsolv, tvd, zscf, zucf
    real cdn ! drag coefficient for neutral conditions

    REAL ri
    ! nombre de Richardson entre la surface et le niveau de reference
    ! zgeop / RG

    !-------------------------------------------------------------------------

    knon = assert_eq([size(speed), size(t), size(q), size(zgeop), size(ts), &
         size(qsurf), size(rugos), size(cdragm), size(cdragh)], "cdrag knon")

    DO i = 1, knon
       du2 = max(cepdu2, speed(i)**2)
       tsolv = ts(i) * (1. + RETV * max(qsurf(i), 0.))
       tvd = (t(i) + zgeop(i) / RCPD / (1. + RVTMP2 * q(i))) &
            * (1. + RETV * q(i))
       ri = zgeop(i) * (tvd - tsolv) / (du2 * tvd)
       cdn = (ckap / log(1. + zgeop(i) / (RG * rugos(i))))**2

       IF (ri < 0.) THEN
          ! Situation instable
          zucf = 1. / (1. + 3. * cb * cc * cdn &
               * SQRT(ABS(ri) * (1. + zgeop(i) / (RG * rugos(i)))))
          cdragm(i) = cdn * max((1. - 2. * cb * ri * zucf), f_ri_cd_min)

          IF (nsrf == is_oce) then
             ! Cf. Miller et al. (1992)
             cdragh(i) = f_cdrag_oce * cdn * (1. + ((0.0016 / (cdn &
                  * SQRT(du2))) &
                  * ABS(tvd - tsolv)**(1. / 3.))**1.25)**(1. / 1.25)
          else
             cdragh(i) = f_cdrag_ter * cdn &
                  * max((1. - 3. * cb * ri * zucf), f_ri_cd_min)
          end IF
       ELSE
          ! Situation stable. Pour \'eviter les incoh\'erences dans
          ! les cas tr\`es stables, on limite ri \`a 20. Cf Hess et
          ! al. (1995).
          ri = min(20., ri)
          zscf = SQRT(1. + cd * ABS(ri))
          cdragm(i) = cdn * max(1. / (1. + 2. * CB * ri / zscf), f_ri_cd_min)
          cdragh(i) = merge(f_cdrag_oce, f_cdrag_ter, nsrf == is_oce) * cdn &
               * max(1. / (1. + 3. * CB * ri * zscf), f_ri_cd_min)
       ENDIF
    END DO

    if (present(pref)) &
         pref = psol * exp(- zgeop / (RD * t * (1. + RETV * max(q, 0.))))

  END SUBROUTINE cdrag

end module cdrag_m
