module cdrag_m

  IMPLICIT NONE

contains

  SUBROUTINE cdrag(nsrf, speed, t, q, zgeop, psol, ts, qsurf, rugos, cdragm, &
       cdragh, pref)

    ! From LMDZ4/libf/phylmd/clcdrag.F90 and
    ! LMDZ4/libf/phylmd/coefcdrag.F90, version 1.1.1.1, 2004/05/19
    ! 12:53:07

    ! Objet : calcul des drag coefficients au sol pour le moment et
    ! les flux de chaleurs sensible et latente et calcul de la
    ! pression au niveau de reference.

    ! Ionela MUSAT, July, 1st, 2002

    ! Louis, J. F., Tiedtke, M. and Geleyn, J. F., 1982. A short
    ! history of the operational PBL parametrization at
    ! ECMWF. Workshop on boundary layer parametrization, November
    ! 1981, ECMWF, Reading, England. Page: 19. Equations in Table 1.

    ! Miller, M. J., A. C. M. Beljaars, T. N. Palmer, 1992. The
    ! sensitivity of the ECMWF model to the parameterization of
    ! evaporation from the tropical oceans. J. Climate, 5:418-434.

    use nr_util, only: assert_eq

    use clesphys, only: f_cdrag_oce, f_cdrag_ter
    use indicesol, only: is_oce
    use SUPHEC_M, only: rcpd, rd, retv, rg
    USE yoethf_m, ONLY: rvtmp2

    INTEGER, intent(in):: nsrf ! indice pour le type de surface

    REAL, intent(in):: speed(:) ! (knon)
    ! norm of the wind at the first model level

    REAL, intent(in):: t(:) ! (knon)
    ! temperature de l'air au 1er niveau du modele

    REAL, intent(in):: q(:) ! (knon) ! humidite de l'air au 1er niveau du modele

    REAL, intent(in):: zgeop(:) ! (knon)
    ! g\'eopotentiel au 1er niveau du mod\`ele

    REAL, intent(in) :: psol(:) ! (knon) pression au sol
    REAL, intent(in):: ts(:) ! (knon) temperature de l'air a la surface
    REAL, intent(in):: qsurf(:) ! (knon) humidite de l'air a la surface
    REAL, intent(in):: rugos(:) ! (knon) rugosit\'e
    REAL, intent(out):: cdragm(:) ! (knon) drag coefficient pour le moment 

    REAL, intent(out):: cdragh(:) ! (knon)
    ! drag coefficient pour les flux de chaleur latente et sensible

    REAL, intent(out), optional:: pref(:) ! (knon) pression au niveau zgeop / RG

    ! Local:

    REAL, PARAMETER:: ckap = 0.4, cb = 5., cc = 5., cd = 5., cepdu2 = 0.1**2
    real, parameter:: f_ri_cd_min = 0.1
    INTEGER i, knon
    REAL zdu2, ztsolv, ztvd, zscf, zucf
    real zcdn ! drag coefficient neutre

    REAL zri
    ! nombre de Richardson entre la surface et le niveau de reference
    ! zgeop / RG

    !-------------------------------------------------------------------------

    knon = assert_eq([size(speed), size(t), size(q), size(zgeop), size(ts), &
         size(qsurf), size(rugos), size(cdragm), size(cdragh)], "cdrag knon")

    DO i = 1, knon
       zdu2 = max(cepdu2, speed(i)**2)
       ztsolv = ts(i) * (1. + RETV * max(qsurf(i), 0.))
       ztvd = (t(i) + zgeop(i) / RCPD / (1. + RVTMP2 * q(i))) &
            * (1. + RETV * q(i))
       zri = zgeop(i) * (ztvd - ztsolv) / (zdu2 * ztvd)
       zcdn = (ckap / log(1. + zgeop(i) / (RG * rugos(i))))**2

       IF (zri < 0.) THEN
          ! situation instable
          zucf = 1. / (1. + 3. * cb * cc * zcdn &
               * SQRT(ABS(zri) * (1. + zgeop(i) / (RG * rugos(i)))))
          cdragm(i) = zcdn * max((1. - 2. * cb * zri * zucf), f_ri_cd_min)

          IF (nsrf == is_oce) then
             ! Cf. Miller et al. (1992).
             cdragh(i) = f_cdrag_oce * zcdn * (1. + ((0.0016 &
                  / (zcdn * SQRT(zdu2))) * ABS(ztvd - ztsolv)**(1. &
                  / 3.))**1.25)**(1. / 1.25)
          else
             cdragh(i) = f_cdrag_ter * zcdn &
                  * max((1. - 3. * cb * zri * zucf), f_ri_cd_min)
          end IF
       ELSE
          ! Situation stable. Pour \'eviter les incoh\'erences dans
          ! les cas tr\`es stables, on limite zri \`a 20. Cf Hess et
          ! al. (1995).
          zri = min(20., zri)
          zscf = SQRT(1. + cd * ABS(zri))
          cdragm(i) = zcdn * max(1. / (1. + 2. * CB * zri / zscf), f_ri_cd_min)
          cdragh(i) = merge(f_cdrag_oce, f_cdrag_ter, nsrf == is_oce) * zcdn &
               * max(1. / (1. + 3. * CB * zri * zscf), f_ri_cd_min)
       ENDIF
    END DO

    if (present(pref)) &
         pref = exp(log(psol) - zgeop / (RD * t * (1. + RETV * max(q, 0.))))

  END SUBROUTINE cdrag

end module cdrag_m
