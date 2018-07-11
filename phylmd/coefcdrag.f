module coefcdrag_m

  IMPLICIT NONE

contains

  SUBROUTINE coefcdrag (nsrf, speed, t, q, zgeop, psol, ts, qsurf, rugos, &
       pcfm, pcfh, pref)

    ! From LMDZ4/libf/phylmd/coefcdrag.F90, version 1.1.1.1, 2004/05/19 12:53:07

    ! Objet : calcul des drag coefficients au sol pour le moment et
    ! les flux de chaleur sensible et latente et calcul de la pression
    ! au niveau de reference.

    ! I. Musat, 01 Jul 2002

    use clesphys, only: f_cdrag_oce, f_cdrag_ter
    use indicesol, only: is_oce
    use nr_util, only: assert_eq
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
    REAL, intent(out):: pcfm(:) ! (knon) drag coefficient pour le moment 

    REAL, intent(out):: pcfh(:) ! (knon)
    ! drag coefficient pour les flux de chaleur latente et sensible

    REAL, intent(out), optional:: pref(:) ! (knon) pression au niveau zgeop/RG

    ! Local:
    REAL, PARAMETER:: ckap=0.40, cb=5.0, cc=5.0, cd=5.0, cepdu2=0.1**2
    INTEGER i, knon
    REAL zdu2, ztsolv, ztvd, zscf, zucf, zcr, friv, frih
    REAL zcfm1, zcfh1, zcfm2, zcfh2
    real zcdn ! drag coefficient neutre

    REAL zri
    ! nb. Richardson entre la surface et la couche zgeop/RG
    ! nombre de Richardson entre la surface et le niveau de reference (zri)

    !-------------------------------------------------------------------------

    knon = assert_eq([size(speed), size(t), size(q), size(zgeop), size(ts), &
         size(qsurf), size(rugos), size(pcfm), size(pcfh), size(pcfm)], &
         "clcdrag knon")
    
    DO i = 1, knon
       zdu2 = max(cepdu2, speed(i)**2)
       ztsolv = ts(i) * (1. + RETV * max(qsurf(i), 0.))
       ztvd = (t(i)+zgeop(i)/RCPD/(1.+RVTMP2*q(i))) *(1.+RETV*q(i))
       zri = zgeop(i)*(ztvd-ztsolv)/(zdu2*ztvd)
       zcdn = (ckap/log(1.+zgeop(i)/(RG*rugos(i))))**2

       IF (zri > 0.) THEN
          ! Situation stable. Pour eviter les inconsistances dans les cas
          ! tres stables on limite zri a 20. cf Hess et al. (1995).
          zri = min(20., zri)
          zscf = SQRT(1.+cd*ABS(zri))
          friv = max(1. / (1.+2.*CB*zri/ zscf), 0.1)
          zcfm1 = zcdn * friv
          frih = max(1./ (1.+3.*CB*zri*zscf), 0.1)
          zcfh1 = f_cdrag_ter * zcdn * frih
          IF (nsrf == is_oce) zcfh1 = f_cdrag_oce * zcdn * frih
          pcfm(i) = zcfm1
          pcfh(i) = zcfh1
       ELSE
          ! situation instable
          zucf = 1./(1.+3.0*cb*cc*zcdn*SQRT(ABS(zri) &
               *(1.0+zgeop(i)/(RG*rugos(i)))))
          zcfm2 = zcdn*max((1.-2.0*cb*zri*zucf), 0.1)
          zcfh2 = f_cdrag_ter * zcdn*max((1.-3.0*cb*zri*zucf), 0.1)
          pcfm(i) = zcfm2
          pcfh(i) = zcfh2

          ! pcfh sur l'ocean cf. Miller et al. (1992)
          zcr = (0.0016/(zcdn*SQRT(zdu2)))*ABS(ztvd-ztsolv)**(1./3.)
          IF (nsrf == is_oce) pcfh(i) = f_cdrag_oce * zcdn &
               * (1. + zcr**1.25)**(1. / 1.25)
       ENDIF
    END DO

    if (present(pref)) &
         pref = exp(log(psol) - zgeop / (RD * t * (1. + RETV * max(q, 0.))))

  END SUBROUTINE coefcdrag

end module coefcdrag_m
