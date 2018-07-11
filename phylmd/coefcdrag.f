module coefcdrag_m

  IMPLICIT NONE

contains

  SUBROUTINE coefcdrag (nsrf, speed, t, q, zgeop, psol, ts, qsurf, rugos, &
       pcfm, pcfh, pref)

    ! From LMDZ4/libf/phylmd/coefcdrag.F90, version 1.1.1.1, 2004/05/19 12:53:07

    ! Objet : calcul des cdrags pour le moment (pcfm) et les flux de
    ! chaleur sensible et latente (pcfh) et de la pression au niveau
    ! de reference (pref).

    ! I. Musat, 01.07.2002

    use clesphys, only: f_cdrag_oce, f_cdrag_ter
    use indicesol, only: is_oce
    use nr_util, only: assert_eq
    use SUPHEC_M, only: rcpd, rd, retv, rg
    use dimphy, only: klon
    USE yoethf_m, ONLY: rvtmp2

    INTEGER, intent(in) :: nsrf
    ! nsrf----input-I- indice pour le type de surface; voir indicesol.inc
    REAL, intent(in) :: speed(:), t(:), q(:), zgeop(:), psol(:) ! (knon)
    ! speed---input-R- module du vent au 1er niveau du modele
    ! t-------input-R- temperature de l'air au 1er niveau du modele
    ! q-------input-R- humidite de l'air au 1er niveau du modele
    ! zgeop---input-R- geopotentiel au 1er niveau du modele
    ! psol----input-R- pression au sol
    REAL, dimension(klon), intent(in) :: ts, qsurf, rugos
    ! ts------input-R- temperature de l'air a la surface
    ! qsurf---input-R- humidite de l'air a la surface
    ! rugos---input-R- rugosite

    REAL, dimension(klon), intent(out) :: pcfm, pcfh
    ! drag coefficients pour le moment et pour les flux de chaleur
    ! latente et sensible
    
    REAL, intent(out), optional:: pref(:) ! (knon) pression au niveau zgeop/RG

    ! Local:
    REAL, parameter :: CKAP=0.40, CB=5.0, CC=5.0, CD=5.0
    INTEGER :: i
    REAL, dimension(klon) :: zdu2, ztsolv, ztvd
    REAL, dimension(klon) :: zscf, friv, frih, zucf, zcr
    REAL, dimension(klon) :: zcfm1, zcfh1
    REAL, dimension(klon) :: zcfm2, zcfh2
    REAL, dimension(klon) :: trm0, trm1
    real cdran(klon) ! drag coefficient neutre
    REAL pref_local(size(speed)) ! (knon) pression au niveau zgeop/RG

    REAL zri1(klon)
    ! nb. Richardson entre la surface et la couche zgeop/RG
    ! nombre de Richardson entre la surface et le niveau de reference (zri1)

    !-------------------------------------------------------------------------

    DO i = 1, size(speed)
       zdu2(i) = speed(i)**2
       pref_local(i) = exp(log(psol(i)) - zgeop(i)/(RD*t(i)* &
            (1.+ RETV * max(q(i), 0.0))))
       ztsolv(i) = ts(i)
       ztvd(i) = (t(i)+zgeop(i)/RCPD/(1.+RVTMP2*q(i))) *(1.+RETV*q(i)) 
       trm0(i) = 1. + RETV * max(qsurf(i), 0.0)
       trm1(i) = 1. + RETV * max(q(i), 0.0)
       ztsolv(i) = ztsolv(i) * trm0(i)
       zri1(i) = zgeop(i)*(ztvd(i)-ztsolv(i))/(zdu2(i)*ztvd(i))
       cdran(i) = (CKAP/log(1.+zgeop(i)/(RG*rugos(i))))**2

       IF (zri1(i) >= 0.) THEN
          ! situation stable : pour eviter les inconsistances dans les cas
          ! tres stables on limite zri1 a 20. cf Hess et al. (1995)
          zri1(i) = min(20., zri1(i))
          zscf(i) = SQRT(1.+CD*ABS(zri1(i)))
          friv(i) = max(1. / (1.+2.*CB*zri1(i)/ zscf(i)), 0.1)
          zcfm1(i) = cdran(i) * friv(i)
          frih(i) = max(1./ (1.+3.*CB*zri1(i)*zscf(i)), 0.1)
          zcfh1(i) = f_cdrag_ter * cdran(i) * frih(i) 
          IF (nsrf == is_oce) zcfh1(i) = f_cdrag_oce * cdran(i) * frih(i)
          pcfm(i) = zcfm1(i)
          pcfh(i) = zcfh1(i)
       ELSE
          ! situation instable
          zucf(i) = 1./(1.+3.0*CB*CC*cdran(i)*SQRT(ABS(zri1(i)) &
               *(1.0+zgeop(i)/(RG*rugos(i)))))
          zcfm2(i) = cdran(i)*max((1.-2.0*CB*zri1(i)*zucf(i)), 0.1)
          zcfh2(i) = f_cdrag_ter * cdran(i)*max((1.-3.0*CB*zri1(i)*zucf(i)), 0.1)
          pcfm(i) = zcfm2(i)
          pcfh(i) = zcfh2(i)

          ! pcfh sur l'ocean cf. Miller et al. (1992)

          zcr(i) = (0.0016/(cdran(i)*SQRT(zdu2(i))))*ABS(ztvd(i)-ztsolv(i)) &
               **(1./3.)
          IF (nsrf == is_oce) pcfh(i) = f_cdrag_oce * cdran(i)*(1.0+zcr(i)**1.25) **(1./1.25)
       ENDIF
    END DO

    if (present(pref)) pref = pref_local

  END SUBROUTINE coefcdrag

end module coefcdrag_m
