module coefcdrag_m

  IMPLICIT NONE

contains

  SUBROUTINE coefcdrag (nsrf, speed, t, q, zgeop, psol, ts, qsurf, rugos, &
       cdram, cdrah, cdran, zri1, pref)

    ! From LMDZ4/libf/phylmd/coefcdrag.F90, version 1.1.1.1, 2004/05/19 12:53:07

    ! Objet : calcul des cdrags pour le moment (cdram) et les flux de
    ! chaleur sensible et latente (cdrah), du drag coefficient neutre
    ! (cdran), du nombre de Richardson entre la surface et le niveau
    ! de reference (zri1) et de la pression au niveau de reference
    ! (pref).

    ! I. Musat, 01.07.2002

    use indicesol, only: is_oce
    use SUPHEC_M, only: rd, retv, rg, rkappa
    use dimphy, only: klon

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

    REAL, dimension(klon), intent(out) :: cdram, cdrah, cdran, zri1, pref
    ! cdram--output-R- drag coefficient pour le moment
    ! cdrah--output-R- drag coefficient pour les flux de chaleur latente et sensible
    ! cdran--output-R- drag coefficient neutre
    ! zri1---output-R- nb. Richardson entre la surface et la couche zgeop/RG
    ! pref---output-R- pression au niveau zgeop/RG

    ! Local:
    REAL, parameter :: RKAR=0.40, CB=5.0, CC=5.0, CD=5.0
    INTEGER :: i
    REAL, dimension(klon) :: zdu2, zdphi, ztsolv, ztvd
    REAL, dimension(klon) :: zscf, friv, frih, zucf, zcr
    REAL, dimension(klon) :: zcfm1, zcfh1
    REAL, dimension(klon) :: zcfm2, zcfh2
    REAL, dimension(klon) :: trm0, trm1

    !-------------------------------------------------------------------------

    DO i = 1, size(speed)
       zdphi(i) = zgeop(i)
       zdu2(i) = speed(i)**2
       pref(i) = exp(log(psol(i)) - zdphi(i)/(RD*t(i)* &
            (1.+ RETV * max(q(i), 0.0))))
       ztsolv(i) = ts(i)
       ztvd(i) = t(i) * (psol(i)/pref(i))**RKAPPA
       trm0(i) = 1. + RETV * max(qsurf(i), 0.0)
       trm1(i) = 1. + RETV * max(q(i), 0.0)
       ztsolv(i) = ztsolv(i) * trm0(i)
       ztvd(i) = ztvd(i) * trm1(i)
       zri1(i) = zdphi(i)*(ztvd(i)-ztsolv(i))/(zdu2(i)*ztvd(i))
       cdran(i) = (RKAR/log(1.+zdphi(i)/(RG*rugos(i))))**2

       IF (zri1(i) >= 0.) THEN
          ! situation stable : pour eviter les inconsistances dans les cas
          ! tres stables on limite zri1 a 20. cf Hess et al. (1995)
          zri1(i) = min(20., zri1(i))
          zscf(i) = SQRT(1.+CD*ABS(zri1(i)))
          friv(i) = max(1. / (1.+2.*CB*zri1(i)/ zscf(i)), 0.1)
          zcfm1(i) = cdran(i) * friv(i)
          frih(i) = max(1./ (1.+3.*CB*zri1(i)*zscf(i)), 0.1)
          zcfh1(i) = cdran(i) * frih(i)
          cdram(i) = zcfm1(i)
          cdrah(i) = zcfh1(i)
       ELSE
          ! situation instable
          zucf(i) = 1./(1.+3.0*CB*CC*cdran(i)*SQRT(ABS(zri1(i)) &
               *(1.0+zdphi(i)/(RG*rugos(i)))))
          zcfm2(i) = cdran(i)*max((1.-2.0*CB*zri1(i)*zucf(i)), 0.1)
          zcfh2(i) = cdran(i)*max((1.-3.0*CB*zri1(i)*zucf(i)), 0.1)
          cdram(i) = zcfm2(i)
          cdrah(i) = zcfh2(i)

          ! cdrah sur l'ocean cf. Miller et al. (1992)

          zcr(i) = (0.0016/(cdran(i)*SQRT(zdu2(i))))*ABS(ztvd(i)-ztsolv(i)) &
               **(1./3.)
          IF (nsrf == is_oce) cdrah(i) = cdran(i)*(1.0+zcr(i)**1.25) &
               **(1./1.25)
       ENDIF
    END DO

  END SUBROUTINE coefcdrag

end module coefcdrag_m
