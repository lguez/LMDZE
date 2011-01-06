!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/coefcdrag.F90,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
!
!
!
      SUBROUTINE coefcdrag (klon, knon, nsrf, zxli, &
                            speed, t, q, zgeop, psol, &
                            ts, qsurf, rugos, okri, ri1, &
                            cdram, cdrah, cdran, zri1, pref)
      use indicesol
      use SUPHEC_M
      use yoethf_m
      IMPLICIT none
!-------------------------------------------------------------------------
! Objet : calcul des cdrags pour le moment (cdram) et les flux de chaleur 
!         sensible et latente (cdrah), du cdrag neutre (cdran), 
!         du nombre de Richardson entre la surface et le niveau de reference 
!         (zri1) et de la pression au niveau de reference (pref).    
!
! I. Musat, 01.07.2002
!-------------------------------------------------------------------------
!
! klon----input-I- dimension de la grille physique (= nb_pts_latitude X nb_pts_longitude)
! knon----input-I- nombre de points pour un type de surface
! nsrf----input-I- indice pour le type de surface; voir indicesol.inc
! zxli----input-L- TRUE si calcul des cdrags selon Laurent Li
! speed---input-R- module du vent au 1er niveau du modele
! t-------input-R- temperature de l'air au 1er niveau du modele
! q-------input-R- humidite de l'air au 1er niveau du modele
! zgeop---input-R- geopotentiel au 1er niveau du modele
! psol----input-R- pression au sol 
! ts------input-R- temperature de l'air a la surface
! qsurf---input-R- humidite de l'air a la surface
! rugos---input-R- rugosite
! okri----input-L- TRUE si on veut tester le nb. Richardson entre la sfce 
!                  et zref par rapport au Ri entre la sfce et la 1ere couche
! ri1-----input-R- nb. Richardson entre la surface et la 1ere couche
!
! cdram--output-R- cdrag pour le moment
! cdrah--output-R- cdrag pour les flux de chaleur latente et sensible
! cdran--output-R- cdrag neutre
! zri1---output-R- nb. Richardson entre la surface et la couche zgeop/RG
! pref---output-R- pression au niveau zgeop/RG
!
      INTEGER, intent(in) :: klon, knon, nsrf
      LOGICAL, intent(in) :: zxli 
      REAL, dimension(klon), intent(in) :: speed, t, q, zgeop, psol
      REAL, dimension(klon), intent(in) :: ts, qsurf, rugos, ri1 
      LOGICAL, intent(in) :: okri    
!
      REAL, dimension(klon), intent(out) :: cdram, cdrah, cdran, zri1, pref
!-------------------------------------------------------------------------
!
! Quelques constantes :
      REAL, parameter :: RKAR=0.40, CB=5.0, CC=5.0, CD=5.0
!
! Variables locales :
      INTEGER :: i
      REAL, dimension(klon) :: zdu2, zdphi, ztsolv, ztvd
      REAL, dimension(klon) :: zscf, friv, frih, zucf, zcr
      REAL, dimension(klon) :: zcfm1, zcfh1
      REAL, dimension(klon) :: zcfm2, zcfh2
      REAL, dimension(klon) :: trm0, trm1
!-------------------------------------------------------------------------
      REAL :: fsta, fins, x
      fsta(x) = 1.0 / (1.0+10.0*x*(1+8.0*x))
      fins(x) = SQRT(1.0-18.0*x)
!-------------------------------------------------------------------------
!
      DO i = 1, knon
!
       zdphi(i) = zgeop(i)
       zdu2(i) = speed(i)**2
       pref(i) = exp(log(psol(i)) - zdphi(i)/(RD*t(i)* &
                 (1.+ RETV * max(q(i),0.0))))
       ztsolv(i) = ts(i)
       ztvd(i) = t(i) * (psol(i)/pref(i))**RKAPPA
       trm0(i) = 1. + RETV * max(qsurf(i),0.0)
       trm1(i) = 1. + RETV * max(q(i),0.0)
       ztsolv(i) = ztsolv(i) * trm0(i)
       ztvd(i) = ztvd(i) * trm1(i)
       zri1(i) = zdphi(i)*(ztvd(i)-ztsolv(i))/(zdu2(i)*ztvd(i))
!
! on teste zri1 par rapport au Richardson de la 1ere couche ri1 
!
!IM +++
       IF(1.EQ.0) THEN
       IF (okri) THEN
         IF (ri1(i).GE.0.0.AND.zri1(i).LT.0.0) THEN
           zri1(i) = ri1(i)
         ELSE IF(ri1(i).LT.0.0.AND.zri1(i).GE.0.0) THEN
           zri1(i) = ri1(i)
         ENDIF 
       ENDIF
       ENDIF
!IM ---
! 
       cdran(i) = (RKAR/log(1.+zdphi(i)/(RG*rugos(i))))**2

       IF (zri1(i) .ge. 0.) THEN 
!
! situation stable : pour eviter les inconsistances dans les cas 
! tres stables on limite zri1 a 20. cf Hess et al. (1995)
!
         zri1(i) = min(20.,zri1(i))
!
         IF (.NOT.zxli) THEN
           zscf(i) = SQRT(1.+CD*ABS(zri1(i)))
           friv(i) = max(1. / (1.+2.*CB*zri1(i)/ zscf(i)), 0.1)
           zcfm1(i) = cdran(i) * friv(i)
           frih(i) = max(1./ (1.+3.*CB*zri1(i)*zscf(i)), 0.1 )
           zcfh1(i) = cdran(i) * frih(i)
           cdram(i) = zcfm1(i)
           cdrah(i) = zcfh1(i)
         ELSE
           cdram(i) = cdran(i)* fsta(zri1(i))
           cdrah(i) = cdran(i)* fsta(zri1(i))
         ENDIF
!
       ELSE
! 
! situation instable
!
         IF (.NOT.zxli) THEN
           zucf(i) = 1./(1.+3.0*CB*CC*cdran(i)*SQRT(ABS(zri1(i)) &
                 *(1.0+zdphi(i)/(RG*rugos(i)))))
           zcfm2(i) = cdran(i)*max((1.-2.0*CB*zri1(i)*zucf(i)),0.1)
           zcfh2(i) = cdran(i)*max((1.-3.0*CB*zri1(i)*zucf(i)),0.1)
           cdram(i) = zcfm2(i)
           cdrah(i) = zcfh2(i)
         ELSE
           cdram(i) = cdran(i)* fins(zri1(i))
           cdrah(i) = cdran(i)* fins(zri1(i))
         ENDIF
!
! cdrah sur l'ocean cf. Miller et al. (1992)
!
         zcr(i) = (0.0016/(cdran(i)*SQRT(zdu2(i))))*ABS(ztvd(i)-ztsolv(i)) &
               **(1./3.)
         IF (nsrf.EQ.is_oce) cdrah(i) = cdran(i)*(1.0+zcr(i)**1.25) &
                  **(1./1.25)
       ENDIF
!
      END DO
      RETURN 
      END SUBROUTINE coefcdrag
