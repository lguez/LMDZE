module clcdrag_m

  IMPLICIT NONE

contains

  SUBROUTINE clcdrag(klon, knon, nsrf, zxli, u, v, t, q, zgeop, ts, qsurf, &
       rugos, pcfm, pcfh)

    ! From LMDZ4/libf/phylmd/clcdrag.F90, version 1.1.1.1 2004/05/19 12:53:07

    USE indicesol, ONLY : is_oce
    USE suphec_m, ONLY : rcpd, retv, rg
    USE yoethf_m, ONLY : rvtmp2

    ! Objet : calcul des cdrags pour le moment (pcfm) et les flux de
    ! chaleur sensible et latente (pcfh).

    ! knon----input-I- nombre de points pour un type de surface
    ! nsrf----input-I- indice pour le type de surface; voir indicesol.inc
    ! zxli----input-L- calcul des cdrags selon Laurent Li
    ! u-------input-R- vent zonal au 1er niveau du modele
    ! v-------input-R- vent meridien au 1er niveau du modele
    ! t-------input-R- temperature de l'air au 1er niveau du modele
    ! q-------input-R- humidite de l'air au 1er niveau du modele
    ! ts------input-R- temperature de l'air a la surface
    ! qsurf---input-R- humidite de l'air a la surface
    ! rugos---input-R- rugosite

    ! pcfm---output-R- cdrag pour le moment 
    ! pcfh---output-R- cdrag pour les flux de chaleur latente et sensible

    INTEGER, intent(in) :: klon
    ! dimension de la grille physique (= nb_pts_latitude X nb_pts_longitude)

    INTEGER, intent(in) :: knon, nsrf

    ! Fonctions thermodynamiques et fonctions d'instabilite
    LOGICAL, intent(in) :: zxli ! utiliser un jeu de fonctions simples

    REAL, intent(in), dimension(klon) :: u, v, t, q
    REAL, intent(in):: zgeop(klon) ! géopotentiel au 1er niveau du modèle
    REAL, intent(in), dimension(klon) :: ts, qsurf
    REAL, intent(in), dimension(klon) :: rugos
    REAL, intent(out):: pcfm(:), pcfh(:) ! (knon)

    ! Quelques constantes et options:
    REAL, PARAMETER :: ckap=0.40, cb=5.0, cc=5.0, cd=5.0, cepdu2=(0.1)**2

    ! Variables locales :
    INTEGER :: i
    REAL :: zdu2, ztsolv, ztvd, zscf
    REAL :: zucf, zcr
    REAL :: friv, frih
    REAL, dimension(klon) :: zcfm1, zcfm2
    REAL, dimension(klon) :: zcfh1, zcfh2
    REAL, dimension(klon) :: zcdn
    REAL, dimension(klon) :: zri

    !--------------------------------------------------------------------

    ! Calculer le frottement au sol (Cdrag)

    DO i = 1, knon
       zdu2 = max(cepdu2,u(i)**2+v(i)**2)
       ztsolv = ts(i) * (1.0+RETV*qsurf(i))
       ztvd = (t(i)+zgeop(i)/RCPD/(1.+RVTMP2*q(i))) &
            *(1.+RETV*q(i))
       zri(i) = zgeop(i)*(ztvd-ztsolv)/(zdu2*ztvd)
       zcdn(i) = (ckap/log(1.+zgeop(i)/(RG*rugos(i))))**2

       IF (zri(i) .gt. 0.) THEN      
          ! situation stable
          zri(i) = min(20.,zri(i))
          IF (.NOT. zxli) THEN
             zscf = SQRT(1.+cd*ABS(zri(i)))
             FRIV = AMAX1(1. / (1.+2.*CB*zri(i)/ZSCF), 0.1)
             zcfm1(i) = zcdn(i) * FRIV
             FRIH = AMAX1(1./ (1.+3.*CB*zri(i)*ZSCF), 0.1 )
             zcfh1(i) = 0.8 * zcdn(i) * FRIH
             pcfm(i) = zcfm1(i)
             pcfh(i) = zcfh1(i)
          ELSE
             pcfm(i) = zcdn(i)* fsta(zri(i))
             pcfh(i) = zcdn(i)* fsta(zri(i))
          ENDIF
       ELSE                          
          ! situation instable
          IF (.NOT. zxli) THEN
             zucf = 1./(1.+3.0*cb*cc*zcdn(i)*SQRT(ABS(zri(i)) &
                  *(1.0+zgeop(i)/(RG*rugos(i)))))
             zcfm2(i) = zcdn(i)*amax1((1.-2.0*cb*zri(i)*zucf),0.1)
             zcfh2(i) = 0.8 * zcdn(i)*amax1((1.-3.0*cb*zri(i)*zucf),0.1)
             pcfm(i) = zcfm2(i)
             pcfh(i) = zcfh2(i)
          ELSE
             pcfm(i) = zcdn(i)* fins(zri(i))
             pcfh(i) = zcdn(i)* fins(zri(i))
          ENDIF
          zcr = (0.0016/(zcdn(i)*SQRT(zdu2)))*ABS(ztvd-ztsolv)**(1./3.)
          IF(nsrf == is_oce) pcfh(i) = 0.8 * zcdn(i) &
               * (1. + zcr**1.25)**(1. / 1.25)
       ENDIF
    END DO

  contains

    ! Fonctions thermodynamiques et fonctions d'instabilite

    function fsta(x)
      REAL fsta
      real, intent(in):: x
      fsta = 1.0 / (1.0+10.0*x*(1+8.0*x))
    end function fsta

    !*******************************************************

    function fins(x)
      REAL fins
      real, intent(in):: x
      fins = SQRT(1.0-18.0*x)
    end function fins

  END SUBROUTINE clcdrag

end module clcdrag_m
