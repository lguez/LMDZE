module clcdrag_m

  IMPLICIT NONE

contains

  SUBROUTINE clcdrag(nsrf, speed, t, q, zgeop, psol, ts, qsurf, rugos, pcfm, &
       pcfh, pref)

    ! From LMDZ4/libf/phylmd/clcdrag.F90, version 1.1.1.1, 2004/05/19 12:53:07

    ! Objet : calcul des cdrags pour le moment (pcfm) et les flux de
    ! chaleur sensible et latente (pcfh).
    ! Calculer le frottement au sol (drag coefficient)

    USE indicesol, ONLY: is_oce
    use nr_util, only: assert_eq
    USE suphec_m, ONLY: rcpd, retv, rg
    USE yoethf_m, ONLY: rvtmp2

    INTEGER, intent(in):: nsrf ! indice pour le type de surface

    REAL, intent(in):: speed(:) ! (knon)
    ! norm of the wind at the first model layer

    REAL, intent(in):: t(:) ! (knon)
    ! temperature de l'air au 1er niveau du modele

    REAL, intent(in):: q(:) ! (knon) ! humidite de l'air au 1er niveau du modele
    REAL, intent(in):: zgeop(:) ! (knon) géopotentiel au 1er niveau du modèle
    REAL, intent(in) :: psol(:) ! (knon) pression au sol
    REAL, intent(in):: ts(:) ! (knon) temperature de l'air a la surface
    REAL, intent(in):: qsurf(:) ! (knon) humidite de l'air a la surface
    REAL, intent(in):: rugos(:) ! (knon) rugosit\'e
    REAL, intent(out):: pcfm(:) ! (knon) drag coefficient pour le moment 

    REAL, intent(out):: pcfh(:) ! (knon)
    ! drag coefficient pour les flux de chaleur latente et sensible

    REAL, intent(out), optional:: pref(:) ! (knon) pression au niveau zgeop/RG

    ! Local:

    ! Quelques constantes et options:
    REAL, PARAMETER:: ckap=0.40, cb=5.0, cc=5.0, cd=5.0, cepdu2=0.1**2

    INTEGER:: i, knon
    REAL:: zdu2, ztsolv, ztvd, zscf
    REAL:: zucf, zcr
    REAL:: friv, frih
    REAL, dimension(size(speed)):: zcfm1, zcfm2
    REAL, dimension(size(speed)):: zcfh1, zcfh2
    REAL, dimension(size(speed)):: zcdn
    REAL, dimension(size(speed)):: zri

    !--------------------------------------------------------------------

    knon = assert_eq([size(speed), size(t), size(q), size(zgeop), size(ts), &
         size(qsurf), size(rugos), size(pcfm), size(pcfh), size(pcfm)], &
         "clcdrag knon")
    
    DO i = 1, knon
       zdu2 = max(cepdu2,speed(i)**2)
       ztsolv = ts(i) * (1.0+RETV*qsurf(i))
       ztvd = (t(i)+zgeop(i)/RCPD/(1.+RVTMP2*q(i))) &
            *(1.+RETV*q(i))
       zri(i) = zgeop(i)*(ztvd-ztsolv)/(zdu2*ztvd)
       zcdn(i) = (ckap/log(1.+zgeop(i)/(RG*rugos(i))))**2

       IF (zri(i) .gt. 0.) THEN      
          ! situation stable
          zri(i) = min(20.,zri(i))
          zscf = SQRT(1.+cd*ABS(zri(i)))
          FRIV = AMAX1(1. / (1.+2.*CB*zri(i)/ZSCF), 0.1)
          zcfm1(i) = zcdn(i) * FRIV
          FRIH = AMAX1(1./ (1.+3.*CB*zri(i)*ZSCF), 0.1 )
          zcfh1(i) = 0.8 * zcdn(i) * FRIH
          pcfm(i) = zcfm1(i)
          pcfh(i) = zcfh1(i)
       ELSE                          
          ! situation instable
          zucf = 1./(1.+3.0*cb*cc*zcdn(i)*SQRT(ABS(zri(i)) &
               *(1.0+zgeop(i)/(RG*rugos(i)))))
          zcfm2(i) = zcdn(i)*amax1((1.-2.0*cb*zri(i)*zucf),0.1)
          zcfh2(i) = 0.8 * zcdn(i)*amax1((1.-3.0*cb*zri(i)*zucf),0.1)
          pcfm(i) = zcfm2(i)
          pcfh(i) = zcfh2(i)
          zcr = (0.0016/(zcdn(i)*SQRT(zdu2)))*ABS(ztvd-ztsolv)**(1./3.)
          IF(nsrf == is_oce) pcfh(i) = 0.8 * zcdn(i) &
               * (1. + zcr**1.25)**(1. / 1.25)
       ENDIF
    END DO

  END SUBROUTINE clcdrag

end module clcdrag_m
