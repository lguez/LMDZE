SUBROUTINE coefkz(nsrf, knon, paprs, pplay, ksta, ksta_ter, ts, rugos, u, v, &
     t, q, qsurf, pcfm, pcfh)

  ! Auteur(s) F. Hourdin, M. Forichon, Z.X. Li (LMD/CNRS) date: 1993/09/22
  !           (une version strictement identique a l'ancien modele)
  ! Objet: calculer le coefficient du frottement du sol (Cdrag) et les
  !        coefficients d'echange turbulent dans l'atmosphere.

  USE indicesol, ONLY : is_oce
  USE dimphy, ONLY : klev, klon, max
  USE iniprint, ONLY : prt_level
  USE suphec_m, ONLY : rcpd, rd, retv, rg, rkappa, rlstt, rlvtt, rtt
  USE yoethf_m, ONLY : r2es, r5ies, r5les, rvtmp2
  USE fcttre, ONLY : dqsatl, dqsats, foede, foeew, qsatl, qsats, thermcep
  USE conf_phys_m, ONLY : iflag_pbl

  IMPLICIT none

  ! Arguments:
  ! nsrf-----input-I- indicateur de la nature du sol
  ! knon-----input-I- nombre de points a traiter
  ! paprs----input-R- pression a chaque intercouche (en Pa)
  ! pplay----input-R- pression au milieu de chaque couche (en Pa)
  ! ts-------input-R- temperature du sol (en Kelvin)
  ! rugos----input-R- longeur de rugosite (en m)
  ! u--------input-R- vitesse u
  ! v--------input-R- vitesse v
  ! q--------input-R- vapeur d'eau (kg/kg)

  ! itop-----output-I- numero de couche du sommet de la couche limite
  ! pcfm-----output-R- coefficients a calculer (vitesse)
  ! pcfh-----output-R- coefficients a calculer (chaleur et humidite)

  INTEGER knon, nsrf
  REAL ts(klon)
  REAL paprs(klon, klev+1), pplay(klon, klev)
  REAL u(klon, klev), v(klon, klev), q(klon, klev)
  REAL, intent(in):: t(klon, klev) ! temperature (K)
  REAL rugos(klon)

  REAL pcfm(klon, klev), pcfh(klon, klev)
  INTEGER itop(klon)

  ! Quelques constantes et options:

  REAL cepdu2, ckap, cb, cc, cd, clam
  PARAMETER (cepdu2 =(0.1)**2)
  PARAMETER (CKAP=0.4)
  PARAMETER (cb=5.0)
  PARAMETER (cc=5.0)
  PARAMETER (cd=5.0)
  PARAMETER (clam=160.0)
  REAL ratqs ! largeur de distribution de vapeur d'eau
  PARAMETER (ratqs=0.05)
  LOGICAL richum ! utilise le nombre de Richardson humide
  PARAMETER (richum=.TRUE.)
  REAL ric ! nombre de Richardson critique
  PARAMETER(ric=0.4)
  REAL prandtl
  PARAMETER (prandtl=0.4)
  REAL kstable ! diffusion minimale (situation stable)
  ! GKtest
  ! PARAMETER (kstable=1.0e-10)
  REAL ksta, ksta_ter
  !IM: 261103     REAL kstable_ter, kstable_sinon
  !IM: 211003 cf GK   PARAMETER (kstable_ter = 1.0e-6)
  !IM: 261103     PARAMETER (kstable_ter = 1.0e-8)
  !IM: 261103   PARAMETER (kstable_ter = 1.0e-10)
  !IM: 261103   PARAMETER (kstable_sinon = 1.0e-10)
  ! fin GKtest
  REAL mixlen ! constante controlant longueur de melange
  PARAMETER (mixlen=35.0)
  INTEGER isommet ! le sommet de la couche limite
  PARAMETER (isommet=klev)
  LOGICAL tvirtu ! calculer Ri d'une maniere plus performante
  PARAMETER (tvirtu=.TRUE.)
  LOGICAL opt_ec ! formule du Centre Europeen dans l'atmosphere
  PARAMETER (opt_ec=.FALSE.)

  ! Variables locales:

  INTEGER i, k, kk !IM 120704
  REAL zgeop(klon, klev)
  REAL zmgeom(klon)
  REAL zri(klon)
  REAL zl2(klon)

  REAL u1(klon), v1(klon), t1(klon), q1(klon), z1(klon)
  REAL pcfm1(klon), pcfh1(klon)

  REAL zdphi, zdu2, ztvd, ztvu, zcdn
  REAL zscf
  REAL zt, zq, zdelta, zcvm5, zcor, zqs, zfr, zdqs
  REAL z2geomf, zalh2, zalm2, zscfh, zscfm
  REAL t_coup
  PARAMETER (t_coup=273.15)
  !IM
  LOGICAL check
  PARAMETER (check=.false.)

  ! contre-gradient pour la chaleur sensible: Kelvin/metre
  REAL gamt(2:klev)
  real qsurf(klon) 

  LOGICAL appel1er
  SAVE appel1er

  ! Fonctions thermodynamiques et fonctions d'instabilite
  REAL fsta, fins, x
  LOGICAL zxli ! utiliser un jeu de fonctions simples
  PARAMETER (zxli=.FALSE.)

  fsta(x) = 1.0 / (1.0+10.0*x*(1+8.0*x))
  fins(x) = SQRT(1.0-18.0*x)

  DATA appel1er /.TRUE./

  !--------------------------------------------------------------------

  IF (appel1er) THEN
     if (prt_level > 9) THEN
        print *, 'coefkz, opt_ec:', opt_ec
        print *, 'coefkz, richum:', richum
        IF (richum) print *, 'coefkz, ratqs:', ratqs
        print *, 'coefkz, isommet:', isommet
        print *, 'coefkz, tvirtu:', tvirtu
        appel1er = .FALSE.
     endif
  ENDIF

  ! Initialiser les sorties

  DO k = 1, klev
     DO i = 1, knon
        pcfm(i, k) = 0.0
        pcfh(i, k) = 0.0
     ENDDO
  ENDDO
  DO i = 1, knon
     itop(i) = 0
  ENDDO

  ! Prescrire la valeur de contre-gradient

  if (iflag_pbl.eq.1) then
     DO k = 3, klev
        gamt(k) = -1.0E-03
     ENDDO
     gamt(2) = -2.5E-03
  else
     DO k = 2, klev
        gamt(k) = 0.0
     ENDDO
  ENDIF

  IF ( nsrf .NE. is_oce ) THEN
     kstable = ksta_ter
  ELSE
     kstable = ksta
  ENDIF

  ! Calculer les géopotentiels de chaque couche

  DO i = 1, knon
     zgeop(i, 1) = RD * t(i, 1) / (0.5 * (paprs(i, 1) + pplay(i, 1))) &
          * (paprs(i, 1) - pplay(i, 1))
  ENDDO
  DO k = 2, klev
     DO i = 1, knon
        zgeop(i, k) = zgeop(i, k-1) &
             + RD * 0.5*(t(i, k-1)+t(i, k)) / paprs(i, k) &
             * (pplay(i, k-1)-pplay(i, k))
     ENDDO
  ENDDO

  ! Calculer le frottement au sol (Cdrag)

  DO i = 1, knon
     u1(i) = u(i, 1)
     v1(i) = v(i, 1)
     t1(i) = t(i, 1)
     q1(i) = q(i, 1)
     z1(i) = zgeop(i, 1)
  ENDDO

  CALL clcdrag(klon, knon, nsrf, zxli,  &
       u1, v1, t1, q1, z1, &
       ts, qsurf, rugos, &
       pcfm1, pcfh1) 
  !IM  $             ts, qsurf, rugos,

  DO i = 1, knon
     pcfm(i, 1)=pcfm1(i)
     pcfh(i, 1)=pcfh1(i)
  ENDDO

  ! Calculer les coefficients turbulents dans l'atmosphere

  DO i = 1, knon
     itop(i) = isommet
  ENDDO

  DO k = 2, isommet
     DO i = 1, knon
        zdu2=MAX(cepdu2, (u(i, k)-u(i, k-1))**2 &
             +(v(i, k)-v(i, k-1))**2)
        zmgeom(i)=zgeop(i, k)-zgeop(i, k-1)
        zdphi =zmgeom(i) / 2.0
        zt = (t(i, k)+t(i, k-1)) * 0.5
        zq = (q(i, k)+q(i, k-1)) * 0.5

        !           calculer Qs et dQs/dT:

        IF (thermcep) THEN
           zdelta = MAX(0., SIGN(1., RTT-zt))
           zcvm5 = R5LES*RLVTT/RCPD/(1.0+RVTMP2*zq)*(1.-zdelta)  &
                + R5IES*RLSTT/RCPD/(1.0+RVTMP2*zq)*zdelta
           zqs = R2ES * FOEEW(zt, zdelta) / pplay(i, k)
           zqs = MIN(0.5, zqs)
           zcor = 1./(1.-RETV*zqs)
           zqs = zqs*zcor
           zdqs = FOEDE(zt, zdelta, zcvm5, zqs, zcor)
        ELSE
           IF (zt .LT. t_coup) THEN
              zqs = qsats(zt) / pplay(i, k)
              zdqs = dqsats(zt, zqs)
           ELSE
              zqs = qsatl(zt) / pplay(i, k)
              zdqs = dqsatl(zt, zqs)
           ENDIF
        ENDIF

        !           calculer la fraction nuageuse (processus humide):

        zfr = (zq+ratqs*zq-zqs) / (2.0*ratqs*zq)
        zfr = MAX(0.0, MIN(1.0, zfr))
        IF (.NOT.richum) zfr = 0.0

        !           calculer le nombre de Richardson:

        IF (tvirtu) THEN
           ztvd =( t(i, k) &
                + zdphi/RCPD/(1.+RVTMP2*zq) &
                *( (1.-zfr) + zfr*(1.+RLVTT*zqs/RD/zt)/(1.+zdqs) ) &
                )*(1.+RETV*q(i, k))
           ztvu =( t(i, k-1) &
                - zdphi/RCPD/(1.+RVTMP2*zq) &
                *( (1.-zfr) + zfr*(1.+RLVTT*zqs/RD/zt)/(1.+zdqs) ) &
                )*(1.+RETV*q(i, k-1))
           zri(i) =zmgeom(i)*(ztvd-ztvu)/(zdu2*0.5*(ztvd+ztvu))
           zri(i) = zri(i) &
                + zmgeom(i)*zmgeom(i)/RG*gamt(k) &
                *(paprs(i, k)/101325.0)**RKAPPA &
                /(zdu2*0.5*(ztvd+ztvu))

        ELSE ! calcul de Ridchardson compatible LMD5

           zri(i) =(RCPD*(t(i, k)-t(i, k-1)) &
                -RD*0.5*(t(i, k)+t(i, k-1))/paprs(i, k) &
                *(pplay(i, k)-pplay(i, k-1)) &
                )*zmgeom(i)/(zdu2*0.5*RCPD*(t(i, k-1)+t(i, k)))
           zri(i) = zri(i) + &
                zmgeom(i)*zmgeom(i)*gamt(k)/RG &
                *(paprs(i, k)/101325.0)**RKAPPA &
                /(zdu2*0.5*(t(i, k-1)+t(i, k)))
        ENDIF

        !           finalement, les coefficients d'echange sont obtenus:

        zcdn=SQRT(zdu2) / zmgeom(i) * RG

        IF (opt_ec) THEN
           z2geomf=zgeop(i, k-1)+zgeop(i, k)
           zalm2=(0.5*ckap/RG*z2geomf &
                /(1.+0.5*ckap/rg/clam*z2geomf))**2
           zalh2=(0.5*ckap/rg*z2geomf &
                /(1.+0.5*ckap/RG/(clam*SQRT(1.5*cd))*z2geomf))**2
           IF (zri(i).LT.0.0) THEN  ! situation instable
              zscf = ((zgeop(i, k)/zgeop(i, k-1))**(1./3.)-1.)**3 &
                   / (zmgeom(i)/RG)**3 / (zgeop(i, k-1)/RG)
              zscf = SQRT(-zri(i)*zscf)
              zscfm = 1.0 / (1.0+3.0*cb*cc*zalm2*zscf)
              zscfh = 1.0 / (1.0+3.0*cb*cc*zalh2*zscf)
              pcfm(i, k)=zcdn*zalm2*(1.-2.0*cb*zri(i)*zscfm)
              pcfh(i, k)=zcdn*zalh2*(1.-3.0*cb*zri(i)*zscfh)
           ELSE ! situation stable
              zscf=SQRT(1.+cd*zri(i))
              pcfm(i, k)=zcdn*zalm2/(1.+2.0*cb*zri(i)/zscf)
              pcfh(i, k)=zcdn*zalh2/(1.+3.0*cb*zri(i)*zscf)
           ENDIF
        ELSE
           zl2(i)=(mixlen*MAX(0.0, (paprs(i, k)-paprs(i, itop(i)+1)) &
                /(paprs(i, 2)-paprs(i, itop(i)+1)) ))**2
           pcfm(i, k)=sqrt(max(zcdn*zcdn*(ric-zri(i))/ric, kstable))
           pcfm(i, k)= zl2(i)* pcfm(i, k)
           pcfh(i, k) = pcfm(i, k) /prandtl ! h et m different
        ENDIF
     ENDDO
  ENDDO

  ! Au-dela du sommet, pas de diffusion turbulente:

  DO i = 1, knon
     IF (itop(i)+1 .LE. klev) THEN
        DO k = itop(i)+1, klev
           pcfh(i, k) = 0.0
           pcfm(i, k) = 0.0
        ENDDO
     ENDIF
  ENDDO

END SUBROUTINE coefkz
