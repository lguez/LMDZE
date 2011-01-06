module inigeom_m

  IMPLICIT NONE

contains

  SUBROUTINE inigeom

    ! Auteur : P. Le Van
    ! Version du 01/04/2001

    ! Calcul des �longations cuij1, ..., cuij4, cvij1, ..., cvij4 aux m�mes
    ! endroits que les aires aireij1_2d, ..., aireij4_2d.

    ! Choix entre une fonction "f(y)" � d�riv�e sinuso�dale ou � d�riv�e
    ! tangente hyperbolique
    ! calcul des coefficients (cu_2d, cv_2d, 1./cu_2d**2, 1./cv_2d**2)

    ! les coef. (cu_2d, cv_2d) permettent de passer des vitesses naturelles
    ! aux vitesses covariantes et contravariantes, ou vice-versa

    ! on a :
    ! u (covariant) = cu_2d * u (naturel), u(contrav)= u(nat)/cu_2d
    ! v (covariant) = cv_2d * v (naturel), v(contrav)= v(nat)/cv_2d

    ! on en tire : 
    ! u(covariant) = cu_2d * cu_2d * u(contravariant)
    ! v(covariant) = cv_2d * cv_2d * v(contravariant)

    ! on a l'application (x(X), y(Y)) avec - im/2 +1 <= X <= im/2
    ! et - jm/2 <= Y <= jm/2

    ! x est la longitude du point en radians.
    ! y est la latitude du point en radians.
    ! 
    ! on a : cu_2d(i, j) = rad * cos(y) * dx/dX
    ! cv(j) = rad * dy/dY
    ! aire_2d(i, j) = cu_2d(i, j) * cv(j)
    ! 
    ! y, dx/dX, dy/dY calcules aux points concernes
    ! cv, bien que dependant de j uniquement, sera ici indice aussi en i
    ! pour un adressage plus facile en ij.

    ! aux points u et v, 
    ! xprimu et xprimv sont respectivement les valeurs de dx/dX
    ! yprimu et yprimv sont respectivement les valeurs de dy/dY
    ! rlatu et rlatv sont respectivement les valeurs de la latitude
    ! cvu et cv_2d sont respectivement les valeurs de cv_2d

    ! aux points u, v, scalaires, et z 
    ! cu_2d, cuv, cuscal, cuz sont respectivement les valeurs de cu_2d
    ! Cf. "inigeom.txt".

    USE dimens_m, ONLY : iim, jjm
    USE paramet_m, ONLY : iip1, jjp1
    USE comconst, ONLY : g, omeg, pi, rad
    USE comdissnew, ONLY : coefdis, nitergdiv, nitergrot, niterh
    USE logic, ONLY : fxyhypb, ysinus
    USE comgeom, ONLY : airesurg_2d, aireu_2d, airev_2d, aire_2d, &
         alpha1p2_2d, alpha1p4_2d, alpha1_2d, &
         alpha2p3_2d, alpha2_2d, alpha3p4_2d, alpha3_2d, alpha4_2d, apoln, &
         apols, constang_2d, cuscvugam_2d, cusurcvu_2d, cuvscvgam1_2d, &
         cuvscvgam2_2d, cuvsurcv_2d, cu_2d, cvscuvgam_2d, cvsurcuv_2d, &
         cvuscugam1_2d, cvuscugam2_2d, cvusurcu_2d, cv_2d, fext_2d, rlatu, &
         rlatv, rlonu, rlonv, unsairez_2d, unsaire_2d, unsairz_gam_2d, &
         unsair_gam1_2d, unsair_gam2_2d, unsapolnga1, unsapolnga2, &
         unsapolsga1, unsapolsga2, unscu2_2d, unscv2_2d, xprimu, xprimv
    USE serre, ONLY : alphax, alphay, clat, clon, dzoomx, dzoomy, grossismx, &
         grossismy, pxo, pyo, taux, tauy, transx, transy

    ! Variables locales

    INTEGER i, j, itmax, itmay, iter
    REAL cvu(iip1, jjp1), cuv(iip1, jjm)
    REAL ai14, ai23, airez, rlatp, rlatm, xprm, xprp, un4rad2, yprp, yprm
    REAL eps, x1, xo1, f, df, xdm, y1, yo1, ydm
    REAL coslatm, coslatp, radclatm, radclatp
    REAL cuij1(iip1, jjp1), cuij2(iip1, jjp1), cuij3(iip1, jjp1), &
         cuij4(iip1, jjp1)
    REAL cvij1(iip1, jjp1), cvij2(iip1, jjp1), cvij3(iip1, jjp1), &
         cvij4(iip1, jjp1)
    REAL rlonvv(iip1), rlatuu(jjp1)
    REAL rlatu1(jjm), yprimu1(jjm), rlatu2(jjm), yprimu2(jjm), yprimv(jjm), &
         yprimu(jjp1)
    REAL gamdi_gdiv, gamdi_grot, gamdi_h

    REAL rlonm025(iip1), xprimm025(iip1), rlonp025(iip1), xprimp025(iip1)
    SAVE rlatu1, yprimu1, rlatu2, yprimu2, yprimv, yprimu
    SAVE rlonm025, xprimm025, rlonp025, xprimp025

    real aireij1_2d(iim + 1, jjm + 1)
    real aireij2_2d(iim + 1, jjm + 1)
    real aireij3_2d(iim + 1, jjm + 1), aireij4_2d(iim + 1, jjm + 1) 
    real airuscv2_2d(iim + 1, jjm) 
    real airvscu2_2d(iim + 1, jjm), aiuscv2gam_2d(iim + 1, jjm) 
    real aivscu2gam_2d(iim + 1, jjm)

    !------------------------------------------------------------------

    PRINT *, 'Call sequence information: inigeom'

    IF (nitergdiv/=2) THEN
       gamdi_gdiv = coefdis/(real(nitergdiv)-2.)
    ELSE
       gamdi_gdiv = 0.
    END IF
    IF (nitergrot/=2) THEN
       gamdi_grot = coefdis/(real(nitergrot)-2.)
    ELSE
       gamdi_grot = 0.
    END IF
    IF (niterh/=2) THEN
       gamdi_h = coefdis/(real(niterh)-2.)
    ELSE
       gamdi_h = 0.
    END IF

    print *, 'gamdi_gdiv = ', gamdi_gdiv
    print *, "gamdi_grot = ", gamdi_grot
    print *, "gamdi_h = ", gamdi_h

    WRITE (6, 990)

    IF (.NOT. fxyhypb) THEN
       IF (ysinus) THEN
          print *, ' Inigeom, Y = Sinus (Latitude) '
          ! utilisation de f(x, y) avec y = sinus de la latitude 
          CALL fxysinus(rlatu, yprimu, rlatv, yprimv, rlatu1, yprimu1, &
               rlatu2, yprimu2, rlonu, xprimu, rlonv, xprimv, rlonm025, &
               xprimm025, rlonp025, xprimp025)
       ELSE
          print *, 'Inigeom, Y = Latitude, der. sinusoid .'
          ! utilisation de f(x, y) a tangente sinusoidale, y etant la latit

          pxo = clon*pi/180.
          pyo = 2.*clat*pi/180.

          ! determination de transx (pour le zoom) par Newton-Raphson

          itmax = 10
          eps = .1E-7

          xo1 = 0.
          DO iter = 1, itmax
             x1 = xo1
             f = x1 + alphax*sin(x1-pxo)
             df = 1. + alphax*cos(x1-pxo)
             x1 = x1 - f/df
             xdm = abs(x1-xo1)
             IF (xdm<=eps) EXIT
             xo1 = x1
          END DO

          transx = xo1

          itmay = 10
          eps = .1E-7

          yo1 = 0.
          DO iter = 1, itmay
             y1 = yo1
             f = y1 + alphay*sin(y1-pyo)
             df = 1. + alphay*cos(y1-pyo)
             y1 = y1 - f/df
             ydm = abs(y1-yo1)
             IF (ydm<=eps) EXIT
             yo1 = y1
          END DO

          transy = yo1

          CALL fxy(rlatu, yprimu, rlatv, yprimv, rlatu1, yprimu1, rlatu2, &
               yprimu2, rlonu, xprimu, rlonv, xprimv, rlonm025, xprimm025, &
               rlonp025, xprimp025)
       END IF
    ELSE
       ! Utilisation de fxyhyper, f(x, y) � d�riv�e tangente hyperbolique
       print *, 'Inigeom, Y = Latitude, d�riv�e tangente hyperbolique'
       CALL fxyhyper(clat, grossismy, dzoomy, tauy, clon, grossismx, dzoomx, &
            taux, rlatu, yprimu, rlatv, yprimv, rlatu1, yprimu1, rlatu2, &
            yprimu2, rlonu, xprimu, rlonv, xprimv, rlonm025, xprimm025, &
            rlonp025, xprimp025)
    END IF

    rlatu(1) = asin(1.)
    rlatu(jjp1) = -rlatu(1)

    ! calcul aux poles 

    yprimu(1) = 0.
    yprimu(jjp1) = 0.

    un4rad2 = 0.25*rad*rad

    ! calcul des aires (aire_2d, aireu_2d, airev_2d, 1./aire_2d, 1./airez)
    ! - et de fext_2d, force de coriolis extensive

    ! A 1 point scalaire P (i, j) de la grille, reguliere en (X, Y), sont
    ! affectees 4 aires entourant P, calculees respectivement aux points
    ! (i + 1/4, j - 1/4) : aireij1_2d (i, j)
    ! (i + 1/4, j + 1/4) : aireij2_2d (i, j)
    ! (i - 1/4, j + 1/4) : aireij3_2d (i, j)
    ! (i - 1/4, j - 1/4) : aireij4_2d (i, j)

    !,
    ! Les cotes de chacun de ces 4 carres etant egaux a 1/2 suivant (X, Y).
    ! Chaque aire centree en 1 point scalaire P(i, j) est egale a la somme
    ! des 4 aires aireij1_2d, aireij2_2d, aireij3_2d, aireij4_2d qui sont
    ! affectees au
    ! point (i, j).
    ! On definit en outre les coefficients alpha comme etant egaux a
    ! (aireij / aire_2d), c.a.d par exp.
    ! alpha1_2d(i, j)=aireij1_2d(i, j)/aire_2d(i, j)

    ! De meme, toute aire centree en 1 point U est egale a la somme des
    ! 4 aires aireij1_2d, aireij2_2d, aireij3_2d, aireij4_2d entourant
    ! le point U.
    ! Idem pour airev_2d, airez.

    ! On a, pour chaque maille : dX = dY = 1

    ! V

    ! aireij4_2d . . aireij1_2d

    ! U . . P . U

    ! aireij3_2d . . aireij2_2d

    ! V

    ! Calcul des 4 aires elementaires aireij1_2d, aireij2_2d,
    ! aireij3_2d, aireij4_2d
    ! qui entourent chaque aire_2d(i, j), ainsi que les 4 elongations
    ! elementaires
    ! cuij et les 4 elongat. cvij qui sont calculees aux memes
    ! endroits que les aireij.

    ! do 35 : boucle sur les jjm + 1 latitudes 

    DO j = 1, jjp1

       IF (j==1) THEN

          yprm = yprimu1(j)
          rlatm = rlatu1(j)

          coslatm = cos(rlatm)
          radclatm = 0.5*rad*coslatm

          DO i = 1, iim
             xprp = xprimp025(i)
             xprm = xprimm025(i)
             aireij2_2d(i, 1) = un4rad2*coslatm*xprp*yprm
             aireij3_2d(i, 1) = un4rad2*coslatm*xprm*yprm
             cuij2(i, 1) = radclatm*xprp
             cuij3(i, 1) = radclatm*xprm
             cvij2(i, 1) = 0.5*rad*yprm
             cvij3(i, 1) = cvij2(i, 1)
          END DO

          DO i = 1, iim
             aireij1_2d(i, 1) = 0.
             aireij4_2d(i, 1) = 0.
             cuij1(i, 1) = 0.
             cuij4(i, 1) = 0.
             cvij1(i, 1) = 0.
             cvij4(i, 1) = 0.
          END DO

       END IF

       IF (j==jjp1) THEN
          yprp = yprimu2(j-1)
          rlatp = rlatu2(j-1)

          coslatp = cos(rlatp)
          radclatp = 0.5*rad*coslatp

          DO i = 1, iim
             xprp = xprimp025(i)
             xprm = xprimm025(i)
             aireij1_2d(i, jjp1) = un4rad2*coslatp*xprp*yprp
             aireij4_2d(i, jjp1) = un4rad2*coslatp*xprm*yprp
             cuij1(i, jjp1) = radclatp*xprp
             cuij4(i, jjp1) = radclatp*xprm
             cvij1(i, jjp1) = 0.5*rad*yprp
             cvij4(i, jjp1) = cvij1(i, jjp1)
          END DO

          DO i = 1, iim
             aireij2_2d(i, jjp1) = 0.
             aireij3_2d(i, jjp1) = 0.
             cvij2(i, jjp1) = 0.
             cvij3(i, jjp1) = 0.
             cuij2(i, jjp1) = 0.
             cuij3(i, jjp1) = 0.
          END DO

       END IF

       IF (j>1 .AND. j<jjp1) THEN

          rlatp = rlatu2(j-1)
          yprp = yprimu2(j-1)
          rlatm = rlatu1(j)
          yprm = yprimu1(j)

          coslatm = cos(rlatm)
          coslatp = cos(rlatp)
          radclatp = 0.5*rad*coslatp
          radclatm = 0.5*rad*coslatm

          DO i = 1, iim
             xprp = xprimp025(i)
             xprm = xprimm025(i)

             ai14 = un4rad2*coslatp*yprp
             ai23 = un4rad2*coslatm*yprm
             aireij1_2d(i, j) = ai14*xprp
             aireij2_2d(i, j) = ai23*xprp
             aireij3_2d(i, j) = ai23*xprm
             aireij4_2d(i, j) = ai14*xprm
             cuij1(i, j) = radclatp*xprp
             cuij2(i, j) = radclatm*xprp
             cuij3(i, j) = radclatm*xprm
             cuij4(i, j) = radclatp*xprm
             cvij1(i, j) = 0.5*rad*yprp
             cvij2(i, j) = 0.5*rad*yprm
             cvij3(i, j) = cvij2(i, j)
             cvij4(i, j) = cvij1(i, j)
          END DO

       END IF

       ! periodicite 

       cvij1(iip1, j) = cvij1(1, j)
       cvij2(iip1, j) = cvij2(1, j)
       cvij3(iip1, j) = cvij3(1, j)
       cvij4(iip1, j) = cvij4(1, j)
       cuij1(iip1, j) = cuij1(1, j)
       cuij2(iip1, j) = cuij2(1, j)
       cuij3(iip1, j) = cuij3(1, j)
       cuij4(iip1, j) = cuij4(1, j)
       aireij1_2d(iip1, j) = aireij1_2d(1, j)
       aireij2_2d(iip1, j) = aireij2_2d(1, j)
       aireij3_2d(iip1, j) = aireij3_2d(1, j)
       aireij4_2d(iip1, j) = aireij4_2d(1, j)

    END DO

    DO j = 1, jjp1
       DO i = 1, iim
          aire_2d(i, j) = aireij1_2d(i, j) + aireij2_2d(i, j) &
               + aireij3_2d(i, j) + aireij4_2d(i, j)
          alpha1_2d(i, j) = aireij1_2d(i, j)/aire_2d(i, j)
          alpha2_2d(i, j) = aireij2_2d(i, j)/aire_2d(i, j)
          alpha3_2d(i, j) = aireij3_2d(i, j)/aire_2d(i, j)
          alpha4_2d(i, j) = aireij4_2d(i, j)/aire_2d(i, j)
          alpha1p2_2d(i, j) = alpha1_2d(i, j) + alpha2_2d(i, j)
          alpha1p4_2d(i, j) = alpha1_2d(i, j) + alpha4_2d(i, j)
          alpha2p3_2d(i, j) = alpha2_2d(i, j) + alpha3_2d(i, j)
          alpha3p4_2d(i, j) = alpha3_2d(i, j) + alpha4_2d(i, j)
       END DO

       aire_2d(iip1, j) = aire_2d(1, j)
       alpha1_2d(iip1, j) = alpha1_2d(1, j)
       alpha2_2d(iip1, j) = alpha2_2d(1, j)
       alpha3_2d(iip1, j) = alpha3_2d(1, j)
       alpha4_2d(iip1, j) = alpha4_2d(1, j)
       alpha1p2_2d(iip1, j) = alpha1p2_2d(1, j)
       alpha1p4_2d(iip1, j) = alpha1p4_2d(1, j)
       alpha2p3_2d(iip1, j) = alpha2p3_2d(1, j)
       alpha3p4_2d(iip1, j) = alpha3p4_2d(1, j)
    END DO

    DO j = 1, jjp1
       DO i = 1, iim
          aireu_2d(i, j) = aireij1_2d(i, j) + aireij2_2d(i, j) + &
               aireij4_2d(i+1, j) + aireij3_2d(i+1, j)
          unsaire_2d(i, j) = 1./aire_2d(i, j)
          unsair_gam1_2d(i, j) = unsaire_2d(i, j)**(-gamdi_gdiv)
          unsair_gam2_2d(i, j) = unsaire_2d(i, j)**(-gamdi_h)
          airesurg_2d(i, j) = aire_2d(i, j)/g
       END DO
       aireu_2d(iip1, j) = aireu_2d(1, j)
       unsaire_2d(iip1, j) = unsaire_2d(1, j)
       unsair_gam1_2d(iip1, j) = unsair_gam1_2d(1, j)
       unsair_gam2_2d(iip1, j) = unsair_gam2_2d(1, j)
       airesurg_2d(iip1, j) = airesurg_2d(1, j)
    END DO

    DO j = 1, jjm

       DO i = 1, iim
          airev_2d(i, j) = aireij2_2d(i, j) + aireij3_2d(i, j) + &
               aireij1_2d(i, j+1) + aireij4_2d(i, j+1)
       END DO
       DO i = 1, iim
          airez = aireij2_2d(i, j) + aireij1_2d(i, j+1) + aireij3_2d(i+1, j) &
               + aireij4_2d(i+1, j+1)
          unsairez_2d(i, j) = 1./airez
          unsairz_gam_2d(i, j) = unsairez_2d(i, j)**(-gamdi_grot)
          fext_2d(i, j) = airez*sin(rlatv(j))*2.*omeg
       END DO
       airev_2d(iip1, j) = airev_2d(1, j)
       unsairez_2d(iip1, j) = unsairez_2d(1, j)
       fext_2d(iip1, j) = fext_2d(1, j)
       unsairz_gam_2d(iip1, j) = unsairz_gam_2d(1, j)

    END DO

    ! Calcul des elongations cu_2d, cv_2d, cvu 

    DO j = 1, jjm
       DO i = 1, iim
          cv_2d(i, j) = 0.5 * &
               (cvij2(i, j) + cvij3(i, j) + cvij1(i, j+1) + cvij4(i, j+1))
          cvu(i, j) = 0.5*(cvij1(i, j)+cvij4(i, j)+cvij2(i, j)+cvij3(i, j))
          cuv(i, j) = 0.5*(cuij2(i, j)+cuij3(i, j)+cuij1(i, j+1)+cuij4(i, j+1))
          unscv2_2d(i, j) = 1./(cv_2d(i, j)*cv_2d(i, j))
       END DO
       DO i = 1, iim
          cuvsurcv_2d(i, j) = airev_2d(i, j)*unscv2_2d(i, j)
          cvsurcuv_2d(i, j) = 1./cuvsurcv_2d(i, j)
          cuvscvgam1_2d(i, j) = cuvsurcv_2d(i, j)**(-gamdi_gdiv)
          cuvscvgam2_2d(i, j) = cuvsurcv_2d(i, j)**(-gamdi_h)
          cvscuvgam_2d(i, j) = cvsurcuv_2d(i, j)**(-gamdi_grot)
       END DO
       cv_2d(iip1, j) = cv_2d(1, j)
       cvu(iip1, j) = cvu(1, j)
       unscv2_2d(iip1, j) = unscv2_2d(1, j)
       cuv(iip1, j) = cuv(1, j)
       cuvsurcv_2d(iip1, j) = cuvsurcv_2d(1, j)
       cvsurcuv_2d(iip1, j) = cvsurcuv_2d(1, j)
       cuvscvgam1_2d(iip1, j) = cuvscvgam1_2d(1, j)
       cuvscvgam2_2d(iip1, j) = cuvscvgam2_2d(1, j)
       cvscuvgam_2d(iip1, j) = cvscuvgam_2d(1, j)
    END DO

    DO j = 2, jjm
       DO i = 1, iim
          cu_2d(i, j) = 0.5 * (cuij1(i, j) + cuij4(i+1, j) + cuij2(i, j) &
               + cuij3(i+1, j))
          unscu2_2d(i, j) = 1./(cu_2d(i, j)*cu_2d(i, j))
          cvusurcu_2d(i, j) = aireu_2d(i, j)*unscu2_2d(i, j)
          cusurcvu_2d(i, j) = 1./cvusurcu_2d(i, j)
          cvuscugam1_2d(i, j) = cvusurcu_2d(i, j)**(-gamdi_gdiv)
          cvuscugam2_2d(i, j) = cvusurcu_2d(i, j)**(-gamdi_h)
          cuscvugam_2d(i, j) = cusurcvu_2d(i, j)**(-gamdi_grot)
       END DO
       cu_2d(iip1, j) = cu_2d(1, j)
       unscu2_2d(iip1, j) = unscu2_2d(1, j)
       cvusurcu_2d(iip1, j) = cvusurcu_2d(1, j)
       cusurcvu_2d(iip1, j) = cusurcvu_2d(1, j)
       cvuscugam1_2d(iip1, j) = cvuscugam1_2d(1, j)
       cvuscugam2_2d(iip1, j) = cvuscugam2_2d(1, j)
       cuscvugam_2d(iip1, j) = cuscvugam_2d(1, j)
    END DO

    ! calcul aux poles 

    DO i = 1, iip1
       cu_2d(i, 1) = 0.
       unscu2_2d(i, 1) = 0.
       cvu(i, 1) = 0.

       cu_2d(i, jjp1) = 0.
       unscu2_2d(i, jjp1) = 0.
       cvu(i, jjp1) = 0.
    END DO

    DO j = 1, jjm
       DO i = 1, iim
          airvscu2_2d(i, j) = airev_2d(i, j)/(cuv(i, j)*cuv(i, j))
          aivscu2gam_2d(i, j) = airvscu2_2d(i, j)**(-gamdi_grot)
       END DO
       airvscu2_2d(iip1, j) = airvscu2_2d(1, j)
       aivscu2gam_2d(iip1, j) = aivscu2gam_2d(1, j)
    END DO

    DO j = 2, jjm
       DO i = 1, iim
          airuscv2_2d(i, j) = aireu_2d(i, j)/(cvu(i, j)*cvu(i, j))
          aiuscv2gam_2d(i, j) = airuscv2_2d(i, j)**(-gamdi_grot)
       END DO
       airuscv2_2d(iip1, j) = airuscv2_2d(1, j)
       aiuscv2gam_2d(iip1, j) = aiuscv2gam_2d(1, j)
    END DO

    ! calcul des aires aux poles :

    apoln = sum(aire_2d(:iim, 1))
    apols = sum(aire_2d(:iim, jjp1))
    unsapolnga1 = 1./(apoln**(-gamdi_gdiv))
    unsapolsga1 = 1./(apols**(-gamdi_gdiv))
    unsapolnga2 = 1./(apoln**(-gamdi_h))
    unsapolsga2 = 1./(apols**(-gamdi_h))

    ! changement F. Hourdin calcul conservatif pour fext_2d
    ! constang_2d contient le produit a * cos (latitude) * omega

    DO i = 1, iim
       constang_2d(i, 1) = 0.
    END DO
    DO j = 1, jjm - 1
       DO i = 1, iim
          constang_2d(i, j+1) = rad*omeg*cu_2d(i, j+1)*cos(rlatu(j+1))
       END DO
    END DO
    DO i = 1, iim
       constang_2d(i, jjp1) = 0.
    END DO

    ! periodicite en longitude

    DO j = 1, jjm
       fext_2d(iip1, j) = fext_2d(1, j)
    END DO
    DO j = 1, jjp1
       constang_2d(iip1, j) = constang_2d(1, j)
    END DO

    ! fin du changement

    print *, ' Coordonnees de la grille '
    WRITE (6, 995)

    print *, ' LONGITUDES aux pts. V (degres) '
    WRITE (6, 995)
    DO i = 1, iip1
       rlonvv(i) = rlonv(i)*180./pi
    END DO
    WRITE (6, 400) rlonvv

    WRITE (6, 995)
    print *, ' LATITUDES aux pts. V (degres) '
    WRITE (6, 995)
    DO i = 1, jjm
       rlatuu(i) = rlatv(i)*180./pi
    END DO
    WRITE (6, 400) (rlatuu(i), i=1, jjm)

    DO i = 1, iip1
       rlonvv(i) = rlonu(i)*180./pi
    END DO
    WRITE (6, 995)
    print *, ' LONGITUDES aux pts. U (degres) '
    WRITE (6, 995)
    WRITE (6, 400) rlonvv
    WRITE (6, 995)

    print *, ' LATITUDES aux pts. U (degres) '
    WRITE (6, 995)
    DO i = 1, jjp1
       rlatuu(i) = rlatu(i)*180./pi
    END DO
    WRITE (6, 400) (rlatuu(i), i=1, jjp1)
    WRITE (6, 995)

400 FORMAT (1X, 8F8.2)
990 FORMAT (//)
995 FORMAT (/)

  END SUBROUTINE inigeom

end module inigeom_m
