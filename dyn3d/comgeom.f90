module comgeom

  implicit none

  real, pointer:: cu_2d(:, :) ! (iim + 1, jjm + 1) in m
  real, pointer:: cv_2d(:, :) ! (iim + 1, jjm) in m
  real, allocatable, target:: cu(:) ! ((iim + 1) * (jjm + 1)) in m
  real, allocatable, target:: cv(:) ! ((iim + 1) * jjm) in m

  real, allocatable:: unscu2_2d(:, :) ! (iim + 1, jjm + 1) in m-2
  real, allocatable:: unscv2_2d(:, :) ! (iim + 1, jjm) in m-2

  real, allocatable, target:: aire(:) ! ((iim + 1) * (jjm + 1)) in m2
  real, pointer:: aire_2d(:, :) ! (iim + 1, jjm + 1) in m2

  real, allocatable:: airesurg_2d(:, :) ! (iim + 1, jjm + 1)
  
  real, pointer:: aireu_2d(:, :) ! (iim + 1, jjm + 1) in m2
  real, allocatable, target:: aireu(:) ! ((iim + 1) * (jjm + 1)) in m2

  real, allocatable, target:: airev(:) ! ((iim + 1) * jjm)
  real, pointer:: airev_2d(:, :) ! (iim + 1, jjm) in m2

  real, allocatable:: unsaire_2d(:, :) ! (iim + 1, jjm + 1) in m-2
  real apoln, apols ! in m2
  real, allocatable:: unsairez_2d(:, :) ! (iim + 1, jjm)
  
  real, pointer:: alpha1_2d(:, :) ! (iim + 1, jjm + 1)
  real, allocatable, target:: alpha1(:) ! ((iim + 1) * (jjm + 1))

  real, pointer:: alpha2_2d(:, :) ! (iim + 1, jjm + 1) 
  real, allocatable, target:: alpha2(:) ! ((iim + 1) * (jjm + 1))

  real, pointer:: alpha3_2d(:, :), alpha4_2d(:, :) ! (iim + 1, jjm + 1)
  real, allocatable, target:: alpha3(:), alpha4(:) ! ((iim + 1) * (jjm + 1))

  real, pointer:: alpha1p2_2d(:, :) ! (iim + 1, jjm + 1) 
  real, allocatable, target:: alpha1p2(:) ! ((iim + 1) * (jjm + 1))

  real, pointer:: alpha1p4_2d(:, :), alpha2p3_2d(:, :) ! (iim + 1, jjm + 1)
  real, allocatable, target:: alpha1p4(:), alpha2p3(:) ! ((iim + 1) * (jjm + 1))

  real, allocatable, target:: alpha3p4(:) ! ((iim + 1) * (jjm + 1))
  real, pointer:: alpha3p4_2d(:, :) ! (iim + 1, jjm + 1) 

  real, pointer:: fext_2d(:, :) ! (iim + 1, jjm)
  real, pointer:: constang_2d(:, :) ! (iim + 1, jjm + 1)
  real, allocatable, target:: fext(:) ! ((iim + 1) * jjm)
  real, allocatable, target:: constang(:) ! ((iim + 1) * (jjm + 1))

  real, allocatable:: cuvsurcv_2d(:, :) ! (iim + 1, jjm) no dimension
  
  real, pointer:: cvsurcuv_2d(:, :) ! (iim + 1, jjm) no dimension
  real, allocatable, target:: cvsurcuv(:) ! ((iim + 1) * jjm) no dimension

  real, allocatable:: cvusurcu_2d(:, :) ! (iim + 1, jjm + 1) no dimension
  
  real, pointer:: cusurcvu_2d(:, :) ! (iim + 1, jjm + 1) no dimension
  real, allocatable, target:: cusurcvu(:) ! ((iim + 1) * (jjm + 1)) no dimension

  real, pointer:: cuvscvgam1_2d(:, :) ! (iim + 1, jjm)
  real, allocatable, target:: cuvscvgam1(:) ! ((iim + 1) * jjm)

  real, pointer:: cuvscvgam2_2d(:, :) ! (iim + 1, jjm)
  real, pointer:: cvuscugam1_2d(:, :) ! (iim + 1, jjm + 1)
  real, allocatable, target:: cuvscvgam2(:) ! ((iim + 1) * jjm)
  real, allocatable, target:: cvuscugam1(:) ! ((iim + 1) * (jjm + 1))

  real, pointer:: cvuscugam2_2d(:, :) ! (iim + 1, jjm + 1)
  real, pointer:: cvscuvgam_2d(:, :) ! (iim + 1, jjm)
  real, allocatable, target:: cvuscugam2(:) ! ((iim + 1) * (jjm + 1))
  real, allocatable, target:: cvscuvgam(:) ! ((iim + 1) * jjm)

  real, allocatable, target:: cuscvugam(:) ! ((iim + 1) * (jjm + 1))
  real, pointer:: cuscvugam_2d(:, :) ! (iim + 1, jjm + 1) 

  real unsapolnga1, unsapolnga2, unsapolsga1, unsapolsga2 

  real, pointer:: unsair_gam1_2d(:, :), unsair_gam2_2d(:, :)
  ! (iim + 1, jjm + 1)

  real, allocatable, target:: unsair_gam1(:), unsair_gam2(:)
  ! ((iim + 1) * (jjm + 1))

  real, pointer:: unsairz_gam_2d(:, :) ! (iim + 1, jjm)
  real, allocatable, target:: unsairz_gam(:) ! ((iim + 1) * jjm)

  save

contains

  SUBROUTINE inigeom

    ! Auteur : P. Le Van

    ! Calcul des élongations cuij1, ..., cuij4, cvij1, ..., cvij4 aux mêmes
    ! endroits que les aires aireij1_2d, ..., aireij4_2d.

    ! Calcul des coefficients cu_2d, cv_2d, 1. / cu_2d**2, 1. /
    ! cv_2d**2. Les coefficients cu_2d et cv_2d permettent de passer
    ! des vitesses naturelles aux vitesses covariantes et
    ! contravariantes, ou vice-versa.

    ! On a :
    ! u(covariant) = cu_2d * u(naturel), u(contravariant) = u(naturel) / cu_2d
    ! v(covariant) = cv_2d * v(naturel), v(contravariant) = v(naturel) / cv_2d

    ! On en tire : 
    ! u(covariant) = cu_2d * cu_2d * u(contravariant)
    ! v(covariant) = cv_2d * cv_2d * v(contravariant)

    ! x est la longitude du point en radians.
    ! y est la latitude du point en radians.
    ! 
    ! On a : cu_2d(i, j) = ra * cos(y) * dx / dX
    ! cv(j) = ra * dy / dY
    ! aire_2d(i, j) = cu_2d(i, j) * cv(j)
    ! 
    ! y, dx / dX, dy / dY calculés aux points concernés. cv, bien que
    ! dépendant de j uniquement, sera ici indicé aussi en i pour un
    ! adressage plus facile en ij.

    ! cv_2d est aux points v. cu_2d est aux points u. Cf. "inigeom.txt".

    USE comconst, ONLY : g, omeg, ra
    USE comdissnew, ONLY : coefdis, nitergdiv, nitergrot, niterh
    use dimensions, only: iim, jjm
    use dynetat0_m, only: xprimp025, xprimm025, rlatu1, rlatu2, rlatu, rlatv, &
         yprimu1, yprimu2
    USE paramet_m, ONLY : iip1, jjp1

    ! Local:
    INTEGER i, j
    REAL ai14, ai23, airez, un4rad2
    REAL coslatm, coslatp, radclatm, radclatp
    REAL, dimension(iip1, jjp1):: cuij1, cuij2, cuij3, cuij4 ! in m
    REAL, dimension(iip1, jjp1):: cvij1, cvij2, cvij3, cvij4 ! in m
    REAL gamdi_gdiv, gamdi_grot, gamdi_h
    real, dimension(iim + 1, jjm + 1):: aireij1_2d, aireij2_2d, aireij3_2d, &
         aireij4_2d ! in m2

    !------------------------------------------------------------------

    PRINT *, 'Call sequence information: inigeom'

    allocate(cu((iim + 1) * (jjm + 1)), cv((iim + 1) * jjm))
    allocate(unscu2_2d(iim + 1, jjm + 1))
    allocate(unscv2_2d(iim + 1, jjm))
    allocate(aire((iim + 1) * (jjm + 1)))
    allocate(aireu((iim + 1) * (jjm + 1)))
    allocate(airev((iim + 1) * jjm))
    allocate(alpha1((iim + 1) * (jjm + 1)))
    allocate(alpha2((iim + 1) * (jjm + 1)))
    allocate(alpha3((iim + 1) * (jjm + 1)), alpha4((iim + 1) * (jjm + 1)))
    allocate(alpha1p2((iim + 1) * (jjm + 1)))
    allocate(alpha1p4((iim + 1) * (jjm + 1)), alpha2p3((iim + 1) * (jjm + 1)))
    allocate(alpha3p4((iim + 1) * (jjm + 1)))
    allocate(fext((iim + 1) * jjm), constang((iim + 1) * (jjm + 1)))
    allocate(cvsurcuv((iim + 1) * jjm))
    allocate(cusurcvu((iim + 1) * (jjm + 1)))
    allocate(cuvscvgam1((iim + 1) * jjm))
    allocate(cuvscvgam2((iim + 1) * jjm), cvuscugam1((iim + 1) * (jjm + 1)))
    allocate(cvuscugam2((iim + 1) * (jjm + 1)), cvscuvgam((iim + 1) * jjm))
    allocate(cuscvugam((iim + 1) * (jjm + 1)))
    allocate(unsair_gam1((iim + 1) * (jjm + 1)), &
         unsair_gam2((iim + 1) * (jjm + 1)))
    allocate(unsairz_gam((iim + 1) * jjm))

    allocate(airesurg_2d(iim + 1, jjm + 1))
    allocate(unsaire_2d(iim + 1, jjm + 1))
    allocate(unsairez_2d(iim + 1, jjm))
    allocate(cuvsurcv_2d(iim + 1, jjm))
    allocate(cvusurcu_2d(iim + 1, jjm + 1))

    cu_2d(1:iim + 1, 1:jjm + 1) => cu
    cv_2d(1:iim + 1, 1:jjm) => cv
    aire_2d(1:iim + 1, 1:jjm + 1) => aire
    aireu_2d(1:iim + 1, 1:jjm + 1) => aireu
    airev_2d(1:iim + 1, 1:jjm) => airev
    alpha1_2d(1:iim + 1, 1:jjm + 1) => alpha1
    alpha2_2d(1:iim + 1, 1:jjm + 1) => alpha2
    alpha3_2d(1:iim + 1, 1:jjm + 1) => alpha3
    alpha4_2d(1:iim + 1, 1:jjm + 1) => alpha4
    alpha1p2_2d(1:iim + 1, 1:jjm + 1) => alpha1p2
    alpha1p4_2d(1:iim + 1, 1:jjm + 1) => alpha1p4
    alpha2p3_2d(1:iim + 1, 1:jjm + 1) => alpha2p3
    alpha3p4_2d(1:iim + 1, 1:jjm + 1) => alpha3p4
    fext_2d(1:iim + 1, 1:jjm) => fext
    constang_2d(1:iim + 1, 1:jjm + 1) => constang
    cvsurcuv_2d(1:iim + 1, 1:jjm) => cvsurcuv
    cusurcvu_2d(1:iim + 1, 1:jjm + 1) => cusurcvu
    cuvscvgam1_2d(1:iim + 1, 1:jjm) => cuvscvgam1
    cuvscvgam2_2d(1:iim + 1, 1:jjm) => cuvscvgam2
    cvuscugam1_2d(1:iim + 1, 1:jjm + 1) => cvuscugam1
    cvuscugam2_2d(1:iim + 1, 1:jjm + 1) => cvuscugam2
    cvscuvgam_2d(1:iim + 1, 1:jjm) => cvscuvgam
    cuscvugam_2d(1:iim + 1, 1:jjm + 1) => cuscvugam
    unsair_gam1_2d(1:iim + 1, 1:jjm + 1) => unsair_gam1
    unsair_gam2_2d(1:iim + 1, 1:jjm + 1) => unsair_gam2
    unsairz_gam_2d(1:iim + 1, 1:jjm) => unsairz_gam

    IF (nitergdiv /= 2) THEN
       gamdi_gdiv = coefdis / (nitergdiv - 2)
    ELSE
       gamdi_gdiv = 0.
    END IF

    IF (nitergrot /= 2) THEN
       gamdi_grot = coefdis / (nitergrot - 2)
    ELSE
       gamdi_grot = 0.
    END IF

    IF (niterh /= 2) THEN
       gamdi_h = coefdis / (niterh - 2)
    ELSE
       gamdi_h = 0.
    END IF

    print *, 'gamdi_gdiv = ', gamdi_gdiv
    print *, "gamdi_grot = ", gamdi_grot
    print *, "gamdi_h = ", gamdi_h

    un4rad2 = 0.25 * ra * ra

    ! Cf. "inigeom.txt". Calcul des quatre aires élémentaires
    ! aireij1_2d, aireij2_2d, aireij3_2d, aireij4_2d qui entourent
    ! chaque aire_2d(i, j), ainsi que les quatre élongations
    ! élémentaires cuij et les quatre élongations cvij qui sont
    ! calculées aux mêmes endroits que les aireij.

    coslatm = cos(rlatu1(1))
    radclatm = 0.5 * ra * coslatm

    aireij1_2d(:iim, 1) = 0.
    aireij2_2d(:iim, 1) = un4rad2 * coslatm * xprimp025(:iim) * yprimu1(1)
    aireij3_2d(:iim, 1) = un4rad2 * coslatm * xprimm025(:iim) * yprimu1(1)
    aireij4_2d(:iim, 1) = 0.

    cuij1(:iim, 1) = 0.
    cuij2(:iim, 1) = radclatm * xprimp025(:iim)
    cuij3(:iim, 1) = radclatm * xprimm025(:iim)
    cuij4(:iim, 1) = 0.

    cvij1(:iim, 1) = 0.
    cvij2(:iim, 1) = 0.5 * ra * yprimu1(1)
    cvij3(:iim, 1) = cvij2(:iim, 1)
    cvij4(:iim, 1) = 0.

    do j = 2, jjm
       coslatm = cos(rlatu1(j))
       coslatp = cos(rlatu2(j-1))
       radclatp = 0.5 * ra * coslatp
       radclatm = 0.5 * ra * coslatm
       ai14 = un4rad2 * coslatp * yprimu2(j-1)
       ai23 = un4rad2 * coslatm * yprimu1(j)

       aireij1_2d(:iim, j) = ai14 * xprimp025(:iim)
       aireij2_2d(:iim, j) = ai23 * xprimp025(:iim)
       aireij3_2d(:iim, j) = ai23 * xprimm025(:iim)
       aireij4_2d(:iim, j) = ai14 * xprimm025(:iim)
       cuij1(:iim, j) = radclatp * xprimp025(:iim)
       cuij2(:iim, j) = radclatm * xprimp025(:iim)
       cuij3(:iim, j) = radclatm * xprimm025(:iim)
       cuij4(:iim, j) = radclatp * xprimm025(:iim)
       cvij1(:iim, j) = 0.5 * ra * yprimu2(j-1)
       cvij2(:iim, j) = 0.5 * ra * yprimu1(j)
       cvij3(:iim, j) = cvij2(:iim, j)
       cvij4(:iim, j) = cvij1(:iim, j)
    end do

    coslatp = cos(rlatu2(jjm))
    radclatp = 0.5 * ra * coslatp

    aireij1_2d(:iim, jjp1) = un4rad2 * coslatp * xprimp025(:iim) * yprimu2(jjm)
    aireij2_2d(:iim, jjp1) = 0.
    aireij3_2d(:iim, jjp1) = 0.
    aireij4_2d(:iim, jjp1) = un4rad2 * coslatp * xprimm025(:iim) * yprimu2(jjm)

    cuij1(:iim, jjp1) = radclatp * xprimp025(:iim)
    cuij2(:iim, jjp1) = 0.
    cuij3(:iim, jjp1) = 0.
    cuij4(:iim, jjp1) = radclatp * xprimm025(:iim)

    cvij1(:iim, jjp1) = 0.5 * ra * yprimu2(jjm)
    cvij2(:iim, jjp1) = 0.
    cvij3(:iim, jjp1) = 0.
    cvij4(:iim, jjp1) = cvij1(:iim, jjp1)

    ! Périodicité :

    cvij1(iip1, :) = cvij1(1, :)
    cvij2(iip1, :) = cvij2(1, :)
    cvij3(iip1, :) = cvij3(1, :)
    cvij4(iip1, :) = cvij4(1, :)

    cuij1(iip1, :) = cuij1(1, :)
    cuij2(iip1, :) = cuij2(1, :)
    cuij3(iip1, :) = cuij3(1, :)
    cuij4(iip1, :) = cuij4(1, :)

    aireij1_2d(iip1, :) = aireij1_2d(1, :)
    aireij2_2d(iip1, :) = aireij2_2d(1, :)
    aireij3_2d(iip1, :) = aireij3_2d(1, :)
    aireij4_2d(iip1, :) = aireij4_2d(1, :)

    DO j = 1, jjp1
       DO i = 1, iim
          aire_2d(i, j) = aireij1_2d(i, j) + aireij2_2d(i, j) &
               + aireij3_2d(i, j) + aireij4_2d(i, j)
          alpha1_2d(i, j) = aireij1_2d(i, j) / aire_2d(i, j)
          alpha2_2d(i, j) = aireij2_2d(i, j) / aire_2d(i, j)
          alpha3_2d(i, j) = aireij3_2d(i, j) / aire_2d(i, j)
          alpha4_2d(i, j) = aireij4_2d(i, j) / aire_2d(i, j)
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
               aireij4_2d(i + 1, j) + aireij3_2d(i + 1, j)
          unsaire_2d(i, j) = 1. / aire_2d(i, j)
          unsair_gam1_2d(i, j) = unsaire_2d(i, j)**(-gamdi_gdiv)
          unsair_gam2_2d(i, j) = unsaire_2d(i, j)**(-gamdi_h)
          airesurg_2d(i, j) = aire_2d(i, j) / g
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
               aireij1_2d(i, j + 1) + aireij4_2d(i, j + 1)
       END DO
       DO i = 1, iim
          airez = aireij2_2d(i, j) + aireij1_2d(i, j + 1) &
               + aireij3_2d(i + 1, j) + aireij4_2d(i + 1, j + 1)
          unsairez_2d(i, j) = 1. / airez
          unsairz_gam_2d(i, j) = unsairez_2d(i, j)**(-gamdi_grot)
          fext_2d(i, j) = airez * sin(rlatv(j)) * 2. * omeg
       END DO
       airev_2d(iip1, j) = airev_2d(1, j)
       unsairez_2d(iip1, j) = unsairez_2d(1, j)
       fext_2d(iip1, j) = fext_2d(1, j)
       unsairz_gam_2d(iip1, j) = unsairz_gam_2d(1, j)
    END DO

    ! Calcul des élongations cu_2d, cv_2d

    DO j = 1, jjm
       DO i = 1, iim
          cv_2d(i, j) = 0.5 * &
               (cvij2(i, j) + cvij3(i, j) + cvij1(i, j + 1) + cvij4(i, j + 1))
          unscv2_2d(i, j) = 1. / cv_2d(i, j)**2
       END DO
       DO i = 1, iim
          cuvsurcv_2d(i, j) = airev_2d(i, j) * unscv2_2d(i, j)
          cvsurcuv_2d(i, j) = 1. / cuvsurcv_2d(i, j)
          cuvscvgam1_2d(i, j) = cuvsurcv_2d(i, j)**(-gamdi_gdiv)
          cuvscvgam2_2d(i, j) = cuvsurcv_2d(i, j)**(-gamdi_h)
          cvscuvgam_2d(i, j) = cvsurcuv_2d(i, j)**(-gamdi_grot)
       END DO
       cv_2d(iip1, j) = cv_2d(1, j)
       unscv2_2d(iip1, j) = unscv2_2d(1, j)
       cuvsurcv_2d(iip1, j) = cuvsurcv_2d(1, j)
       cvsurcuv_2d(iip1, j) = cvsurcuv_2d(1, j)
       cuvscvgam1_2d(iip1, j) = cuvscvgam1_2d(1, j)
       cuvscvgam2_2d(iip1, j) = cuvscvgam2_2d(1, j)
       cvscuvgam_2d(iip1, j) = cvscuvgam_2d(1, j)
    END DO

    DO j = 2, jjm
       DO i = 1, iim
          cu_2d(i, j) = 0.5 * (cuij1(i, j) + cuij4(i + 1, j) + cuij2(i, j) &
               + cuij3(i + 1, j))
          unscu2_2d(i, j) = 1. / cu_2d(i, j)**2
          cvusurcu_2d(i, j) = aireu_2d(i, j) * unscu2_2d(i, j)
          cusurcvu_2d(i, j) = 1. / cvusurcu_2d(i, j)
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

    ! Calcul aux pôles

    cu_2d(:, 1) = 0.
    unscu2_2d(:, 1) = 0.

    cu_2d(:, jjp1) = 0.
    unscu2_2d(:, jjp1) = 0.

    ! Calcul des aires aux pôles :

    apoln = sum(aire_2d(:iim, 1))
    apols = sum(aire_2d(:iim, jjp1))
    unsapolnga1 = 1. / (apoln**(-gamdi_gdiv))
    unsapolsga1 = 1. / (apols**(-gamdi_gdiv))
    unsapolnga2 = 1. / (apoln**(-gamdi_h))
    unsapolsga2 = 1. / (apols**(-gamdi_h))

    ! Changement F. Hourdin calcul conservatif pour fext_2d
    ! constang_2d contient le produit a * cos (latitude) * omega

    DO i = 1, iim
       constang_2d(i, 1) = 0.
    END DO
    DO j = 1, jjm - 1
       DO i = 1, iim
          constang_2d(i, j + 1) = ra * omeg * cu_2d(i, j + 1) &
               * cos(rlatu(j + 1))
       END DO
    END DO
    DO i = 1, iim
       constang_2d(i, jjp1) = 0.
    END DO

    ! Périodicité en longitude
    DO j = 1, jjp1
       constang_2d(iip1, j) = constang_2d(1, j)
    END DO

  END SUBROUTINE inigeom

end module comgeom
