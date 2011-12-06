module inifilr_m

  use dimens_m, only: iim

  IMPLICIT NONE

  INTEGER jfiltnu, jfiltsu, jfiltnv, jfiltsv
  INTEGER, PARAMETER:: nfilun=3, nfilus=2, nfilvn=2, nfilvs=2

  real matriceun(iim,iim,nfilun), matriceus(iim,iim,nfilus)
  real matricevn(iim,iim,nfilvn), matricevs(iim,iim,nfilvs)
  real matrinvn(iim,iim,nfilun), matrinvs(iim,iim,nfilus)

  private iim, nfilun, nfilus, nfilvn, nfilvs

contains

  SUBROUTINE inifilr

    ! From filtrez/inifilr.F, version 1.1.1.1 2004/05/19 12:53:09
    ! H. Upadhyaya, O. Sharma

    ! This routine computes the eigenfunctions of the laplacian on the
    ! stretched grid, and the filtering coefficients.
    ! We designate:
    ! eignfn eigenfunctions of the discrete laplacian
    ! eigenvl eigenvalues
    ! jfiltn index of the last scalar line filtered in NH
    ! jfilts index of the first line filtered in SH
    ! modfrst index of the mode from where modes are filtered
    ! modemax maximum number of modes (im)
    ! coefil filtering coefficients (lamda_max * cos(rlat) / lamda)
    ! sdd SQRT(dx)

    ! The modes are filtered from modfrst to modemax.

    USE dimens_m, ONLY : iim, jjm
    USE logic, ONLY : fxyhypb, ysinus
    USE comgeom, ONLY : rlatu, rlatv, xprimu
    use nr_util, only: pi
    USE serre, ONLY : alphax
    USE coefils, ONLY : coefilu, coefilu2, coefilv, coefilv2, eignfnu, &
         eignfnv, modfrstu, modfrstv

    ! Local:
    REAL dlonu(iim), dlatu(jjm)
    REAL rlamda(2: iim), eignvl(iim)

    REAL lamdamax, cof
    INTEGER i, j, modemax, imx, k, kf
    REAL dymin, dxmin, colat0
    REAL eignft(iim, iim), coff
    EXTERNAL inifgn

    !-----------------------------------------------------------

    print *, "Call sequence information: inifilr"

    DO i = 1, iim
       dlonu(i) = xprimu(i)
    END DO

    CALL inifgn(eignvl)

    PRINT *, 'EIGNVL '
    PRINT "(1X, 5E13.6)", eignvl

    ! compute eigenvalues and eigenfunctions
    ! compute the filtering coefficients for scalar lines and
    ! meridional wind v-lines
    ! we filter all those latitude lines where coefil < 1
    ! NO FILTERING AT POLES
    ! colat0 is to be used when alpha (stretching coefficient)
    ! is set equal to zero for the regular grid case

    ! Calcul de colat0

    DO j = 1, jjm
       dlatu(j) = rlatu(j) - rlatu(j+1)
    END DO

    dxmin = dlonu(1)
    DO i = 2, iim
       dxmin = min(dxmin, dlonu(i))
    END DO
    dymin = dlatu(1)
    DO j = 2, jjm
       dymin = min(dymin, dlatu(j))
    END DO

    colat0 = min(0.5, dymin/dxmin)

    IF (.NOT. fxyhypb .AND. ysinus) THEN
       colat0 = 0.6
       ! À revoir pour ysinus
       alphax = 0.
    END IF

    PRINT *, 'colat0 = ', colat0
    PRINT *, 'alphax = ', alphax

    IF (alphax == 1.) THEN
       PRINT *, 'alphax doit etre < a 1. Corriger '
       STOP 1
    END IF

    lamdamax = iim / (pi * colat0 * (1. - alphax))
    rlamda = lamdamax / sqrt(abs(eignvl(2: iim)))

    DO j = 1, jjm
       DO i = 1, iim
          coefilu(i, j) = 0.
          coefilv(i, j) = 0.
          coefilu2(i, j) = 0.
          coefilv2(i, j) = 0.
       end DO
    END DO

    ! Determination de jfiltnu, jfiltnv, jfiltsu, jfiltsv

    modemax = iim
    imx = iim

    PRINT *, 'TRUNCATION AT ', imx

    DO j = 2, jjm / 2 + 1
       IF (cos(rlatu(j)) / colat0 < 1. &
            .and. rlamda(imx) * cos(rlatu(j)) < 1.) jfiltnu = j

       IF (cos(rlatu(jjm - j + 2)) / colat0 < 1. &
            .and. rlamda(imx) * cos(rlatu(jjm - j + 2)) < 1.) &
            jfiltsu = jjm - j + 2
    END DO

    DO j = 1, jjm/2
       cof = cos(rlatv(j))/colat0
       IF (cof < 1.) THEN
          IF (rlamda(imx)*cos(rlatv(j)) < 1.) jfiltnv = j
       END IF

       cof = cos(rlatv(jjm-j+1))/colat0
       IF (cof < 1.) THEN
          IF (rlamda(imx)*cos(rlatv(jjm-j+1)) < 1.) jfiltsv = jjm - j + 1
       END IF
    END DO

    IF (jfiltnu <= 0) jfiltnu = 1
    IF (jfiltnu > jjm/2+1) THEN
       PRINT *, 'jfiltnu en dehors des valeurs acceptables ', jfiltnu
       STOP 1
    END IF

    IF (jfiltsu <= 0) jfiltsu = 1
    IF (jfiltsu > jjm + 1) THEN
       PRINT *, 'jfiltsu en dehors des valeurs acceptables ', jfiltsu
       STOP 1
    END IF

    IF (jfiltnv <= 0) jfiltnv = 1
    IF (jfiltnv > jjm/2) THEN
       PRINT *, 'jfiltnv en dehors des valeurs acceptables ', jfiltnv
       STOP 1
    END IF

    IF (jfiltsv <= 0) jfiltsv = 1
    IF (jfiltsv > jjm) THEN
       PRINT *, 'jfiltsv en dehors des valeurs acceptables ', jfiltsv
       STOP 1
    END IF

    PRINT *, 'jfiltnv jfiltsv jfiltnu jfiltsu ', jfiltnv, jfiltsv, jfiltnu, &
         jfiltsu

    ! Determination de coefilu, coefilv, n=modfrstu, modfrstv

    DO j = 1, jjm
       modfrstu(j) = iim
       modfrstv(j) = iim
    END DO

    DO j = 2, jfiltnu
       DO k = 2, modemax
          cof = rlamda(k) * cos(rlatu(j))
          IF (cof < 1.) exit
       end DO
       if (k == modemax + 1) cycle
       modfrstu(j) = k

       kf = modfrstu(j)
       DO k = kf, modemax
          cof = rlamda(k)*cos(rlatu(j))
          coefilu(k, j) = cof - 1.
          coefilu2(k, j) = cof*cof - 1.
       end DO
    END DO

    DO j = 1, jfiltnv
       DO k = 2, modemax
          cof = rlamda(k)*cos(rlatv(j))
          IF (cof < 1.) exit
       end DO
       if (k == modemax + 1) cycle
       modfrstv(j) = k

       kf = modfrstv(j)
       DO k = kf, modemax
          cof = rlamda(k)*cos(rlatv(j))
          coefilv(k, j) = cof - 1.
          coefilv2(k, j) = cof*cof - 1.
       end DO
    end DO

    DO j = jfiltsu, jjm
       DO k = 2, modemax
          cof = rlamda(k)*cos(rlatu(j))
          IF (cof < 1.) exit
       end DO
       if (k == modemax + 1) cycle
       modfrstu(j) = k

       kf = modfrstu(j)
       DO k = kf, modemax
          cof = rlamda(k)*cos(rlatu(j))
          coefilu(k, j) = cof - 1.
          coefilu2(k, j) = cof*cof - 1.
       end DO
    end DO

    DO j = jfiltsv, jjm
       DO k = 2, modemax
          cof = rlamda(k)*cos(rlatv(j))
          IF (cof < 1.) exit
       end DO
       if (k == modemax + 1) cycle
       modfrstv(j) = k

       kf = modfrstv(j)
       DO k = kf, modemax
          cof = rlamda(k)*cos(rlatv(j))
          coefilv(k, j) = cof - 1.
          coefilv2(k, j) = cof*cof - 1.
       end DO
    END DO

    IF (jfiltnv>=jjm/2 .OR. jfiltnu>=jjm/2) THEN
       IF (jfiltnv == jfiltsv) jfiltsv = 1 + jfiltnv
       IF (jfiltnu == jfiltsu) jfiltsu = 1 + jfiltnu

       PRINT *, 'jfiltnv jfiltsv jfiltnu jfiltsu', jfiltnv, jfiltsv, jfiltnu, &
            jfiltsu
    END IF

    PRINT *, 'Modes premiers v '
    PRINT 334, modfrstv
    PRINT *, 'Modes premiers u '
    PRINT 334, modfrstu

    IF (nfilun < jfiltnu) THEN
       PRINT *, 'le parametre nfilun utilise pour la matrice ', &
            'matriceun est trop petit ! '
       PRINT *, 'Le changer dans parafilt.h et le mettre a ', jfiltnu
       PRINT *, 'Pour information, nfilun, nfilus, nfilvn, nfilvs ', &
            'doivent etre egaux successivement a ', jfiltnu, &
            jjm - jfiltsu + 1, jfiltnv, jjm - jfiltsv + 1
       STOP 1
    END IF
    IF (nfilun > jfiltnu+2) THEN
       PRINT *, 'le parametre nfilun utilise pour la matrice ', &
            'matriceun est trop grand ! Gachis de memoire ! '
       PRINT *, 'Le changer dans parafilt.h et le mettre a ', jfiltnu
       PRINT *, 'Pour information, nfilun, nfilus, nfilvn, nfilvs ', &
            'doivent etre egaux successivement a ', &
            jfiltnu, jjm - jfiltsu + 1, jfiltnv, jjm - jfiltsv + 1
    END IF
    IF (nfilus < jjm-jfiltsu+1) THEN
       PRINT *, 'le parametre nfilus utilise pour la matrice ', &
            'matriceus est trop petit ! '
       PRINT *, 'Le changer dans parafilt.h et le mettre a ', &
            jjm - jfiltsu + 1
       PRINT *, 'Pour information , nfilun, nfilus, nfilvn, nfilvs ', &
            'doivent etre egaux successivement a ', &
            jfiltnu, jjm - jfiltsu + 1, jfiltnv, jjm - jfiltsv + 1
       STOP 1
    END IF
    IF (nfilus > jjm-jfiltsu+3) THEN
       PRINT *, 'le parametre nfilus utilise pour la matrice ', &
            'matriceus est trop grand ! '
       PRINT *, 'Le changer dans parafilt.h et le mettre a ', &
            jjm - jfiltsu + 1
       PRINT *, 'Pour information , nfilun, nfilus, nfilvn, nfilvs ', &
            'doivent etre egaux successivement a ', &
            jfiltnu, jjm - jfiltsu + 1, jfiltnv, jjm - jfiltsv + 1
    END IF
    IF (nfilvn < jfiltnv) THEN
       PRINT *, 'le parametre nfilvn utilise pour la matrice ', &
            'matricevn est trop petit ! '
       PRINT *, 'Le changer dans parafilt.h et le mettre a ', jfiltnv
       PRINT *, 'Pour information , nfilun, nfilus, nfilvn, nfilvs ', &
            'doivent etre egaux successivement a ', &
            jfiltnu, jjm - jfiltsu + 1, jfiltnv, jjm - jfiltsv + 1
       STOP 1
    END IF
    IF (nfilvn > jfiltnv+2) THEN
       PRINT *, 'le parametre nfilvn utilise pour la matrice ', &
            'matricevn est trop grand ! Gachis de memoire ! '
       PRINT *, 'Le changer dans parafilt.h et le mettre a ', jfiltnv
       PRINT *, 'Pour information , nfilun, nfilus, nfilvn, nfilvs ', &
            'doivent etre egaux successivement a ', &
            jfiltnu, jjm - jfiltsu + 1, jfiltnv, jjm - jfiltsv + 1
    END IF
    IF (nfilvs < jjm-jfiltsv+1) THEN
       PRINT *, 'le parametre nfilvs utilise pour la matrice ', &
            'matricevs est trop petit ! Le changer dans parafilt.h '
       PRINT *, 'Le changer dans parafilt.h et le mettre a ', &
            jjm - jfiltsv + 1
       PRINT *, 'Pour information , nfilun, nfilus, nfilvn, nfilvs ', &
            'doivent etre egaux successivement a ', &
            jfiltnu, jjm - jfiltsu + 1, jfiltnv, jjm - jfiltsv + 1
       STOP 1
    END IF
    IF (nfilvs > jjm-jfiltsv+3) THEN
       PRINT *, 'le parametre nfilvs utilise pour la matrice ', &
            'matricevs est trop grand ! Gachis de memoire ! '
       PRINT *, 'Le changer dans parafilt.h et le mettre a ', &
            jjm - jfiltsv + 1
       PRINT *, 'Pour information , nfilun, nfilus, nfilvn, nfilvs ', &
            'doivent etre egaux successivement a ', &
            jfiltnu, jjm - jfiltsu + 1, jfiltnv, jjm - jfiltsv + 1
    END IF

    ! Calcul de la matrice filtre 'matriceu' pour les champs situes
    ! sur la grille scalaire

    DO j = 2, jfiltnu
       DO i = 1, iim
          IF (i < modfrstu(j)) then
             coff = 0.
          else
             coff = coefilu(i, j)
          end IF
          eignft(i, :) = eignfnv(:, i)*coff
       END DO
       matriceun(:, :, j) = matmul(eignfnv, eignft)
    END DO

    DO j = jfiltsu, jjm
       DO i = 1, iim
          IF (i < modfrstu(j)) then
             coff = 0.
          else
             coff = coefilu(i, j)
          end IF
          eignft(i, :) = eignfnv(:, i) * coff
       END DO
       matriceus(:, :, j - jfiltsu + 1) = matmul(eignfnv, eignft)
    END DO

    ! Calcul de la matrice filtre 'matricev' pour les champs situes
    ! sur la grille de V ou de Z

    DO j = 1, jfiltnv
       DO i = 1, iim
          IF (i < modfrstv(j)) then
             coff = 0.
          else
             coff = coefilv(i, j)
          end IF
          eignft(i, :) = eignfnu(:, i)*coff
       END DO
       matricevn(:, :, j) = matmul(eignfnu, eignft)
    END DO

    DO j = jfiltsv, jjm
       DO i = 1, iim
          IF (i < modfrstv(j)) then
             coff = 0.
          else
             coff = coefilv(i, j)
          end IF
          eignft(i, :) = eignfnu(:, i)*coff
       END DO
       matricevs(:, :, j-jfiltsv+1) = matmul(eignfnu, eignft)
    END DO

    ! Calcul de la matrice filtre 'matrinv' pour les champs situes
    ! sur la grille scalaire , pour le filtre inverse

    DO j = 2, jfiltnu
       DO i = 1, iim
          IF (i < modfrstu(j)) then
             coff = 0.
          else
             coff = coefilu(i, j)/(1.+coefilu(i, j))
          end IF
          eignft(i, :) = eignfnv(:, i)*coff
       END DO
       matrinvn(:, :, j) = matmul(eignfnv, eignft)
    END DO

    DO j = jfiltsu, jjm
       DO i = 1, iim
          IF (i < modfrstu(j)) then
             coff = 0.
          else
             coff = coefilu(i, j)/(1.+coefilu(i, j))
          end IF
          eignft(i, :) = eignfnv(:, i)*coff
       END DO
       matrinvs(:, :, j-jfiltsu+1) = matmul(eignfnv, eignft)
    END DO

334 FORMAT (1X, 24I3)

  END SUBROUTINE inifilr

end module inifilr_m
