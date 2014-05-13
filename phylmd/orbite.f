MODULE orbite_m

  ! From phylmd/orbite.F, v 1.1.1.1 2004/05/19 12:53:08

  IMPLICIT NONE

CONTAINS

  SUBROUTINE orbite(xjour, longi, dist)

    USE yomcst, ONLY : r_ecc, r_peri
    use nr_util, only: pi

    ! Auteur(s): Z.X. Li (LMD/CNRS)
    ! Date: 1993/08/18
    ! Pour un jour donné, calcule la longitude vraie de la Terre (par
    ! rapport au point vernal, 21 mars) dans son orbite solaire. Calcule aussi
    ! la distance Terre-Soleil, c'est-à-dire l'unité astronomique.

    REAL, INTENT (IN):: xjour ! jour de l'année à compter du premier janvier

    REAL, INTENT (OUT):: longi
    ! longitude vraie de la Terre dans son orbite solaire, par
    ! rapport au point vernal (21 mars), en degrés

    REAL, INTENT (OUT), OPTIONAL:: dist
    ! distance terre-soleil (par rapport a la moyenne)

    ! Variables locales
    REAL pir, xl, xllp, xee, xse, xlam, anm, ranm, ranv

    !----------------------------------------------------------------------

    pir = pi / 180.
    xl = r_peri + 180.
    xllp = xl * pir
    xee = r_ecc * r_ecc
    xse = sqrt(1. - xee)
    xlam = (r_ecc / 2 + r_ecc * xee / 8.) * (1. + xse) * sin(xllp) &
         - xee / 4. * (0.5 + xse) * sin(2.*xllp) &
         + r_ecc * xee / 8. * (1. / 3. + xse) * sin(3.*xllp)
    xlam = 2. * xlam / pir
    anm = xlam + (xjour - 81.) - xl
    ranm = anm * pir
    xee = xee * r_ecc
    ranv = ranm + (2. * r_ecc - xee / 4.) * sin(ranm) + &
         5. / 4. * r_ecc * r_ecc * sin(2 * ranm) &
         + 13. / 12. * xee * sin(3.*ranm)

    longi = ranv / pir + xl

    IF (present(dist)) then
       dist = (1 - r_ecc*r_ecc) &
            / (1 + r_ecc*cos(pir*(longi - (r_peri + 180.))))
    end IF

  END SUBROUTINE orbite

END MODULE orbite_m
