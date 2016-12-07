MODULE orbite_m

  IMPLICIT NONE

CONTAINS

  pure SUBROUTINE orbite(xjour, longi, dist)

    ! From phylmd/orbite.F, version 1.1.1.1 2004/05/19 12:53:08

    ! Author: Z.X. Li (LMD/CNRS)
    ! Date: 1993/08/18

    ! Pour un jour donn\'e, calcule la longitude vraie de la Terre et la
    ! distance Terre-Soleil, c'est-\`a-dire l'unit\'e astronomique.

    use nr_util, only: pi
    USE yomcst, ONLY: r_ecc, r_peri

    REAL, INTENT (IN):: xjour ! jour de l'ann\'ee \`a compter du premier janvier

    REAL, INTENT (OUT):: longi
    ! longitude vraie de la Terre dans son orbite solaire, par rapport
    ! au point vernal (21 mars), en degr\'es

    REAL, INTENT (OUT), OPTIONAL:: dist
    ! distance terre-soleil (par rapport \`a la moyenne)

    ! Local:
    REAL pir, xl, xllp, xee, xse, ranm

    !----------------------------------------------------------------------

    pir = pi / 180.
    xl = r_peri + 180.
    xllp = xl * pir
    xee = r_ecc**2
    xse = sqrt(1. - xee)
    ranm = 2. * ((r_ecc / 2 + r_ecc * xee / 8.) * (1. + xse) * sin(xllp) &
         - xee / 4. * (0.5 + xse) * sin(2.*xllp) + r_ecc * xee / 8. &
         * (1. / 3. + xse) * sin(3. * xllp)) + (xjour - 81.) * pir - xllp
    xee = xee * r_ecc
    longi = (ranm + (2. * r_ecc - xee / 4.) * sin(ranm) + 5. / 4. * r_ecc**2 &
         * sin(2 * ranm) + 13. / 12. * xee * sin(3. * ranm)) / pir + xl

    IF (present(dist)) dist = (1 - r_ecc * r_ecc) &
         / (1 + r_ecc * cos(pir * (longi - (r_peri + 180.))))

  END SUBROUTINE orbite

END MODULE orbite_m
