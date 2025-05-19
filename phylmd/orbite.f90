MODULE orbite_m

  IMPLICIT NONE

CONTAINS

  pure SUBROUTINE orbite(xjour, longi, dist)

    ! From phylmd/orbite.F, version 1.1.1.1, 2004/05/19 12:53:08

    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 1993/08/18

    ! Pour un jour donn\'e, calcule la longitude vraie de la Terre et
    ! la distance Terre-Soleil.

    use jumble, only: deg_to_rad, rad_to_deg

    USE yomcst, ONLY: ecc, peri

    REAL, INTENT (IN):: xjour ! jour de l'ann\'ee \`a compter du premier janvier

    REAL, INTENT (OUT):: longi
    ! longitude vraie de la Terre dans son orbite solaire, par rapport
    ! au point vernal (21 mars), en degr\'es

    REAL, INTENT (OUT), OPTIONAL:: dist ! distance terre-soleil, en ua

    ! Local:
    REAL xl, xllp, xee, xse, ranm

    !----------------------------------------------------------------------

    xl = peri + 180.
    xllp = xl * deg_to_rad
    xee = ecc**2
    xse = sqrt(1. - xee)
    ranm = 2. * ((ecc / 2 + ecc * xee / 8.) * (1. + xse) * sin(xllp) &
         - xee / 4. * (0.5 + xse) * sin(2.*xllp) + ecc * xee / 8. &
         * (1. / 3. + xse) * sin(3. * xllp)) + (xjour - 81.) * deg_to_rad - xllp
    xee = xee * ecc
    longi = ranm + (2. * ecc - xee / 4.) * sin(ranm) &
         + 5. / 4. * ecc**2 * sin(2 * ranm) + 13. / 12. * xee * sin(3. * ranm)
    IF (present(dist)) dist = (1 - ecc**2) / (1 + ecc * cos(longi))
    longi = longi * rad_to_deg + xl

  END SUBROUTINE orbite

END MODULE orbite_m
