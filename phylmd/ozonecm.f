module ozonecm_m

  IMPLICIT NONE

contains

  function ozonecm(rjour, paprs)

    ! From phylmd/ozonecm.F, version 1.3 2005/06/06 13:16:33

    ! The ozone climatology is based on an analytic formula which fits the
    ! Krueger and Mintzner (1976) profile, as well as the variations with
    ! altitude and latitude of the maximum ozone concentrations and the total
    ! column ozone concentration of Keating and Young (1986). The analytic
    ! formula have been established by J.-F. Royer (CRNM, Meteo France), who
    ! also provided us the code.

    ! A. J. Krueger and R. A. Minzner, A Mid-Latitude Ozone Model for the
    ! 1976 U. S. Standard Atmosphere, J. Geophys. Res., 81, 4477 (1976).

    ! Keating, G. M. and D. F. Young, 1985: Interim reference ozone
    ! models for the middle atmosphere, Handbook for MAP, vol. 16,
    ! 205-229.

    use nr_util, only: assert, twopi, deg_to_rad

    use dimensions, only: llm
    use phyetat0_m, only: rlat

    REAL, INTENT (IN) :: rjour

    REAL, INTENT (IN) :: paprs(:, :) ! (klon, llm + 1)
    ! pression pour chaque inter-couche, en Pa

    REAL ozonecm(size(paprs, 1), llm) ! (klon, llm)
    ! "ozonecm(j, k)" is the column-density of ozone in cell "(j, k)", that is
    ! between interface "k" and interface "k + 1", in kDU.

    ! Local:

    REAL tozon ! equivalent pressure of ozone above interface "k", in Pa
    INTEGER i, k

    REAL field(llm + 1)
    ! "field(k)" is the column-density of ozone between interface
    ! "k" and the top of the atmosphere (interface "llm + 1"), in kDU.

    real, PARAMETER:: ps = 101325.
    REAL, parameter:: an = 360., zo3q3 = 4E-8
    real, parameter:: zo3a3 = zo3q3 / ps / 2.
    REAL, parameter:: dobson_unit = 2.1415E-5 ! in kg m-2
    REAL gms, slat, slat2, sint, cost, ppm, a
    REAL asec, bsec, aprim

    !----------------------------------------------------------

    call assert(size(paprs, 2) == llm + 1, "ozonecm")

    sint = sin(twopi * (rjour + 15.) / an)
    cost = cos(twopi * (rjour + 15.) / an)
    field(llm + 1) = 0.

    DO i = 1, size(paprs, 1)
       slat = sin(deg_to_rad * rlat(i))
       slat2 = slat * slat
       gms = 0.0531 + sint * (- 0.001595 + 0.009443 * slat) + cost &
            * (- 0.001344 - 0.00346 * slat) + slat2 * (0.056222 + slat2 &
            * (- 0.037609 + 0.012248 * sint + 0.00521 * cost + 0.00889 &
            * slat)) - zo3q3 * ps
       ppm = 800. - 500. * slat2 - 150. * cost * slat
       bsec = 2650. + 5000. * slat2
       a = 4. * bsec**1.5 * ppm**1.5 * (1. + (bsec / ps)**1.5) &
            / (bsec**1.5 + ppm**1.5)**2
       aprim = max(0., (2.666666 * (1.74E-5 - 7.5E-6 * slat2 &
            - 1.7E-6 * cost * slat) * ppm - a * gms) / (1. - a))
       asec = max(0., (gms - aprim) * (1. + (bsec / ps)**1.5))
       aprim = gms - asec / (1. + (bsec / ps)**1.5)

       DO k = 1, llm
          tozon = aprim / (1. + 3. * (ppm / paprs(i, k))**2) &
               + asec / (1. + (bsec / paprs(i, k))**1.5) &
               + zo3a3 * paprs(i, k) * paprs(i, k)
          ! Convert from Pa to kDU:
          field(k) = tozon / 9.81 / dobson_unit / 1e3
       END DO

       forall (k = 1: llm) ozonecm(i, k) = field(k) - field(k + 1)
    END DO

    ozonecm = max(ozonecm, 1e-12)

  END function ozonecm

end module ozonecm_m
