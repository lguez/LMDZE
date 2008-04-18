module regr_pr_comb_coefoz_m

  ! This module is clean: no C preprocessor directive, no include line.

  use dimens_m, only: llm
  use dimphy, only: klon

  implicit none

  ! The five module variables declared here are on the "physics" grid.
  ! The value of each variable for index "(i, k)" is at longitude
  ! "rlon(i)", latitude "rlat(i)" and middle of layer "k".

  real, save:: c_Mob(klon, llm)
  ! (sum of Mobidic terms in the net mass production rate of ozone
  ! by chemistry, per unit mass of air, in s-1)

  real, save:: a2(klon, llm)
  ! (derivative of mass production rate of ozone per unit mass of
  ! air with respect to ozone mass fraction, in s-1)

  real, save:: a4_mass(klon, llm)
  ! (derivative of mass production rate of ozone per unit mass of
  ! air with respect to temperature, in s-1 K-1)

  real, save:: a6_mass(klon, llm)
  ! (derivative of mass production rate of ozone per unit mass of
  ! air with respect to mass column-density of ozone above, in m2 s-1 kg-1)

  real, save:: r_het_interm(klon, llm)
  ! (net mass production rate by heterogeneous chemistry, per unit
  ! mass of ozone, corrected for chlorine content and latitude, but
  ! not for temperature and sun direction, in s-1)

  private klon, llm

contains

  subroutine regr_pr_comb_coefoz(julien)

    ! "regr_pr_comb_coefoz" stands for "regrid pressure combine
    ! coefficients ozone".

    ! This subroutine :
    ! -- reads from a file all eight coefficients for ozone chemistry,
    !    at the current day ;
    ! -- regrids the coefficients in pressure to the LMDZ vertical grid ;
    ! -- packs the coefficients to the "physics" horizontal grid ;
    ! -- combines the eight coefficients to define the five module variables.

    ! We assume that, in "coefoz_LMDZ.nc", the pressure levels are in hPa
    ! and strictly increasing.

    use netcdf95, only: nf95_open, nf95_close, nf95_get_coord
    use netcdf, only: nf90_nowrite
    use regr_pr_coefoz, only: regr_pr_av_coefoz, regr_pr_int_coefoz
    use phyetat0_m, only: rlat

    integer, intent(in):: julien ! jour julien, 1 <= julien <= 360

    ! Variables local to the procedure:
    integer ncid ! for NetCDF

    real, pointer:: plev(:)
    ! (pressure level of input data, converted to Pa, in strictly
    ! increasing order)

    integer n_plev ! number of pressure levels in the input data

    real, allocatable:: press_in_edg(:)
    ! (edges of pressure intervals for input data, in Pa, in strictly
    ! increasing order)

    real coefoz(klon, llm)
    ! (temporary storage for an ozone coefficient)
    ! (On the "physics" grid.
    ! "coefoz(i, k)" is at longitude "rlon(i)", latitude "rlat(i)",
    ! middle of layer "k".)

    real a6(klon, llm)
    ! (derivative of "P_net_Mob" with respect to column-density of ozone
    ! above, in cm2 s-1)
    ! (On the "physics" grid.
    ! "a6(i, k)" is at longitude "rlon(i)", latitude "rlat(i)",
    ! middle of layer "k".)

    real, parameter:: amu = 1.6605402e-27 ! atomic mass unit, in kg

    real, parameter:: Clx = 3.8e-9
    ! (total chlorine content in the upper stratosphere)

    integer k

    !------------------------------------

    print *, "Call sequence information: regr_pr_comb_coefoz"

    call nf95_open("coefoz_LMDZ.nc", nf90_nowrite, ncid)

    call nf95_get_coord(ncid, "plev", plev)
    ! Convert from hPa to Pa because "regr_pr_av" and "regr_pr_int" require so:
    plev = plev * 100.
    n_plev = size(plev)

    ! Compute edges of pressure intervals:
    allocate(press_in_edg(n_plev + 1))
    press_in_edg(1) = 0.
    ! We choose edges halfway in logarithm:
    forall (k = 2:n_plev) press_in_edg(k) = sqrt(plev(k - 1) * plev(k))
    press_in_edg(n_plev + 1) = huge(0.)
    ! (infinity, but any value guaranteed to be greater than the
    ! surface pressure would do)

    call regr_pr_av_coefoz(ncid, "a2", julien, press_in_edg, a2)

    call regr_pr_av_coefoz(ncid, "a4", julien, press_in_edg, a4_mass)
    a4_mass = a4_mass * 48. / 29.

    call regr_pr_av_coefoz(ncid, "a6", julien, press_in_edg, a6)

    ! Compute "a6_mass" avoiding underflow, do not divide by 1e4
    ! before dividing by molecular mass:
    a6_mass = a6 / (1e4 * 29. * amu)
    ! (factor 1e4: conversion from cm2 to m2)

    ! Combine coefficients to get "c_Mob":
    ! (We use as few local variables as possible, in order to spare
    ! main memory.)

    call regr_pr_av_coefoz(ncid, "P_net_Mob", julien, press_in_edg, c_Mob)

    call regr_pr_av_coefoz(ncid, "r_Mob", julien, press_in_edg, coefoz)
    c_mob = c_mob - a2 * coefoz

    call regr_pr_int_coefoz(ncid, "Sigma_Mob", julien, plev, top_value=0., &
         v3=coefoz)
    c_mob = (c_mob - a6 * coefoz) * 48. / 29.

    call regr_pr_av_coefoz(ncid, "temp_Mob", julien, press_in_edg, coefoz)
    c_mob = c_mob - a4_mass * coefoz

    call regr_pr_av_coefoz(ncid, "R_Het", julien, press_in_edg, r_het_interm)
    ! Heterogeneous chemistry is only at high latitudes:
    forall (k = 1: llm)
       where (abs(rlat) <= 45.) r_het_interm(:, k) = 0.
    end forall
    where (r_het_interm  /= 0.) r_het_interm = r_het_interm * (Clx / 3.8e-9)**2

    deallocate(plev) ! pointer
    call nf95_close(ncid)

  end subroutine regr_pr_comb_coefoz

end module regr_pr_comb_coefoz_m
