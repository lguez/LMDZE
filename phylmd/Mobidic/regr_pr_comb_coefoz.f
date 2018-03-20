module regr_pr_comb_coefoz_m

  use dimensions, only: llm
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

  subroutine regr_pr_comb_coefoz(julien, paprs, pplay)

    ! "regr_pr_comb_coefoz" stands for "regrid pressure combine
    ! coefficients ozone".

    ! This subroutine :
    ! -- reads from a file all eight coefficients for ozone chemistry,
    !    at the current day ;
    ! -- regrids the coefficients in pressure to the LMDZ vertical grid ;
    ! -- packs the coefficients to the "physics" horizontal grid ;
    ! -- combines the eight coefficients to define the five module variables.

    use netcdf, only: nf90_nowrite
    use netcdf95, only: nf95_open, nf95_close
    use phyetat0_m, only: rlat
    use regr_pr_av_m, only: regr_pr_av
    use regr_pr_int_m, only: regr_pr_int

    integer, intent(in):: julien ! jour julien, 1 <= julien <= 360

    real, intent(in):: paprs(:, :) ! (klon, llm + 1)
    ! (pression pour chaque inter-couche, en Pa)

    real, intent(in):: pplay(:, :) ! (klon, llm)
    ! (pression pour le mileu de chaque couche, en Pa)

    ! Variables local to the procedure:

    integer ncid ! for NetCDF

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

    call regr_pr_av(ncid, "a2", julien, paprs, a2)

    call regr_pr_av(ncid, "a4", julien, paprs, a4_mass)
    a4_mass = a4_mass * 48. / 29.

    call regr_pr_av(ncid, "a6", julien, paprs, a6)

    ! Compute "a6_mass" avoiding underflow, do not divide by 1e4
    ! before dividing by molecular mass:
    a6_mass = a6 / (1e4 * 29. * amu)
    ! (factor 1e4: conversion from cm2 to m2)

    ! Combine coefficients to get "c_Mob":
    ! (We use as few local variables as possible, in order to spare
    ! main memory.)

    call regr_pr_av(ncid, "P_net_Mob", julien, paprs, c_Mob)

    call regr_pr_av(ncid, "r_Mob", julien, paprs, coefoz)
    c_mob = c_mob - a2 * coefoz

    call regr_pr_int(ncid, "Sigma_Mob", julien, pplay, top_value=0., v3=coefoz)
    c_mob = (c_mob - a6 * coefoz) * 48. / 29.

    call regr_pr_av(ncid, "temp_Mob", julien, paprs, coefoz)
    c_mob = c_mob - a4_mass * coefoz

    call regr_pr_av(ncid, "R_Het", julien, paprs, r_het_interm)
    ! Heterogeneous chemistry is only at high latitudes:
    forall (k = 1: llm)
       where (abs(rlat) <= 45.) r_het_interm(:, k) = 0.
    end forall
    r_het_interm = r_het_interm * (Clx / 3.8e-9)**2

    call nf95_close(ncid)

  end subroutine regr_pr_comb_coefoz

end module regr_pr_comb_coefoz_m
