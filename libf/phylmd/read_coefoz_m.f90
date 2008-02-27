module read_coefoz_m

  ! This module is clean: no C preprocessor directive, no include line.

  use dimens_m, only: llm
  use dimphy, only: klon

  implicit none

  real, save:: c_Mob(klon, llm, 12)
  ! (sum of Mobidic terms in the net mass production rate of ozone
  ! by chemistry, per unit mass of air, in s-1)
  ! (On the "physics" grid.
  ! Third dimension is the number of the month in the year.
  ! "c_Mob(i, k, month)" is at longitude "rlon(i)", latitude "rlat(i)",
  ! middle of layer "k".)

  real, save:: a2(klon, llm, 12)
  ! (derivative of mass production rate of ozone per unit mass of
  ! air with respect to ozone mass fraction, in s-1)
  ! (On the "physics" grid.
  ! Third dimension is the number of the month in the year.
  ! "a2(i, k, month)" is at longitude "rlon(i)", latitude "rlat(i)",
  ! middle of layer "k".)

  real, save:: a4_mass(klon, llm, 12)
  ! (derivative of mass production rate of ozone per unit mass of
  ! air with respect to temperature, in s-1 K-1)
  ! (On the "physics" grid.
  ! Third dimension is the number of the month in the year.
  ! "a4_mass(i, k, month)" is at longitude "rlon(i)", latitude "rlat(i)",
  ! middle of layer "k".)

  real, save:: a6_mass(klon, llm, 12)
  ! (derivative of mass production rate of ozone per unit mass of
  ! air with respect to mass column-density of ozone above, in m2 s-1 kg-1)
  ! (On the "physics" grid.
  ! Third dimension is the number of the month in the year.
  ! "a6_mass(i, k, month)" is at longitude "rlon(i)", latitude "rlat(i)",
  ! middle of layer "k".)

  real, save:: r_het_interm(klon, llm, 12)
  ! (net mass production rate by heterogeneous chemistry, per unit
  ! mass of ozone, corrected for chlorine content and latitude, but
  ! not for temperature and sun direction, in s-1)
  ! (On the "physics" grid.
  ! Third dimension is the number of the month in the year.
  ! "r_het_interm(i, k, month)" is at longitude "rlon(i)", latitude "rlat(i)",
  ! middle of layer "k".)

  private klon, llm

contains

  subroutine read_coefoz

    ! This subroutine reads from a file all eight parameters for ozone
    ! chemistry, at all months.
    ! The parameters are packed to the "physics" grid.
    ! Finally, the eight parameters are combined to define the five
    ! module variables.

    use netcdf95, only: nf95_open, nf95_close, nf90_nowrite
    use o3_Mob_ph_m, only: o3_Mob_ph
    use phyetat0_m, only: rlat

    ! Variables local to the procedure:
    integer ncid ! for NetCDF

    real a6(klon, llm, 12)
    ! (derivative of P_net_Mob with respect to column-density of ozone
    ! above, in cm2 s-1)
    ! (On the "physics" grid.
    ! Third dimension is the number of the month in the year.
    ! "a6(i, k, month)" is at longitude "rlon(i)", latitude "rlat(i)",
    ! middle of layer "k".)

    real, parameter:: amu = 1.6605402e-27 ! atomic mass unit, in kg

    real, parameter:: Clx = 3.8e-9
    ! (total chlorine content in the upper stratosphere)

    integer k, month

    !------------------------------------

    print *, "Call sequence information: read_coefoz"

    call nf95_open("coefoz_LMDZ.nc", nf90_nowrite, ncid)

    a2 = o3_Mob_ph(ncid, "a2")
    a4_mass = o3_Mob_ph(ncid, "a4") * 48. / 29.
    a6 = o3_Mob_ph(ncid, "a6")

    ! Compute "a6_mass" avoiding underflow, do not divide by 1e4
    ! before dividing by molecular mass:
    a6_mass = a6 / (1e4 * 29. * amu)
    ! (factor 1e4: conversion from cm2 to m2)

    c_Mob = 48. / 29. * (o3_Mob_ph(ncid, "P_net_Mob") &
         - a2 * o3_Mob_ph(ncid, "r_Mob") - a6 * o3_Mob_ph(ncid, "Sigma_Mob")) &
         - a4_mass * o3_Mob_ph(ncid, "temp_Mob") 

    r_het_interm = o3_Mob_ph(ncid, "R_Het") * (Clx / 3.8e-9)**2
    ! Heterogeneous chemistry is only at high latitudes:
    forall (k = 1: llm, month = 1: 12)
       where (abs(rlat) <= 45.) r_het_interm(:, k, month) = 0.
    end forall

    call nf95_close(ncid)

  end subroutine read_coefoz

end module read_coefoz_m
