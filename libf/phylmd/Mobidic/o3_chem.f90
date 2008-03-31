module o3_chem_m

  ! This module is clean: no C preprocessor directive, no include line.

  IMPLICIT none

  private o3_prod

contains

  subroutine o3_chem(julien, gmtime, t_seri, zmasse, pdtphys, q)

    ! This procedure evolves the ozone mass fraction through a time
    ! step taking only chemistry into account.

    use nrutil, only: assert
    use dimphy, only: klon
    use dimens_m, only: llm
    use regr_pr_comb_coefoz_m, only: c_Mob, a4_mass, a2, r_het_interm
    use orbite_m, only: orbite, zenang
    use nrtype, only: pi

    integer, intent(in):: julien ! jour julien, 1 <= julien <= 360
    real, intent(in):: gmtime ! heure de la journée en fraction de jour
    real, intent(in):: t_seri(:, :) ! temperature,  in K

    real, intent(in):: zmasse(:, :)
    ! (column-density of mass of air in a cell, in kg m-2)
    ! (On the "physics" grid.
    ! "zmasse(i, k)" is at longitude "rlon(i)", latitude "rlat(i)", for
    ! layer "k".)

    real, intent(in):: pdtphys ! time step for physics, in s

    real, intent(inout):: q(:, :) ! mass fraction of ozone
    ! (On the "physics" grid.
    ! "q(i, k)" is at longitude "rlon(i)", latitude "rlat(i)", middle of
    ! layer "k".)

    ! Variables local to the procedure:
    integer k

    real c(klon, llm)
    ! (constant term during a time step in the net mass production
    ! rate of ozone by chemistry, per unit mass of air, in s-1)
    ! (On the "physics" grid.
    ! "c(i, k)" is at longitude "rlon(i)", latitude "rlat(i)", middle of
    ! layer "k".)

    real b(klon, llm)
    ! (coefficient of "q" in the net mass production
    ! rate of ozone by chemistry, per unit mass of air, in s-1)
    ! (On the "physics" grid.
    ! "b(i, k)" is at longitude "rlon(i)", latitude "rlat(i)", middle of
    ! layer "k".)

    real dq_o3_chem(klon, llm)
    ! (variation of ozone mass fraction due to chemistry during a time step)
    ! (On the "physics" grid.
    ! "dq_o3_chem(i, k)" is at longitude "rlon(i)", latitude
    ! "rlat(i)", middle of layer "k".)

    real earth_long
    ! (longitude vraie de la Terre dans son orbite solaire, par
    ! rapport au point vernal (21 mars), en degrés)

    real pmu0(klon) ! mean of cosine of solar zenith angle during "pdtphys"

    !-------------------------------------------------------------

    call assert(klon == (/size(q, 1), size(t_seri, 1), size(zmasse, 1)/), &
         "o3_chem klon")
    call assert(llm == (/size(q, 2), size(t_seri, 2), size(zmasse, 2)/), &
         "o3_chem llm")

    c = c_Mob + a4_mass * t_seri

    ! Compute coefficient "b":

    ! Heterogeneous chemistry is only at low temperature:
    where (t_seri < 195.)
       b = r_het_interm
    elsewhere
       b = 0.
    end where

    ! Heterogeneous chemistry is only during daytime:
    call orbite(real(julien), earth_long)
    call zenang(earth_long, gmtime, pdtphys, pmu0)
    forall (k = 1: llm)
       where (pmu0 <= cos(87. / 180. * pi)) b(:, k) = 0.
    end forall

    b = b + a2

    ! Midpoint method:

    ! Trial step to the midpoint:
    dq_o3_chem = o3_prod(q, zmasse, c, b) * pdtphys  / 2
    ! "Real" step across the whole interval:
    dq_o3_chem = o3_prod(q + dq_o3_chem, zmasse, c, b) * pdtphys
    q = q + dq_o3_chem

    ! Confine the mass fraction:
    q = min(max(q, 0.), .01)

  end subroutine o3_chem

  !*************************************************

  function o3_prod(q, zmasse, c, b)

    ! This function computes the production rate of ozone by chemistry.

    use regr_pr_comb_coefoz_m, only: a6_mass
    use nrutil, only: assert
    use dimens_m, only: llm
    use dimphy, only: klon

    real, intent(in):: q(:, :) ! mass fraction of ozone
    ! (On the "physics" grid.
    ! "q(i, k)" is at longitude "rlon(i)", latitude "rlat(i)", middle of
    ! layer "k".)

    real, intent(in):: zmasse(:, :)
    ! (column-density of mass of air in a layer, in kg m-2)
    ! (On the "physics" grid.
    ! "zmasse(i, k)" is at longitude "rlon(i)", latitude "rlat(i)", middle of
    ! layer "k".)

    real, intent(in):: c(:, :)
    ! (constant term during a time step in the net mass production
    ! rate of ozone by chemistry, per unit mass of air, in s-1)
    ! (On the "physics" grid.
    ! "c(i, k)" is at longitude "rlon(i)", latitude "rlat(i)", middle of
    ! layer "k".)

    real, intent(in):: b(:, :)
    ! (coefficient of "q" in the net mass production
    ! rate of ozone by chemistry, per unit mass of air, in s-1)
    ! (On the "physics" grid.
    ! "b(i, k)" is at longitude "rlon(i)", latitude "rlat(i)", middle of
    ! layer "k".)

    real o3_prod(klon, llm)
    ! (net mass production rate of ozone by chemistry, per unit mass
    ! of air, in s-1)
    ! (On the "physics" grid.
    ! "o3_prod(i, k)" is at longitude "rlon(i)", latitude "rlat(i)", middle of
    ! layer "k".)

    ! Variables local to the procedure:

    real sigma_mass(klon, llm)
    ! (mass column-density of ozone above point, in kg m-2)
    ! (On the "physics" grid.
    ! "sigma_mass(i, k)" is at longitude "rlon(i)", latitude
    ! "rlat(i)", middle of layer "k".)

    integer k

    !-------------------------------------------------------------------

    call assert(klon == (/size(q, 1), size(zmasse, 1), size(c, 1), &
         size(b, 1)/), "o3_prod 1")
    call assert(llm == (/size(q, 2), size(zmasse, 2), size(c, 2), &
         size(b, 2)/), "o3_prod 2")

    ! Compute the column-density above the base of layer
    ! "k", and, as a first approximation, take it as column-density
    ! above the middle of layer "k":
    sigma_mass(:, llm) = zmasse(:, llm) * q(:, llm) ! top layer
    do k =  llm - 1, 1, -1
       sigma_mass(:, k) = sigma_mass(:, k+1) + zmasse(:, k) * q(:, k)
    end do

    o3_prod = c + b * q + a6_mass * sigma_mass

  end function o3_prod

end module o3_chem_m
