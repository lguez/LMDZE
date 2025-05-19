module o3_chem_m

  IMPLICIT none

  private o3_prod

contains

  subroutine o3_chem(julien, gmtime, t_seri, zmasse, pdtphys, q)

    ! This procedure evolves the ozone mass fraction through a time
    ! step taking only chemistry into account.

    ! All the 2-dimensional arrays are on the "physics" grid.
    ! Their shape is "(/klon, llm/)".
    ! Index "(i, :)" is for longitude "rlon(i)", latitude "rlat(i)".

    use jumble, only: assert, pi
    use dimphy, only: klon
    use dimensions, only: llm
    use regr_pr_comb_coefoz_m, only: c_Mob, a4_mass, a2, r_het_interm
    use orbite_m, only: orbite
    use zenang_m, only: zenang

    integer, intent(in):: julien ! jour julien, 1 <= julien <= 360
    real, intent(in):: gmtime ! heure de la journée en fraction de jour
    real, intent(in):: t_seri(:, :) ! (klon, llm) temperature, in K

    real, intent(in):: zmasse(:, :) ! (klon, llm)
    ! (column-density of mass of air in a cell, in kg m-2)
    ! "zmasse(:, k)" is for layer "k".)

    real, intent(in):: pdtphys ! time step for physics, in s

    real, intent(inout):: q(:, :) ! (klon, llm) mass fraction of ozone
    ! "q(:, k)" is at middle of layer "k".)

    ! Variables local to the procedure:
    integer k

    real c(klon, llm)
    ! (constant term during a time step in the net mass production
    ! rate of ozone by chemistry, per unit mass of air, in s-1)
    ! "c(:, k)" is at middle of layer "k".)

    real b(klon, llm)
    ! (coefficient of "q" in the net mass production
    ! rate of ozone by chemistry, per unit mass of air, in s-1)
    ! "b(:, k)" is at middle of layer "k".)

    real dq_o3_chem(klon, llm)
    ! (variation of ozone mass fraction due to chemistry during a time step)
    ! "dq_o3_chem(:, k)" is at middle of layer "k".)

    real earth_long
    ! (longitude vraie de la Terre dans son orbite solaire, par
    ! rapport au point vernal (21 mars), en degrés)

    real mu0(klon) ! mean of cosine of solar zenith angle during "pdtphys"

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
    call zenang(earth_long, gmtime, pdtphys, mu0)
    forall (k = 1: llm)
       where (mu0 <= cos(87. / 180. * pi)) b(:, k) = 0.
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

    ! All the 2-dimensional arrays are on the "physics" grid.
    ! Their shape is "(/klon, llm/)".
    ! Index "(i, :)" is for longitude "rlon(i)", latitude "rlat(i)".

    use jumble, only: assert

    use regr_pr_comb_coefoz_m, only: a6_mass
    use dimensions, only: llm
    use dimphy, only: klon

    real, intent(in):: q(:, :) ! mass fraction of ozone
    ! "q(:, k)" is at middle of layer "k".)

    real, intent(in):: zmasse(:, :)
    ! (column-density of mass of air in a layer, in kg m-2)
    ! ("zmasse(:, k)" is for layer "k".)

    real, intent(in):: c(:, :)
    ! (constant term during a time step in the net mass production
    ! rate of ozone by chemistry, per unit mass of air, in s-1)
    ! "c(:, k)" is at middle of layer "k".)

    real, intent(in):: b(:, :)
    ! (coefficient of "q" in the net mass production rate of ozone by
    ! chemistry, per unit mass of air, in s-1)
    ! ("b(:, k)" is at middle of layer "k".)

    real o3_prod(klon, llm)
    ! (net mass production rate of ozone by chemistry, per unit mass
    ! of air, in s-1)
    ! ("o3_prod(:, k)" is at middle of layer "k".)

    ! Variables local to the procedure:

    real sigma_mass(klon, llm)
    ! (mass column-density of ozone above point, in kg m-2)
    ! ("sigma_mass(:, k)" is at middle of layer "k".)

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
