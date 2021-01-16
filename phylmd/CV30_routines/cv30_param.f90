module cv30_param_m

  ! From LMDZ4/libf/phylmd/cvparam3.h, version 1.1.1.1, 2004/05/19 12:53:09

  ! Set parameters for Emanuel convection scheme. Includes
  ! microphysical parameters and parameters that control the rate of
  ! approach to quasi-equilibrium.

  implicit none

  integer, parameter:: minorig = 1 
  ! first level of convection (limit levels for convection)

  integer, protected:: nl
  ! Limit for convection. The maximum number of levels to which
  ! convection can penetrate, plus 1. We should have:
  ! 6 <= nl <= KLEV - 1
  ! (because of locate in cv30_feed).

  real, parameter:: sigd = 0.01 
  ! fractional area covered by unsaturated downdraft
  ! \sigma_d in Emanuel (1991 928)

  real, parameter:: dtcrit = - 2.
  ! Critical buoyancy (K) used to adjust the approach to
  ! quasi-equilibrium. It must be < 0.

  real, protected:: beta, alpha
  ! Parameters that control the rate of approach to quasi-equilibrium
  ! (their standard values are 1. and 0.96, respectively) (beta must
  ! be less than or equal to 1).

contains

  SUBROUTINE cv30_param

    ! From LMDZ4/libf/phylmd/cv3_routines.F, version 1.5, 2005/07/11 15:20:02
    ! SET CONSTANTS AND PARAMETERS

    use conf_gcm_m, only: dtphys
    USE dimphy, ONLY: klev

    ! Local:

    real, parameter:: tau = 8000.
    ! characteristic timescale used to compute alpha and beta

    !------------------------------------------------------------

    nl = klev - 1
    beta = 1. - dtphys / tau

    alpha = 1.5E-3 * dtphys / tau * 1.5
    ! factor 1.5: increase alpha to compensate W decrease

  end SUBROUTINE cv30_param

end module cv30_param_m
