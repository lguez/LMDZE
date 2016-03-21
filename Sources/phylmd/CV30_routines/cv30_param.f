module cv30_param_m

  ! From LMDZ4/libf/phylmd/cvparam3.h, version 1.1.1.1, 2004/05/19 12:53:09

  ! Parameters for Emanuel convection scheme:
  ! - microphysical parameters
  ! - parameters that control the rate of approach to quasi-equilibrium

  USE dimphy, ONLY: klev

  implicit none

  integer minorig ! first level of convection

  integer, parameter:: nl = klev - 1
  ! Limit for convection. The maximum number of levels to which
  ! convection can penetrate, plus 1.  nl must be <= KLEV-1.

  real sigd ! FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT 
  real spfac ! FRACTION OF PRECIPITATION FALLING OUTSIDE OF CLOUD 

  real pbcrit
  !  CRITICAL CLOUD DEPTH (MB) BENEATH WHICH THE PRECIPITATION
  ! EFFICIENCY IS ASSUMED TO BE ZERO

  real ptcrit
  ! CLOUD DEPTH (MB) ABOVE WHICH THE PRECIPitation EFFICIENCY IS
  ! ASSUMED TO BE UNITY

  real omtrain
  real dtovsh, dpbase, dttrig

  real dtcrit
  ! CRITICAL BUOYANCY (K) USED TO ADJUST THE APPROACH TO
  ! QUASI-EQUILIBRIUM. IT MUST BE LESS THAN 0.

  real beta, alpha
  ! PARAMETERS THAT CONTROL THE RATE OF APPROACH TO QUASI-EQUILIBRIUM
  ! (THEIR STANDARD VALUES ARE 1.0 AND 0.96, RESPECTIVELY) (BETA MUST
  ! BE LESS THAN OR EQUAL TO 1)

  real tau ! CHARACTERISTIC TIMESCALE USED TO COMPUTE ALPHA and BETA
  real delta
  real betad

  private klev

contains

  SUBROUTINE cv30_param(delt)

    ! From LMDZ4/libf/phylmd/cv3_routines.F, version 1.5, 2005/07/11 15:20:02

    ! Set parameters for Emanuel convection scheme

    real, intent(in):: delt ! timestep (seconds)

    !------------------------------------------------------------

    ! Limit levels for convection:
    minorig = 1

    ! "Microphysical" parameters:

    sigd = 0.01
    spfac = 0.15
    pbcrit = 150.0
    ptcrit = 500.0
    ! cf. FH epmax = 0.993

    omtrain = 45.0 ! used also for snow (no distinction rain/snow)

    ! Misc:
    dtovsh = -0.2 ! dT for overshoot
    dpbase = -40. ! definition cloud base (400m above LCL)
    dttrig = 5. ! (loose) condition for triggering

    ! Rate of approach to quasi-equilibrium:
    dtcrit = -2.0
    tau = 8000.
    beta = 1.0 - delt/tau
    alpha = 1.5E-3 * delt/tau
    ! Increase alpha to compensate W decrease:
    alpha = alpha*1.5

    ! Interface cloud parameterization:
    delta=0.01 ! cld

    ! Interface with boundary-layer (gust factor): (sb)
    betad=10.0 ! original value (from convect 4.3)

  end SUBROUTINE cv30_param

end module cv30_param_m
