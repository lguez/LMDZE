!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/cv3_routines.F,v 1.5 2005/07/11 15:20:02 lmdzadmin Exp $
!
!
!
      SUBROUTINE cv3_param(nd,delt)
      use conema3_m
            use cvparam3
      implicit none

!------------------------------------------------------------
! Set parameters for convectL for iflag_con = 3
!------------------------------------------------------------

!
!   ***  PBCRIT IS THE CRITICAL CLOUD DEPTH (MB) BENEATH WHICH THE ***
!   ***      PRECIPITATION EFFICIENCY IS ASSUMED TO BE ZERO     ***
!   ***  PTCRIT IS THE CLOUD DEPTH (MB) ABOVE WHICH THE PRECIP. ***
!   ***            EFFICIENCY IS ASSUMED TO BE UNITY            ***
!   ***  SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT  ***
!   ***  SPFAC IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE ***
!   ***                        OF CLOUD                         ***
!
! [TAU: CHARACTERISTIC TIMESCALE USED TO COMPUTE ALPHA & BETA]
!   ***    ALPHA AND BETA ARE PARAMETERS THAT CONTROL THE RATE OF ***
!   ***                 APPROACH TO QUASI-EQUILIBRIUM           ***
!   ***    (THEIR STANDARD VALUES ARE 1.0 AND 0.96, RESPECTIVELY) ***
!   ***           (BETA MUST BE LESS THAN OR EQUAL TO 1)        ***
!
!   ***    DTCRIT IS THE CRITICAL BUOYANCY (K) USED TO ADJUST THE ***
!   ***                 APPROACH TO QUASI-EQUILIBRIUM           ***
!   ***                     IT MUST BE LESS THAN 0              ***


      integer nd
      real, intent(in):: delt ! timestep (seconds)

! noff: integer limit for convection (nd-noff)
! minorig: First level of convection

! -- limit levels for convection:

      noff    = 1
      minorig = 1
      nl=nd-noff
      nlp=nl+1
      nlm=nl-1

! -- "microphysical" parameters:

      sigd   = 0.01
      spfac  = 0.15
      pbcrit = 150.0
      ptcrit = 500.0
!IM cf. FH     epmax  = 0.993

      omtrain = 45.0 ! used also for snow (no disctinction rain/snow)

! -- misc:

      dtovsh = -0.2 ! dT for overshoot
      dpbase = -40. ! definition cloud base (400m above LCL)
      dttrig = 5.   ! (loose) condition for triggering

! -- rate of approach to quasi-equilibrium:

      dtcrit = -2.0
      tau    = 8000.
      beta   = 1.0 - delt/tau
      alpha  = 1.5E-3 * delt/tau
! increase alpha to compensate W decrease:
      alpha  = alpha*1.5

! -- interface cloud parameterization:

      delta=0.01  ! cld

! -- interface with boundary-layer (gust factor): (sb)

      betad=10.0   ! original value (from convect 4.3)

      return
      end
