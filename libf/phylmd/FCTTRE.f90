module FCTTRE

  ! From phylmd/FCTTRE.inc,v 1.2 2004/06/22 11:45:32

  !      This COMDECK includes the Thermodynamical functions for the cy39
  !       ECMWF Physics package.
  !       Consistent with YOMCST Basic physics constants, assuming the
  !       partial pressure of water vapour is given by a first order
  !       Taylor expansion of Qs(T) w.r.t. to Temperature, using constants
  !       in YOETHF

  implicit none

  LOGICAL, PARAMETER:: thermcep=.TRUE.

contains

  REAL function FOEEW ( PTARG,PDELARG )

    use yoethf, only: R3LES, R3IES, R4LES, R4IES
    use YOMCST, only: rtt

    REAL, intent(in):: PTARG, PDELARG

    !-----------------------

    FOEEW = EXP ((R3LES*(1.-PDELARG)+R3IES*PDELARG) * (PTARG-RTT) &
         / (PTARG-(R4LES*(1.-PDELARG)+R4IES*PDELARG)) )

  end function FOEEW

  !******************************************

  REAL function FOEDE(PTARG,PDELARG,P5ARG,PQSARG,PCOARG)

    use yoethf, only: R4LES, R4IES

    REAL, intent(in):: PTARG, PDELARG
    real, intent(in):: P5ARG, PQSARG, PCOARG

    !-----------------------

    FOEDE = PQSARG*PCOARG*P5ARG / (PTARG-(R4LES*(1.-PDELARG)+R4IES*PDELARG))**2

  end function FOEDE

  !******************************************

  REAL function qsats(ptarg)

    REAL, intent(in):: PTARG

    !-----------------------

    qsats = 100.0 * 0.622 * 10.0 &
         ** (2.07023 - 0.00320991 * ptarg &
         - 2484.896 / ptarg + 3.56654 * LOG10(ptarg))

  end function qsats

  !******************************************

  REAL function qsatl(ptarg)

    REAL, intent(in):: PTARG

    !-----------------------

    qsatl = 100.0 * 0.622 * 10.0 &
         ** (23.8319 - 2948.964 / ptarg &
         - 5.028 * LOG10(ptarg) &
         - 29810.16 * EXP( - 0.0699382 * ptarg) &
         + 25.21935 * EXP( - 2999.924 / ptarg))

  end function qsatl

  !******************************************

  REAL function dqsats(ptarg,pqsarg)

    use YOMCST, only: RLVTT, rcpd

    REAL, intent(in):: PTARG, pqsarg

    !-----------------------

    dqsats = RLVTT/RCPD*pqsarg * (3.56654/ptarg &
         +2484.896*LOG(10.)/ptarg**2 &
         -0.00320991*LOG(10.))

  end function dqsats

  !******************************************

  REAL function dqsatl(ptarg,pqsarg)

    use YOMCST, only: RLVTT, rcpd

    REAL, intent(in):: PTARG, pqsarg

    !-----------------------

    dqsatl = RLVTT/RCPD*pqsarg*LOG(10.)* &
         (2948.964/ptarg**2-5.028/LOG(10.)/ptarg &
         +25.21935*2999.924/ptarg**2*EXP(-2999.924/ptarg) &
         +29810.16*0.0699382*EXP(-0.0699382*ptarg))

  end function dqsatl

end module FCTTRE
