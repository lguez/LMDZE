SUBROUTINE TLIFT(P,T,RR,RS,GZ,PLCL,ICB,NK, TVP,TPK,CLW,ND,NL, DTVPDT1,DTVPDQ1)

  ! From phylmd/tlift.F, v 1.1.1.1 2004/05/19 12:53:08

  ! Argument NK ajoute (jyg) = Niveau de depart de la convection

  use YOMCST, only: rcpd, RCPV, rcw, rcs, rv, rd, rlvtt, RLMLT

  implicit none

  integer, PARAMETER:: NA=60
  integer nd, icb, nk, nl
  real plcl
  REAL GZ(ND),TPK(ND),CLW(ND)
  REAL T(ND),RR(ND),RS(ND),TVP(ND),P(ND)

  REAL DTVPDT1(ND),DTVPDQ1(ND)
  ! Derivatives of parcel virtual temperature with regard to T1 and Q1

  REAL QI(NA)

  REAL DTPDT1(NA),DTPDQ1(NA)
  ! Derivatives of parcel temperature with regard to T1 and Q1

  LOGICAL ICE_CONV
  real gravity, cpd, cpv, cl, ci, CPVMCL, CLMCI, EPS, alv0, alf0, CPP, cpinv
  real ah0, alv, alf, tg, s, ahg, tc, denom, es, esi, qsat_new, snew
  integer icb1, i, IMIN, j

  !--------------------------------------------------------------

  ! ***   ASSIGN VALUES OF THERMODYNAMIC CONSTANTS     ***
  ! on utilise les constantes thermo du Centre Europeen: (SB)

  CPD = RCPD
  CPV = RCPV
  CL = RCW
  CI = RCS
  CPVMCL = CL-CPV
  CLMCI = CL-CI
  EPS = RD/RV
  ALV0 = RLVTT
  ALF0 = RLMLT ! (ALF0 = RLSTT-RLVTT)

  ! ***  CALCULATE CERTAIN PARCEL QUANTITIES, INCLUDING STATIC ENERGY   ***

  ICB1=MIN(ICB,NL)

  CPP=CPD*(1.-RR(NK))+RR(NK)*CPV
  CPINV=1./CPP

  !         ICB may be below condensation level
  DO I=1,ICB1
     CLW(I)=0.0
  end DO

  DO I=NK,ICB1
     TPK(I)=T(NK)-(GZ(I) - GZ(NK))*CPINV
     TVP(I)=TPK(I)*(1.+RR(NK)/EPS-RR(NK))
     DTVPDT1(I) = 1.+RR(NK)/EPS-RR(NK)
     DTVPDQ1(I) = TPK(I)*(1./EPS-1.)
  end DO

  !    ***  FIND LIFTED PARCEL TEMPERATURE AND MIXING RATIO    ***

  AH0=(CPD*(1.-RR(NK))+CL*RR(NK))*T(NK) &
       +RR(NK)*(ALV0-CPVMCL*(T(NK)-273.15)) + GZ(NK)
  IMIN = ICB1
  !         If ICB is below LCL, start loop at ICB+1
  IF (PLCL .LT. P(ICB1)) IMIN = MIN(IMIN+1,NL)

  DO  I=IMIN,NL
     ALV=ALV0-CPVMCL*(T(I)-273.15)
     ALF=ALF0+CLMCI*(T(I)-273.15)

     GRAVITY=RS(I)
     TG=T(I)
     S=CPD*(1.-RR(NK))+CL*RR(NK)+ALV*ALV*GRAVITY/(RV*T(I)*T(I))

     S=1./S

     DO  J=1,2
        AHG=CPD*TG+(CL-CPD)*RR(NK)*TG+ALV*GRAVITY+GZ(I)

        TG=TG+S*(AH0-AHG)
        TC=TG-273.15
        DENOM=243.5+TC
        DENOM=MAX(DENOM,1.0)

        !       FORMULE DE BOLTON POUR PSAT
        ES=6.112*EXP(17.67*TC/DENOM)
        GRAVITY=EPS*ES/(P(I)-ES*(1.-EPS))
     end DO

     TPK(I)=(AH0-GZ(I)-ALV*GRAVITY)/(CPD+(CL-CPD)*RR(NK))
     CLW(I)=RR(NK)-GRAVITY
     CLW(I)=MAX(0.0,CLW(I))
     TVP(I)=TPK(I)*(1.+GRAVITY/EPS-RR(NK))

     !jyg1       Derivatives
     DTPDT1(I) = CPD*S
     DTPDQ1(I) = ALV*S

     DTVPDT1(I) = DTPDT1(I)*(1. + GRAVITY/EPS - &
          RR(NK) + ALV*GRAVITY/(RD*TPK(I)) )
     DTVPDQ1(I) = DTPDQ1(I)*(1. + GRAVITY/EPS - &
          RR(NK) + ALV*GRAVITY/(RD*TPK(I)) ) - TPK(I)
  end DO

  ICE_CONV = .FALSE.

  IF (ICE_CONV) THEN
     DO  I=ICB1,NL
        IF (T(I).LT.263.15) THEN
           TG=TPK(I)
           TC=TPK(I)-273.15
           DENOM=243.5+TC
           ES=6.112*EXP(17.67*TC/DENOM)
           ALV=ALV0-CPVMCL*(T(I)-273.15)
           ALF=ALF0+CLMCI*(T(I)-273.15)

           DO J=1,4
              ESI=EXP(23.33086-(6111.72784/TPK(I))+0.15215*LOG(TPK(I)))
              QSAT_NEW=EPS*ESI/(P(I)-ESI*(1.-EPS))

              SNEW= CPD*(1.-RR(NK))+CL*RR(NK) &
                   +ALV*ALV*QSAT_NEW/(RV*TPK(I)*TPK(I))

              SNEW=1./SNEW
              TPK(I)=TG+(ALF*QI(I)+ALV*GRAVITY*(1.-(ESI/ES)))*SNEW
           ENDDO
           CLW(I)=RR(NK)-QSAT_NEW
           CLW(I)=MAX(0.0,CLW(I))
           TVP(I)=TPK(I)*(1.+QSAT_NEW/EPS-RR(NK))
        ENDIF
     end DO
  ENDIF

  !* BK :  RAJOUT DE LA TEMPERATURE DES ASCENDANCES
  !*   NON DILUES AU  NIVEAU KLEV = ND
  !*   POSONS LE ENVIRON EGAL A CELUI DE KLEV-1
  TPK(NL+1)=TPK(NL)

END SUBROUTINE TLIFT
