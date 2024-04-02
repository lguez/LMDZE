module lwb_m

  IMPLICIT none

contains

  SUBROUTINE LWB(PDT0,PTAVE,PTL &
       , PB,PBINT,PBSUIN,PBSUR,PBTOP,PDBSL &
       , PGA,PGB,PGASUR,PGBSUR,PGATOP,PGBTOP)
    use dimensions
    use dimphy
    use raddimlw
    !
    !-----------------------------------------------------------------------
    !     PURPOSE.
    !     --------
    !           COMPUTES PLANCK FUNCTIONS
    !
    !        EXPLICIT ARGUMENTS :
    !        --------------------
    !     ==== INPUTS ===
    ! PDT0   : (KLON)             ; SURFACE TEMPERATURE DISCONTINUITY
    ! PTAVE  : (KLON,LLM)       ; TEMPERATURE
    ! PTL    : (KLON,0:LLM)     ; HALF LEVEL TEMPERATURE
    !     ==== OUTPUTS ===
    ! PB     : (KLON,Ninter,LLM+1); SPECTRAL HALF LEVEL PLANCK FUNCTION
    ! PBINT  : (KLON,LLM+1)     ; HALF LEVEL PLANCK FUNCTION
    ! PBSUIN : (KLON)             ; SURFACE PLANCK FUNCTION
    ! PBSUR  : (KLON,Ninter)        ; SURFACE SPECTRAL PLANCK FUNCTION
    ! PBTOP  : (KLON,Ninter)        ; TOP SPECTRAL PLANCK FUNCTION
    ! PDBSL  : (KLON,Ninter,LLM*2); SUB-LAYER PLANCK FUNCTION GRADIENT
    ! PGA    : (KLON,8,2,LLM); dB/dT-weighted LAYER PADE APPROXIMANTS
    ! PGB    : (KLON,8,2,LLM); dB/dT-weighted LAYER PADE APPROXIMANTS
    ! PGASUR, PGBSUR (KLON,8,2)   ; SURFACE PADE APPROXIMANTS
    ! PGATOP, PGBTOP (KLON,8,2)   ; T.O.A. PADE APPROXIMANTS
    !
    !        IMPLICIT ARGUMENTS :   NONE
    !        --------------------
    !
    !     METHOD.
    !     -------
    !
    !          1. COMPUTES THE PLANCK FUNCTION ON ALL LEVELS AND HALF LEVELS
    !     FROM A POLYNOMIAL DEVELOPMENT OF PLANCK FUNCTION
    !
    !     REFERENCE.
    !     ----------
    !
    !        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
    !        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS           "
    !
    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*
    !
    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 89-07-14
    !
    !-----------------------------------------------------------------------
    !
    ! ARGUMENTS:
    !
    DOUBLE PRECISION PDT0(KLON)
    DOUBLE PRECISION PTAVE(KLON,LLM)
    DOUBLE PRECISION PTL(KLON,LLM+1)
    !
    DOUBLE PRECISION PB(KLON,Ninter,LLM+1) ! SPECTRAL HALF LEVEL PLANCK FUNCTION
    DOUBLE PRECISION PBINT(KLON,LLM+1) ! HALF LEVEL PLANCK FUNCTION
    DOUBLE PRECISION PBSUIN(KLON) ! SURFACE PLANCK FUNCTION
    DOUBLE PRECISION PBSUR(KLON,Ninter) ! SURFACE SPECTRAL PLANCK FUNCTION
    DOUBLE PRECISION PBTOP(KLON,Ninter) ! TOP SPECTRAL PLANCK FUNCTION
    DOUBLE PRECISION PDBSL(KLON,Ninter,LLM*2) ! SUB-LAYER PLANCK FUNCTION GRADIENT
    DOUBLE PRECISION PGA(KLON,8,2,LLM) ! dB/dT-weighted LAYER PADE APPROXIMANTS
    DOUBLE PRECISION PGB(KLON,8,2,LLM) ! dB/dT-weighted LAYER PADE APPROXIMANTS
    DOUBLE PRECISION PGASUR(KLON,8,2) ! SURFACE PADE APPROXIMANTS
    DOUBLE PRECISION PGBSUR(KLON,8,2) ! SURFACE PADE APPROXIMANTS
    DOUBLE PRECISION PGATOP(KLON,8,2) ! T.O.A. PADE APPROXIMANTS
    DOUBLE PRECISION PGBTOP(KLON,8,2) ! T.O.A. PADE APPROXIMANTS
    !
    !-------------------------------------------------------------------------
    !*  LOCAL VARIABLES:
    INTEGER INDB(KLON),INDS(KLON)
    DOUBLE PRECISION ZBLAY(KLON,LLM),ZBLEV(KLON,LLM+1)
    DOUBLE PRECISION ZRES(KLON),ZRES2(KLON),ZTI(KLON),ZTI2(KLON)
    !
    INTEGER jk, jl, ic, jnu, jf, jg
    INTEGER jk1, jk2
    INTEGER k, j, ixtox, indto, ixtx, indt
    INTEGER indsu, indtp
    DOUBLE PRECISION zdsto1, zdstox, zdst1, zdstx
    !
    !* Quelques parametres:
    DOUBLE PRECISION TSTAND
    PARAMETER (TSTAND=250.0)
    DOUBLE PRECISION TSTP
    PARAMETER (TSTP=12.5)
    INTEGER MXIXT
    PARAMETER (MXIXT=10)
    !
    !* Used Data Block:
    DOUBLE PRECISION TINTP(11)
    SAVE TINTP
    DOUBLE PRECISION GA(11,16,3), GB(11,16,3)
    SAVE GA, GB
    DOUBLE PRECISION XP(6,6)
    SAVE XP
    !
    DATA TINTP / 187.5d0, 200.d0, 212.5d0, 225.d0, 237.5d0, 250.d0, &
         262.5d0, 275.d0, 287.5d0, 300.d0, 312.5d0 /
    !-----------------------------------------------------------------------
    !-- WATER VAPOR -- INT.1 -- 0- 500 CM-1 -- FROM ABS225 ----------------
    !
    !
    !
    !
    !-- R.D. -- G = - 0.2 SLA
    !
    !
    !----- INTERVAL = 1 ----- T =  187.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA( 1, 1,IC),IC=1,3) / &
         0.63499072D-02,-0.99506586D-03, 0.00000000D+00/
    DATA (GB( 1, 1,IC),IC=1,3) / &
         0.63499072D-02, 0.97222852D-01, 0.10000000D+01/
    DATA (GA( 1, 2,IC),IC=1,3) / &
         0.77266491D-02,-0.11661515D-02, 0.00000000D+00/
    DATA (GB( 1, 2,IC),IC=1,3) / &
         0.77266491D-02, 0.10681591D+00, 0.10000000D+01/
    !
    !----- INTERVAL = 1 ----- T =  200.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA( 2, 1,IC),IC=1,3) / &
         0.65566348D-02,-0.10184169D-02, 0.00000000D+00/
    DATA (GB( 2, 1,IC),IC=1,3) / &
         0.65566348D-02, 0.98862238D-01, 0.10000000D+01/
    DATA (GA( 2, 2,IC),IC=1,3) / &
         0.81323287D-02,-0.11886130D-02, 0.00000000D+00/
    DATA (GB( 2, 2,IC),IC=1,3) / &
         0.81323287D-02, 0.10921298D+00, 0.10000000D+01/
    !
    !----- INTERVAL = 1 ----- T =  212.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA( 3, 1,IC),IC=1,3) / &
         0.67849730D-02,-0.10404730D-02, 0.00000000D+00/
    DATA (GB( 3, 1,IC),IC=1,3) / &
         0.67849730D-02, 0.10061504D+00, 0.10000000D+01/
    DATA (GA( 3, 2,IC),IC=1,3) / &
         0.86507620D-02,-0.12139929D-02, 0.00000000D+00/
    DATA (GB( 3, 2,IC),IC=1,3) / &
         0.86507620D-02, 0.11198225D+00, 0.10000000D+01/
    !
    !----- INTERVAL = 1 ----- T =  225.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA( 4, 1,IC),IC=1,3) / &
         0.70481947D-02,-0.10621792D-02, 0.00000000D+00/
    DATA (GB( 4, 1,IC),IC=1,3) / &
         0.70481947D-02, 0.10256222D+00, 0.10000000D+01/
    DATA (GA( 4, 2,IC),IC=1,3) / &
         0.92776391D-02,-0.12445811D-02, 0.00000000D+00/
    DATA (GB( 4, 2,IC),IC=1,3) / &
         0.92776391D-02, 0.11487826D+00, 0.10000000D+01/
    !
    !----- INTERVAL = 1 ----- T =  237.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA( 5, 1,IC),IC=1,3) / &
         0.73585943D-02,-0.10847662D-02, 0.00000000D+00/
    DATA (GB( 5, 1,IC),IC=1,3) / &
         0.73585943D-02, 0.10475952D+00, 0.10000000D+01/
    DATA (GA( 5, 2,IC),IC=1,3) / &
         0.99806312D-02,-0.12807672D-02, 0.00000000D+00/
    DATA (GB( 5, 2,IC),IC=1,3) / &
         0.99806312D-02, 0.11751113D+00, 0.10000000D+01/
    !
    !----- INTERVAL = 1 ----- T =  250.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA( 6, 1,IC),IC=1,3) / &
         0.77242818D-02,-0.11094726D-02, 0.00000000D+00/
    DATA (GB( 6, 1,IC),IC=1,3) / &
         0.77242818D-02, 0.10720986D+00, 0.10000000D+01/
    DATA (GA( 6, 2,IC),IC=1,3) / &
         0.10709803D-01,-0.13208251D-02, 0.00000000D+00/
    DATA (GB( 6, 2,IC),IC=1,3) / &
         0.10709803D-01, 0.11951535D+00, 0.10000000D+01/
    !
    !----- INTERVAL = 1 ----- T =  262.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA( 7, 1,IC),IC=1,3) / &
         0.81472693D-02,-0.11372949D-02, 0.00000000D+00/
    DATA (GB( 7, 1,IC),IC=1,3) / &
         0.81472693D-02, 0.10985370D+00, 0.10000000D+01/
    DATA (GA( 7, 2,IC),IC=1,3) / &
         0.11414739D-01,-0.13619034D-02, 0.00000000D+00/
    DATA (GB( 7, 2,IC),IC=1,3) / &
         0.11414739D-01, 0.12069945D+00, 0.10000000D+01/
    !
    !----- INTERVAL = 1 ----- T =  275.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA( 8, 1,IC),IC=1,3) / &
         0.86227527D-02,-0.11687683D-02, 0.00000000D+00/
    DATA (GB( 8, 1,IC),IC=1,3) / &
         0.86227527D-02, 0.11257633D+00, 0.10000000D+01/
    DATA (GA( 8, 2,IC),IC=1,3) / &
         0.12058772D-01,-0.14014165D-02, 0.00000000D+00/
    DATA (GB( 8, 2,IC),IC=1,3) / &
         0.12058772D-01, 0.12108524D+00, 0.10000000D+01/
    !
    !----- INTERVAL = 1 ----- T =  287.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA( 9, 1,IC),IC=1,3) / &
         0.91396814D-02,-0.12038314D-02, 0.00000000D+00/
    DATA (GB( 9, 1,IC),IC=1,3) / &
         0.91396814D-02, 0.11522980D+00, 0.10000000D+01/
    DATA (GA( 9, 2,IC),IC=1,3) / &
         0.12623992D-01,-0.14378639D-02, 0.00000000D+00/
    DATA (GB( 9, 2,IC),IC=1,3) / &
         0.12623992D-01, 0.12084229D+00, 0.10000000D+01/
    !
    !----- INTERVAL = 1 ----- T =  300.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA(10, 1,IC),IC=1,3) / &
         0.96825438D-02,-0.12418367D-02, 0.00000000D+00/
    DATA (GB(10, 1,IC),IC=1,3) / &
         0.96825438D-02, 0.11766343D+00, 0.10000000D+01/
    DATA (GA(10, 2,IC),IC=1,3) / &
         0.13108146D-01,-0.14708488D-02, 0.00000000D+00/
    DATA (GB(10, 2,IC),IC=1,3) / &
         0.13108146D-01, 0.12019005D+00, 0.10000000D+01/
    !
    !----- INTERVAL = 1 ----- T =  312.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA(11, 1,IC),IC=1,3) / &
         0.10233955D-01,-0.12817135D-02, 0.00000000D+00/
    DATA (GB(11, 1,IC),IC=1,3) / &
         0.10233955D-01, 0.11975320D+00, 0.10000000D+01/
    DATA (GA(11, 2,IC),IC=1,3) / &
         0.13518390D-01,-0.15006791D-02, 0.00000000D+00/
    DATA (GB(11, 2,IC),IC=1,3) / &
         0.13518390D-01, 0.11932684D+00, 0.10000000D+01/
    !
    !
    !
    !--- WATER VAPOR --- INTERVAL 2 -- 500-800 CM-1--- FROM ABS225 ---------
    !
    !
    !
    !
    !--- R.D.  ---  G = 0.02 + 0.50 / ( 1 + 4.5 U )
    !
    !
    !----- INTERVAL = 2 ----- T =  187.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 1, 3,IC),IC=1,3) / &
         0.11644593d+01, 0.41243390d+00, 0.00000000d+00/
    DATA (GB( 1, 3,IC),IC=1,3) / &
         0.11644593d+01, 0.10346097d+01, 0.10000000d+01/
    DATA (GA( 1, 4,IC),IC=1,3) / &
         0.12006968d+01, 0.48318936d+00, 0.00000000d+00/
    DATA (GB( 1, 4,IC),IC=1,3) / &
         0.12006968d+01, 0.10626130d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  200.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 2, 3,IC),IC=1,3) / &
         0.11747203d+01, 0.43407282d+00, 0.00000000d+00/
    DATA (GB( 2, 3,IC),IC=1,3) / &
         0.11747203d+01, 0.10433655d+01, 0.10000000d+01/
    DATA (GA( 2, 4,IC),IC=1,3) / &
         0.12108196d+01, 0.50501827d+00, 0.00000000d+00/
    DATA (GB( 2, 4,IC),IC=1,3) / &
         0.12108196d+01, 0.10716026d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  212.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 3, 3,IC),IC=1,3) / &
         0.11837872d+01, 0.45331413d+00, 0.00000000d+00/
    DATA (GB( 3, 3,IC),IC=1,3) / &
         0.11837872d+01, 0.10511933d+01, 0.10000000d+01/
    DATA (GA( 3, 4,IC),IC=1,3) / &
         0.12196717d+01, 0.52409502d+00, 0.00000000d+00/
    DATA (GB( 3, 4,IC),IC=1,3) / &
         0.12196717d+01, 0.10795108d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  225.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 4, 3,IC),IC=1,3) / &
         0.11918561d+01, 0.47048604d+00, 0.00000000d+00/
    DATA (GB( 4, 3,IC),IC=1,3) / &
         0.11918561d+01, 0.10582150d+01, 0.10000000d+01/
    DATA (GA( 4, 4,IC),IC=1,3) / &
         0.12274493d+01, 0.54085277d+00, 0.00000000d+00/
    DATA (GB( 4, 4,IC),IC=1,3) / &
         0.12274493d+01, 0.10865006d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  237.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 5, 3,IC),IC=1,3) / &
         0.11990757d+01, 0.48586286d+00, 0.00000000d+00/
    DATA (GB( 5, 3,IC),IC=1,3) / &
         0.11990757d+01, 0.10645317d+01, 0.10000000d+01/
    DATA (GA( 5, 4,IC),IC=1,3) / &
         0.12343189d+01, 0.55565422d+00, 0.00000000d+00/
    DATA (GB( 5, 4,IC),IC=1,3) / &
         0.12343189d+01, 0.10927103d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  250.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 6, 3,IC),IC=1,3) / &
         0.12055643d+01, 0.49968044d+00, 0.00000000d+00/
    DATA (GB( 6, 3,IC),IC=1,3) / &
         0.12055643d+01, 0.10702313d+01, 0.10000000d+01/
    DATA (GA( 6, 4,IC),IC=1,3) / &
         0.12404147d+01, 0.56878618d+00, 0.00000000d+00/
    DATA (GB( 6, 4,IC),IC=1,3) / &
         0.12404147d+01, 0.10982489d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  262.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 7, 3,IC),IC=1,3) / &
         0.12114186d+01, 0.51214132d+00, 0.00000000d+00/
    DATA (GB( 7, 3,IC),IC=1,3) / &
         0.12114186d+01, 0.10753907d+01, 0.10000000d+01/
    DATA (GA( 7, 4,IC),IC=1,3) / &
         0.12458431d+01, 0.58047395d+00, 0.00000000d+00/
    DATA (GB( 7, 4,IC),IC=1,3) / &
         0.12458431d+01, 0.11032019d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  275.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 8, 3,IC),IC=1,3) / &
         0.12167192d+01, 0.52341830d+00, 0.00000000d+00/
    DATA (GB( 8, 3,IC),IC=1,3) / &
         0.12167192d+01, 0.10800762d+01, 0.10000000d+01/
    DATA (GA( 8, 4,IC),IC=1,3) / &
         0.12506907d+01, 0.59089894d+00, 0.00000000d+00/
    DATA (GB( 8, 4,IC),IC=1,3) / &
         0.12506907d+01, 0.11076379d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  287.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 9, 3,IC),IC=1,3) / &
         0.12215344d+01, 0.53365803d+00, 0.00000000d+00/
    DATA (GB( 9, 3,IC),IC=1,3) / &
         0.12215344d+01, 0.10843446d+01, 0.10000000d+01/
    DATA (GA( 9, 4,IC),IC=1,3) / &
         0.12550299d+01, 0.60021475d+00, 0.00000000d+00/
    DATA (GB( 9, 4,IC),IC=1,3) / &
         0.12550299d+01, 0.11116160d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  300.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA(10, 3,IC),IC=1,3) / &
         0.12259226d+01, 0.54298448d+00, 0.00000000d+00/
    DATA (GB(10, 3,IC),IC=1,3) / &
         0.12259226d+01, 0.10882439d+01, 0.10000000d+01/
    DATA (GA(10, 4,IC),IC=1,3) / &
         0.12589256d+01, 0.60856112d+00, 0.00000000d+00/
    DATA (GB(10, 4,IC),IC=1,3) / &
         0.12589256d+01, 0.11151910d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  312.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA(11, 3,IC),IC=1,3) / &
         0.12299344d+01, 0.55150227d+00, 0.00000000d+00/
    DATA (GB(11, 3,IC),IC=1,3) / &
         0.12299344d+01, 0.10918144d+01, 0.10000000d+01/
    DATA (GA(11, 4,IC),IC=1,3) / &
         0.12624402d+01, 0.61607594d+00, 0.00000000d+00/
    DATA (GB(11, 4,IC),IC=1,3) / &
         0.12624402d+01, 0.11184188d+01, 0.10000000d+01/
    !
    !
    !
    !
    !
    !
    !- WATER VAPOR - INT. 3 -- 800-970 + 1110-1250 CM-1 -- FIT FROM 215 IS -
    !
    !
    !-- WATER VAPOR LINES IN THE WINDOW REGION (800-1250 CM-1)
    !
    !
    !
    !--- G = 3.875d-03 ---------------
    !
    !----- INTERVAL = 3 ----- T =  187.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 1, 7,IC),IC=1,3) / &
         0.10192131d+02, 0.80737799d+01, 0.00000000d+00/
    DATA (GB( 1, 7,IC),IC=1,3) / &
         0.10192131d+02, 0.82623280d+01, 0.10000000d+01/
    DATA (GA( 1, 8,IC),IC=1,3) / &
         0.92439050d+01, 0.77425778d+01, 0.00000000d+00/
    DATA (GB( 1, 8,IC),IC=1,3) / &
         0.92439050d+01, 0.79342219d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 3 ----- T =  200.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 2, 7,IC),IC=1,3) / &
         0.97258602d+01, 0.79171158d+01, 0.00000000d+00/
    DATA (GB( 2, 7,IC),IC=1,3) / &
         0.97258602d+01, 0.81072291d+01, 0.10000000d+01/
    DATA (GA( 2, 8,IC),IC=1,3) / &
         0.87567422d+01, 0.75443460d+01, 0.00000000d+00/
    DATA (GB( 2, 8,IC),IC=1,3) / &
         0.87567422d+01, 0.77373458d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 3 ----- T =  212.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 3, 7,IC),IC=1,3) / &
         0.92992890d+01, 0.77609605d+01, 0.00000000d+00/
    DATA (GB( 3, 7,IC),IC=1,3) / &
         0.92992890d+01, 0.79523834d+01, 0.10000000d+01/
    DATA (GA( 3, 8,IC),IC=1,3) / &
         0.83270144d+01, 0.73526151d+01, 0.00000000d+00/
    DATA (GB( 3, 8,IC),IC=1,3) / &
         0.83270144d+01, 0.75467334d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 3 ----- T =  225.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 4, 7,IC),IC=1,3) / &
         0.89154021d+01, 0.76087371d+01, 0.00000000d+00/
    DATA (GB( 4, 7,IC),IC=1,3) / &
         0.89154021d+01, 0.78012527d+01, 0.10000000d+01/
    DATA (GA( 4, 8,IC),IC=1,3) / &
         0.79528337d+01, 0.71711188d+01, 0.00000000d+00/
    DATA (GB( 4, 8,IC),IC=1,3) / &
         0.79528337d+01, 0.73661786d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 3 ----- T =  237.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 5, 7,IC),IC=1,3) / &
         0.85730084d+01, 0.74627112d+01, 0.00000000d+00/
    DATA (GB( 5, 7,IC),IC=1,3) / &
         0.85730084d+01, 0.76561458d+01, 0.10000000d+01/
    DATA (GA( 5, 8,IC),IC=1,3) / &
         0.76286839d+01, 0.70015571d+01, 0.00000000d+00/
    DATA (GB( 5, 8,IC),IC=1,3) / &
         0.76286839d+01, 0.71974319d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 3 ----- T =  250.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 6, 7,IC),IC=1,3) / &
         0.82685838d+01, 0.73239981d+01, 0.00000000d+00/
    DATA (GB( 6, 7,IC),IC=1,3) / &
         0.82685838d+01, 0.75182174d+01, 0.10000000d+01/
    DATA (GA( 6, 8,IC),IC=1,3) / &
         0.73477879d+01, 0.68442532d+01, 0.00000000d+00/
    DATA (GB( 6, 8,IC),IC=1,3) / &
         0.73477879d+01, 0.70408543d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 3 ----- T =  262.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 7, 7,IC),IC=1,3) / &
         0.79978921d+01, 0.71929934d+01, 0.00000000d+00/
    DATA (GB( 7, 7,IC),IC=1,3) / &
         0.79978921d+01, 0.73878952d+01, 0.10000000d+01/
    DATA (GA( 7, 8,IC),IC=1,3) / &
         0.71035818d+01, 0.66987996d+01, 0.00000000d+00/
    DATA (GB( 7, 8,IC),IC=1,3) / &
         0.71035818d+01, 0.68960649d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 3 ----- T =  275.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 8, 7,IC),IC=1,3) / &
         0.77568055d+01, 0.70697065d+01, 0.00000000d+00/
    DATA (GB( 8, 7,IC),IC=1,3) / &
         0.77568055d+01, 0.72652133d+01, 0.10000000d+01/
    DATA (GA( 8, 8,IC),IC=1,3) / &
         0.68903312d+01, 0.65644820d+01, 0.00000000d+00/
    DATA (GB( 8, 8,IC),IC=1,3) / &
         0.68903312d+01, 0.67623672d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 3 ----- T =  287.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 9, 7,IC),IC=1,3) / &
         0.75416266d+01, 0.69539626d+01, 0.00000000d+00/
    DATA (GB( 9, 7,IC),IC=1,3) / &
         0.75416266d+01, 0.71500151d+01, 0.10000000d+01/
    DATA (GA( 9, 8,IC),IC=1,3) / &
         0.67032875d+01, 0.64405267d+01, 0.00000000d+00/
    DATA (GB( 9, 8,IC),IC=1,3) / &
         0.67032875d+01, 0.66389989d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 3 ----- T =  300.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA(10, 7,IC),IC=1,3) / &
         0.73491694d+01, 0.68455144d+01, 0.00000000d+00/
    DATA (GB(10, 7,IC),IC=1,3) / &
         0.73491694d+01, 0.70420667d+01, 0.10000000d+01/
    DATA (GA(10, 8,IC),IC=1,3) / &
         0.65386461d+01, 0.63262376d+01, 0.00000000d+00/
    DATA (GB(10, 8,IC),IC=1,3) / &
         0.65386461d+01, 0.65252707d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 3 ----- T =  312.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA(11, 7,IC),IC=1,3) / &
         0.71767400d+01, 0.67441020d+01, 0.00000000d+00/
    DATA (GB(11, 7,IC),IC=1,3) / &
         0.71767400d+01, 0.69411177d+01, 0.10000000d+01/
    DATA (GA(11, 8,IC),IC=1,3) / &
         0.63934377d+01, 0.62210701d+01, 0.00000000d+00/
    DATA (GB(11, 8,IC),IC=1,3) / &
         0.63934377d+01, 0.64206412d+01, 0.10000000d+01/
    !
    !
    !-- WATER VAPOR -- 970-1110 CM-1 ----------------------------------------
    !
    !-- G = 3.6d-03
    !
    !----- INTERVAL = 4 ----- T =  187.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 1, 9,IC),IC=1,3) / &
         0.24870635d+02, 0.10542131d+02, 0.00000000d+00/
    DATA (GB( 1, 9,IC),IC=1,3) / &
         0.24870635d+02, 0.10656640d+02, 0.10000000d+01/
    DATA (GA( 1,10,IC),IC=1,3) / &
         0.24586283d+02, 0.10490353d+02, 0.00000000d+00/
    DATA (GB( 1,10,IC),IC=1,3) / &
         0.24586283d+02, 0.10605856d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  200.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 2, 9,IC),IC=1,3) / &
         0.24725591d+02, 0.10515895d+02, 0.00000000d+00/
    DATA (GB( 2, 9,IC),IC=1,3) / &
         0.24725591d+02, 0.10630910d+02, 0.10000000d+01/
    DATA (GA( 2,10,IC),IC=1,3) / &
         0.24441465d+02, 0.10463512d+02, 0.00000000d+00/
    DATA (GB( 2,10,IC),IC=1,3) / &
         0.24441465d+02, 0.10579514d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  212.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 3, 9,IC),IC=1,3) / &
         0.24600320d+02, 0.10492949d+02, 0.00000000d+00/
    DATA (GB( 3, 9,IC),IC=1,3) / &
         0.24600320d+02, 0.10608399d+02, 0.10000000d+01/
    DATA (GA( 3,10,IC),IC=1,3) / &
         0.24311657d+02, 0.10439183d+02, 0.00000000d+00/
    DATA (GB( 3,10,IC),IC=1,3) / &
         0.24311657d+02, 0.10555632d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  225.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 4, 9,IC),IC=1,3) / &
         0.24487300d+02, 0.10472049d+02, 0.00000000d+00/
    DATA (GB( 4, 9,IC),IC=1,3) / &
         0.24487300d+02, 0.10587891d+02, 0.10000000d+01/
    DATA (GA( 4,10,IC),IC=1,3) / &
         0.24196167d+02, 0.10417324d+02, 0.00000000d+00/
    DATA (GB( 4,10,IC),IC=1,3) / &
         0.24196167d+02, 0.10534169d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  237.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 5, 9,IC),IC=1,3) / &
         0.24384935d+02, 0.10452961d+02, 0.00000000d+00/
    DATA (GB( 5, 9,IC),IC=1,3) / &
         0.24384935d+02, 0.10569156d+02, 0.10000000d+01/
    DATA (GA( 5,10,IC),IC=1,3) / &
         0.24093406d+02, 0.10397704d+02, 0.00000000d+00/
    DATA (GB( 5,10,IC),IC=1,3) / &
         0.24093406d+02, 0.10514900d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  250.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 6, 9,IC),IC=1,3) / &
         0.24292341d+02, 0.10435562d+02, 0.00000000d+00/
    DATA (GB( 6, 9,IC),IC=1,3) / &
         0.24292341d+02, 0.10552075d+02, 0.10000000d+01/
    DATA (GA( 6,10,IC),IC=1,3) / &
         0.24001597d+02, 0.10380038d+02, 0.00000000d+00/
    DATA (GB( 6,10,IC),IC=1,3) / &
         0.24001597d+02, 0.10497547d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  262.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 7, 9,IC),IC=1,3) / &
         0.24208572d+02, 0.10419710d+02, 0.00000000d+00/
    DATA (GB( 7, 9,IC),IC=1,3) / &
         0.24208572d+02, 0.10536510d+02, 0.10000000d+01/
    DATA (GA( 7,10,IC),IC=1,3) / &
         0.23919098d+02, 0.10364052d+02, 0.00000000d+00/
    DATA (GB( 7,10,IC),IC=1,3) / &
         0.23919098d+02, 0.10481842d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  275.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 8, 9,IC),IC=1,3) / &
         0.24132642d+02, 0.10405247d+02, 0.00000000d+00/
    DATA (GB( 8, 9,IC),IC=1,3) / &
         0.24132642d+02, 0.10522307d+02, 0.10000000d+01/
    DATA (GA( 8,10,IC),IC=1,3) / &
         0.23844511d+02, 0.10349509d+02, 0.00000000d+00/
    DATA (GB( 8,10,IC),IC=1,3) / &
         0.23844511d+02, 0.10467553d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  287.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA( 9, 9,IC),IC=1,3) / &
         0.24063614d+02, 0.10392022d+02, 0.00000000d+00/
    DATA (GB( 9, 9,IC),IC=1,3) / &
         0.24063614d+02, 0.10509317d+02, 0.10000000d+01/
    DATA (GA( 9,10,IC),IC=1,3) / &
         0.23776708d+02, 0.10336215d+02, 0.00000000d+00/
    DATA (GB( 9,10,IC),IC=1,3) / &
         0.23776708d+02, 0.10454488d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  300.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA(10, 9,IC),IC=1,3) / &
         0.24000649d+02, 0.10379892d+02, 0.00000000d+00/
    DATA (GB(10, 9,IC),IC=1,3) / &
         0.24000649d+02, 0.10497402d+02, 0.10000000d+01/
    DATA (GA(10,10,IC),IC=1,3) / &
         0.23714816d+02, 0.10324018d+02, 0.00000000d+00/
    DATA (GB(10,10,IC),IC=1,3) / &
         0.23714816d+02, 0.10442501d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  312.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   28   37   45
    DATA (GA(11, 9,IC),IC=1,3) / &
         0.23943021d+02, 0.10368736d+02, 0.00000000d+00/
    DATA (GB(11, 9,IC),IC=1,3) / &
         0.23943021d+02, 0.10486443d+02, 0.10000000d+01/
    DATA (GA(11,10,IC),IC=1,3) / &
         0.23658197d+02, 0.10312808d+02, 0.00000000d+00/
    DATA (GB(11,10,IC),IC=1,3) / &
         0.23658197d+02, 0.10431483d+02, 0.10000000d+01/
    !
    !
    !
    !-- H2O -- WEAKER PARTS OF THE STRONG BANDS  -- FROM ABS225 ----
    !
    !-- WATER VAPOR --- 350 - 500 CM-1
    !
    !-- G = - 0.2*SLA, 0.0 +0.5/(1+0.5U)
    !
    !----- INTERVAL = 5 ----- T =  187.5
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA( 1, 5,IC),IC=1,3) / &
         0.15750172d+00,-0.22159303d-01, 0.00000000d+00/
    DATA (GB( 1, 5,IC),IC=1,3) / &
         0.15750172d+00, 0.38103212d+00, 0.10000000d+01/
    DATA (GA( 1, 6,IC),IC=1,3) / &
         0.17770551d+00,-0.24972399d-01, 0.00000000d+00/
    DATA (GB( 1, 6,IC),IC=1,3) / &
         0.17770551d+00, 0.41646579d+00, 0.10000000d+01/
    !
    !----- INTERVAL = 5 ----- T =  200.0
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA( 2, 5,IC),IC=1,3) / &
         0.16174076d+00,-0.22748917d-01, 0.00000000d+00/
    DATA (GB( 2, 5,IC),IC=1,3) / &
         0.16174076d+00, 0.38913800d+00, 0.10000000d+01/
    DATA (GA( 2, 6,IC),IC=1,3) / &
         0.18176757d+00,-0.25537247d-01, 0.00000000d+00/
    DATA (GB( 2, 6,IC),IC=1,3) / &
         0.18176757d+00, 0.42345095d+00, 0.10000000d+01/
    !
    !----- INTERVAL = 5 ----- T =  212.5
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA( 3, 5,IC),IC=1,3) / &
         0.16548628d+00,-0.23269898d-01, 0.00000000d+00/
    DATA (GB( 3, 5,IC),IC=1,3) / &
         0.16548628d+00, 0.39613651d+00, 0.10000000d+01/
    DATA (GA( 3, 6,IC),IC=1,3) / &
         0.18527967d+00,-0.26025624d-01, 0.00000000d+00/
    DATA (GB( 3, 6,IC),IC=1,3) / &
         0.18527967d+00, 0.42937476d+00, 0.10000000d+01/
    !
    !----- INTERVAL = 5 ----- T =  225.0
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA( 4, 5,IC),IC=1,3) / &
         0.16881124d+00,-0.23732392d-01, 0.00000000d+00/
    DATA (GB( 4, 5,IC),IC=1,3) / &
         0.16881124d+00, 0.40222421d+00, 0.10000000d+01/
    DATA (GA( 4, 6,IC),IC=1,3) / &
         0.18833348d+00,-0.26450280d-01, 0.00000000d+00/
    DATA (GB( 4, 6,IC),IC=1,3) / &
         0.18833348d+00, 0.43444062d+00, 0.10000000d+01/
    !
    !----- INTERVAL = 5 ----- T =  237.5
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA( 5, 5,IC),IC=1,3) / &
         0.17177839d+00,-0.24145123d-01, 0.00000000d+00/
    DATA (GB( 5, 5,IC),IC=1,3) / &
         0.17177839d+00, 0.40756010d+00, 0.10000000d+01/
    DATA (GA( 5, 6,IC),IC=1,3) / &
         0.19100108d+00,-0.26821236d-01, 0.00000000d+00/
    DATA (GB( 5, 6,IC),IC=1,3) / &
         0.19100108d+00, 0.43880316d+00, 0.10000000d+01/
    !
    !----- INTERVAL = 5 ----- T =  250.0
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA( 6, 5,IC),IC=1,3) / &
         0.17443933d+00,-0.24515269d-01, 0.00000000d+00/
    DATA (GB( 6, 5,IC),IC=1,3) / &
         0.17443933d+00, 0.41226954d+00, 0.10000000d+01/
    DATA (GA( 6, 6,IC),IC=1,3) / &
         0.19334122d+00,-0.27146657d-01, 0.00000000d+00/
    DATA (GB( 6, 6,IC),IC=1,3) / &
         0.19334122d+00, 0.44258354d+00, 0.10000000d+01/
    !
    !----- INTERVAL = 5 ----- T =  262.5
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA( 7, 5,IC),IC=1,3) / &
         0.17683622d+00,-0.24848690d-01, 0.00000000d+00/
    DATA (GB( 7, 5,IC),IC=1,3) / &
         0.17683622d+00, 0.41645142d+00, 0.10000000d+01/
    DATA (GA( 7, 6,IC),IC=1,3) / &
         0.19540288d+00,-0.27433354d-01, 0.00000000d+00/
    DATA (GB( 7, 6,IC),IC=1,3) / &
         0.19540288d+00, 0.44587882d+00, 0.10000000d+01/
    !
    !----- INTERVAL = 5 ----- T =  275.0
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA( 8, 5,IC),IC=1,3) / &
         0.17900375d+00,-0.25150210d-01, 0.00000000d+00/
    DATA (GB( 8, 5,IC),IC=1,3) / &
         0.17900375d+00, 0.42018474d+00, 0.10000000d+01/
    DATA (GA( 8, 6,IC),IC=1,3) / &
         0.19722732d+00,-0.27687065d-01, 0.00000000d+00/
    DATA (GB( 8, 6,IC),IC=1,3) / &
         0.19722732d+00, 0.44876776d+00, 0.10000000d+01/
    !
    !----- INTERVAL = 5 ----- T =  287.5
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA( 9, 5,IC),IC=1,3) / &
         0.18097099d+00,-0.25423873d-01, 0.00000000d+00/
    DATA (GB( 9, 5,IC),IC=1,3) / &
         0.18097099d+00, 0.42353379d+00, 0.10000000d+01/
    DATA (GA( 9, 6,IC),IC=1,3) / &
         0.19884918d+00,-0.27912608d-01, 0.00000000d+00/
    DATA (GB( 9, 6,IC),IC=1,3) / &
         0.19884918d+00, 0.45131451d+00, 0.10000000d+01/
    !
    !----- INTERVAL = 5 ----- T =  300.0
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA(10, 5,IC),IC=1,3) / &
         0.18276283d+00,-0.25673139d-01, 0.00000000d+00/
    DATA (GB(10, 5,IC),IC=1,3) / &
         0.18276283d+00, 0.42655211d+00, 0.10000000d+01/
    DATA (GA(10, 6,IC),IC=1,3) / &
         0.20029696d+00,-0.28113944d-01, 0.00000000d+00/
    DATA (GB(10, 6,IC),IC=1,3) / &
         0.20029696d+00, 0.45357095d+00, 0.10000000d+01/
    !
    !----- INTERVAL = 5 ----- T =  312.5
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA(11, 5,IC),IC=1,3) / &
         0.18440117d+00,-0.25901055d-01, 0.00000000d+00/
    DATA (GB(11, 5,IC),IC=1,3) / &
         0.18440117d+00, 0.42928533d+00, 0.10000000d+01/
    DATA (GA(11, 6,IC),IC=1,3) / &
         0.20159300d+00,-0.28294180d-01, 0.00000000d+00/
    DATA (GB(11, 6,IC),IC=1,3) / &
         0.20159300d+00, 0.45557797d+00, 0.10000000d+01/
    !
    !
    !
    !
    !- WATER VAPOR - WINGS OF VIBRATION-ROTATION BAND - 1250-1450+1880-2820 -
    !--- G = 0.0
    !
    !
    !----- INTERVAL = 6 ----- T =  187.5
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA( 1,11,IC),IC=1,3) / &
         0.11990218d+02,-0.12823142d+01, 0.00000000d+00/
    DATA (GB( 1,11,IC),IC=1,3) / &
         0.11990218d+02, 0.26681588d+02, 0.10000000d+01/
    DATA (GA( 1,12,IC),IC=1,3) / &
         0.79709806d+01,-0.74805226d+00, 0.00000000d+00/
    DATA (GB( 1,12,IC),IC=1,3) / &
         0.79709806d+01, 0.18377807d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 6 ----- T =  200.0
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA( 2,11,IC),IC=1,3) / &
         0.10904073d+02,-0.10571588d+01, 0.00000000d+00/
    DATA (GB( 2,11,IC),IC=1,3) / &
         0.10904073d+02, 0.24728346d+02, 0.10000000d+01/
    DATA (GA( 2,12,IC),IC=1,3) / &
         0.75400737d+01,-0.56252739d+00, 0.00000000d+00/
    DATA (GB( 2,12,IC),IC=1,3) / &
         0.75400737d+01, 0.17643148d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 6 ----- T =  212.5
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA( 3,11,IC),IC=1,3) / &
         0.89126838d+01,-0.74864953d+00, 0.00000000d+00/
    DATA (GB( 3,11,IC),IC=1,3) / &
         0.89126838d+01, 0.20551342d+02, 0.10000000d+01/
    DATA (GA( 3,12,IC),IC=1,3) / &
         0.81804377d+01,-0.46188072d+00, 0.00000000d+00/
    DATA (GB( 3,12,IC),IC=1,3) / &
         0.81804377d+01, 0.19296161d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 6 ----- T =  225.0
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA( 4,11,IC),IC=1,3) / &
         0.85622405d+01,-0.58705980d+00, 0.00000000d+00/
    DATA (GB( 4,11,IC),IC=1,3) / &
         0.85622405d+01, 0.19955244d+02, 0.10000000d+01/
    DATA (GA( 4,12,IC),IC=1,3) / &
         0.10564339d+02,-0.40712065d+00, 0.00000000d+00/
    DATA (GB( 4,12,IC),IC=1,3) / &
         0.10564339d+02, 0.24951120d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 6 ----- T =  237.5
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA( 5,11,IC),IC=1,3) / &
         0.94892164d+01,-0.49305772d+00, 0.00000000d+00/
    DATA (GB( 5,11,IC),IC=1,3) / &
         0.94892164d+01, 0.22227100d+02, 0.10000000d+01/
    DATA (GA( 5,12,IC),IC=1,3) / &
         0.46896789d+02,-0.15295996d+01, 0.00000000d+00/
    DATA (GB( 5,12,IC),IC=1,3) / &
         0.46896789d+02, 0.10957372d+03, 0.10000000d+01/
    !
    !----- INTERVAL = 6 ----- T =  250.0
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA( 6,11,IC),IC=1,3) / &
         0.13580937d+02,-0.51461431d+00, 0.00000000d+00/
    DATA (GB( 6,11,IC),IC=1,3) / &
         0.13580937d+02, 0.31770288d+02, 0.10000000d+01/
    DATA (GA( 6,12,IC),IC=1,3) / &
         -0.30926524d+01, 0.43555255d+00, 0.00000000d+00/
    DATA (GB( 6,12,IC),IC=1,3) / &
         -0.30926524d+01,-0.67432659d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 6 ----- T =  262.5
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA( 7,11,IC),IC=1,3) / &
         -0.32050918d+03, 0.12373350d+02, 0.00000000d+00/
    DATA (GB( 7,11,IC),IC=1,3) / &
         -0.32050918d+03,-0.74061287d+03, 0.10000000d+01/
    DATA (GA( 7,12,IC),IC=1,3) / &
         0.85742941d+00, 0.50380874d+00, 0.00000000d+00/
    DATA (GB( 7,12,IC),IC=1,3) / &
         0.85742941d+00, 0.24550746d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 6 ----- T =  275.0
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA( 8,11,IC),IC=1,3) / &
         -0.37133165d+01, 0.44809588d+00, 0.00000000d+00/
    DATA (GB( 8,11,IC),IC=1,3) / &
         -0.37133165d+01,-0.81329826d+01, 0.10000000d+01/
    DATA (GA( 8,12,IC),IC=1,3) / &
         0.19164038d+01, 0.68537352d+00, 0.00000000d+00/
    DATA (GB( 8,12,IC),IC=1,3) / &
         0.19164038d+01, 0.49089917d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 6 ----- T =  287.5
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA( 9,11,IC),IC=1,3) / &
         0.18890836d+00, 0.46548918d+00, 0.00000000d+00/
    DATA (GB( 9,11,IC),IC=1,3) / &
         0.18890836d+00, 0.90279822d+00, 0.10000000d+01/
    DATA (GA( 9,12,IC),IC=1,3) / &
         0.23513199d+01, 0.89437630d+00, 0.00000000d+00/
    DATA (GB( 9,12,IC),IC=1,3) / &
         0.23513199d+01, 0.59008712d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 6 ----- T =  300.0
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA(10,11,IC),IC=1,3) / &
         0.14209226d+01, 0.59121475d+00, 0.00000000d+00/
    DATA (GB(10,11,IC),IC=1,3) / &
         0.14209226d+01, 0.37532746d+01, 0.10000000d+01/
    DATA (GA(10,12,IC),IC=1,3) / &
         0.25566644d+01, 0.11127003d+01, 0.00000000d+00/
    DATA (GB(10,12,IC),IC=1,3) / &
         0.25566644d+01, 0.63532616d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 6 ----- T =  312.5
    !
    !-- INDICES FOR PADE APPROXIMATION   1 35 40 45
    DATA (GA(11,11,IC),IC=1,3) / &
         0.19817679d+01, 0.74676119d+00, 0.00000000d+00/
    DATA (GB(11,11,IC),IC=1,3) / &
         0.19817679d+01, 0.50437916d+01, 0.10000000d+01/
    DATA (GA(11,12,IC),IC=1,3) / &
         0.26555181d+01, 0.13329782d+01, 0.00000000d+00/
    DATA (GB(11,12,IC),IC=1,3) / &
         0.26555181d+01, 0.65558627d+01, 0.10000000d+01/
    !
    !
    !
    !
    !
    !-- END WATER VAPOR
    !
    !
    !-- CO2 -- INT.2 -- 500-800 CM-1 --- FROM ABS225 ----------------------
    !
    !
    !
    !-- FIU = 0.8 + MAX(0.35,(7-IU)*0.9)  , X/T,  9
    !
    !----- INTERVAL = 2 ----- T =  187.5
    !
    !-- INDICES FOR PADE APPROXIMATION   1 30 38 45
    DATA (GA( 1,13,IC),IC=1,3) / &
         0.87668459d-01, 0.13845511d+01, 0.00000000d+00/
    DATA (GB( 1,13,IC),IC=1,3) / &
         0.87668459d-01, 0.23203798d+01, 0.10000000d+01/
    DATA (GA( 1,14,IC),IC=1,3) / &
         0.74878820d-01, 0.11718758d+01, 0.00000000d+00/
    DATA (GB( 1,14,IC),IC=1,3) / &
         0.74878820d-01, 0.20206726d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  200.0
    !
    !-- INDICES FOR PADE APPROXIMATION   1 30 38 45
    DATA (GA( 2,13,IC),IC=1,3) / &
         0.83754276d-01, 0.13187042d+01, 0.00000000d+00/
    DATA (GB( 2,13,IC),IC=1,3) / &
         0.83754276d-01, 0.22288925d+01, 0.10000000d+01/
    DATA (GA( 2,14,IC),IC=1,3) / &
         0.71650966d-01, 0.11216131d+01, 0.00000000d+00/
    DATA (GB( 2,14,IC),IC=1,3) / &
         0.71650966d-01, 0.19441824d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  212.5
    !
    !-- INDICES FOR PADE APPROXIMATION   1 30 38 45
    DATA (GA( 3,13,IC),IC=1,3) / &
         0.80460283d-01, 0.12644396d+01, 0.00000000d+00/
    DATA (GB( 3,13,IC),IC=1,3) / &
         0.80460283d-01, 0.21515593d+01, 0.10000000d+01/
    DATA (GA( 3,14,IC),IC=1,3) / &
         0.68979615d-01, 0.10809473d+01, 0.00000000d+00/
    DATA (GB( 3,14,IC),IC=1,3) / &
         0.68979615d-01, 0.18807257d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  225.0
    !
    !-- INDICES FOR PADE APPROXIMATION   1 30 38 45
    DATA (GA( 4,13,IC),IC=1,3) / &
         0.77659686d-01, 0.12191543d+01, 0.00000000d+00/
    DATA (GB( 4,13,IC),IC=1,3) / &
         0.77659686d-01, 0.20855896d+01, 0.10000000d+01/
    DATA (GA( 4,14,IC),IC=1,3) / &
         0.66745345d-01, 0.10476396d+01, 0.00000000d+00/
    DATA (GB( 4,14,IC),IC=1,3) / &
         0.66745345d-01, 0.18275618d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  237.5
    !
    !-- INDICES FOR PADE APPROXIMATION   1 30 38 45
    DATA (GA( 5,13,IC),IC=1,3) / &
         0.75257056d-01, 0.11809511d+01, 0.00000000d+00/
    DATA (GB( 5,13,IC),IC=1,3) / &
         0.75257056d-01, 0.20288489d+01, 0.10000000d+01/
    DATA (GA( 5,14,IC),IC=1,3) / &
         0.64857571d-01, 0.10200373d+01, 0.00000000d+00/
    DATA (GB( 5,14,IC),IC=1,3) / &
         0.64857571d-01, 0.17825910d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  250.0
    !
    !-- INDICES FOR PADE APPROXIMATION   1 30 38 45
    DATA (GA( 6,13,IC),IC=1,3) / &
         0.73179175d-01, 0.11484154d+01, 0.00000000d+00/
    DATA (GB( 6,13,IC),IC=1,3) / &
         0.73179175d-01, 0.19796791d+01, 0.10000000d+01/
    DATA (GA( 6,14,IC),IC=1,3) / &
         0.63248495d-01, 0.99692726d+00, 0.00000000d+00/
    DATA (GB( 6,14,IC),IC=1,3) / &
         0.63248495d-01, 0.17442308d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  262.5
    !
    !-- INDICES FOR PADE APPROXIMATION   1 30 38 45
    DATA (GA( 7,13,IC),IC=1,3) / &
         0.71369063d-01, 0.11204723d+01, 0.00000000d+00/
    DATA (GB( 7,13,IC),IC=1,3) / &
         0.71369063d-01, 0.19367778d+01, 0.10000000d+01/
    DATA (GA( 7,14,IC),IC=1,3) / &
         0.61866970d-01, 0.97740923d+00, 0.00000000d+00/
    DATA (GB( 7,14,IC),IC=1,3) / &
         0.61866970d-01, 0.17112809d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  275.0
    !
    !-- INDICES FOR PADE APPROXIMATION   1 30 38 45
    DATA (GA( 8,13,IC),IC=1,3) / &
         0.69781812d-01, 0.10962918d+01, 0.00000000d+00/
    DATA (GB( 8,13,IC),IC=1,3) / &
         0.69781812d-01, 0.18991112d+01, 0.10000000d+01/
    DATA (GA( 8,14,IC),IC=1,3) / &
         0.60673632d-01, 0.96080188d+00, 0.00000000d+00/
    DATA (GB( 8,14,IC),IC=1,3) / &
         0.60673632d-01, 0.16828137d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  287.5
    !
    !-- INDICES FOR PADE APPROXIMATION   1 30 38 45
    DATA (GA( 9,13,IC),IC=1,3) / &
         0.68381606d-01, 0.10752229d+01, 0.00000000d+00/
    DATA (GB( 9,13,IC),IC=1,3) / &
         0.68381606d-01, 0.18658501d+01, 0.10000000d+01/
    DATA (GA( 9,14,IC),IC=1,3) / &
         0.59637277d-01, 0.94657562d+00, 0.00000000d+00/
    DATA (GB( 9,14,IC),IC=1,3) / &
         0.59637277d-01, 0.16580908d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  300.0
    !
    !-- INDICES FOR PADE APPROXIMATION   1 30 38 45
    DATA (GA(10,13,IC),IC=1,3) / &
         0.67139539d-01, 0.10567474d+01, 0.00000000d+00/
    DATA (GB(10,13,IC),IC=1,3) / &
         0.67139539d-01, 0.18363226d+01, 0.10000000d+01/
    DATA (GA(10,14,IC),IC=1,3) / &
         0.58732178d-01, 0.93430511d+00, 0.00000000d+00/
    DATA (GB(10,14,IC),IC=1,3) / &
         0.58732178d-01, 0.16365014d+01, 0.10000000d+01/
    !
    !----- INTERVAL = 2 ----- T =  312.5
    !
    !-- INDICES FOR PADE APPROXIMATION   1 30 38 45
    DATA (GA(11,13,IC),IC=1,3) / &
         0.66032012d-01, 0.10404465d+01, 0.00000000d+00/
    DATA (GB(11,13,IC),IC=1,3) / &
         0.66032012d-01, 0.18099779d+01, 0.10000000d+01/
    DATA (GA(11,14,IC),IC=1,3) / &
         0.57936092d-01, 0.92363528d+00, 0.00000000d+00/
    DATA (GB(11,14,IC),IC=1,3) / &
         0.57936092d-01, 0.16175164d+01, 0.10000000d+01/
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !-- CARBON DIOXIDE LINES IN THE WINDOW REGION (800-1250 CM-1)
    !
    !
    !-- G = 0.0
    !
    !
    !----- INTERVAL = 4 ----- T =  187.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA( 1,15,IC),IC=1,3) / &
         0.13230067d+02, 0.22042132d+02, 0.00000000d+00/
    DATA (GB( 1,15,IC),IC=1,3) / &
         0.13230067d+02, 0.22051750d+02, 0.10000000d+01/
    DATA (GA( 1,16,IC),IC=1,3) / &
         0.13183816d+02, 0.22169501d+02, 0.00000000d+00/
    DATA (GB( 1,16,IC),IC=1,3) / &
         0.13183816d+02, 0.22178972d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  200.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA( 2,15,IC),IC=1,3) / &
         0.13213564d+02, 0.22107298d+02, 0.00000000d+00/
    DATA (GB( 2,15,IC),IC=1,3) / &
         0.13213564d+02, 0.22116850d+02, 0.10000000d+01/
    DATA (GA( 2,16,IC),IC=1,3) / &
         0.13189991d+02, 0.22270075d+02, 0.00000000d+00/
    DATA (GB( 2,16,IC),IC=1,3) / &
         0.13189991d+02, 0.22279484d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  212.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA( 3,15,IC),IC=1,3) / &
         0.13209140d+02, 0.22180915d+02, 0.00000000d+00/
    DATA (GB( 3,15,IC),IC=1,3) / &
         0.13209140d+02, 0.22190410d+02, 0.10000000d+01/
    DATA (GA( 3,16,IC),IC=1,3) / &
         0.13209485d+02, 0.22379193d+02, 0.00000000d+00/
    DATA (GB( 3,16,IC),IC=1,3) / &
         0.13209485d+02, 0.22388551d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  225.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA( 4,15,IC),IC=1,3) / &
         0.13213894d+02, 0.22259478d+02, 0.00000000d+00/
    DATA (GB( 4,15,IC),IC=1,3) / &
         0.13213894d+02, 0.22268925d+02, 0.10000000d+01/
    DATA (GA( 4,16,IC),IC=1,3) / &
         0.13238789d+02, 0.22492992d+02, 0.00000000d+00/
    DATA (GB( 4,16,IC),IC=1,3) / &
         0.13238789d+02, 0.22502309d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  237.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA( 5,15,IC),IC=1,3) / &
         0.13225963d+02, 0.22341039d+02, 0.00000000d+00/
    DATA (GB( 5,15,IC),IC=1,3) / &
         0.13225963d+02, 0.22350445d+02, 0.10000000d+01/
    DATA (GA( 5,16,IC),IC=1,3) / &
         0.13275017d+02, 0.22608508d+02, 0.00000000d+00/
    DATA (GB( 5,16,IC),IC=1,3) / &
         0.13275017d+02, 0.22617792d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  250.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA( 6,15,IC),IC=1,3) / &
         0.13243806d+02, 0.22424247d+02, 0.00000000d+00/
    DATA (GB( 6,15,IC),IC=1,3) / &
         0.13243806d+02, 0.22433617d+02, 0.10000000d+01/
    DATA (GA( 6,16,IC),IC=1,3) / &
         0.13316096d+02, 0.22723843d+02, 0.00000000d+00/
    DATA (GB( 6,16,IC),IC=1,3) / &
         0.13316096d+02, 0.22733099d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  262.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA( 7,15,IC),IC=1,3) / &
         0.13266104d+02, 0.22508089d+02, 0.00000000d+00/
    DATA (GB( 7,15,IC),IC=1,3) / &
         0.13266104d+02, 0.22517429d+02, 0.10000000d+01/
    DATA (GA( 7,16,IC),IC=1,3) / &
         0.13360555d+02, 0.22837837d+02, 0.00000000d+00/
    DATA (GB( 7,16,IC),IC=1,3) / &
         0.13360555d+02, 0.22847071d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  275.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA( 8,15,IC),IC=1,3) / &
         0.13291782d+02, 0.22591771d+02, 0.00000000d+00/
    DATA (GB( 8,15,IC),IC=1,3) / &
         0.13291782d+02, 0.22601086d+02, 0.10000000d+01/
    DATA (GA( 8,16,IC),IC=1,3) / &
         0.13407324d+02, 0.22949751d+02, 0.00000000d+00/
    DATA (GB( 8,16,IC),IC=1,3) / &
         0.13407324d+02, 0.22958967d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  287.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA( 9,15,IC),IC=1,3) / &
         0.13319961d+02, 0.22674661d+02, 0.00000000d+00/
    DATA (GB( 9,15,IC),IC=1,3) / &
         0.13319961d+02, 0.22683956d+02, 0.10000000d+01/
    DATA (GA( 9,16,IC),IC=1,3) / &
         0.13455544d+02, 0.23059032d+02, 0.00000000d+00/
    DATA (GB( 9,16,IC),IC=1,3) / &
         0.13455544d+02, 0.23068234d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  300.0
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA(10,15,IC),IC=1,3) / &
         0.13349927d+02, 0.22756246d+02, 0.00000000d+00/
    DATA (GB(10,15,IC),IC=1,3) / &
         0.13349927d+02, 0.22765522d+02, 0.10000000d+01/
    DATA (GA(10,16,IC),IC=1,3) / &
         0.13504450d+02, 0.23165146d+02, 0.00000000d+00/
    DATA (GB(10,16,IC),IC=1,3) / &
         0.13504450d+02, 0.23174336d+02, 0.10000000d+01/
    !
    !----- INTERVAL = 4 ----- T =  312.5
    !
    !-- INDICES FOR PADE APPROXIMATION     1   15   29   45
    DATA (GA(11,15,IC),IC=1,3) / &
         0.13381108d+02, 0.22836093d+02, 0.00000000d+00/
    DATA (GB(11,15,IC),IC=1,3) / &
         0.13381108d+02, 0.22845354d+02, 0.10000000d+01/
    DATA (GA(11,16,IC),IC=1,3) / &
         0.13553282d+02, 0.23267456d+02, 0.00000000d+00/
    DATA (GB(11,16,IC),IC=1,3) / &
         0.13553282d+02, 0.23276638d+02, 0.10000000d+01/

    !     ------------------------------------------------------------------
    DATA (( XP(  J,K),J=1,6),       K=1,6) / &
         0.46430621d+02, 0.12928299d+03, 0.20732648d+03, &
         0.31398411d+03, 0.18373177d+03,-0.11412303d+03, &
         0.73604774d+02, 0.27887914d+03, 0.27076947d+03, &
         -0.57322111d+02,-0.64742459d+02, 0.87238280d+02, &
         0.37050866d+02, 0.20498759d+03, 0.37558029d+03, &
         0.17401171d+03,-0.13350302d+03,-0.37651795d+02, &
         0.14930141d+02, 0.89161160d+02, 0.17793062d+03, &
         0.93433860d+02,-0.70646020d+02,-0.26373150d+02, &
         0.40386780d+02, 0.10855270d+03, 0.50755010d+02, &
         -0.31496190d+02, 0.12791300d+00, 0.18017770d+01, &
         0.90811926d+01, 0.75073923d+02, 0.24654438d+03, &
         0.39332612d+03, 0.29385281d+03, 0.89107921d+02 /
    !
    !
    !*         1.0     PLANCK FUNCTIONS AND GRADIENTS
    !                  ------------------------------
    !
    DO JK = 1 , LLM+1
       DO JL = 1, KLON
          PBINT(JL,JK) = 0.
       end DO
    end DO
    DO JL = 1, KLON
       PBSUIN(JL) = 0.
    end DO
    !
    DO JNU=1,Ninter
       !
       !
       !*         1.1   LEVELS FROM SURFACE TO LLM
       !                ----------------------------
       !
       DO JK = 1 , LLM
          DO JL = 1, KLON
             ZTI(JL)=(PTL(JL,JK)-TSTAND)/TSTAND
             ZRES(JL) = XP(1,JNU)+ZTI(JL)*(XP(2,JNU)+ZTI(JL)*(XP(3,JNU) &
                  +ZTI(JL)*(XP(4,JNU)+ZTI(JL)*(XP(5,JNU)+ZTI(JL)*(XP(6,JNU) &
                  )))))
             PBINT(JL,JK)=PBINT(JL,JK)+ZRES(JL)
             PB(JL,JNU,JK)= ZRES(JL)
             ZBLEV(JL,JK) = ZRES(JL)
             ZTI2(JL)=(PTAVE(JL,JK)-TSTAND)/TSTAND
             ZRES2(JL)=XP(1,JNU)+ZTI2(JL)*(XP(2,JNU)+ZTI2(JL)*(XP(3,JNU) &
                  +ZTI2(JL)*(XP(4,JNU)+ZTI2(JL)*(XP(5,JNU)+ZTI2(JL)*(XP(6,JNU) &
                  )))))
             ZBLAY(JL,JK) = ZRES2(JL)
          end DO
       end DO
       !
       !
       !*         1.2   TOP OF THE ATMOSPHERE AND SURFACE
       !                ---------------------------------
       !
       DO JL = 1, KLON
          ZTI(JL)=(PTL(JL,LLM+1)-TSTAND)/TSTAND
          ZTI2(JL) = (PTL(JL,1) + PDT0(JL) - TSTAND) / TSTAND
          ZRES(JL) = XP(1,JNU)+ZTI(JL)*(XP(2,JNU)+ZTI(JL)*(XP(3,JNU) &
               +ZTI(JL)*(XP(4,JNU)+ZTI(JL)*(XP(5,JNU)+ZTI(JL)*(XP(6,JNU) &
               )))))
          ZRES2(JL) = XP(1,JNU)+ZTI2(JL)*(XP(2,JNU)+ZTI2(JL)*(XP(3,JNU) &
               +ZTI2(JL)*(XP(4,JNU)+ZTI2(JL)*(XP(5,JNU)+ZTI2(JL)*(XP(6,JNU) &
               )))))
          PBINT(JL,LLM+1) = PBINT(JL,LLM+1)+ZRES(JL)
          PB(JL,JNU,LLM+1)= ZRES(JL)
          ZBLEV(JL,LLM+1) = ZRES(JL)
          PBTOP(JL,JNU) = ZRES(JL)
          PBSUR(JL,JNU) = ZRES2(JL)
          PBSUIN(JL) = PBSUIN(JL) + ZRES2(JL)
       end DO
       !
       !
       !*         1.3   GRADIENTS IN SUB-LAYERS
       !                -----------------------
       !
       DO JK = 1 , LLM
          JK2 = 2 * JK
          JK1 = JK2 - 1
          DO JL = 1, KLON
             PDBSL(JL,JNU,JK1) = ZBLAY(JL,JK  ) - ZBLEV(JL,JK)
             PDBSL(JL,JNU,JK2) = ZBLEV(JL,JK+1) - ZBLAY(JL,JK)
          end DO
       end DO
       !
    end DO
    !
    !*         2.0   CHOOSE THE RELEVANT SETS OF PADE APPROXIMANTS
    !                ---------------------------------------------
    !
    DO JL=1, KLON
       ZDSTO1 = (PTL(JL,LLM+1)-TINTP(1)) / TSTP
       IXTOX = MAX( 1, MIN( MXIXT, INT( ZDSTO1 + 1. ) ) )
       ZDSTOX = (PTL(JL,LLM+1)-TINTP(IXTOX))/TSTP
       IF (ZDSTOX.LT.0.5) THEN
          INDTO=IXTOX
       ELSE
          INDTO=IXTOX+1
       END IF
       INDB(JL)=INDTO
       ZDST1 = (PTL(JL,1)-TINTP(1)) / TSTP
       IXTX = MAX( 1, MIN( MXIXT, INT( ZDST1 + 1. ) ) )
       ZDSTX = (PTL(JL,1)-TINTP(IXTX))/TSTP
       IF (ZDSTX.LT.0.5) THEN
          INDT=IXTX
       ELSE
          INDT=IXTX+1
       END IF
       INDS(JL)=INDT
    end DO
    !
    DO JF=1,2
       DO JG=1, 8
          DO JL=1, KLON
             INDSU=INDS(JL)
             PGASUR(JL,JG,JF)=GA(INDSU,2*JG-1,JF)
             PGBSUR(JL,JG,JF)=GB(INDSU,2*JG-1,JF)
             INDTP=INDB(JL)
             PGATOP(JL,JG,JF)=GA(INDTP,2*JG-1,JF)
             PGBTOP(JL,JG,JF)=GB(INDTP,2*JG-1,JF)
          end DO
       end DO
    end DO
    !
    DO JK=1,LLM
       DO JL=1, KLON
          ZDST1 = (PTAVE(JL,JK)-TINTP(1)) / TSTP
          IXTX = MAX( 1, MIN( MXIXT, INT( ZDST1 + 1. ) ) )
          ZDSTX = (PTAVE(JL,JK)-TINTP(IXTX))/TSTP
          IF (ZDSTX.LT.0.5) THEN
             INDT=IXTX
          ELSE
             INDT=IXTX+1
          END IF
          INDB(JL)=INDT
       end DO
       !
       DO JF=1,2
          DO JG=1, 8
             DO JL=1, KLON
                INDT=INDB(JL)
                PGA(JL,JG,JF,JK)=GA(INDT,2*JG,JF)
                PGB(JL,JG,JF,JK)=GB(INDT,2*JG,JF)
             end DO
          end DO
       end DO
    end DO
    !
    !     ------------------------------------------------------------------
    !
    RETURN
  END SUBROUTINE LWB

end module lwb_m
