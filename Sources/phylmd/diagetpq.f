module diagetpq_m

  IMPLICIT NONE

contains

  SUBROUTINE diagetpq(airephy, tit, iprt, idiag, idiag2, dtime, t, q, ql, &
       u, v, paprs, d_h_vcol, d_qt, d_ec)

    ! From LMDZ4/libf/phylmd/diagphy.F, version 1.1.1.1, 2004/05/19 12:53:08

    ! Purpose:

    ! Calcule la différence d'enthalpie et de masse d'eau entre deux
    ! appels et calcule le flux de chaleur et le flux d'eau
    ! nécessaires à ces changements. Ces valeurs sont moyennées sur la
    ! surface de tout le globe et sont exprimées en W/m2 et
    ! kg/s/m2. Outil pour diagnostiquer la conservation de l'énergie
    ! et de la masse dans la physique. Suppose que les niveaux de
    ! pression entre les couches ne varient pas entre deux appels.

    ! Plusieurs de ces diagnostics peuvent être faits en parallèle :
    ! les bilans sont sauvegardés dans des tableaux indices. On
    ! parlera "d'indice de diagnostic".

    ! Jean-Louis Dufresne, July 2002

    USE dimphy, ONLY: klev, klon
    USE suphec_m, ONLY: rcpd, rcpv, rcs, rcw, rg, rlstt, rlvtt

    ! Arguments: 

    ! Input variables
    real, intent(in):: airephy(klon) ! grid area
    CHARACTER(len=*), intent(in):: tit ! comment added in PRINT
    INTEGER, intent(in):: iprt ! PRINT level ( <=1 : no PRINT)

    INTEGER, intent(in):: idiag
    ! indice dans lequel seront rangés les nouveaux bilans d'enthalpie et
    ! de masse

    INTEGER, intent(in):: idiag2
    ! Les nouveaux bilans d'enthalpie et de masse sont comparés au
    ! bilan de d'enthalpie de masse de l'indice numéro idiag2. Cas
    ! particulier : si idiag2=0, pas de comparaison, on sort
    ! directement les bilans d'enthalpie et de masse.

    REAL, intent(in):: dtime ! time step (s)
    REAL, intent(in):: t(klon, klev) ! temperature (K)
    REAL, intent(in):: q(klon, klev) ! vapeur d'eau (kg/kg)
    REAL, intent(in):: ql(klon, klev) ! liquid water (kg/kg)
    REAL, intent(in):: u(klon, klev), v(klon, klev) ! vitesse
    REAL, intent(in):: paprs(klon, klev+1) ! pression a intercouche (Pa)

    ! The following total values are computed per UNIT of Earth surface

    REAL, intent(out):: d_h_vcol
    ! heat flux (W/m2) defined as the enthalpy change (J/m2) during
    ! one time step (dtime) for the whole atmosphere (air, water
    ! vapour, liquid and solid)

    REAL, intent(out):: d_qt
    ! total water mass flux (kg/m2/s) defined as the 
    ! total water (kg/m2) change during one time step (dtime)

    REAL, intent(out):: d_ec
    ! kinetic energy budget (W/m2) for vertical air column

    ! Local:

    REAL d_qw
    ! water vapour mass flux (kg/m2/s) defined as the water vapour
    ! (kg/m2) change during one time step (dtime)

    REAL d_ql ! same, for the liquid water only (kg/m2/s)

    REAL h_vcol_tot
    ! total enthalpy of vertical air column (air with water vapour,
    ! liquid and solid) (J/m2)

    REAL h_dair_tot ! total enthalpy of dry air (J/m2)
    REAL h_qw_tot ! total enthalpy of water vapour (J/m2)
    REAL h_ql_tot ! total enthalpy of liquid water (J/m2)
    REAL qw_tot ! total mass of water vapour (kg/m2)
    REAL ql_tot ! total mass of liquid water (kg/m2)
    real ec_tot ! total kinetic energy (kg/m2)
    REAL zairm(klon, klev) ! layer air mass (kg/m2)
    REAL zqw_col(klon)
    REAL zql_col(klon)
    REAL zec_col(klon)
    REAL zh_dair_col(klon)
    REAL zh_qw_col(klon), zh_ql_col(klon)
    REAL d_h_dair, d_h_qw, d_h_ql
    REAL airetot, zcpvap, zcwat, zcice
    INTEGER i, k
    INTEGER, PARAMETER:: ndiag = 10 ! max number of diagnostic in parallel
    integer:: pas(ndiag) = 0
    REAL, save:: h_vcol_pre(ndiag), h_dair_pre(ndiag), h_qw_pre(ndiag)
    REAL, save:: h_ql_pre(ndiag), qw_pre(ndiag), ql_pre(ndiag)
    REAL, save:: ec_pre(ndiag)

    !-------------------------------------------------------------

    DO k = 1, klev
       DO i = 1, klon
          zairm(i, k) = (paprs(i, k)-paprs(i, k+1))/RG
       ENDDO
    END DO

    ! Reset variables
    DO i = 1, klon
       zqw_col(i)=0.
       zql_col(i)=0.
       zec_col(i) = 0.
       zh_dair_col(i) = 0.
       zh_qw_col(i) = 0.
       zh_ql_col(i) = 0.
    ENDDO

    zcpvap=RCPV
    zcwat=RCW
    zcice=RCS

    ! Compute vertical sum for each atmospheric column
    DO k = 1, klev
       DO i = 1, klon
          ! Water mass
          zqw_col(i) = zqw_col(i) + q(i, k)*zairm(i, k)
          zql_col(i) = zql_col(i) + ql(i, k)*zairm(i, k)
          ! Kinetic Energy
          zec_col(i) = zec_col(i) +0.5*(u(i, k)**2+v(i, k)**2)*zairm(i, k)
          ! Air enthalpy
          zh_dair_col(i) = zh_dair_col(i) &
               + RCPD*(1.-q(i, k)-ql(i, k))*zairm(i, k)*t(i, k)
          zh_qw_col(i) = zh_qw_col(i) + zcpvap*q(i, k)*zairm(i, k)*t(i, k) 
          zh_ql_col(i) = zh_ql_col(i) &
               + zcwat*ql(i, k)*zairm(i, k)*t(i, k) &
               - RLVTT*ql(i, k)*zairm(i, k)
       END DO
    ENDDO

    ! Mean over the planet surface
    qw_tot = 0.
    ql_tot = 0.
    ec_tot = 0.
    h_vcol_tot = 0.
    h_dair_tot = 0.
    h_qw_tot = 0.
    h_ql_tot = 0.
    airetot=0.

    do i=1, klon
       qw_tot = qw_tot + zqw_col(i)*airephy(i)
       ql_tot = ql_tot + zql_col(i)*airephy(i)
       ec_tot = ec_tot + zec_col(i)*airephy(i)
       h_dair_tot = h_dair_tot + zh_dair_col(i)*airephy(i)
       h_qw_tot = h_qw_tot + zh_qw_col(i)*airephy(i)
       h_ql_tot = h_ql_tot + zh_ql_col(i)*airephy(i)
       airetot=airetot+airephy(i)
    END DO

    qw_tot = qw_tot/airetot
    ql_tot = ql_tot/airetot
    ec_tot = ec_tot/airetot
    h_dair_tot = h_dair_tot/airetot
    h_qw_tot = h_qw_tot/airetot
    h_ql_tot = h_ql_tot/airetot

    h_vcol_tot = h_dair_tot+h_qw_tot+h_ql_tot

    ! Compute the change of the atmospheric state compared to the one
    ! stored in "idiag2", and convert it in flux. This computation is
    ! performed if idiag2 /= 0 and if it is not the first call for
    ! "idiag".

    IF (idiag2 > 0 .and. pas(idiag2) /= 0) THEN
       d_h_vcol = (h_vcol_tot - h_vcol_pre(idiag2) )/dtime
       d_h_dair = (h_dair_tot- h_dair_pre(idiag2))/dtime
       d_h_qw = (h_qw_tot - h_qw_pre(idiag2) )/dtime
       d_h_ql = (h_ql_tot - h_ql_pre(idiag2) )/dtime 
       d_qw = (qw_tot - qw_pre(idiag2) )/dtime
       d_ql = (ql_tot - ql_pre(idiag2) )/dtime
       d_ec = (ec_tot - ec_pre(idiag2) )/dtime
       d_qt = d_qw + d_ql
    ELSE 
       d_h_vcol = 0.
       d_h_dair = 0.
       d_h_qw = 0.
       d_h_ql = 0.
       d_qw = 0.
       d_ql = 0.
       d_ec = 0.
       d_qt = 0.
    ENDIF

    IF (iprt >= 2) THEN
       print 9000, tit, pas(idiag), d_qt, d_qw, d_ql
       print 9001, tit, pas(idiag), d_h_vcol
       print 9002, tit, pas(idiag), d_ec
    END IF

    ! Store the new atmospheric state in "idiag"
    pas(idiag)=pas(idiag)+1
    h_vcol_pre(idiag) = h_vcol_tot
    h_dair_pre(idiag) = h_dair_tot
    h_qw_pre(idiag) = h_qw_tot
    h_ql_pre(idiag) = h_ql_tot
    qw_pre(idiag) = qw_tot
    ql_pre(idiag) = ql_tot
    ec_pre (idiag) = ec_tot

9000 format('Physics water mass budget (kg/m2/s)', A15, 1i6, 10(1pE14.6))
9001 format('Physics enthalpy budget (W/m2) ', A15, 1i6, 10(F8.2))
9002 format('Physics kinetic energy budget (W/m2) ', A15, 1i6, 10(F8.2))

  END SUBROUTINE diagetpq

end module diagetpq_m
