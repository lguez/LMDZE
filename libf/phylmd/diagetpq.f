      SUBROUTINE diagetpq(airephy,tit,iprt,idiag,idiag2,dtime
     e  ,t,q,ql,qs,u,v,paprs
     s  , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/diagphy.F,v 1.1.1.1 2004/05/19 12:53:08 lmdzadmin Exp $
!
C======================================================================
C
C Purpose:
C    Calcul la difference d'enthalpie et de masse d'eau entre 2 appels,
C    et calcul le flux de chaleur et le flux d'eau necessaire a ces 
C    changements. Ces valeurs sont moyennees sur la surface de tout
C    le globe et sont exprime en W/2 et kg/s/m2
C    Outil pour diagnostiquer la conservation de l'energie
C    et de la masse dans la physique. Suppose que les niveau de
c    pression entre couche ne varie pas entre 2 appels.
C
C Plusieurs de ces diagnostics peuvent etre fait en parallele: les
c bilans sont sauvegardes dans des tableaux indices. On parlera
C "d'indice de diagnostic"
c 
C
c======================================================================
C Arguments: 
C airephy-------input-R-  grid area
C tit-----imput-A15- Comment added in PRINT (CHARACTER*15)
C iprt----input-I-  PRINT level ( <=1 : no PRINT)
C idiag---input-I- indice dans lequel sera range les nouveaux
C                  bilans d' entalpie et de masse
C idiag2--input-I-les nouveaux bilans d'entalpie et de masse 
C                 sont compare au bilan de d'enthalpie de masse de
C                 l'indice numero idiag2 
C                 Cas parriculier : si idiag2=0, pas de comparaison, on
c                 sort directement les bilans d'enthalpie et de masse 
C dtime----input-R- time step (s)
c t--------input-R- temperature (K)
c q--------input-R- vapeur d'eau (kg/kg)
c ql-------input-R- liquid watter (kg/kg)
c qs-------input-R- solid watter (kg/kg)
c u--------input-R- vitesse u
c v--------input-R- vitesse v
c paprs----input-R- pression a intercouche (Pa)
c
C the following total value are computed by UNIT of earth surface
C
C d_h_vcol--output-R- Heat flux (W/m2) define as the Enthalpy 
c            change (J/m2) during one time step (dtime) for the whole 
C            atmosphere (air, watter vapour, liquid and solid)
C d_qt------output-R- total water mass flux (kg/m2/s) defined as the 
C           total watter (kg/m2) change during one time step (dtime),
C d_qw------output-R- same, for the watter vapour only (kg/m2/s)
C d_ql------output-R- same, for the liquid watter only (kg/m2/s)
C d_qs------output-R- same, for the solid watter only (kg/m2/s)
C d_ec------output-R- Cinetic Energy Budget (W/m2) for vertical air column
C
C     other (COMMON...)
C     RCPD, RCPV, ....
C
C J.L. Dufresne, July 2002
c======================================================================
 
      use dimens_m
      use dimphy
      use SUPHEC_M
      use yoethf_m
      IMPLICIT NONE
C
C
c     Input variables
      real airephy(klon)
      CHARACTER*15 tit
      INTEGER iprt,idiag, idiag2
      REAL, intent(in):: dtime
      REAL t(klon,klev), q(klon,klev), ql(klon,klev), qs(klon,klev)
      REAL u(klon,klev), v(klon,klev)
      REAL, intent(in):: paprs(klon,klev+1)
c     Output variables
      REAL d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec
C
C     Local variables
c
      REAL h_vcol_tot, h_dair_tot, h_qw_tot, h_ql_tot
     .  , h_qs_tot, qw_tot, ql_tot, qs_tot , ec_tot
c h_vcol_tot--  total enthalpy of vertical air column 
C            (air with watter vapour, liquid and solid) (J/m2)
c h_dair_tot-- total enthalpy of dry air (J/m2)
c h_qw_tot----  total enthalpy of watter vapour (J/m2)
c h_ql_tot----  total enthalpy of liquid watter (J/m2)
c h_qs_tot----  total enthalpy of solid watter  (J/m2)
c qw_tot------  total mass of watter vapour (kg/m2)
c ql_tot------  total mass of liquid watter (kg/m2)
c qs_tot------  total mass of solid watter (kg/m2)
c ec_tot------  total cinetic energy (kg/m2)
C
      REAL zairm(klon,klev) ! layer air mass (kg/m2)
      REAL  zqw_col(klon)
      REAL  zql_col(klon)
      REAL  zqs_col(klon)
      REAL  zec_col(klon)
      REAL  zh_dair_col(klon)
      REAL  zh_qw_col(klon), zh_ql_col(klon), zh_qs_col(klon)
C
      REAL      d_h_dair, d_h_qw, d_h_ql, d_h_qs
C
      REAL airetot, zcpvap, zcwat, zcice
C
      INTEGER i, k
C
      INTEGER ndiag     ! max number of diagnostic in parallel
      PARAMETER (ndiag=10)
      integer pas(ndiag)
      save pas
      data pas/ndiag*0/
C     
      REAL      h_vcol_pre(ndiag), h_dair_pre(ndiag), h_qw_pre(ndiag)
     $    , h_ql_pre(ndiag), h_qs_pre(ndiag), qw_pre(ndiag)
     $    , ql_pre(ndiag), qs_pre(ndiag) , ec_pre(ndiag)
      SAVE      h_vcol_pre, h_dair_pre, h_qw_pre, h_ql_pre
     $        , h_qs_pre, qw_pre, ql_pre, qs_pre , ec_pre

c======================================================================
C
      DO k = 1, klev
        DO i = 1, klon
C         layer air mass
          zairm(i,k) = (paprs(i,k)-paprs(i,k+1))/RG
        ENDDO
      END DO
C
C     Reset variables
      DO i = 1, klon
        zqw_col(i)=0.
        zql_col(i)=0.
        zqs_col(i)=0.
        zec_col(i) = 0.
        zh_dair_col(i) = 0.
        zh_qw_col(i) = 0.
        zh_ql_col(i) = 0.
        zh_qs_col(i) = 0.
      ENDDO
C
      zcpvap=RCPV
      zcwat=RCW
      zcice=RCS
C
C     Compute vertical sum for each atmospheric column
C     ================================================
      DO k = 1, klev
        DO i = 1, klon
C         Watter mass
          zqw_col(i) = zqw_col(i) + q(i,k)*zairm(i,k)
          zql_col(i) = zql_col(i) + ql(i,k)*zairm(i,k)
          zqs_col(i) = zqs_col(i) + qs(i,k)*zairm(i,k)
C         Cinetic Energy
          zec_col(i) =  zec_col(i)
     $        +0.5*(u(i,k)**2+v(i,k)**2)*zairm(i,k)
C         Air enthalpy
          zh_dair_col(i) = zh_dair_col(i) 
     $        + RCPD*(1.-q(i,k)-ql(i,k)-qs(i,k))*zairm(i,k)*t(i,k)
          zh_qw_col(i) = zh_qw_col(i)
     $        + zcpvap*q(i,k)*zairm(i,k)*t(i,k) 
          zh_ql_col(i) = zh_ql_col(i)
     $        + zcwat*ql(i,k)*zairm(i,k)*t(i,k) 
     $        - RLVTT*ql(i,k)*zairm(i,k)
          zh_qs_col(i) = zh_qs_col(i)
     $        + zcice*qs(i,k)*zairm(i,k)*t(i,k) 
     $        - RLSTT*qs(i,k)*zairm(i,k)

        END DO
      ENDDO
C
C     Mean over the planete surface
C     =============================
      qw_tot = 0.
      ql_tot = 0.
      qs_tot = 0.
      ec_tot = 0.
      h_vcol_tot = 0.
      h_dair_tot = 0.
      h_qw_tot = 0.
      h_ql_tot = 0.
      h_qs_tot = 0.
      airetot=0.
C
      do i=1,klon
        qw_tot = qw_tot + zqw_col(i)*airephy(i)
        ql_tot = ql_tot + zql_col(i)*airephy(i)
        qs_tot = qs_tot + zqs_col(i)*airephy(i)
        ec_tot = ec_tot + zec_col(i)*airephy(i)
        h_dair_tot = h_dair_tot + zh_dair_col(i)*airephy(i)
        h_qw_tot = h_qw_tot + zh_qw_col(i)*airephy(i)
        h_ql_tot = h_ql_tot + zh_ql_col(i)*airephy(i)
        h_qs_tot = h_qs_tot + zh_qs_col(i)*airephy(i)
        airetot=airetot+airephy(i)
      END DO
C
      qw_tot = qw_tot/airetot
      ql_tot = ql_tot/airetot
      qs_tot = qs_tot/airetot
      ec_tot = ec_tot/airetot
      h_dair_tot = h_dair_tot/airetot
      h_qw_tot = h_qw_tot/airetot
      h_ql_tot = h_ql_tot/airetot
      h_qs_tot = h_qs_tot/airetot
C
      h_vcol_tot = h_dair_tot+h_qw_tot+h_ql_tot+h_qs_tot
C
C     Compute the change of the atmospheric state compare to the one 
C     stored in "idiag2", and convert it in flux. THis computation
C     is performed IF idiag2 /= 0 and IF it is not the first CALL
c     for "idiag"
C     ===================================
C
      IF ( (idiag2.gt.0) .and. (pas(idiag2) .ne. 0) ) THEN
        d_h_vcol  = (h_vcol_tot - h_vcol_pre(idiag2) )/dtime
        d_h_dair = (h_dair_tot- h_dair_pre(idiag2))/dtime
        d_h_qw   = (h_qw_tot  - h_qw_pre(idiag2)  )/dtime
        d_h_ql   = (h_ql_tot  - h_ql_pre(idiag2)  )/dtime 
        d_h_qs   = (h_qs_tot  - h_qs_pre(idiag2)  )/dtime 
        d_qw     = (qw_tot    - qw_pre(idiag2)    )/dtime
        d_ql     = (ql_tot    - ql_pre(idiag2)    )/dtime
        d_qs     = (qs_tot    - qs_pre(idiag2)    )/dtime
        d_ec     = (ec_tot    - ec_pre(idiag2)    )/dtime
        d_qt = d_qw + d_ql + d_qs
      ELSE 
        d_h_vcol = 0.
        d_h_dair = 0.
        d_h_qw   = 0.
        d_h_ql   = 0.
        d_h_qs   = 0. 
        d_qw     = 0.
        d_ql     = 0.
        d_qs     = 0.
        d_ec     = 0.
        d_qt     = 0.
      ENDIF 
C
      IF (iprt.ge.2) THEN
        WRITE(6,9000) tit,pas(idiag),d_qt,d_qw,d_ql,d_qs
 9000   format('Phys. Watter Mass Budget (kg/m2/s)',A15
     $      ,1i6,10(1pE14.6))
        WRITE(6,9001) tit,pas(idiag), d_h_vcol
 9001   format('Phys. Enthalpy Budget (W/m2) ',A15,1i6,10(F8.2))
        WRITE(6,9002) tit,pas(idiag), d_ec
 9002   format('Phys. Cinetic Energy Budget (W/m2) ',A15,1i6,10(F8.2))
      END IF 
C
C     Store the new atmospheric state in "idiag"
C
      pas(idiag)=pas(idiag)+1
      h_vcol_pre(idiag)  = h_vcol_tot
      h_dair_pre(idiag) = h_dair_tot
      h_qw_pre(idiag)   = h_qw_tot
      h_ql_pre(idiag)   = h_ql_tot
      h_qs_pre(idiag)   = h_qs_tot
      qw_pre(idiag)     = qw_tot
      ql_pre(idiag)     = ql_tot
      qs_pre(idiag)     = qs_tot
      ec_pre (idiag)    = ec_tot
C
      RETURN 
      END
