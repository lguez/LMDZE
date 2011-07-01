!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/diagphy.F,v 1.1.1.1 2004/05/19 12:53:08 lmdzadmin Exp $
!
      SUBROUTINE diagphy(airephy,tit,iprt
     $    , tops, topl, sols, soll, sens
     $    , evap, rain_fall, snow_fall, ts
     $    , d_etp_tot, d_qt_tot, d_ec_tot
     $    , fs_bound, fq_bound)
C======================================================================
C
C Purpose:
C    Compute the thermal flux and the watter mass flux at the atmosphere
c    boundaries. Print them and also the atmospheric enthalpy change and
C    the  atmospheric mass change.
C
C Arguments: 
C airephy-------input-R-  grid area
C tit---------input-A15- Comment to be added in PRINT (CHARACTER*15)
C iprt--------input-I-  PRINT level ( <=0 : no PRINT)
C tops(klon)--input-R-  SW rad. at TOA (W/m2), positive up.
C topl(klon)--input-R-  LW rad. at TOA (W/m2), positive down
C sols(klon)--input-R-  Net SW flux above surface (W/m2), positive up 
C                   (i.e. -1 * flux absorbed by the surface)
C soll(klon)--input-R-  Net LW flux above surface (W/m2), positive up 
C                   (i.e. flux emited - flux absorbed by the surface)
C sens(klon)--input-R-  Sensible Flux at surface  (W/m2), positive down
C evap(klon)--input-R-  Evaporation + sublimation watter vapour mass flux
C                   (kg/m2/s), positive up
C rain_fall(klon)
C           --input-R- Liquid  watter mass flux (kg/m2/s), positive down
C snow_fall(klon)
C           --input-R- Solid  watter mass flux (kg/m2/s), positive down
C ts(klon)----input-R- Surface temperature (K)
C d_etp_tot---input-R- Heat flux equivalent to atmospheric enthalpy 
C                    change (W/m2)
C d_qt_tot----input-R- Mass flux equivalent to atmospheric watter mass 
C                    change (kg/m2/s)
C d_ec_tot----input-R- Flux equivalent to atmospheric cinetic energy
C                    change (W/m2)
C
C fs_bound---output-R- Thermal flux at the atmosphere boundaries (W/m2)
C fq_bound---output-R- Watter mass flux at the atmosphere boundaries (kg/m2/s)
C
C J.L. Dufresne, July 2002
C======================================================================
C 
      use dimens_m
      use dimphy
      use SUPHEC_M
      use yoethf_m
      implicit none

C
C     Input variables
      real airephy(klon)
      CHARACTER*15 tit
      INTEGER iprt
      real tops(klon),topl(klon),sols(klon),soll(klon)
      real sens(klon),evap(klon),rain_fall(klon),snow_fall(klon)
      REAL ts(klon)
      REAL d_etp_tot, d_qt_tot, d_ec_tot
c     Output variables
      REAL fs_bound, fq_bound
C
C     Local variables
      real stops,stopl,ssols,ssoll
      real ssens,sfront,slat
      real airetot, zcpvap, zcwat, zcice
      REAL rain_fall_tot, snow_fall_tot, evap_tot
C
      integer i
C
      integer pas
      save pas
      data pas/0/
C
      pas=pas+1
      stops=0.
      stopl=0.
      ssols=0.
      ssoll=0.
      ssens=0.
      sfront = 0.
      evap_tot = 0.
      rain_fall_tot = 0.
      snow_fall_tot = 0.
      airetot=0.
C
C     Pour les chaleur specifiques de la vapeur d'eau, de l'eau et de
C     la glace, on travaille par difference a la chaleur specifique de l'
c     air sec. En effet, comme on travaille a niveau de pression donne,
C     toute variation de la masse d'un constituant est totalement
c     compense par une variation de masse d'air.
C
      zcpvap=RCPV-RCPD
      zcwat=RCW-RCPD
      zcice=RCS-RCPD
C
      do i=1,klon
           stops=stops+tops(i)*airephy(i)
           stopl=stopl+topl(i)*airephy(i)
           ssols=ssols+sols(i)*airephy(i)
           ssoll=ssoll+soll(i)*airephy(i)
           ssens=ssens+sens(i)*airephy(i)
           sfront = sfront
     $         + ( evap(i)*zcpvap-rain_fall(i)*zcwat-snow_fall(i)*zcice
     $           ) *ts(i) *airephy(i)
           evap_tot = evap_tot + evap(i)*airephy(i)
           rain_fall_tot = rain_fall_tot + rain_fall(i)*airephy(i)
           snow_fall_tot = snow_fall_tot + snow_fall(i)*airephy(i)
           airetot=airetot+airephy(i)
      enddo
      stops=stops/airetot
      stopl=stopl/airetot
      ssols=ssols/airetot
      ssoll=ssoll/airetot
      ssens=ssens/airetot
      sfront = sfront/airetot
      evap_tot = evap_tot /airetot
      rain_fall_tot = rain_fall_tot/airetot
      snow_fall_tot = snow_fall_tot/airetot
C
      slat = RLVTT * rain_fall_tot + RLSTT * snow_fall_tot
C     Heat flux at atm. boundaries
      fs_bound = stops-stopl - (ssols+ssoll)+ssens+sfront
     $    + slat
C     Watter flux at atm. boundaries
      fq_bound = evap_tot - rain_fall_tot -snow_fall_tot
C
      IF (iprt.ge.1) write(6,6666) 
     $    tit, pas, fs_bound, d_etp_tot, fq_bound, d_qt_tot
C
      IF (iprt.ge.1) write(6,6668) 
     $    tit, pas, d_etp_tot+d_ec_tot-fs_bound, d_qt_tot-fq_bound
C
      IF (iprt.ge.2) write(6,6667) 
     $    tit, pas, stops,stopl,ssols,ssoll,ssens,slat,evap_tot
     $    ,rain_fall_tot+snow_fall_tot

      return

 6666 format('Phys. Flux Budget ',a15,1i6,2f8.2,2(1pE13.5))
 6667 format('Phys. Boundary Flux ',a15,1i6,6f8.2,2(1pE13.5))
 6668 format('Phys. Total Budget ',a15,1i6,f8.2,2(1pE13.5))

      end
