module diagphy_m

  implicit none

contains

  SUBROUTINE diagphy(airephy, tit, iprt, tops, topl, sols, soll, sens, evap, &
       rain_fall, snow_fall, ts, d_etp_tot, d_qt_tot, d_ec_tot, fs_bound, &
       fq_bound)

    ! From LMDZ4/libf/phylmd/diagphy.F, version 1.1.1.1 2004/05/19 12:53:08

    ! Purpose: compute the thermal flux and the water mass flux at
    ! the atmospheric boundaries. Print them and print the atmospheric
    ! enthalpy change and the atmospheric mass change.

    ! J.-L. Dufresne, July 2002

    USE dimphy, ONLY: klon
    USE suphec_m, ONLY: rcpd, rcpv, rcs, rcw, rlstt, rlvtt

    ! Arguments: 

    ! Input variables
    real airephy(klon)
    ! airephy-------input-R- grid area
    CHARACTER(len=15) tit
    ! tit---------input-A15- Comment to be added in PRINT (CHARACTER*15)
    INTEGER iprt
    ! iprt--------input-I- PRINT level (<=0 : no PRINT)
    real tops(klon), sols(klon)
    ! tops(klon)--input-R- SW rad. at TOA (W/m2), positive up.
    ! sols(klon)--input-R- Net SW flux above surface (W/m2), positive up 
    ! (i.e. -1 * flux absorbed by the surface)

    real, intent(in):: soll(klon)
    ! net longwave flux above surface (W/m2), positive up (i. e. flux emited
    ! - flux absorbed by the surface)

    real, intent(in):: topl(klon) !LW rad. at TOA (W/m2), positive down
    real sens(klon)
    ! sens(klon)--input-R- Sensible Flux at surface (W/m2), positive down
    real evap(klon)
    ! evap(klon)--input-R- Evaporation + sublimation water vapour mass flux
    ! (kg/m2/s), positive up

    real, intent(in):: rain_fall(klon)
    ! liquid water mass flux (kg/m2/s), positive down

    real snow_fall(klon)
    ! snow_fall(klon)
    ! --input-R- Solid water mass flux (kg/m2/s), positive down
    REAL ts(klon)
    ! ts(klon)----input-R- Surface temperature (K)
    REAL d_etp_tot, d_qt_tot, d_ec_tot
    ! d_etp_tot---input-R- Heat flux equivalent to atmospheric enthalpy 
    ! change (W/m2)
    ! d_qt_tot----input-R- Mass flux equivalent to atmospheric water mass 
    ! change (kg/m2/s)
    ! d_ec_tot----input-R- Flux equivalent to atmospheric cinetic energy
    ! change (W/m2)

    ! Output variables
    REAL fs_bound
    ! fs_bound---output-R- Thermal flux at the atmosphere boundaries (W/m2)
    real fq_bound
    ! fq_bound---output-R- Water mass flux at the atmosphere
    ! boundaries (kg/m2/s)

    ! Local variables:

    real stops, stopl, ssols, ssoll
    real ssens, sfront, slat
    real airetot, zcpvap, zcwat, zcice
    REAL rain_fall_tot, snow_fall_tot, evap_tot

    integer i
    integer:: pas = 0

    !------------------------------------------------------------------

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

    ! Pour les chaleur specifiques de la vapeur d'eau, de l'eau et de
    ! la glace, on travaille par difference a la chaleur specifique de
    ! l' air sec. En effet, comme on travaille a niveau de pression
    ! donne, toute variation de la masse d'un constituant est
    ! totalement compense par une variation de masse d'air.

    zcpvap=RCPV-RCPD
    zcwat=RCW-RCPD
    zcice=RCS-RCPD

    do i=1, klon
       stops=stops+tops(i)*airephy(i)
       stopl=stopl+topl(i)*airephy(i)
       ssols=ssols+sols(i)*airephy(i)
       ssoll=ssoll+soll(i)*airephy(i)
       ssens=ssens+sens(i)*airephy(i)
       sfront = sfront &
            + (evap(i)*zcpvap-rain_fall(i)*zcwat-snow_fall(i)*zcice) * ts(i) &
            * airephy(i)
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

    slat = RLVTT * rain_fall_tot + RLSTT * snow_fall_tot
    ! Heat flux at atm. boundaries
    fs_bound = stops-stopl - (ssols+ssoll)+ssens+sfront + slat
    ! Water flux at atm. boundaries
    fq_bound = evap_tot - rain_fall_tot -snow_fall_tot

    IF (iprt >= 1) print 6666, tit, pas, fs_bound, d_etp_tot, fq_bound, d_qt_tot

    IF (iprt >= 1) print 6668, tit, pas, d_etp_tot+d_ec_tot-fs_bound, &
         d_qt_tot-fq_bound

    IF (iprt >= 2) print 6667, tit, pas, stops, stopl, ssols, ssoll, ssens, &
         slat, evap_tot, rain_fall_tot+snow_fall_tot

6666 format('Phys. Flux Budget ', a15, 1i6, 2f8.2, 2(1pE13.5))
6667 format('Phys. Boundary Flux ', a15, 1i6, 6f8.2, 2(1pE13.5))
6668 format('Phys. Total Budget ', a15, 1i6, f8.2, 2(1pE13.5))

  end SUBROUTINE diagphy

end module diagphy_m
