module cv30_undilute2_m

  implicit none

contains

  SUBROUTINE cv30_undilute2(ncum, icb, icbs, nk, tnk, qnk, gznk, t, qs, gz, &
       p, h, tv, lv, pbase, buoybase, plcl, inb, tp, tvp, clw, hp, ep, sigp, &
       buoy)

    ! Undilute (adiabatic) updraft, second part. Purpose: find the
    ! rest of the lifted parcel temperatures; compute the
    ! precipitation efficiencies and the fraction of precipitation
    ! falling outside of cloud; find the level of neutral buoyancy.

    ! Vertical profile of buoyancy computed here (use of buoybase).

    use conema3_m, only: epmax
    use cv30_param_m, only: dtovsh, minorig, nl, pbcrit, ptcrit, spfac
    use cv_thermo_m, only: cl, clmcpv, cpd, cpv, eps, lv0, rrv
    USE dimphy, ONLY: klon, klev

    integer, intent(in):: ncum

    integer, intent(in):: icb(klon), icbs(klon)
    ! icbs is the first level above LCL (may differ from icb)

    integer, intent(in):: nk(klon)
    real, intent(in):: tnk(klon), qnk(klon), gznk(klon)
    real, intent(in):: t(klon, klev), qs(klon, klev), gz(klon, klev)
    real, intent(in):: p(klon, klev), h(klon, klev)
    real, intent(in):: tv(klon, klev), lv(klon, klev)
    real, intent(in):: pbase(klon), buoybase(klon), plcl(klon)

    ! outputs:
    integer, intent(out):: inb(:) ! (ncum)
    ! first model level above the level of neutral buoyancy of the
    ! parcel (<= nl - 1)

    real tp(klon, klev), tvp(klon, klev), clw(klon, klev)
    ! condensed water not removed from tvp
    real hp(klon, klev), ep(klon, klev), sigp(klon, klev)
    real buoy(klon, klev)

    ! Local:
    integer i, k
    real tg, qg, ahg, alv, s, tc, es, denom
    real pden
    real ah0(klon)

    !---------------------------------------------------------------------

    ! SOME INITIALIZATIONS

    do k = 1, nl
       do i = 1, ncum
          ep(i, k) = 0.0
          sigp(i, k) = spfac
       end do
    end do

    ! FIND THE REST OF THE LIFTED PARCEL TEMPERATURES

    ! The procedure is to solve the equation.
    ! cp * tp + L * qp + phi = cp * tnk + L * qnk + gznk.

    ! Calculate certain parcel quantities, including static energy

    do i = 1, ncum
       ah0(i) = (cpd * (1. - qnk(i)) + cl * qnk(i)) * tnk(i) &
            + qnk(i) * (lv0 - clmcpv * (tnk(i) - 273.15)) + gznk(i)
    end do

    ! Find lifted parcel quantities above cloud base

    do k = minorig + 1, nl
       do i = 1, ncum
          if (k >= (icbs(i) + 1)) then
             tg = t(i, k)
             qg = qs(i, k)
             alv = lv0 - clmcpv * (t(i, k) - 273.15)

             ! First iteration.

             s = cpd * (1. - qnk(i)) + cl * qnk(i) &
                  + alv * alv * qg / (rrv * t(i, k) * t(i, k))
             s = 1. / s

             ahg = cpd * tg + (cl - cpd) * qnk(i) * tg + alv * qg + gz(i, k)
             tg = tg + s * (ah0(i) - ahg)

             tc = tg - 273.15
             denom = 243.5 + tc
             denom = MAX(denom, 1.0)

             es = 6.112 * exp(17.67 * tc / denom)

             qg = eps * es / (p(i, k) - es * (1. - eps))

             ! Second iteration.

             ahg = cpd * tg + (cl - cpd) * qnk(i) * tg + alv * qg + gz(i, k)
             tg = tg + s * (ah0(i) - ahg)

             tc = tg - 273.15
             denom = 243.5 + tc
             denom = MAX(denom, 1.0)

             es = 6.112 * exp(17.67 * tc / denom)

             qg = eps * es / (p(i, k) - es * (1. - eps))

             alv = lv0 - clmcpv * (t(i, k) - 273.15)

             ! no approximation:
             tp(i, k) = (ah0(i) - gz(i, k) - alv * qg) &
                  / (cpd + (cl - cpd) * qnk(i))

             clw(i, k) = qnk(i) - qg
             clw(i, k) = max(0.0, clw(i, k))
             ! qg utilise au lieu du vrai mixing ratio rg:
             tvp(i, k) = tp(i, k) * (1. + qg / eps - qnk(i)) ! whole thing
          endif
       end do
    end do

    ! SET THE PRECIPITATION EFFICIENCIES AND THE FRACTION OF
    ! PRECIPITATION FALLING OUTSIDE OF CLOUD
    ! THESE MAY BE FUNCTIONS OF TP(I), P(I) AND CLW(I)
    do k = 1, nl
       do i = 1, ncum
          pden = ptcrit - pbcrit
          ep(i, k) = (plcl(i) - p(i, k) - pbcrit) / pden * epmax
          ep(i, k) = amax1(ep(i, k), 0.0)
          ep(i, k) = amin1(ep(i, k), epmax)
          sigp(i, k) = spfac
       end do
    end do

    ! CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL
    ! VIRTUAL TEMPERATURE

    ! tvp est calcule en une seule fois, et sans retirer
    ! l'eau condensee (~> reversible CAPE)
    do i = 1, ncum
       tp(i, nl + 1) = tp(i, nl)
    end do

    ! EFFECTIVE VERTICAL PROFILE OF BUOYANCY:

    ! first estimate of buoyancy:
    do i = 1, ncum
       do k = 1, nl
          buoy(i, k) = tvp(i, k) - tv(i, k)
       end do
    end do

    ! set buoyancy = buoybase for all levels below base
    ! for safety, set buoy(icb) = buoybase
    do i = 1, ncum
       do k = 1, nl
          if ((k >= icb(i)) .and. (k <= nl) .and. (p(i, k) >= pbase(i))) then
             buoy(i, k) = buoybase(i)
          endif
       end do
       buoy(icb(i), k) = buoybase(i)
    end do

    ! Compute inb:

    inb = nl - 1

    do i = 1, ncum
       do k = 1, nl - 1
          if ((k >= icb(i)) .and. (buoy(i, k) < dtovsh)) then
             inb(i) = MIN(inb(i), k)
          endif
       end do
    end do

    ! CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL

    do k = 1, nl + 1
       do i = 1, ncum
          hp(i, k) = h(i, k)
       enddo
    enddo

    do k = minorig + 1, nl
       do i = 1, ncum
          if (k >= icb(i) .and. k <= inb(i)) hp(i, k) = h(i, nk(i)) &
               + (lv(i, k) + (cpd - cpv) * t(i, k)) * ep(i, k) * clw(i, k)
       end do
    end do

  end SUBROUTINE cv30_undilute2

end module cv30_undilute2_m
