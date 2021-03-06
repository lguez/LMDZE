module cv30_undilute2_m

  implicit none

contains

  SUBROUTINE cv30_undilute2(icb, icbs, tnk, qnk, gznk, t, qs, gz, p, h, tv, &
       lv, pbase, buoybase, plcl, inb, tp, tvp, clw, hp, ep, buoy)

    ! Undilute (adiabatic) updraft, second part. Purpose: find the
    ! rest of the lifted parcel temperatures; compute the
    ! precipitation efficiencies and the fraction of precipitation
    ! falling outside of cloud; find the level of neutral buoyancy.

    ! Vertical profile of buoyancy computed here (use of buoybase).

    use conf_phys_m, only: epmax
    use cv30_param_m, only: minorig, nl
    use cv_thermo, only: clmcpv, eps
    USE dimphy, ONLY: klon, klev
    use SUPHEC_M, only: rcw, rlvtt, rcpd, rcpv, rv

    integer, intent(in):: icb(:) ! (ncum) {2 <= icb <= nl - 3}

    integer, intent(in):: icbs(:) ! (ncum)
    ! icbs is the first level above LCL (may differ from icb)

    real, intent(in):: tnk(:), qnk(:), gznk(:) ! (klon)
    real, intent(in):: t(klon, klev), qs(klon, klev), gz(klon, klev)
    real, intent(in):: p(klon, klev), h(klon, klev)
    real, intent(in):: tv(klon, klev)
    real, intent(in):: lv(:, :) ! (ncum, nl)
    real, intent(in):: pbase(:), buoybase(:), plcl(:) ! (ncum)

    ! outputs:
    integer, intent(out):: inb(:) ! (ncum)
    ! first model level above the level of neutral buoyancy of the
    ! parcel (1 <= inb <= nl - 1)

    real tp(klon, klev), tvp(klon, klev), clw(klon, klev)
    ! condensed water not removed from tvp
    real hp(klon, klev), ep(klon, klev)
    real buoy(klon, klev)

    ! Local:

    integer ncum

    real, parameter:: pbcrit = 150.
    ! critical cloud depth (mbar) beneath which the precipitation
    ! efficiency is assumed to be zero

    real, parameter:: ptcrit = 500.
    ! cloud depth (mbar) above which the precipitation efficiency is
    ! assumed to be unity

    real, parameter:: dtovsh = - 0.2 ! dT for overshoot

    integer i, k
    real tg, qg, ahg, alv, s, tc, es, denom
    real pden
    real ah0(klon)

    !---------------------------------------------------------------------

    ncum = size(icb)

    ! SOME INITIALIZATIONS

    do k = 1, nl
       do i = 1, ncum
          ep(i, k) = 0.
       end do
    end do

    ! FIND THE REST OF THE LIFTED PARCEL TEMPERATURES

    ! The procedure is to solve the equation.
    ! cp * tp + L * qp + phi = cp * tnk + L * qnk + gznk.

    ! Calculate certain parcel quantities, including static energy

    do i = 1, ncum
       ah0(i) = (rcpd * (1. - qnk(i)) + rcw * qnk(i)) * tnk(i) &
            + qnk(i) * (rlvtt - clmcpv * (tnk(i) - 273.15)) + gznk(i)
    end do

    ! Find lifted parcel quantities above cloud base

    do k = minorig + 1, nl
       do i = 1, ncum
          if (k >= (icbs(i) + 1)) then
             tg = t(i, k)
             qg = qs(i, k)
             alv = rlvtt - clmcpv * (t(i, k) - 273.15)

             ! First iteration.

             s = rcpd * (1. - qnk(i)) + rcw * qnk(i) &
                  + alv * alv * qg / (rv * t(i, k) * t(i, k))
             s = 1. / s

             ahg = rcpd * tg + (rcw - rcpd) * qnk(i) * tg + alv * qg + gz(i, k)
             tg = tg + s * (ah0(i) - ahg)

             tc = tg - 273.15
             denom = 243.5 + tc
             denom = MAX(denom, 1.)

             es = 6.112 * exp(17.67 * tc / denom)

             qg = eps * es / (p(i, k) - es * (1. - eps))

             ! Second iteration.

             ahg = rcpd * tg + (rcw - rcpd) * qnk(i) * tg + alv * qg + gz(i, k)
             tg = tg + s * (ah0(i) - ahg)

             tc = tg - 273.15
             denom = 243.5 + tc
             denom = MAX(denom, 1.)

             es = 6.112 * exp(17.67 * tc / denom)

             qg = eps * es / (p(i, k) - es * (1. - eps))

             alv = rlvtt - clmcpv * (t(i, k) - 273.15)

             ! no approximation:
             tp(i, k) = (ah0(i) - gz(i, k) - alv * qg) &
                  / (rcpd + (rcw - rcpd) * qnk(i))

             clw(i, k) = qnk(i) - qg
             clw(i, k) = max(0., clw(i, k))
             ! qg utilise au lieu du vrai mixing ratio rg:
             tvp(i, k) = tp(i, k) * (1. + qg / eps - qnk(i)) ! whole thing
          endif
       end do
    end do

    ! SET THE PRECIPITATION EFFICIENCIES
    ! It MAY BE a FUNCTION OF TP(I), P(I) AND CLW(I)
    do k = 1, nl
       do i = 1, ncum
          pden = ptcrit - pbcrit
          ep(i, k) = (plcl(i) - p(i, k) - pbcrit) / pden * epmax
          ep(i, k) = max(ep(i, k), 0.)
          ep(i, k) = min(ep(i, k), epmax)
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
          if (k >= icb(i) .and. k <= inb(i)) hp(i, k) = h(i, minorig) &
               + (lv(i, k) + (rcpd - rcpv) * t(i, k)) * ep(i, k) * clw(i, k)
       end do
    end do

  end SUBROUTINE cv30_undilute2

end module cv30_undilute2_m
