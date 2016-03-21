module cv30_undilute2_m

  implicit none

contains

  SUBROUTINE cv30_undilute2(nloc, ncum, nd, icb, icbs, nk, tnk, qnk, gznk, t, &
       qs, gz, p, h, tv, lv, pbase, buoybase, plcl, inb, tp, tvp, clw, hp, &
       ep, sigp, buoy)

    ! Purpose: find the rest of the lifted parcel temperatures;
    ! compute the precipitation efficiencies and the fraction of
    ! precipitation falling outside of cloud; find the level of
    ! neutral buoyancy.

    ! Vertical profile of buoyancy computed here (use of buoybase).

    use conema3_m, only: epmax
    use cv30_param_m, only: dtovsh, minorig, nl, pbcrit, ptcrit, spfac
    use cvthermo, only: cl, clmcpv, cpd, cpv, eps, lv0, rrv

    ! inputs:
    integer, intent(in):: nloc, ncum, nd
    integer icb(nloc), icbs(nloc), nk(nloc)
    ! icbs (input) is the first level above LCL (may differ from icb)
    real tnk(nloc), qnk(nloc), gznk(nloc)
    real t(nloc, nd), qs(nloc, nd), gz(nloc, nd)
    real p(nloc, nd), h(nloc, nd)
    real tv(nloc, nd), lv(nloc, nd)
    real pbase(nloc), buoybase(nloc), plcl(nloc)

    ! outputs:
    integer inb(nloc)
    real tp(nloc, nd), tvp(nloc, nd), clw(nloc, nd)
    ! condensed water not removed from tvp
    real hp(nloc, nd), ep(nloc, nd), sigp(nloc, nd)
    real buoy(nloc, nd)

    ! Local:
    integer i, k
    real tg, qg, ahg, alv, s, tc, es, denom
    real pden
    real ah0(nloc)

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

    ! FIND THE FIRST MODEL LEVEL (INB) ABOVE THE PARCEL'S
    ! LEVEL OF NEUTRAL BUOYANCY

    do i = 1, ncum
       inb(i) = nl - 1
    end do

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
