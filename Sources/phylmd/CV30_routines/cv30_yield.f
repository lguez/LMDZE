module cv30_yield_m

  implicit none

contains

  SUBROUTINE cv30_yield(icb, inb, delt, t, rr, u, v, gz, p, ph, h, hp, lv, &
       cpn, th, ep, clw, m, tp, mp, rp, up, vp, wt, water, evap, b, ment, &
       qent, uent, vent, nent, elij, sig, tv, tvp, iflag, precip, VPrecip, &
       ft, fr, fu, fv, upwd, dnwd, dnwd0, ma, mike, tls, tps, qcondc)

    ! Tendencies, precipitation, variables of interface with other
    ! processes, etc.

    use conema3_m, only: iflag_clw
    use cv30_param_m, only: minorig, nl, sigd
    use cv_thermo_m, only: cl, cpd, cpv, rowl, rrd, rrv
    USE dimphy, ONLY: klev, klon
    use SUPHEC_M, only: rg

    ! inputs:
    integer, intent(in):: icb(:), inb(:) ! (ncum)
    real, intent(in):: delt
    real, intent(in):: t(klon, klev), rr(klon, klev)
    real, intent(in):: u(klon, klev), v(klon, klev)
    real gz(klon, klev)
    real p(klon, klev)
    real ph(klon, klev + 1), h(klon, klev), hp(klon, klev)
    real lv(klon, klev), cpn(klon, klev)
    real th(klon, klev)
    real ep(klon, klev), clw(klon, klev)
    real m(klon, klev)
    real tp(klon, klev)
    real mp(klon, klev), rp(klon, klev), up(klon, klev)
    real, intent(in):: vp(:, 2:) ! (ncum, 2:nl)
    real, intent(in):: wt(:, :) ! (ncum, nl - 1)
    real, intent(in):: water(:, :), evap(:, :) ! (ncum, nl)
    real, intent(in):: b(:, :) ! (ncum, nl - 1)
    real ment(klon, klev, klev), qent(klon, klev, klev), uent(klon, klev, klev)
    real vent(klon, klev, klev)
    integer nent(klon, klev)
    real elij(klon, klev, klev)
    real sig(klon, klev)
    real tv(klon, klev), tvp(klon, klev)

    ! outputs:
    integer, intent(out):: iflag(:) ! (ncum)
    real precip(klon)
    real VPrecip(klon, klev + 1)
    real ft(klon, klev), fr(klon, klev), fu(klon, klev), fv(klon, klev)
    real upwd(klon, klev), dnwd(klon, klev)
    real dnwd0(klon, klev)
    real ma(klon, klev)
    real mike(klon, klev)
    real tls(klon, klev), tps(klon, klev)
    real qcondc(klon, klev)

    ! Local:
    real, parameter:: delta = 0.01 ! interface cloud parameterization
    integer ncum
    integer i, k, il, n, j
    real awat, delti
    real ax, bx, cx, dx
    real cpinv, rdcp, dpinv
    real lvcp(klon, klev)
    real am(klon), work(klon), ad(klon), amp1(klon)
    real up1(klon, klev, klev), dn1(klon, klev, klev)
    real asum(klon), bsum(klon), csum(klon), dsum(klon)
    real qcond(klon, klev), nqcond(klon, klev), wa(klon, klev)
    real siga(klon, klev), sax(klon, klev), mac(klon, klev)

    !-------------------------------------------------------------

    ncum = size(icb)
    iflag = 0

    ! initialization:

    delti = 1.0 / delt

    do il = 1, ncum
       precip(il) = 0.0
       VPrecip(il, klev + 1) = 0.
    enddo

    do i = 1, klev
       do il = 1, ncum
          VPrecip(il, i) = 0.0
          ft(il, i) = 0.0
          fr(il, i) = 0.0
          fu(il, i) = 0.0
          fv(il, i) = 0.0
          qcondc(il, i) = 0.0
          qcond(il, i) = 0.0
          nqcond(il, i) = 0.0
       enddo
    enddo

    do i = 1, nl
       do il = 1, ncum
          lvcp(il, i) = lv(il, i) / cpn(il, i)
       enddo
    enddo

    ! calculate surface precipitation in mm / day

    do il = 1, ncum
       if (ep(il, inb(il)) >= 1e-4) precip(il) = wt(il, 1) * sigd &
            * water(il, 1) * 86400. * 1000. / (rowl * rg)
    enddo

    ! CALCULATE VERTICAL PROFILE OF PRECIPITATIONs IN kg / m2 / s ===

    ! MAF rajout pour lessivage
    do k = 1, nl - 1
       do il = 1, ncum
          if (k <= inb(il)) VPrecip(il, k) = wt(il, k) * sigd * water(il, k) &
               / rg
       end do
    end do

    ! calculate tendencies of lowest level potential temperature
    ! and mixing ratio

    do il = 1, ncum
       work(il) = 1.0 / (ph(il, 1) - ph(il, 2))
       am(il) = 0.0
    enddo

    do k = 2, nl
       do il = 1, ncum
          if (k <= inb(il)) am(il) = am(il) + m(il, k)
       enddo
    enddo

    do il = 1, ncum
       if (0.01 * rg * work(il) * am(il) >= delti) iflag(il) = 1

       ft(il, 1) = 0.01 * rg * work(il) * am(il) * (t(il, 2) - t(il, 1) &
            + (gz(il, 2) - gz(il, 1)) / cpn(il, 1)) - 0.5 * lvcp(il, 1) &
            * sigd * (evap(il, 1) + evap(il, 2)) - 0.009 * rg * sigd &
            * mp(il, 2) * t(il, 1) * b(il, 1) * work(il) + 0.01 * sigd &
            * wt(il, 1) * (cl - cpd) * water(il, 2) * (t(il, 2) - t(il, 1)) &
            * work(il) / cpn(il, 1)

       !jyg1 Correction pour mieux conserver l'eau (conformite avec CONVECT4.3)
       ! (sb: pour l'instant, on ne fait que le chgt concernant rg, pas evap)
       fr(il, 1) = 0.01 * rg * mp(il, 2) * (rp(il, 2) - rr(il, 1)) &
            * work(il) + sigd * 0.5 * (evap(il, 1) + evap(il, 2))
       ! + tard : + sigd * evap(il, 1)

       fr(il, 1) = fr(il, 1) + 0.01 * rg * am(il) * (rr(il, 2) - rr(il, 1)) &
            * work(il)

       fu(il, 1) = fu(il, 1) + 0.01 * rg * work(il) * (mp(il, 2) &
            * (up(il, 2) - u(il, 1)) + am(il) * (u(il, 2) - u(il, 1)))
       fv(il, 1) = fv(il, 1) + 0.01 * rg * work(il) * (mp(il, 2) &
            * (vp(il, 2) - v(il, 1)) + am(il) * (v(il, 2) - v(il, 1)))
    enddo

    do j = 2, nl
       do il = 1, ncum
          if (j <= inb(il)) then
             fr(il, 1) = fr(il, 1) + 0.01 * rg * work(il) * ment(il, j, 1) &
                  * (qent(il, j, 1) - rr(il, 1))
             fu(il, 1) = fu(il, 1) + 0.01 * rg * work(il) * ment(il, j, 1) &
                  * (uent(il, j, 1) - u(il, 1))
             fv(il, 1) = fv(il, 1) + 0.01 * rg * work(il) * ment(il, j, 1) &
                  * (vent(il, j, 1) - v(il, 1))
          endif
       enddo
    enddo

    ! calculate tendencies of potential temperature and mixing ratio
    ! at levels above the lowest level

    ! first find the net saturated updraft and downdraft mass fluxes
    ! through each level

    loop_i: do i = 2, nl - 1
       if (any(inb >= i)) then
          amp1(:ncum) = 0.
          ad(:ncum) = 0.

          do k = i + 1, nl + 1
             do il = 1, ncum
                if (i <= inb(il) .and. k <= (inb(il) + 1)) then
                   amp1(il) = amp1(il) + m(il, k)
                endif
             end do
          end do

          do k = 1, i
             do j = i + 1, nl + 1
                do il = 1, ncum
                   if (i <= inb(il) .and. j <= (inb(il) + 1)) then
                      amp1(il) = amp1(il) + ment(il, k, j)
                   endif
                end do
             end do
          end do

          do k = 1, i - 1
             do j = i, nl + 1 ! newvecto: nl au lieu nl + 1?
                do il = 1, ncum
                   if (i <= inb(il) .and. j <= inb(il)) then
                      ad(il) = ad(il) + ment(il, j, k)
                   endif
                end do
             end do
          end do

          do il = 1, ncum
             if (i <= inb(il)) then
                dpinv = 1.0 / (ph(il, i) - ph(il, i + 1))
                cpinv = 1.0 / cpn(il, i)

                if (0.01 * rg * dpinv * amp1(il) >= delti) iflag(il) = 1

                ft(il, i) = 0.01 * rg * dpinv * (amp1(il) * (t(il, i + 1) &
                     - t(il, i) + (gz(il, i + 1) - gz(il, i)) * cpinv) &
                     - ad(il) * (t(il, i) - t(il, i - 1) + (gz(il, i) &
                     - gz(il, i - 1)) * cpinv)) - 0.5 * sigd * lvcp(il, i) &
                     * (evap(il, i) + evap(il, i + 1)) - 0.009 * rg * sigd &
                     * (mp(il, i + 1) * t(il, i) * b(il, i) - mp(il, i) &
                     * t(il, i - 1) * cpn(il, i - 1) * cpinv * b(il, i - 1)) &
                     * dpinv + 0.01 * rg * dpinv * ment(il, i, i) &
                     * (hp(il, i) - h(il, i) + t(il, i) * (cpv - cpd) &
                     * (rr(il, i) - qent(il, i, i))) * cpinv + 0.01 * sigd &
                     * wt(il, i) * (cl - cpd) * water(il, i + 1) &
                     * (t(il, i + 1) - t(il, i)) * dpinv * cpinv
                fr(il, i) = 0.01 * rg * dpinv * (amp1(il) * (rr(il, i + 1) &
                     - rr(il, i)) - ad(il) * (rr(il, i) - rr(il, i - 1)))
                fu(il, i) = fu(il, i) + 0.01 * rg * dpinv * (amp1(il) &
                     * (u(il, i + 1) - u(il, i)) - ad(il) * (u(il, i) &
                     - u(il, i - 1)))
                fv(il, i) = fv(il, i) + 0.01 * rg * dpinv * (amp1(il) &
                     * (v(il, i + 1) - v(il, i)) - ad(il) * (v(il, i) &
                     - v(il, i - 1)))
             endif
          end do

          do k = 1, i - 1
             do il = 1, ncum
                if (i <= inb(il)) then
                   dpinv = 1.0 / (ph(il, i) - ph(il, i + 1))
                   cpinv = 1.0 / cpn(il, i)

                   awat = elij(il, k, i) - (1. - ep(il, i)) * clw(il, i)
                   awat = amax1(awat, 0.0)

                   fr(il, i) = fr(il, i) + 0.01 * rg * dpinv &
                        * ment(il, k, i) * (qent(il, k, i) - awat - rr(il, i))
                   fu(il, i) = fu(il, i) + 0.01 * rg * dpinv &
                        * ment(il, k, i) * (uent(il, k, i) - u(il, i))
                   fv(il, i) = fv(il, i) + 0.01 * rg * dpinv &
                        * ment(il, k, i) * (vent(il, k, i) - v(il, i))

                   ! (saturated updrafts resulting from mixing)
                   qcond(il, i) = qcond(il, i) + (elij(il, k, i) - awat)
                   nqcond(il, i) = nqcond(il, i) + 1.
                endif ! i
             end do
          end do

          do k = i, nl + 1
             do il = 1, ncum
                if (i <= inb(il) .and. k <= inb(il)) then
                   dpinv = 1.0 / (ph(il, i) - ph(il, i + 1))
                   cpinv = 1.0 / cpn(il, i)

                   fr(il, i) = fr(il, i) + 0.01 * rg * dpinv &
                        * ment(il, k, i) * (qent(il, k, i) - rr(il, i))
                   fu(il, i) = fu(il, i) + 0.01 * rg * dpinv &
                        * ment(il, k, i) * (uent(il, k, i) - u(il, i))
                   fv(il, i) = fv(il, i) + 0.01 * rg * dpinv &
                        * ment(il, k, i) * (vent(il, k, i) - v(il, i))
                endif
             end do
          end do

          do il = 1, ncum
             if (i <= inb(il)) then
                dpinv = 1.0 / (ph(il, i) - ph(il, i + 1))
                cpinv = 1.0 / cpn(il, i)

                ! sb: on ne fait pas encore la correction permettant de mieux
                ! conserver l'eau:
                fr(il, i) = fr(il, i) + 0.5 * sigd * (evap(il, i) &
                     + evap(il, i + 1)) + 0.01 * rg * (mp(il, i + 1) &
                     * (rp(il, i + 1) - rr(il, i)) - mp(il, i) * (rp(il, i) &
                     - rr(il, i - 1))) * dpinv

                fu(il, i) = fu(il, i) + 0.01 * rg * (mp(il, i + 1) &
                     * (up(il, i + 1) - u(il, i)) - mp(il, i) * (up(il, i) &
                     - u(il, i - 1))) * dpinv
                fv(il, i) = fv(il, i) + 0.01 * rg * (mp(il, i + 1) &
                     * (vp(il, i + 1) - v(il, i)) - mp(il, i) * (vp(il, i) &
                     - v(il, i - 1))) * dpinv
             endif
          end do

          ! sb: interface with the cloud parameterization:

          do k = i + 1, nl
             do il = 1, ncum
                if (k <= inb(il) .and. i <= inb(il)) then
                   ! (saturated downdrafts resulting from mixing)
                   qcond(il, i) = qcond(il, i) + elij(il, k, i)
                   nqcond(il, i) = nqcond(il, i) + 1.
                endif
             enddo
          enddo

          ! (particular case: no detraining level is found)
          do il = 1, ncum
             if (i <= inb(il) .and. nent(il, i) == 0) then
                qcond(il, i) = qcond(il, i) + (1. - ep(il, i)) * clw(il, i)
                nqcond(il, i) = nqcond(il, i) + 1.
             endif
          enddo

          do il = 1, ncum
             if (i <= inb(il) .and. nqcond(il, i) /= 0.) then
                qcond(il, i) = qcond(il, i) / nqcond(il, i)
             endif
          enddo
       end if
    end do loop_i

    ! move the detrainment at level inb down to level inb - 1
    ! in such a way as to preserve the vertically
    ! integrated enthalpy and water tendencies

    do il = 1, ncum
       ax = 0.1 * ment(il, inb(il), inb(il)) * (hp(il, inb(il)) &
            - h(il, inb(il)) + t(il, inb(il)) * (cpv - cpd) &
            * (rr(il, inb(il)) - qent(il, inb(il), inb(il)))) &
            / (cpn(il, inb(il)) * (ph(il, inb(il)) - ph(il, inb(il) + 1)))
       ft(il, inb(il)) = ft(il, inb(il)) - ax
       ft(il, inb(il) - 1) = ft(il, inb(il) - 1) + ax * cpn(il, inb(il)) &
            * (ph(il, inb(il)) - ph(il, inb(il) + 1)) / (cpn(il, inb(il) - 1) &
            * (ph(il, inb(il) - 1) - ph(il, inb(il))))

       bx = 0.1 * ment(il, inb(il), inb(il)) * (qent(il, inb(il), inb(il)) &
            - rr(il, inb(il))) / (ph(il, inb(il)) - ph(il, inb(il) + 1))
       fr(il, inb(il)) = fr(il, inb(il)) - bx
       fr(il, inb(il) - 1) = fr(il, inb(il) - 1) &
            + bx * (ph(il, inb(il)) - ph(il, inb(il) + 1)) &
            / (ph(il, inb(il) - 1) - ph(il, inb(il)))

       cx = 0.1 * ment(il, inb(il), inb(il)) * (uent(il, inb(il), inb(il)) &
            - u(il, inb(il))) / (ph(il, inb(il)) - ph(il, inb(il) + 1))
       fu(il, inb(il)) = fu(il, inb(il)) - cx
       fu(il, inb(il) - 1) = fu(il, inb(il) - 1) &
            + cx * (ph(il, inb(il)) - ph(il, inb(il) + 1)) &
            / (ph(il, inb(il) - 1) - ph(il, inb(il)))

       dx = 0.1 * ment(il, inb(il), inb(il)) * (vent(il, inb(il), inb(il)) &
            - v(il, inb(il))) / (ph(il, inb(il)) - ph(il, inb(il) + 1))
       fv(il, inb(il)) = fv(il, inb(il)) - dx
       fv(il, inb(il) - 1) = fv(il, inb(il) - 1) &
            + dx * (ph(il, inb(il)) - ph(il, inb(il) + 1)) &
            / (ph(il, inb(il) - 1) - ph(il, inb(il)))

    end do

    ! homoginize tendencies below cloud base

    do il = 1, ncum
       asum(il) = 0.0
       bsum(il) = 0.0
       csum(il) = 0.0
       dsum(il) = 0.0
    enddo

    do i = 1, nl
       do il = 1, ncum
          if (i <= (icb(il) - 1)) then
             asum(il) = asum(il) + ft(il, i) * (ph(il, i) - ph(il, i + 1))
             bsum(il) = bsum(il) + fr(il, i) * (lv(il, i) + (cl - cpd) &
                  * (t(il, i) - t(il, 1))) * (ph(il, i) - ph(il, i + 1))
             csum(il) = csum(il) + (lv(il, i) + (cl - cpd) * (t(il, i) &
                  - t(il, 1))) * (ph(il, i) - ph(il, i + 1))
             dsum(il) = dsum(il) + t(il, i) * (ph(il, i) - ph(il, i + 1)) &
                  / th(il, i)
          endif
       enddo
    enddo

    do i = 1, nl
       do il = 1, ncum
          if (i <= (icb(il) - 1)) then
             ft(il, i) = asum(il) * t(il, i) / (th(il, i) * dsum(il))
             fr(il, i) = bsum(il) / csum(il)
          endif
       enddo
    enddo

    ! reset counter and return

    do il = 1, ncum
       sig(il, klev) = 2.0
    enddo

    do i = 1, klev
       do il = 1, ncum
          upwd(il, i) = 0.0
          dnwd(il, i) = 0.0
       enddo
    enddo

    do i = 1, nl
       do il = 1, ncum
          dnwd0(il, i) = - mp(il, i)
       enddo
    enddo
    do i = nl + 1, klev
       do il = 1, ncum
          dnwd0(il, i) = 0.
       enddo
    enddo

    do i = 1, nl
       do il = 1, ncum
          if (i >= icb(il) .and. i <= inb(il)) then
             upwd(il, i) = 0.0
             dnwd(il, i) = 0.0
          endif
       enddo
    enddo

    do i = 1, nl
       do k = 1, nl
          do il = 1, ncum
             up1(il, k, i) = 0.0
             dn1(il, k, i) = 0.0
          enddo
       enddo
    enddo

    do i = 1, nl
       do k = i, nl
          do n = 1, i - 1
             do il = 1, ncum
                if (i >= icb(il).and.i <= inb(il).and.k <= inb(il)) then
                   up1(il, k, i) = up1(il, k, i) + ment(il, n, k)
                   dn1(il, k, i) = dn1(il, k, i) - ment(il, k, n)
                endif
             enddo
          enddo
       enddo
    enddo

    do i = 2, nl
       do k = i, nl
          do il = 1, ncum
             if (i <= inb(il).and.k <= inb(il)) then
                upwd(il, i) = upwd(il, i) + m(il, k) + up1(il, k, i)
                dnwd(il, i) = dnwd(il, i) + dn1(il, k, i)
             endif
          enddo
       enddo
    enddo

    ! D\'etermination de la variation de flux ascendant entre
    ! deux niveaux non dilu\'es mike

    do i = 1, nl
       do il = 1, ncum
          mike(il, i) = m(il, i)
       enddo
    enddo

    do i = nl + 1, klev
       do il = 1, ncum
          mike(il, i) = 0.
       enddo
    enddo

    do i = 1, klev
       do il = 1, ncum
          ma(il, i) = 0
       enddo
    enddo

    do i = 1, nl
       do j = i, nl
          do il = 1, ncum
             ma(il, i) = ma(il, i) + m(il, j)
          enddo
       enddo
    enddo

    do i = nl + 1, klev
       do il = 1, ncum
          ma(il, i) = 0.
       enddo
    enddo

    do i = 1, nl
       do il = 1, ncum
          if (i <= (icb(il) - 1)) then
             ma(il, i) = 0
          endif
       enddo
    enddo

    ! icb repr\'esente le niveau o\`u se trouve la base du nuage, et
    ! inb le sommet du nuage

    do i = 1, klev
       DO il = 1, ncum
          rdcp = (rrd * (1. - rr(il, i)) - rr(il, i) * rrv) &
               / (cpd * (1. - rr(il, i)) + rr(il, i) * cpv)
          tls(il, i) = t(il, i) * (1000.0 / p(il, i))**rdcp
          tps(il, i) = tp(il, i)
       end DO
    enddo

    ! Diagnose the in-cloud mixing ratio of condensed water

    do i = 1, klev
       do il = 1, ncum
          mac(il, i) = 0.0
          wa(il, i) = 0.0
          siga(il, i) = 0.0
          sax(il, i) = 0.0
       enddo
    enddo

    do i = minorig, nl
       do k = i + 1, nl + 1
          do il = 1, ncum
             if (i <= inb(il) .and. k <= (inb(il) + 1)) then
                mac(il, i) = mac(il, i) + m(il, k)
             endif
          enddo
       enddo
    enddo

    do i = 1, nl
       do j = 1, i
          do il = 1, ncum
             if (i >= icb(il) .and. i <= (inb(il) - 1) &
                  .and. j >= icb(il)) then
                sax(il, i) = sax(il, i) + rrd * (tvp(il, j) - tv(il, j)) &
                     * (ph(il, j) - ph(il, j + 1)) / p(il, j)
             endif
          enddo
       enddo
    enddo

    do i = 1, nl
       do il = 1, ncum
          if (i >= icb(il) .and. i <= (inb(il) - 1) &
               .and. sax(il, i) > 0.0) then
             wa(il, i) = sqrt(2. * sax(il, i))
          endif
       enddo
    enddo

    do i = 1, nl
       do il = 1, ncum
          if (wa(il, i) > 0.0) siga(il, i) = mac(il, i) / wa(il, i) * rrd &
               * tvp(il, i) / p(il, i) / 100. / delta
          siga(il, i) = min(siga(il, i), 1.0)

          if (iflag_clw == 0) then
             qcondc(il, i) = siga(il, i) * clw(il, i) * (1. - ep(il, i)) &
                  + (1. - siga(il, i)) * qcond(il, i)
          else if (iflag_clw == 1) then
             qcondc(il, i) = qcond(il, i)
          endif
       enddo
    enddo

  end SUBROUTINE cv30_yield

end module cv30_yield_m
