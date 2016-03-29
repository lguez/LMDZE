module cv30_unsat_m

  implicit none

contains

  SUBROUTINE cv30_unsat(icb, inb, t, q, qs, gz, u, v, p, ph, th, tv, lv, cpn, &
       ep, sigp, clw, m, ment, elij, delt, plcl, mp, qp, up, vp, wt, water, &
       evap, b)

    use cv30_param_m, only: nl, sigd
    use cvthermo, only: cpd, ginv, grav
    USE dimphy, ONLY: klon, klev

    ! inputs:
    integer, intent(in):: icb(:), inb(:) ! (ncum)
    real, intent(in):: t(:, :), q(:, :), qs(:, :) ! (klon, klev)
    real, intent(in):: gz(:, :) ! (klon, klev)
    real, intent(in):: u(:, :), v(:, :) ! (klon, klev)
    real p(klon, klev), ph(klon, klev + 1)
    real th(klon, klev)
    real tv(klon, klev)
    real lv(klon, klev)
    real cpn(klon, klev)
    real, intent(in):: ep(klon, klev), sigp(klon, klev), clw(klon, klev)
    real m(klon, klev), ment(klon, klev, klev), elij(klon, klev, klev)
    real, intent(in):: delt
    real plcl(klon)

    ! outputs:
    real, intent(out):: mp(klon, klev)
    real, intent(out):: qp(:, :), up(:, :), vp(:, :) ! (ncum, nl)
    real wt(klon, klev), water(klon, klev), evap(klon, klev)
    real, intent(out):: b(:, :) ! (ncum, nl - 1)

    ! Local:
    integer ncum
    integer i, j, il, num1
    real tinv, delti
    real awat, afac, afac1, afac2, bfac
    real pr1, pr2, sigt, b6, c6, revap, tevap, delth
    real amfac, amp2, xf, tf, fac2, ur, sru, fac, d, af, bf
    real ampmax
    real lvcp(klon, klev)
    real wdtrain(size(icb))
    logical lwork(size(icb))

    !------------------------------------------------------

    ncum = size(icb)
    delti = 1. / delt
    tinv = 1. / 3.
    mp = 0.
    b = 0.

    do i = 1, nl
       do il = 1, ncum
          qp(il, i) = q(il, i)
          up(il, i) = u(il, i)
          vp(il, i) = v(il, i)
          wt(il, i) = 0.001
          water(il, i) = 0.
          evap(il, i) = 0.
          lvcp(il, i) = lv(il, i) / cpn(il, i)
       enddo
    enddo

    ! check whether ep(inb) = 0, if so, skip precipitating 
    ! downdraft calculation 
    forall (il = 1:ncum) lwork(il) = ep(il, inb(il)) >= 1e-4

    wdtrain = 0.

    downdraft_loop: DO i = nl - 1, 1, - 1
       num1 = 0

       do il = 1, ncum
          if (i <= inb(il) .and. lwork(il)) num1 = num1 + 1
       enddo

       if (num1 > 0) then
          ! integrate liquid water equation to find condensed water 
          ! and condensed water flux 

          ! calculate detrained precipitation 

          do il = 1, ncum
             if (i <= inb(il) .and. lwork(il)) then
                wdtrain(il) = grav * ep(il, i) * m(il, i) * clw(il, i)
             endif
          enddo

          if (i > 1) then
             do j = 1, i - 1
                do il = 1, ncum
                   if (i <= inb(il) .and. lwork(il)) then
                      awat = elij(il, j, i) - (1. - ep(il, i)) * clw(il, i)
                      awat = max(awat, 0.)
                      wdtrain(il) = wdtrain(il) + grav * awat * ment(il, j, i)
                   endif
                enddo
             end do
          endif

          ! find rain water and evaporation using provisional 
          ! estimates of qp(i) and qp(i - 1) 

          do il = 1, ncum
             if (i <= inb(il) .and. lwork(il)) then
                wt(il, i) = 45.

                if (i < inb(il)) then
                   qp(il, i) = qp(il, i + 1) + (cpd * (t(il, i + 1) &
                        - t(il, i)) + gz(il, i + 1) - gz(il, i)) / lv(il, i)
                   qp(il, i) = 0.5 * (qp(il, i) + q(il, i))
                endif

                qp(il, i) = max(qp(il, i), 0.)
                qp(il, i) = min(qp(il, i), qs(il, i))
                qp(il, inb(il)) = q(il, inb(il))

                if (i == 1) then
                   afac = p(il, 1) * (qs(il, 1) - qp(il, 1)) &
                        / (1e4 + 2000. * p(il, 1) * qs(il, 1))
                else
                   qp(il, i - 1) = qp(il, i) + (cpd * (t(il, i) &
                        - t(il, i - 1)) + gz(il, i) - gz(il, i - 1)) / lv(il, i)
                   qp(il, i - 1) = 0.5 * (qp(il, i - 1) + q(il, i - 1))
                   qp(il, i - 1) = min(qp(il, i - 1), qs(il, i - 1))
                   qp(il, i - 1) = max(qp(il, i - 1), 0.)
                   afac1 = p(il, i) * (qs(il, i) - qp(il, i)) &
                        / (1e4 + 2000. * p(il, i) * qs(il, i))
                   afac2 = p(il, i - 1) * (qs(il, i - 1) - qp(il, i - 1)) &
                        / (1e4 + 2000. * p(il, i - 1) * qs(il, i - 1))
                   afac = 0.5 * (afac1 + afac2)
                endif

                if (i == inb(il)) afac = 0.
                afac = max(afac, 0.)
                bfac = 1. / (sigd * wt(il, i))

                ! prise en compte de la variation progressive de sigt dans
                ! les couches icb et icb - 1:
                ! pour plcl < ph(i + 1), pr1 = 0 & pr2 = 1
                ! pour plcl > ph(i), pr1 = 1 & pr2 = 0
                ! pour ph(i + 1) < plcl < ph(i), pr1 est la proportion a cheval
                ! sur le nuage, et pr2 est la proportion sous la base du
                ! nuage.
                pr1 = (plcl(il) - ph(il, i + 1)) / (ph(il, i) - ph(il, i + 1))
                pr1 = max(0., min(1., pr1))
                pr2 = (ph(il, i) - plcl(il)) / (ph(il, i) - ph(il, i + 1))
                pr2 = max(0., min(1., pr2))
                sigt = sigp(il, i) * pr1 + pr2

                b6 = bfac * 50. * sigd * (ph(il, i) - ph(il, i + 1)) * sigt &
                     * afac
                c6 = water(il, i + 1) + bfac * wdtrain(il) - 50. * sigd * bfac &
                     * (ph(il, i) - ph(il, i + 1)) * evap(il, i + 1)
                if (c6 > 0.) then
                   revap = 0.5 * (- b6 + sqrt(b6 * b6 + 4. * c6))
                   evap(il, i) = sigt * afac * revap
                   water(il, i) = revap * revap
                else
                   evap(il, i) = - evap(il, i + 1) + 0.02 * (wdtrain(il) &
                        + sigd * wt(il, i) * water(il, i + 1)) &
                        / (sigd * (ph(il, i) - ph(il, i + 1)))
                end if

                ! calculate precipitating downdraft mass flux under 
                ! hydrostatic approximation 

                if (i /= 1) then
                   tevap = max(0., evap(il, i))
                   delth = max(0.001, (th(il, i) - th(il, i - 1)))
                   mp(il, i) = 100. * ginv * lvcp(il, i) * sigd * tevap &
                        * (p(il, i - 1) - p(il, i)) / delth

                   ! If hydrostatic assumption fails, solve cubic
                   ! difference equation for downdraft theta and mass
                   ! flux from two simultaneous differential equations

                   amfac = sigd * sigd * 70. * ph(il, i) &
                        * (p(il, i - 1) - p(il, i)) &
                        * (th(il, i) - th(il, i - 1)) / (tv(il, i) * th(il, i))
                   amp2 = abs(mp(il, i + 1) * mp(il, i + 1) - mp(il, i) &
                        * mp(il, i))

                   if (amp2 > 0.1 * amfac) then
                      xf = 100. * sigd * sigd * sigd * (ph(il, i) &
                           - ph(il, i + 1))
                      tf = b(il, i) - 5. * (th(il, i) - th(il, i - 1)) &
                           * t(il, i) / (lvcp(il, i) * sigd * th(il, i))
                      af = xf * tf + mp(il, i + 1) * mp(il, i + 1) * tinv
                      bf = 2. * (tinv * mp(il, i + 1))**3 + tinv &
                           * mp(il, i + 1) * xf * tf + 50. * (p(il, i - 1) &
                           - p(il, i)) * xf * tevap
                      fac2 = 1.
                      if (bf < 0.) fac2 = - 1.
                      bf = abs(bf)
                      ur = 0.25 * bf * bf - af * af * af * tinv * tinv * tinv

                      if (ur >= 0.) then
                         sru = sqrt(ur)
                         fac = 1.
                         if ((0.5 * bf - sru) < 0.) fac = - 1.
                         mp(il, i) = mp(il, i + 1) * tinv &
                              + (0.5 * bf + sru)**tinv &
                              + fac * (abs(0.5 * bf - sru))**tinv
                      else
                         d = atan(2. * sqrt(- ur) / (bf + 1e-28))
                         if (fac2 < 0.)d = 3.14159 - d
                         mp(il, i) = mp(il, i + 1) * tinv + 2. &
                              * sqrt(af * tinv) * cos(d * tinv)
                      endif

                      mp(il, i) = max(0., mp(il, i))

                      ! Il y a vraisemblablement une erreur dans la
                      ! ligne suivante : il faut diviser par (mp(il,
                      ! i) * sigd * grav) et non par (mp(il, i) + sigd
                      ! * 0.1).  Et il faut bien revoir les facteurs
                      ! 100.
                      b(il, i - 1) = b(il, i) + 100. * (p(il, i - 1) &
                           - p(il, i)) * tevap / (mp(il, i) + sigd * 0.1) &
                           - 10. * (th(il, i) - th(il, i - 1)) * t(il, i) &
                           / (lvcp(il, i) * sigd * th(il, i))
                      b(il, i - 1) = max(b(il, i - 1), 0.)
                   endif

                   ! limit magnitude of mp(i) to meet cfl condition 
                   ampmax = 2. * (ph(il, i) - ph(il, i + 1)) * delti
                   amp2 = 2. * (ph(il, i - 1) - ph(il, i)) * delti
                   ampmax = min(ampmax, amp2)
                   mp(il, i) = min(mp(il, i), ampmax)

                   ! force mp to decrease linearly to zero 
                   ! between cloud base and the surface 
                   if (p(il, i) > p(il, icb(il))) mp(il, i) = mp(il, icb(il)) &
                        * (p(il, 1) - p(il, i)) / (p(il, 1) - p(il, icb(il)))
                endif ! i == 1

                ! find mixing ratio of precipitating downdraft 

                if (i /= inb(il)) then
                   qp(il, i) = q(il, i)

                   if (mp(il, i) > mp(il, i + 1)) then
                      qp(il, i) = qp(il, i + 1) * mp(il, i + 1) + q(il, i) &
                           * (mp(il, i) - mp(il, i + 1)) + 100. * ginv &
                           * 0.5 * sigd * (ph(il, i) - ph(il, i + 1)) &
                           * (evap(il, i + 1) + evap(il, i))
                      qp(il, i) = qp(il, i) / mp(il, i)
                      up(il, i) = up(il, i + 1) * mp(il, i + 1) + u(il, i) &
                           * (mp(il, i) - mp(il, i + 1))
                      up(il, i) = up(il, i) / mp(il, i)
                      vp(il, i) = vp(il, i + 1) * mp(il, i + 1) + v(il, i) &
                           * (mp(il, i) - mp(il, i + 1))
                      vp(il, i) = vp(il, i) / mp(il, i)
                   else
                      if (mp(il, i + 1) > 1e-16) then
                         qp(il, i) = qp(il, i + 1) + 100. * ginv * 0.5 * sigd &
                              * (ph(il, i) - ph(il, i + 1)) &
                              * (evap(il, i + 1) + evap(il, i)) / mp(il, i + 1)
                         up(il, i) = up(il, i + 1)
                         vp(il, i) = vp(il, i + 1)
                      endif
                   endif

                   qp(il, i) = min(qp(il, i), qs(il, i))
                   qp(il, i) = max(qp(il, i), 0.)
                endif
             endif
          end do
       end if
    end DO downdraft_loop

  end SUBROUTINE cv30_unsat

end module cv30_unsat_m
