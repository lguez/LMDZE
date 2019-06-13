module cv30_unsat_m

  implicit none

contains

  SUBROUTINE cv30_unsat(icb, inb, t, q, qs, gz, u, v, p, ph, th, tv, lv, cpn, &
       ep, clw, m, ment, elij, delt, plcl, mp, qp, up, vp, wt, water, evap, b)

    ! Unsaturated (precipitating) downdrafts

    use cv30_param_m, only: nl, sigd
    use cv_thermo, only: ginv
    use SUPHEC_M, only: rg, rcpd

    integer, intent(in):: icb(:) ! (ncum)
    ! {2 <= icb <= nl - 3}

    integer, intent(in):: inb(:) ! (ncum)
    ! first model level above the level of neutral buoyancy of the
    ! parcel (1 <= inb <= nl - 1)

    real, intent(in):: t(:, :) ! (ncum, nl) temperature (K)
    real, intent(in):: q(:, :), qs(:, :) ! (ncum, nl)
    real, intent(in):: gz(:, :) ! (klon, klev)
    real, intent(in):: u(:, :), v(:, :) ! (ncum, nl)
    real, intent(in):: p(:, :) ! (klon, klev) pressure at full level, in hPa
    real, intent(in):: ph(:, :) ! (ncum, klev + 1)
    real, intent(in):: th(:, :) ! (ncum, nl - 1) potential temperature, in K
    real, intent(in):: tv(:, :) ! (klon, klev)

    real, intent(in):: lv(:, :) ! (ncum, nl)
    ! specific latent heat of vaporization of water, in J kg-1

    real, intent(in):: cpn(:, :) ! (ncum, nl)
    ! specific heat capacity at constant pressure of humid air, in J K-1 kg-1

    real, intent(in):: ep(:, :) ! (ncum, klev)
    real, intent(in):: clw(:, :) ! (ncum, klev)
    real, intent(in):: m(:, :) ! (ncum, klev)
    real, intent(in):: ment(:, :, :) ! (ncum, klev, klev)
    real, intent(in):: elij(:, :, :) ! (ncum, klev, klev)
    real, intent(in):: delt
    real, intent(in):: plcl(:) ! (ncum)

    real, intent(out):: mp(:, :)
    ! (ncum, nl) Mass flux of the unsaturated downdraft, defined
    ! positive downward, in kg m-2 s-1. M_p in Emanuel (1991 928).

    real, intent(out):: qp(:, :), up(:, :), vp(:, :) ! (ncum, nl)
    real, intent(out):: wt(:, :) ! (ncum, nl)

    real, intent(out):: water(:, :) ! (ncum, nl)
    ! precipitation mixing ratio, l_p in Emanuel (1991 928)

    real, intent(out):: evap(:, :) ! (ncum, nl)
    ! sigt * rate of evaporation of precipitation, in s-1
    ! \sigma_s E in Emanuel (1991 928)

    real, intent(out):: b(:, :) ! (ncum, nl - 1)

    ! Local:

    real, parameter:: sigp = 0.15
    ! fraction of precipitation falling outside of cloud, \sig_s in
    ! Emanuel (1991 928)

    integer ncum
    integer i, il, imax
    real, parameter:: tinv = 1. / 3.
    real  delti
    real afac, afac1, afac2, bfac
    real pr1, sigt, b6, c6, revap, tevap
    real xf, tf, fac2, ur, sru, fac, d, af, bf
    real ampmax
    real lvcp(size(icb), nl) ! (ncum, nl) L_v / C_p, in K
    real wdtrain(size(icb)) ! (ncum)
    logical lwork(size(icb)) ! (ncum)

    !------------------------------------------------------

    ncum = size(icb)
    delti = 1. / delt
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

    ! Check whether ep(inb) = 0. If so, skip precipitating downdraft
    ! calculation.

    forall (il = 1:ncum) lwork(il) = ep(il, inb(il)) >= 1e-4

    imax = nl - 1
    do while (.not. any(inb >= imax .and. lwork) .and. imax >= 1)
       imax = imax - 1
    end do

    downdraft_loop: DO i = imax, 1, - 1
       ! Integrate liquid water equation to find condensed water 
       ! and condensed water flux 

       ! Calculate detrained precipitation 
       forall (il = 1:ncum, inb(il) >= i .and. lwork(il)) wdtrain(il) = rg &
            * (ep(il, i) * m(il, i) * clw(il, i) &
            + sum(max(elij(il, :i - 1, i) - (1. - ep(il, i)) * clw(il, i), 0.) &
            * ment(il, :i - 1, i)))

       ! Find rain water and evaporation using provisional 
       ! estimates of qp(i) and qp(i - 1) 

       loop_horizontal: do il = 1, ncum
          if (i <= inb(il) .and. lwork(il)) then
             wt(il, i) = 45.

             if (i < inb(il)) then
                qp(il, i) = qp(il, i + 1) + (rcpd * (t(il, i + 1) - t(il, i)) &
                     + gz(il, i + 1) - gz(il, i)) / lv(il, i)
                qp(il, i) = 0.5 * (qp(il, i) + q(il, i))
             endif

             qp(il, i) = max(qp(il, i), 0.)
             qp(il, i) = min(qp(il, i), qs(il, i))
             qp(il, inb(il)) = q(il, inb(il))

             if (i == 1) then
                afac = p(il, 1) * (qs(il, 1) - qp(il, 1)) &
                     / (1e4 + 2000. * p(il, 1) * qs(il, 1))
             else
                qp(il, i - 1) = qp(il, i) + (rcpd * (t(il, i) - t(il, i - 1)) &
                     + gz(il, i) - gz(il, i - 1)) / lv(il, i)
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

             if (i <= icb(il)) then
                ! Prise en compte de la variation progressive de sigt dans
                ! les couches icb et icb - 1 :
                ! pour plcl <= ph(i + 1), pr1 = 0
                ! pour plcl >= ph(i), pr1 = 1
                ! pour ph(i + 1) < plcl < ph(i), pr1 est la proportion
                ! \`a cheval sur le nuage.
                pr1 = max(0., min(1., &
                     (plcl(il) - ph(il, i + 1)) / (ph(il, i) - ph(il, i + 1))))
                sigt = sigp * pr1 + 1. - pr1
             else
                ! {i >= icb(il) + 1}
                sigt = sigp
             end if

             b6 = bfac * 50. * sigd * (ph(il, i) - ph(il, i + 1)) * sigt * afac
             c6 = water(il, i + 1) + bfac * wdtrain(il) - 50. * sigd * bfac &
                  * (ph(il, i) - ph(il, i + 1)) * evap(il, i + 1)

             if (c6 > 0.) then
                revap = 0.5 * (- b6 + sqrt(b6 * b6 + 4. * c6))
                evap(il, i) = sigt * afac * revap
                water(il, i) = revap * revap
             else
                evap(il, i) = - evap(il, i + 1) + 0.02 * (wdtrain(il) + sigd &
                     * wt(il, i) * water(il, i + 1)) / (sigd * (ph(il, i) &
                     - ph(il, i + 1)))
             end if

             ! Calculate precipitating downdraft mass flux under 
             ! hydrostatic approximation 

             test_above_surface: if (i /= 1) then
                tevap = max(0., evap(il, i))
                mp(il, i) = 100. * ginv * lvcp(il, i) * sigd * tevap &
                     * (p(il, i - 1) - p(il, i)) &
                     / max(0.001, th(il, i) - th(il, i - 1))

                ! If hydrostatic assumption fails, solve cubic
                ! difference equation for downdraft theta and mass
                ! flux from two simultaneous differential equations
                if (abs(mp(il, i + 1)**2 - mp(il, i)**2) > 0.1 * sigd**2 &
                     * 70. * ph(il, i) * (p(il, i - 1) - p(il, i)) &
                     * (th(il, i) - th(il, i - 1)) / (tv(il, i) * th(il, i))) &
                     then
                   xf = 100. * sigd**3 * (ph(il, i) - ph(il, i + 1))
                   tf = b(il, i) - 5. * (th(il, i) - th(il, i - 1)) &
                        * t(il, i) / (lvcp(il, i) * sigd * th(il, i))
                   af = xf * tf + mp(il, i + 1)**2 * tinv
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
                      mp(il, i) = mp(il, i + 1) * tinv + (0.5 * bf &
                           + sru)**tinv + fac * (abs(0.5 * bf - sru))**tinv
                   else
                      d = atan(2. * sqrt(- ur) / (bf + 1e-28))
                      if (fac2 < 0.) d = 3.14159 - d
                      mp(il, i) = mp(il, i + 1) * tinv + 2. * sqrt(af * tinv) &
                           * cos(d * tinv)
                   endif

                   mp(il, i) = max(0., mp(il, i))

                   ! Il y a vraisemblablement une erreur dans la ligne
                   ! suivante : il faut diviser par (mp(il, i) * sigd
                   ! * rg) et non par (mp(il, i) + sigd * 0.1).  Et il
                   ! faut bien revoir les facteurs 100.
                   b(il, i - 1) = max(b(il, i) + 100. * (p(il, i - 1) &
                        - p(il, i)) * tevap / (mp(il, i) + sigd * 0.1) - 10. &
                        * (th(il, i) - th(il, i - 1)) * t(il, i) &
                        / (lvcp(il, i) * sigd * th(il, i)), 0.)
                endif

                ! Limit magnitude of mp to meet CFL condition:
                ampmax = 2. * (ph(il, i) - ph(il, i + 1)) * delti
                ampmax = min(ampmax, 2. * (ph(il, i - 1) - ph(il, i)) * delti)
                mp(il, i) = min(mp(il, i), ampmax)

                ! Force mp to decrease linearly to zero between cloud
                ! base and the surface:
                if (p(il, i) > p(il, icb(il))) mp(il, i) = mp(il, icb(il)) &
                     * (p(il, 1) - p(il, i)) / (p(il, 1) - p(il, icb(il)))
             endif test_above_surface

             ! Find mixing ratio of precipitating downdraft 

             if (i /= inb(il)) then
                qp(il, i) = q(il, i)

                if (mp(il, i) > mp(il, i + 1)) then
                   qp(il, i) = qp(il, i + 1) * mp(il, i + 1) + q(il, i) &
                        * (mp(il, i) - mp(il, i + 1)) + 100. * ginv * 0.5 &
                        * sigd * (ph(il, i) - ph(il, i + 1)) &
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
                           * (ph(il, i) - ph(il, i + 1)) * (evap(il, i + 1) &
                           + evap(il, i)) / mp(il, i + 1)
                      up(il, i) = up(il, i + 1)
                      vp(il, i) = vp(il, i + 1)
                   endif
                endif

                qp(il, i) = min(qp(il, i), qs(il, i))
                qp(il, i) = max(qp(il, i), 0.)
             endif
          endif
       end do loop_horizontal
    end DO downdraft_loop

  end SUBROUTINE cv30_unsat

end module cv30_unsat_m
