module aaam_bud_m

  implicit none

contains

  subroutine aaam_bud(rg, ome, phis, dragu, liftu, phyu, dragv, liftv, phyv, &
       p, u, v, aam, torsfc)

    ! Author: F. Lott (LMD/CNRS). Date: 2003/10/20. Object: Compute
    ! different terms of the axial AAAM budget and mountain torque.
    ! Only valid for regular rectangular grids. Should be called after
    ! "lift_noro".

    USE dimens_m, ONLY : iim, jjm
    use nr_util, only: assert_eq, assert, pi
    use phyetat0_m, only: rlat, rlon
    USE suphec_m, ONLY: ra

    real, intent(in):: rg ! gravity constant
    real, intent(in):: ome ! Earth rotation rate
    real, intent(in):: phis(:) ! (nlon) Geopotential at the ground
    REAL, intent(in):: dragu(:) ! (nlon) orodrag stress (zonal)
    REAL, intent(in):: liftu(:) ! (nlon) orolift stress (zonal)
    REAL, intent(in):: phyu(:) ! (nlon) Stress total de la physique (zonal)
    REAL, intent(in):: dragv(:) ! (nlon) orodrag stress (Meridional)
    REAL, intent(in):: liftv(:) ! (nlon) orolift stress (Meridional)
    REAL, intent(in):: phyv(:) ! (nlon) Stress total de la physique (Meridional)

    REAL, intent(in):: p(:, :) 
    ! (nlon, nlev + 1) pressure (Pa) at model half levels

    real, intent(in):: u(:, :), v(:, :) ! (nlon, nlev) horizontal wind (m/s)
    REAL, intent(out):: aam ! axial component of wind AAM
    REAL, intent(out):: torsfc ! axial component of total surface torque

    ! Local Variables:

    INTEGER nlev ! number of vertical levels
    INTEGER i, j, k, l
    REAL dlat, dlon ! latitude and longitude increments (radians)

    REAL raam(3) ! wind AAM (components 1 & 2: equatorial; component 3: axial)
    REAL oaam(3) ! mass AAM (components 1 & 2: equatorial; component 3: axial)
    REAL tmou(3) ! resolved mountain torque (3 components)
    REAL tsso(3) ! parameterised moutain drag torque (3 components)
    REAL tbls(3) ! parameterised boundary layer torque (3 components)

    REAL ZS(801, 401) ! topographic height
    REAL PS(801, 401) ! surface pressure 
    REAL UB(801, 401), VB(801, 401) ! barotropic wind, zonal and meridional
    REAL SSOU(801, 401), SSOV(801, 401)
    REAL BLSU(801, 401), BLSV(801, 401)
    REAL ZLON(801), ZLAT(401) ! longitude and latitude in radians

    !-------------------------------------------------------------------

    call assert(size(phis) == (/size(dragu), size(liftu), size(phyu), &
         size(dragv), size(liftv), size(phyv), size(p, 1), size(u, 1), &
         size(v, 1)/), "aaam_bud nlon")
    nlev = assert_eq(size(p, 2) - 1, size(u, 2), size(v, 2), "aaam_bud nlev")

    if (iim + 1 > 801 .or. jjm + 1 > 401) then
       print *, ' Problème de dimension dans aaam_bud'
       stop 1
    endif

    dlat = pi / jjm
    dlon = 2 * pi / real(iim) 

    oaam = 0.
    raam = 0.
    tmou = 0.
    tsso = 0.
    tbls = 0.

    ! Mountain height, pressure and barotropic wind:

    ! North pole values (j = 1):

    ub(1, 1) = 0.
    vb(1, 1) = 0.
    do k = 1, nlev
       ub(1, 1) = ub(1, 1) + u(1, k) * (p(1, k) - p(1, k + 1)) / rg
       vb(1, 1) = vb(1, 1) + v(1, k) * (p(1, k) - p(1, k + 1)) / rg
    enddo

    zlat(1) = rlat(1) * pi / 180.

    do i = 1, iim + 1
       zs(i, 1) = phis(1) / rg
       ps(i, 1) = p(1, 1)
       ub(i, 1) = ub(1, 1) 
       vb(i, 1) = vb(1, 1) 
       ssou(i, 1) = dragu(1) + liftu(1)
       ssov(i, 1) = dragv(1) + liftv(1)
       blsu(i, 1) = phyu(1) - dragu(1) - liftu(1)
       blsv(i, 1) = phyv(1) - dragv(1) - liftv(1)
    enddo

    l = 1
    do j = 2, jjm
       ! Values at Greenwich (Periodicity)

       zs(iim + 1, j) = phis(l + 1) / rg
       ps(iim + 1, j) = p(l + 1, 1)
       ssou(iim + 1, j) = dragu(l + 1) + liftu(l + 1)
       ssov(iim + 1, j) = dragv(l + 1) + liftv(l + 1)
       blsu(iim + 1, j) = phyu(l + 1) - dragu(l + 1) - liftu(l + 1)
       blsv(iim + 1, j) = phyv(l + 1) - dragv(l + 1) - liftv(l + 1)
       zlon(iim + 1) = - rlon(l + 1) * pi / 180.
       zlat(j) = rlat(l + 1) * pi / 180.

       ub(iim + 1, j) = 0.
       vb(iim + 1, j) = 0.
       do k = 1, nlev
          ub(iim + 1, j) = ub(iim + 1, j) &
               + u(l + 1, k) * (p(l + 1, k) - p(l + 1, k + 1)) / rg
          vb(iim + 1, j) = vb(iim + 1, j) &
               + v(l + 1, k) * (p(l + 1, k) - p(l + 1, k + 1)) / rg
       enddo

       do i = 1, iim
          l = l + 1
          zs(i, j) = phis(l) / rg
          ps(i, j) = p(l, 1)
          ssou(i, j) = dragu(l) + liftu(l)
          ssov(i, j) = dragv(l) + liftv(l)
          blsu(i, j) = phyu(l) - dragu(l) - liftu(l)
          blsv(i, j) = phyv(l) - dragv(l) - liftv(l)
          zlon(i) = rlon(l) * pi / 180.

          ub(i, j) = 0.
          vb(i, j) = 0.
          do k = 1, nlev
             ub(i, j) = ub(i, j) + u(l, k) * (p(l, k) - p(l, k + 1)) / rg
             vb(i, j) = vb(i, j) + v(l, k) * (p(l, k) - p(l, k + 1)) / rg
          enddo
       enddo
    enddo

    ! South Pole

    l = l + 1
    ub(1, jjm + 1) = 0.
    vb(1, jjm + 1) = 0.
    do k = 1, nlev
       ub(1, jjm + 1) = ub(1, jjm + 1) + u(l, k) * (p(l, k) - p(l, k + 1)) / rg
       vb(1, jjm + 1) = vb(1, jjm + 1) + v(l, k) * (p(l, k) - p(l, k + 1)) / rg
    enddo
    zlat(jjm + 1) = rlat(l) * pi / 180.

    do i = 1, iim + 1
       zs(i, jjm + 1) = phis(l) / rg
       ps(i, jjm + 1) = p(l, 1)
       ssou(i, jjm + 1) = dragu(l) + liftu(l)
       ssov(i, jjm + 1) = dragv(l) + liftv(l)
       blsu(i, jjm + 1) = phyu(l) - dragu(l) - liftu(l)
       blsv(i, jjm + 1) = phyv(l) - dragv(l) - liftv(l)
       ub(i, jjm + 1) = ub(1, jjm + 1) 
       vb(i, jjm + 1) = vb(1, jjm + 1) 
    enddo

    ! Moment angulaire 

    DO j = 1, jjm 
       DO i = 1, iim
          raam(1) = raam(1) - ra**3 * dlon * dlat * 0.5 * (cos(zlon(i )) &
               * sin(zlat(j )) * cos(zlat(j )) * ub(i , j ) + cos(zlon(i )) &
               * sin(zlat(j + 1)) * cos(zlat(j + 1)) * ub(i , j + 1)) &
               + ra**3 * dlon * dlat * 0.5 * (sin(zlon(i )) * cos(zlat(j )) &
               * vb(i , j ) + sin(zlon(i )) * cos(zlat(j + 1)) * vb(i , j + 1))

          oaam(1) = oaam(1) - ome * ra**4 * dlon * dlat / rg * 0.5 &
               * (cos(zlon(i )) * cos(zlat(j ))**2 * sin(zlat(j )) &
               * ps(i , j ) + cos(zlon(i )) * cos(zlat(j + 1))**2 &
               * sin(zlat(j + 1)) * ps(i , j + 1))

          raam(2) = raam(2) - ra**3 * dlon * dlat * 0.5 * (sin(zlon(i )) &
               * sin(zlat(j )) * cos(zlat(j )) * ub(i , j ) + sin(zlon(i )) &
               * sin(zlat(j + 1)) * cos(zlat(j + 1)) * ub(i , j + 1)) &
               - ra**3 * dlon * dlat * 0.5 * (cos(zlon(i )) * cos(zlat(j )) &
               * vb(i , j ) + cos(zlon(i )) * cos(zlat(j + 1)) * vb(i , j + 1))

          oaam(2) = oaam(2) - ome * ra**4 * dlon * dlat / rg * 0.5 &
               * (sin(zlon(i )) * cos(zlat(j ))**2 * sin(zlat(j )) &
               * ps(i , j ) + sin(zlon(i )) * cos(zlat(j + 1))**2 &
               * sin(zlat(j + 1)) * ps(i , j + 1))

          raam(3) = raam(3) + ra**3 * dlon * dlat * 0.5 * (cos(zlat(j))**2 &
               * ub(i, j) + cos(zlat(j + 1))**2 * ub(i, j + 1))

          oaam(3) = oaam(3) + ome * ra**4 * dlon * dlat / rg * 0.5 &
               * (cos(zlat(j))**3 * ps(i, j) + cos(zlat(j + 1))**3 &
               * ps(i, j + 1))
       ENDDO
    ENDDO

    ! Couple des montagnes :

    DO j = 1, jjm
       DO i = 1, iim
          tmou(1) = tmou(1) - ra**2 * dlon * 0.5 * sin(zlon(i)) &
               * (zs(i, j) - zs(i, j + 1)) &
               * (cos(zlat(j + 1)) * ps(i, j + 1) + cos(zlat(j)) * ps(i, j)) 
          tmou(2) = tmou(2) + ra**2 * dlon * 0.5 * cos(zlon(i)) &
               * (zs(i, j) - zs(i, j + 1)) &
               * (cos(zlat(j + 1)) * ps(i, j + 1) + cos(zlat(j)) * ps(i, j)) 
       ENDDO
    ENDDO

    DO j = 2, jjm 
       DO i = 1, iim
          tmou(1) = tmou(1) + ra**2 * dlat * 0.5 * sin(zlat(j)) &
               * (zs(i + 1, j) - zs(i, j)) &
               * (cos(zlon(i + 1)) * ps(i + 1, j) + cos(zlon(i)) * ps(i, j))
          tmou(2) = tmou(2) + ra**2 * dlat * 0.5 * sin(zlat(j)) &
               * (zs(i + 1, j) - zs(i, j)) &
               * (sin(zlon(i + 1)) * ps(i + 1, j) + sin(zlon(i)) * ps(i, j))
          tmou(3) = tmou(3) - ra**2 * dlat * 0.5* cos(zlat(j)) &
               * (zs(i + 1, j) - zs(i, j)) * (ps(i + 1, j) + ps(i, j))
       ENDDO
    ENDDO

    ! Couples des differentes friction au sol :

    DO j = 2, jjm
       DO i = 1, iim
          tsso(1) = tsso(1) - ra**3 * cos(zlat(j)) * dlon * dlat* &
               ssou(i, j) * sin(zlat(j)) * cos(zlon(i)) &
               + ra**3 * cos(zlat(j)) * dlon * dlat* &
               ssov(i, j) * sin(zlon(i))

          tsso(2) = tsso(2) - ra**3 * cos(zlat(j)) * dlon * dlat* &
               ssou(i, j) * sin(zlat(j)) * sin(zlon(i)) &
               - ra**3 * cos(zlat(j)) * dlon * dlat* &
               ssov(i, j) * cos(zlon(i))

          tsso(3) = tsso(3) + ra**3 * cos(zlat(j)) * dlon * dlat* &
               ssou(i, j) * cos(zlat(j))

          tbls(1) = tbls(1) - ra**3 * cos(zlat(j)) * dlon * dlat* &
               blsu(i, j) * sin(zlat(j)) * cos(zlon(i)) &
               + ra**3 * cos(zlat(j)) * dlon * dlat* &
               blsv(i, j) * sin(zlon(i))

          tbls(2) = tbls(2) - ra**3 * cos(zlat(j)) * dlon * dlat* &
               blsu(i, j) * sin(zlat(j)) * sin(zlon(i)) &
               - ra**3 * cos(zlat(j)) * dlon * dlat* &
               blsv(i, j) * cos(zlon(i))

          tbls(3) = tbls(3) + ra**3 * cos(zlat(j)) * dlon * dlat* &
               blsu(i, j) * cos(zlat(j))
       ENDDO
    ENDDO

    aam = raam(3)
    torsfc = tmou(3) + tsso(3) + tbls(3)

  END subroutine aaam_bud

end module aaam_bud_m
