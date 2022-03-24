module grid_noro_m

  implicit none

  REAL, SAVE, allocatable, protected:: mask(:, :) ! (iim + 1, jjm + 1)
  ! interpolated fraction of land

  REAL, SAVE, allocatable, protected:: phis(:, :) ! (iim + 1, jjm + 1)
  ! surface geopotential, not smoothed, in m2 s-2

contains

  SUBROUTINE grid_noro(xdata, ydata, relief, x, y, zmea, zstd, zsig, zgam, &
       zthe, zpic, zval)

    ! From dyn3d/grid_noro.F, version 1.1.1.1 2004/05/19 12:53:06

    ! Authors: Fran\c{}cois Lott, Laurent Li, A. Harzallah and Laurent
    ! Fairhead

    ! Compute the parameters of the sub-grid scale orography scheme as
    ! described in Lott and Miller (1997) and Lott (1999). Compute
    ! all the parameters needed for the gravity wave drag code.

    ! Target points are on a rectangular grid:
    ! jjm + 1 latitudes including North and South Poles;
    ! iim + 1 longitudes, with periodicity: longitude(iim + 1) = longitude(1)
    ! At the poles the field value is repeated iim + 1 times.

    ! The parameters a, b, c, d represent the limit of the target
    ! grid-point region. The mean over this region is calculated from
    ! US Navy data, ponderated by a weight proportional to the surface
    ! occupied by the data inside the model grid-point area. In most
    ! circumstances, this weight is the ratio between the surface of
    ! the US Navy gridpoint area and the surface of the model
    ! grid-point area. See documentation.

    ! Libraries:
    use jumble, only: assert, pi

    use comconst, only: ra
    use dimensions, only: iim, jjm
    use indicesol, only: epsfra
    use mva9_m, only: mva9
    USE yoegwd, only: gtsec

    ! Coordinates of input field:
    REAL, intent(in):: xdata(:) ! (iusn)
    REAL, intent(in):: ydata(:) ! (jusn)

    REAL, intent(in):: relief(:, :) ! (iusn, jusn) input field, in m
    REAL, intent(in):: x(:), y(:) ! coordinates of output field

    ! Correlations of US Navy orography gradients:

    real, intent(out):: zmea(:, :) ! (iim + 1, jjm + 1) smoothed mean orography
    real, intent(out):: zstd(:, :) ! (iim + 1, jjm + 1) standard deviation
    REAL, intent(out):: zsig(:, :) ! (iim + 1, jjm + 1) slope
    real, intent(out):: zgam(:, :) ! (iim + 1, jjm + 1) anisotropy

    real, intent(out):: zthe(:, :) ! (iim + 1, jjm + 1)
    ! orientation of the small axis (direction of maximum variance)

    REAL, intent(out):: zpic(:, :) ! (iim + 1, jjm + 1) maximum altitude
    real, intent(out):: zval(:, :) ! (iim + 1, jjm + 1) minimum altitude

    ! Local:

    integer iusn, jusn, iext
    REAL xusn((size(xdata) * 6) / 5), yusn(size(ydata) + 2)
    REAL zusn((size(xdata) * 6) / 5, size(ydata) + 2) ! in m

    ! Intermediate fields (correlations of orography gradient)
    REAL, dimension(iim + 1, jjm + 1):: ztz, zxtzx, zytzy, zxtzy, weight

    ! Correlations of US Navy orography gradients:
    REAL, dimension((size(xdata) * 6) / 5, size(ydata) + 2):: zxtzxusn, &
         zytzyusn, zxtzyusn

    real, dimension(iim + 1, jjm + 1):: mask_tmp, num_tot, num_lan
    REAL a(iim + 1), b(iim + 1), c(jjm + 1), d(jjm + 1)
    real weighx, weighy, xincr, xk, xp, xm, xw, xq, xl
    real zbordnor, zbordsud, zbordoue, zlenx, zleny, zmeasud

    real zdeltax, zdeltay
    ! utilisés pour calculer des dérivées zonales et méridiennes à
    ! partir de paire de points de grille séparés de deux pas de
    ! grille, et non pas adjacents
    
    real zweinor, zweisud, zmeanor, zbordest
    integer ii, i, jj, j

    !--------------------------------------------------------------------

    print *, "Call sequence information: grid_noro"
    allocate(mask(iim + 1, jjm + 1), phis(iim + 1, jjm + 1))
    iusn = size(xdata)
    jusn = size(ydata)
    iext = iusn / 10
    call assert(iusn == size(relief, 1), "grid_noro iusn")
    call assert(jusn == size(relief, 2), "grid_noro jusn")

    call assert([size(x), size(zmea, 1), size(zstd, 1), size(zsig, 1), &
         size(zgam, 1), size(zthe, 1), size(zpic, 1), &
         size(zval, 1)] == iim + 1, "grid_noro iim")

    call assert([size(y), size(zmea, 2), size(zstd, 2), size(zsig, 2), &
         size(zgam, 2), size(zthe, 2), size(zpic, 2), &
         size(zval, 2)] == jjm + 1, "grid_noro jjm")

    zdeltay = 2. * pi / real(jusn) * ra

    ! Extension of the US Navy database for computations at boundaries:

    DO j = 1, jusn
       yusn(j + 1) = ydata(j)
       DO i = 1, iusn
          zusn(i + iext, j + 1) = relief(i, j)
          xusn(i + iext) = xdata(i)
       ENDDO
       DO i = 1, iext
          zusn(i, j + 1) = relief(iusn - iext + i, j)
          xusn(i) = xdata(iusn - iext + i) - 2. * pi
          zusn(iusn + iext + i, j + 1) = relief(i, j)
          xusn(iusn + iext + i) = xdata(i) + 2. * pi
       ENDDO
    ENDDO

    yusn(1) = ydata(1) + (ydata(1) - ydata(2))
    yusn(jusn + 2) = ydata(jusn) + (ydata(jusn) - ydata(jusn - 1))
    DO i = 1, iusn / 2 + iext
       zusn(i, 1) = zusn(i + iusn / 2, 2)
       zusn(i + iusn / 2 + iext, 1) = zusn(i, 2)
       zusn(i, jusn + 2) = zusn(i + iusn / 2, jusn + 1)
       zusn(i + iusn / 2 + iext, jusn + 2) = zusn(i, jusn + 1)
    ENDDO

    ! Compute limits of model gridpoint area (regular grid)

    a(1) = x(1) - (x(2) - x(1)) / 2.0
    b(1) = (x(1) + x(2)) / 2.0
    DO i = 2, iim
       a(i) = b(i - 1)
       b(i) = (x(i) + x(i + 1)) / 2.0
    ENDDO
    a(iim + 1) = b(iim)
    b(iim + 1) = x(iim + 1) + (x(iim + 1) - x(iim)) / 2.0

    c(1) = y(1) - (y(2) - y(1)) / 2.0
    d(1) = (y(1) + y(2)) / 2.0
    DO j = 2, jjm
       c(j) = d(j - 1)
       d(j) = (y(j) + y(j + 1)) / 2.0
    ENDDO
    c(jjm + 1) = d(jjm)
    d(jjm + 1) = y(jjm + 1) + (y(jjm + 1) - y(jjm)) / 2.0

    ! Initialisations :
    weight = 0.
    zxtzx = 0.
    zytzy = 0.
    zxtzy = 0.
    ztz = 0.
    zmea = 0.
    zpic = - 1E10
    zval = 1E10

    ! Compute slopes correlations on US Navy grid

    zytzyusn = 0.
    zxtzxusn = 0.
    zxtzyusn = 0.

    DO j = 2, jusn + 1 
       zdeltax = zdeltay * cos(yusn(j))
       DO i = 2, iusn + 2 * iext - 1
          zytzyusn(i, j) = (zusn(i, j + 1) - zusn(i, j - 1))**2 / zdeltay**2
          zxtzxusn(i, j) = (zusn(i + 1, j) - zusn(i - 1, j))**2 / zdeltax**2
          zxtzyusn(i, j) = (zusn(i, j + 1) - zusn(i, j - 1)) / zdeltay &
               * (zusn(i + 1, j) - zusn(i - 1, j)) / zdeltax
       ENDDO
    ENDDO

    ! Summation over gridpoint area

    zleny = pi / real(jusn) * ra
    xincr = pi / 2. / real(jusn)
    DO ii = 1, iim + 1
       DO jj = 1, jjm + 1
          num_tot(ii, jj) = 0.
          num_lan(ii, jj) = 0.
          DO j = 2, jusn + 1 
             zlenx = zleny * cos(yusn(j))
             zdeltax = zdeltay * cos(yusn(j))
             zbordnor = (c(jj) - yusn(j) + xincr) * ra
             zbordsud = (yusn(j) - d(jj) + xincr) * ra
             weighy = MAX(0., min(zbordnor, zbordsud, zleny))
             IF (weighy /= 0) THEN
                DO i = 2, iusn + 2 * iext - 1
                   zbordest = (xusn(i) - a(ii) + xincr) * ra * cos(yusn(j))
                   zbordoue = (b(ii) + xincr - xusn(i)) * ra * cos(yusn(j))
                   weighx = MAX(0., min(zbordest, zbordoue, zlenx))
                   IF (weighx /= 0) THEN
                      num_tot(ii, jj) = num_tot(ii, jj) + 1.
                      if (zusn(i, j) >= 1.) &
                           num_lan(ii, jj) = num_lan(ii, jj) + 1.
                      weight(ii, jj) = weight(ii, jj) + weighx * weighy
                      zxtzx(ii, jj) = zxtzx(ii, jj) &
                           + zxtzxusn(i, j) * weighx * weighy
                      zytzy(ii, jj) = zytzy(ii, jj) &
                           + zytzyusn(i, j) * weighx * weighy
                      zxtzy(ii, jj) = zxtzy(ii, jj) &
                           + zxtzyusn(i, j) * weighx * weighy
                      ztz(ii, jj) = ztz(ii, jj) &
                           + zusn(i, j) * zusn(i, j) * weighx * weighy
                      ! mean
                      zmea(ii, jj) = zmea(ii, jj) + zusn(i, j) * weighx * weighy
                      ! peacks
                      zpic(ii, jj) = max(zpic(ii, jj), zusn(i, j))
                      ! valleys
                      zval(ii, jj) = min(zval(ii, jj), zusn(i, j))
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    if (any(weight == 0.)) then
       print *, "zero weight in grid_noro"
       stop 1
    end if

    ! Compute parameters needed by the Lott & Miller (1997) and Lott
    ! (1999) subgrid-scale orographic scheme.

    DO ii = 1, iim + 1
       DO jj = 1, jjm + 1
          mask(ii, jj) = num_lan(ii, jj) / num_tot(ii, jj)
          ! Mean orography:
          zmea (ii, jj) = zmea (ii, jj) / weight(ii, jj)
          zxtzx(ii, jj) = zxtzx(ii, jj) / weight(ii, jj)
          zytzy(ii, jj) = zytzy(ii, jj) / weight(ii, jj)
          zxtzy(ii, jj) = zxtzy(ii, jj) / weight(ii, jj)
          ztz(ii, jj) = ztz(ii, jj) / weight(ii, jj)
          ! Standard deviation:
          zstd(ii, jj) = sqrt(MAX(0., ztz(ii, jj) - zmea(ii, jj)**2))
       ENDDO
    ENDDO

    ! Correct values of horizontal slope near the poles:
    DO ii = 1, iim + 1
       zxtzx(ii, 1) = zxtzx(ii, 2)
       zxtzx(ii, jjm + 1) = zxtzx(ii, jjm)
       zxtzy(ii, 1) = zxtzy(ii, 2)
       zxtzy(ii, jjm + 1) = zxtzy(ii, jjm)
       zytzy(ii, 1) = zytzy(ii, 2)
       zytzy(ii, jjm + 1) = zytzy(ii, jjm)
    ENDDO

    ! Masque prenant en compte maximum de terre. On met un seuil \`a 10
    ! % de terre car en dessous les param\`etres de surface n'ont pas de
    ! sens.
    mask_tmp = merge(1., 0., mask >= 0.1)

    phis(:iim, 2:jjm) = zmea(:iim, 2:jjm) * mask_tmp(:iim, 2:jjm) * 9.81
    ! (zmea is not yet smoothed)

    ! Filters to smooth out fields for input into subgrid-scale
    ! orographic scheme.

    ! First filter, moving average over 9 points.
    CALL MVA9(zmea)
    CALL MVA9(zstd)
    CALL MVA9(zpic)
    CALL MVA9(zval)
    CALL MVA9(zxtzx)
    CALL MVA9(zxtzy) 
    CALL MVA9(zytzy)

    DO ii = 1, iim
       DO jj = 1, jjm + 1
          ! Coefficients K, L et M:
          xk = (zxtzx(ii, jj) + zytzy(ii, jj)) / 2.
          xl = (zxtzx(ii, jj) - zytzy(ii, jj)) / 2.
          xm = zxtzy(ii, jj)
          xp = xk - sqrt(xl**2 + xm**2)
          xq = xk + sqrt(xl**2 + xm**2)
          xw = 1e-8
          if (xp <= xw) xp = 0.
          if (xq <= xw) xq = xw
          if (abs(xm) <= xw) xm = xw * sign(1., xm)
          ! modification pour masque de terre fractionnaire
          ! slope: 
          zsig(ii, jj) = sqrt(xq) * mask_tmp(ii, jj)
          ! isotropy:
          zgam(ii, jj) = xp / xq * mask_tmp(ii, jj)
          ! angle theta:
          zthe(ii, jj) = 57.29577951 * atan2(xm, xl) / 2. * mask_tmp(ii, jj)
       ENDDO
    ENDDO

    zmea(:iim, :) = zmea(:iim, :) * mask_tmp(:iim, :)
    zpic(:iim, :) = zpic(:iim, :) * mask_tmp(:iim, :)
    zval(:iim, :) = zval(:iim, :) * mask_tmp(:iim, :)
    zstd(:iim, :) = zstd(:iim, :) * mask_tmp(:iim, :)

    ! gamma and theta at 1. and 0. at poles
    zmea(iim + 1, :) = zmea(1, :)
    phis(iim + 1, 2:jjm) = phis(1, 2:jjm)
    zpic(iim + 1, :) = zpic(1, :)
    zval(iim + 1, :) = zval(1, :)
    zstd(iim + 1, :) = zstd(1, :)
    zsig(iim + 1, :) = zsig(1, :)
    zgam(iim + 1, :) = zgam(1, :)
    zthe(iim + 1, :) = zthe(1, :)

    zweinor = sum(weight(:iim, 1))
    zweisud = sum(weight(:iim, jjm + 1))
    zmeanor = sum(zmea(:iim, 1) * weight(:iim, 1))
    zmeasud = sum(zmea(:iim, jjm + 1) * weight(:iim, jjm + 1))

    zmea(:, 1) = zmeanor / zweinor
    zmea(:, jjm + 1) = zmeasud / zweisud

    phis(:, 1) = zmeanor / zweinor * 9.81
    phis(:, jjm + 1) = zmeasud / zweisud * 9.81

    zpic(:, 1) = sum(zpic(:iim, 1) * weight(:iim, 1)) / zweinor
    zpic(:, jjm + 1) = sum(zpic(:iim, jjm + 1) * weight(:iim, jjm + 1)) &
         / zweisud

    zval(:, 1) = sum(zval(:iim, 1) * weight(:iim, 1)) / zweinor
    zval(:, jjm + 1) = sum(zval(:iim, jjm + 1) * weight(:iim, jjm + 1)) &
         / zweisud

    zstd(:, 1) = sum(zstd(:iim, 1) * weight(:iim, 1)) / zweinor
    zstd(:, jjm + 1) = sum(zstd(:iim, jjm + 1) * weight(:iim, jjm + 1)) &
         / zweisud

    zsig(:, 1) = sum(zsig(:iim, 1) * weight(:iim, 1)) / zweinor
    zsig(:, jjm + 1) = sum(zsig(:iim, jjm + 1) * weight(:iim, jjm + 1)) &
         / zweisud

    zgam(:, 1) = 1.
    zgam(:, jjm + 1) = 1.
    zgam = max(zgam, gtsec)

    zthe(:, 1) = 0.
    zthe(:, jjm + 1) = 0.

    mask(2:, 1) = mask(1, 1) ! north pole
    mask(2:, jjm + 1) = mask(1, jjm + 1) ! south pole
    mask(iim + 1, 2:jjm) = mask(1, 2:jjm) ! Greenwich

    WHERE (mask < EPSFRA)
       mask = 0.
    elsewhere (1. - mask < EPSFRA)
       mask = 1.
    endwhere

  END SUBROUTINE grid_noro

  !*********************************************************************

  subroutine read_phis(ncid_start)

    use dimensions, only: iim, jjm
    use netcdf95, only: NF95_GET_VAR, nf95_inq_varid

    integer, intent(in):: ncid_start

    ! Local:
    integer varid

    !------------------------------------------------------------------

    ALLOCATE(phis(iim + 1, jjm + 1))
    call NF95_INQ_VARID (ncid_start, "phis", varid)
    call NF95_GET_VAR(ncid_start, varid, phis)

  end subroutine read_phis

end module grid_noro_m
