module grid_noro_m

  ! Clean: no C preprocessor directive, no include line

  implicit none

  private mva9

contains

  SUBROUTINE grid_noro(xdata, ydata, zdata, x, y, zphi, zmea, zstd, zsig, &
       zgam, zthe, zpic, zval, mask)

    ! From dyn3d/grid_noro.F, version 1.1.1.1 2004/05/19 12:53:06

    ! Authors: F. Lott, Z. X. Li, A. Harzallah and L. Fairhead

    ! Compute the parameters of the SSO scheme as described in
    ! Lott and Miller (1997) and Lott (1999).
    ! Target points are on a rectangular grid:
    ! jjm+1 latitudes including North and South Poles;
    ! iim+1 longitudes, with periodicity: longitude(iim+1)=longitude(1)
    ! At the poles the field value is repeated iim+1 times.

    ! The parameters a, b, c, d represent the limite of the target
    ! gridpoint region. The means over this region are calculated
    ! from USN data, ponderated by a weight proportional to the 
    ! surface occupied by the data inside the model gridpoint area.
    ! In most circumstances, this weight is the ratio between the
    ! surface of the USN gridpoint area and the surface of the
    ! model gridpoint area. 

    !           (c)
    !        ----d-----
    !        | . . . .|
    !        |        |
    !     (b)a . * . .b(a)
    !        |        |
    !        | . . . .|
    !        ----c-----
    !           (d)

    use dimens_m, only: iim, jjm
    use comconst, only: pi
    use numer_rec, only: assert

    REAL, intent(in):: xdata(:), ydata(:) ! coordinates of input field
    REAL, intent(in):: zdata(:, :) ! input field
    REAL, intent(in):: x(:), y(:) ! ccordinates output field

    ! Correlations of USN orography gradients:

    REAL zphi(:, :)
    real, intent(out):: zmea(:, :) ! Mean orography
    real, intent(out):: zstd(:, :) ! Standard deviation
    REAL zsig(:, :) ! Slope
    real zgam(:, :) ! Anisotropy
    real zthe(:, :) ! Orientation of the small axis
    REAL, intent(out):: zpic(:, :) ! Maximum altitude
    real, intent(out):: zval(:, :) ! Minimum altitude

    real, intent(out):: mask(:, :) ! fraction of land

    ! Variables local to the procedure:

    ! In this version it is assumed that the input data come from
    ! the US Navy dataset:
    integer, parameter:: iusn=2160, jusn=1080

    integer, parameter:: iext=216
    REAL xusn(iusn+2*iext), yusn(jusn+2)
    REAL zusn(iusn+2*iext, jusn+2)

    ! Intermediate fields  (correlations of orography gradient)

    REAL ztz(iim+1, jjm+1), zxtzx(iim+1, jjm+1)
    REAL zytzy(iim+1, jjm+1), zxtzy(iim+1, jjm+1)
    REAL weight(iim+1, jjm+1)

    ! Correlations of USN orography gradients:

    REAL zxtzxusn(iusn+2*iext, jusn+2), zytzyusn(iusn+2*iext, jusn+2)
    REAL zxtzyusn(iusn+2*iext, jusn+2)

    real mask_tmp(size(x), size(y))
    real num_tot(2200, 1100), num_lan(2200, 1100)

    REAL a(2200), b(2200), c(1100), d(1100)
    real rad, weighx, weighy, xincr, xk, xp, xm, xw, xq, xl
    real zbordnor, zdeltax, zbordsud, zdeltay, zbordoue, zlenx, zleny, zmeasud
    real zllmpic, zllmmea, zllmgam, zllmthe, zllmstd, zllmsig, zllmval
    real zpicnor, zminthe, zsigsud, zstdnor, zstdsud, zvalsud, zvalnor
    real zweinor, zweisud, zsignor, zpicsud, zmeanor, zbordest
    integer ii, i, jj, j

    !-------------------------------

    print *, "Call sequence information: grid_noro"

    call assert((/size(xdata), size(zdata, 1)/) == iusn, "grid_noro iusn")
    call assert((/size(ydata), size(zdata, 2)/) == jusn, "grid_noro jusn")

    call assert((/size(x), size(zphi, 1), size(zmea, 1), size(zstd, 1), &
         size(zsig, 1), size(zgam, 1), size(zthe, 1), size(zpic, 1), &
         size(zval, 1), size(mask, 1)/) == iim + 1, "grid_noro iim")

    call assert((/size(y), size(zphi, 2), size(zmea, 2), size(zstd, 2), &
         size(zsig, 2), size(zgam, 2), size(zthe, 2), size(zpic, 2), &
         size(zval, 2), size(mask, 2)/) == jjm + 1, "grid_noro jjm")

    IF (iim > 2200 .OR. jjm > 1099) THEN
       print *, "iim = ", iim, ", jjm = ", jjm
       stop '"iim" or "jjm" is too big'
    ENDIF

    print *, "Paramètres de l'orographie à l'échelle sous-maille" 
    rad = 6371229.
    zdeltay = 2. * pi / real(jusn) * rad

    ! Extension of the USN database to POCEED computations at boundaries:

    DO j=1, jusn
       yusn(j+1)=ydata(j)
       DO i=1, iusn
          zusn(i+iext, j+1)=zdata(i, j)
          xusn(i+iext)=xdata(i)
       ENDDO
       DO i=1, iext
          zusn(i, j+1)=zdata(iusn-iext+i, j)
          xusn(i)=xdata(iusn-iext+i)-2.*pi
          zusn(iusn+iext+i, j+1)=zdata(i, j)
          xusn(iusn+iext+i)=xdata(i)+2.*pi
       ENDDO
    ENDDO

    yusn(1)=ydata(1)+(ydata(1)-ydata(2))
    yusn(jusn+2)=ydata(jusn)+(ydata(jusn)-ydata(jusn-1))
    DO i=1, iusn/2+iext
       zusn(i, 1)=zusn(i+iusn/2, 2)
       zusn(i+iusn/2+iext, 1)=zusn(i, 2)
       zusn(i, jusn+2)=zusn(i+iusn/2, jusn+1)
       zusn(i+iusn/2+iext, jusn+2)=zusn(i, jusn+1)
    ENDDO

    ! COMPUTE LIMITS OF MODEL GRIDPOINT AREA
    !     ( REGULAR GRID)

    a(1) = x(1) - (x(2)-x(1))/2.0
    b(1) = (x(1)+x(2))/2.0
    DO i = 2, iim
       a(i) = b(i-1)
       b(i) = (x(i)+x(i+1))/2.0
    ENDDO
    a(iim+1) = b(iim)
    b(iim+1) = x(iim+1) + (x(iim+1)-x(iim))/2.0

    c(1) = y(1) - (y(2)-y(1))/2.0
    d(1) = (y(1)+y(2))/2.0
    DO j = 2, jjm
       c(j) = d(j-1)
       d(j) = (y(j)+y(j+1))/2.0
    ENDDO
    c(jjm + 1) = d(jjm)
    d(jjm + 1) = y(jjm + 1) + (y(jjm + 1)-y(jjm))/2.0

    ! Initialisations :
    weight(:, :) = 0.
    zxtzx(:, :)  = 0.
    zytzy(:, :)  = 0.
    zxtzy(:, :)  = 0.
    ztz(:, :)    = 0.
    zmea(:, :)   = 0.
    zpic(:, :)  =-1.E+10
    zval(:, :)  = 1.E+10

    !  COMPUTE SLOPES CORRELATIONS ON USN GRID

    zytzyusn(:, :)=0.
    zxtzxusn(:, :)=0.
    zxtzyusn(:, :)=0.

    DO j = 2, jusn+1 
       zdeltax=zdeltay*cos(yusn(j))
       DO i = 2, iusn+2*iext-1
          zytzyusn(i, j)=(zusn(i, j+1)-zusn(i, j-1))**2/zdeltay**2
          zxtzxusn(i, j)=(zusn(i+1, j)-zusn(i-1, j))**2/zdeltax**2
          zxtzyusn(i, j)=(zusn(i, j+1)-zusn(i, j-1))/zdeltay &
               *(zusn(i+1, j)-zusn(i-1, j))/zdeltax
       ENDDO
    ENDDO

    !  SUMMATION OVER GRIDPOINT AREA

    zleny=pi/real(jusn)*rad
    xincr=pi/2./real(jusn)
    DO ii = 1, iim+1
       DO jj = 1, jjm + 1
          num_tot(ii, jj)=0.
          num_lan(ii, jj)=0.
          DO j = 2, jusn+1 
             zlenx=zleny*cos(yusn(j))
             zdeltax=zdeltay*cos(yusn(j))
             zbordnor=(c(jj)-yusn(j)+xincr)*rad
             zbordsud=(yusn(j)-d(jj)+xincr)*rad
             weighy=AMAX1(0., amin1(zbordnor, zbordsud, zleny))
             IF (weighy /= 0) THEN
                DO i = 2, iusn+2*iext-1
                   zbordest=(xusn(i)-a(ii)+xincr)*rad*cos(yusn(j))
                   zbordoue=(b(ii)+xincr-xusn(i))*rad*cos(yusn(j))
                   weighx=AMAX1(0., amin1(zbordest, zbordoue, zlenx))
                   IF (weighx /= 0) THEN
                      num_tot(ii, jj) = num_tot(ii, jj) + 1.
                      if (zusn(i, j) >= 1.) then
                         num_lan(ii, jj) = num_lan(ii, jj) + 1.
                      end if
                      weight(ii, jj) = weight(ii, jj) + weighx * weighy
                      zxtzx(ii, jj)=zxtzx(ii, jj)+zxtzxusn(i, j)*weighx*weighy
                      zytzy(ii, jj)=zytzy(ii, jj)+zytzyusn(i, j)*weighx*weighy
                      zxtzy(ii, jj)=zxtzy(ii, jj)+zxtzyusn(i, j)*weighx*weighy
                      ztz(ii, jj) = ztz(ii, jj) &
                           + zusn(i, j) * zusn(i, j) * weighx * weighy
                      ! mean
                      zmea(ii, jj) =zmea(ii, jj)+zusn(i, j)*weighx*weighy
                      ! peacks
                      zpic(ii, jj)=amax1(zpic(ii, jj), zusn(i, j))
                      ! valleys
                      zval(ii, jj)=amin1(zval(ii, jj), zusn(i, j))
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    if (any(weight == 0.)) stop "zero weight in grid_noro"

    !  COMPUTE PARAMETERS NEEDED BY THE LOTT & MILLER (1997) AND
    !  LOTT (1999) SSO SCHEME.

    zllmmea=0.
    zllmstd=0.
    zllmsig=0.
    zllmgam=0.
    zllmpic=0.
    zllmval=0.
    zllmthe=0.
    zminthe=0.
    DO ii = 1, iim+1
       DO jj = 1, jjm + 1
          mask(ii, jj) = num_lan(ii, jj)/num_tot(ii, jj)
          !  Mean Orography:
          zmea (ii, jj)=zmea (ii, jj)/weight(ii, jj)
          zxtzx(ii, jj)=zxtzx(ii, jj)/weight(ii, jj)
          zytzy(ii, jj)=zytzy(ii, jj)/weight(ii, jj)
          zxtzy(ii, jj)=zxtzy(ii, jj)/weight(ii, jj)
          ztz(ii, jj)  =ztz(ii, jj)/weight(ii, jj)
          !  Standard deviation:
          zstd(ii, jj)=sqrt(AMAX1(0., ztz(ii, jj)-zmea(ii, jj)**2))
       ENDDO
    ENDDO

    ! CORRECT VALUES OF HORIZONTAL SLOPE NEAR THE POLES:

    DO ii = 1, iim+1
       zxtzx(ii, 1)=zxtzx(ii, 2)
       zxtzx(ii, jjm + 1)=zxtzx(ii, jjm)
       zxtzy(ii, 1)=zxtzy(ii, 2)
       zxtzy(ii, jjm + 1)=zxtzy(ii, jjm)
       zytzy(ii, 1)=zytzy(ii, 2)
       zytzy(ii, jjm + 1)=zytzy(ii, jjm)
    ENDDO

    !  FILTERS TO SMOOTH OUT FIELDS FOR INPUT INTO SSO SCHEME.

    !  FIRST FILTER, MOVING AVERAGE OVER 9 POINTS.

    CALL MVA9(zmea, iim+1, jjm+1)
    CALL MVA9(zstd, iim+1, jjm+1)
    CALL MVA9(zpic, iim+1, jjm+1)
    CALL MVA9(zval, iim+1, jjm+1)
    CALL MVA9(zxtzx, iim+1, jjm+1)
    CALL MVA9(zxtzy, iim+1, jjm+1) 
    CALL MVA9(zytzy, iim+1, jjm+1)

    ! Masque prenant en compte maximum de terre
    ! On seuil a 10% de terre de terre car en dessous les parametres
    ! de surface n'ont pas de sens (PB)
    mask_tmp= 0.
    WHERE (mask >= 0.1) mask_tmp = 1.

    DO ii = 1, iim
       DO jj = 1, jjm + 1
          IF (weight(ii, jj) /= 0.) THEN
             !  Coefficients K, L et M:
             xk=(zxtzx(ii, jj)+zytzy(ii, jj))/2.
             xl=(zxtzx(ii, jj)-zytzy(ii, jj))/2.
             xm=zxtzy(ii, jj)
             xp=xk-sqrt(xl**2+xm**2)
             xq=xk+sqrt(xl**2+xm**2)
             xw=1.e-8
             if(xp.le.xw) xp=0.
             if(xq.le.xw) xq=xw
             if(abs(xm).le.xw) xm=xw*sign(1., xm)
             !$$* PB modif pour maque de terre fractionnaire
             ! slope: 
             zsig(ii, jj)=sqrt(xq)*mask_tmp(ii, jj)
             ! isotropy:
             zgam(ii, jj)=xp/xq*mask_tmp(ii, jj)
             ! angle theta:
             zthe(ii, jj)=57.29577951*atan2(xm, xl)/2.*mask_tmp(ii, jj)
             zphi(ii, jj)=zmea(ii, jj)*mask_tmp(ii, jj)
             zmea(ii, jj)=zmea(ii, jj)*mask_tmp(ii, jj)
             zpic(ii, jj)=zpic(ii, jj)*mask_tmp(ii, jj)
             zval(ii, jj)=zval(ii, jj)*mask_tmp(ii, jj)
             zstd(ii, jj)=zstd(ii, jj)*mask_tmp(ii, jj)
          ENDIF
          zllmmea=AMAX1(zmea(ii, jj), zllmmea)
          zllmstd=AMAX1(zstd(ii, jj), zllmstd)
          zllmsig=AMAX1(zsig(ii, jj), zllmsig)
          zllmgam=AMAX1(zgam(ii, jj), zllmgam)
          zllmthe=AMAX1(zthe(ii, jj), zllmthe)
          zminthe=amin1(zthe(ii, jj), zminthe)
          zllmpic=AMAX1(zpic(ii, jj), zllmpic)
          zllmval=AMAX1(zval(ii, jj), zllmval)
       ENDDO
    ENDDO
    print *, 'MEAN ORO: ', zllmmea
    print *, 'ST. DEV.: ', zllmstd
    print *, 'PENTE: ', zllmsig
    print *, 'ANISOTROP: ', zllmgam
    print *, 'ANGLE: ', zminthe, zllmthe
    print *, 'pic: ', zllmpic
    print *, 'val: ', zllmval

    ! gamma and theta a 1. and 0. at poles
    zmea(iim+1, :)=zmea(1, :)
    zphi(iim+1, :)=zphi(1, :)
    zpic(iim+1, :)=zpic(1, :)
    zval(iim+1, :)=zval(1, :)
    zstd(iim+1, :)=zstd(1, :)
    zsig(iim+1, :)=zsig(1, :)
    zgam(iim+1, :)=zgam(1, :)
    zthe(iim+1, :)=zthe(1, :)

    zmeanor=0.
    zmeasud=0.
    zstdnor=0.
    zstdsud=0.
    zsignor=0.
    zsigsud=0.
    zweinor=0.
    zweisud=0.
    zpicnor=0.
    zpicsud=0.                                   
    zvalnor=0.
    zvalsud=0. 

    DO ii=1, iim
       zweinor=zweinor+              weight(ii,   1)
       zweisud=zweisud+              weight(ii, jjm + 1)
       zmeanor=zmeanor+zmea(ii,   1)*weight(ii,   1)
       zmeasud=zmeasud+zmea(ii, jjm + 1)*weight(ii, jjm + 1)
       zstdnor=zstdnor+zstd(ii,   1)*weight(ii,   1)
       zstdsud=zstdsud+zstd(ii, jjm + 1)*weight(ii, jjm + 1)
       zsignor=zsignor+zsig(ii,   1)*weight(ii,   1)
       zsigsud=zsigsud+zsig(ii, jjm + 1)*weight(ii, jjm + 1)
       zpicnor=zpicnor+zpic(ii,   1)*weight(ii,   1)
       zpicsud=zpicsud+zpic(ii, jjm + 1)*weight(ii, jjm + 1)
       zvalnor=zvalnor+zval(ii,   1)*weight(ii,   1)
       zvalsud=zvalsud+zval(ii, jjm + 1)*weight(ii, jjm + 1)
    ENDDO

    zmea(:,   1)=zmeanor/zweinor
    zmea(:, jjm + 1)=zmeasud/zweisud

    zphi(:,   1)=zmeanor/zweinor
    zphi(:, jjm + 1)=zmeasud/zweisud

    zpic(:,   1)=zpicnor/zweinor
    zpic(:, jjm + 1)=zpicsud/zweisud

    zval(:,   1)=zvalnor/zweinor
    zval(:, jjm + 1)=zvalsud/zweisud

    zstd(:,   1)=zstdnor/zweinor
    zstd(:, jjm + 1)=zstdsud/zweisud

    zsig(:,   1)=zsignor/zweinor
    zsig(:, jjm + 1)=zsigsud/zweisud

    zgam(:,   1)=1.
    zgam(:, jjm + 1)=1.

    zthe(:,   1)=0.
    zthe(:, jjm + 1)=0.

  END SUBROUTINE grid_noro

  !******************************************

  SUBROUTINE MVA9(X, IMAR, JMAR)

    ! From dyn3d/grid_noro.F, v 1.1.1.1 2004/05/19 12:53:06

    ! MAKE A MOVING AVERAGE OVER 9 GRIDPOINTS OF THE X FIELDS

    integer, intent(in):: imar, jmar
    REAL, intent(inout):: X(IMAR, JMAR)

    integer, PARAMETER:: ISMo=300, JSMo=200
    real XF(ISMo, JSMo)
    real WEIGHTpb(-1:1, -1:1)
    real sum
    integer i, is, js, j

    if(imar>ismo) stop 'surdimensionner ismo dans mva9 (grid_noro)'
    if(jmar>jsmo) stop 'surdimensionner jsmo dans mva9 (grid_noro)'

    SUM=0.
    DO IS=-1, 1
       DO JS=-1, 1
          WEIGHTpb(IS, JS)=1./FLOAT((1+IS**2)*(1+JS**2))
          SUM=SUM+WEIGHTpb(IS, JS)
       ENDDO
    ENDDO

    DO IS=-1, 1
       DO JS=-1, 1
          WEIGHTpb(IS, JS)=WEIGHTpb(IS, JS)/SUM
       ENDDO
    ENDDO

    DO J=2, JMAR-1
       DO I=2, IMAR-1
          XF(I, J)=0.
          DO IS=-1, 1
             DO JS=-1, 1
                XF(I, J)=XF(I, J)+X(I+IS, J+JS)*WEIGHTpb(IS, JS)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    DO J=2, JMAR-1
       XF(1, J)=0.
       IS=IMAR-1
       DO JS=-1, 1 
          XF(1, J)=XF(1, J)+X(IS, J+JS)*WEIGHTpb(-1, JS)
       ENDDO
       DO IS=0, 1 
          DO JS=-1, 1 
             XF(1, J)=XF(1, J)+X(1+IS, J+JS)*WEIGHTpb(IS, JS)
          ENDDO
       ENDDO
       XF(IMAR, J)=XF(1, J)
    ENDDO

    DO I=1, IMAR
       XF(I, 1)=XF(I, 2)
       XF(I, JMAR)=XF(I, JMAR-1)
    ENDDO

    DO I=1, IMAR
       DO J=1, JMAR
          X(I, J)=XF(I, J)
       ENDDO
    ENDDO

  END SUBROUTINE MVA9

end module grid_noro_m
