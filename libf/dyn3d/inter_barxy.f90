module inter_barxy_m

  ! From inter_barxy.F, version 1.1.1.1 2004/05/19 12:53:07

  implicit none

  private
  public inter_barxy

contains

  SUBROUTINE inter_barxy(dlonid, dlatid, champ, rlonimod, rlatimod, champint)

    ! Author: P. Le Van

    use nr_util, only: assert_eq, assert
    use dimens_m, only: iim, jjm
    use comgeom, only: aire_2d, apoln, apols

    REAL, intent(in):: dlonid(:)
    ! (longitude from input file, in rad, from -pi to pi)

    REAL, intent(in):: dlatid(:), champ(:, :), rlonimod(:)

    REAL, intent(in):: rlatimod(:)
    ! (latitude angle, in degrees or rad, in strictly decreasing order)

    real, intent(out):: champint(:, :)
    ! Si taille de la seconde dim = jjm + 1, on veut interpoler sur les
    ! jjm+1 latitudes rlatu du modele (latitudes des scalaires et de U)
    ! Si taille de la seconde dim = jjm, on veut interpoler sur les
    ! jjm latitudes rlatv du modèle (latitudes de V) 

    ! Variables local to the procedure:

    REAL champy(iim, size(champ, 2))
    integer j, i, jnterfd, jmods

    REAL yjmod(size(champint, 2))
    ! (angle, in degrees, in strictly increasing order)

    REAL   yjdat(size(dlatid) + 1) ! angle, in degrees, in increasing order
    LOGICAL decrois ! "dlatid" is in decreasing order

    !-----------------------------------

    jnterfd = assert_eq(size(champ, 2) - 1, size(dlatid), &
         "inter_barxy jnterfd")
    jmods = size(champint, 2)
    call assert(size(champ, 1) == size(dlonid), "inter_barxy size(champ, 1)")
    call assert((/size(rlonimod), size(champint, 1)/) == iim, &
         "inter_barxy iim")
    call assert(any(jmods == (/jjm, jjm + 1/)), 'inter_barxy jmods')
    call assert(size(rlatimod) == jjm, "inter_barxy size(rlatimod)")

    ! Check decreasing order for "rlatimod":
    DO i = 2, jjm
       IF (rlatimod(i) >= rlatimod(i-1)) stop &
            '"inter_barxy": "rlatimod" should be strictly decreasing'
    ENDDO

    yjmod(:jjm) = ord_coordm(rlatimod)
    IF (jmods == jjm + 1) THEN
       IF (90. - yjmod(jjm) < 0.01) stop &
            '"inter_barxy": with jmods = jjm + 1, yjmod(jjm) should be < 90.'
    ELSE
       ! jmods = jjm
       IF (ABS(yjmod(jjm) - 90.) > 0.01) stop &
            '"inter_barxy": with jmods = jjm, yjmod(jjm) should be 90.'
    ENDIF

    if (jmods == jjm + 1) yjmod(jjm + 1) = 90.

    DO j = 1, jnterfd + 1
       champy(:, j) = inter_barx(dlonid, champ(:, j), rlonimod)
    ENDDO

    CALL ord_coord(dlatid, yjdat, decrois) 
    IF (decrois) champy(:, :) = champy(:, jnterfd + 1:1:-1)
    DO i = 1, iim
       champint(i, :) = inter_bary(yjdat, champy(i, :), yjmod)
    ENDDO
    champint(:, :) = champint(:, jmods:1:-1)

    IF (jmods == jjm + 1) THEN
       ! Valeurs uniques aux poles
       champint(:, 1) = SUM(aire_2d(:iim,  1) * champint(:, 1)) / apoln
       champint(:, jjm + 1) = SUM(aire_2d(:iim, jjm + 1) &
            * champint(:, jjm + 1)) / apols
    ENDIF

  END SUBROUTINE inter_barxy

  !******************************

  function inter_barx(dlonid, fdat, rlonimod) 

    ! From dyn3d/inter_barx.F, v 1.1.1.1 2004/05/19 12:53:06

    ! Auteurs : Robert Sadourny, P. Le Van

    !        INTERPOLATION BARYCENTRIQUE BASEE SUR LES AIRES
    !            VERSION UNIDIMENSIONNELLE  ,   EN  LONGITUDE .

    !     idat : indice du champ de donnees, de 1 a idatmax
    !     imod : indice du champ du modele,  de 1 a  imodmax
    !     fdat(idat) : champ de donnees (entrees)
    !     inter_barx(imod) : champ du modele (sorties)
    !     dlonid(idat): abscisses des interfaces des mailles donnees
    !     rlonimod(imod): abscisses des interfaces des mailles modele
    !      ( L'indice 1 correspond a l'interface mailLE 1 / maille 2)
    !      ( Les abscisses sont exprimées en degres)

    use nr_util, only: assert_eq

    IMPLICIT NONE

    REAL, intent(in):: dlonid(:)
    real, intent(in):: fdat(:)
    real, intent(in):: rlonimod(:)

    real inter_barx(size(rlonimod))

    !    ...  Variables locales ... 

    INTEGER idatmax, imodmax
    REAL xxid(size(dlonid)+1), xxd(size(dlonid)+1), fdd(size(dlonid)+1)
    REAL  fxd(size(dlonid)+1), xchan(size(dlonid)+1), fdchan(size(dlonid)+1) 
    REAL  xxim(size(rlonimod))

    REAL x0, xim0, dx, dxm
    REAL chmin, chmax, pi

    INTEGER imod, idat, i, ichang, id0, id1, nid, idatmax1

    !-----------------------------------------------------

    idatmax = assert_eq(size(dlonid), size(fdat), "inter_barx idatmax")
    imodmax = size(rlonimod)

    pi = 2. * ASIN(1.)

    !   REDEFINITION DE L'ORIGINE DES ABSCISSES
    !    A L'INTERFACE OUEST DE LA PREMIERE MAILLE DU MODELE  
    DO imod = 1, imodmax
       xxim(imod) = rlonimod(imod)
    ENDDO

    CALL minmax( imodmax, xxim, chmin, chmax)
    IF( chmax.LT.6.50 )   THEN
       DO imod = 1, imodmax
          xxim(imod) = xxim(imod) * 180./pi
       ENDDO
    ENDIF

    xim0 = xxim(imodmax) - 360.

    DO imod = 1, imodmax
       xxim(imod) = xxim(imod) - xim0
    ENDDO

    idatmax1 = idatmax +1

    DO idat = 1, idatmax
       xxd(idat) = dlonid(idat)
    ENDDO

    CALL minmax( idatmax, xxd, chmin, chmax)
    IF( chmax.LT.6.50 )  THEN
       DO idat = 1, idatmax
          xxd(idat) = xxd(idat) * 180./pi
       ENDDO
    ENDIF

    DO idat = 1, idatmax
       xxd(idat) = AMOD( xxd(idat) - xim0, 360. )
       fdd(idat) = fdat (idat)
    ENDDO

    i = 2
    DO while (xxd(i) >= xxd(i-1) .and. i < idatmax)
       i = i + 1
    ENDDO
    IF (xxd(i) < xxd(i-1)) THEN
       ichang = i
       !  ***  reorganisation  des longitudes entre 0. et 360. degres ****
       nid = idatmax - ichang +1
       DO i = 1, nid
          xchan (i) = xxd(i+ichang -1 )
          fdchan(i) = fdd(i+ichang -1 )
       ENDDO
       DO i=1, ichang -1
          xchan (i+ nid) = xxd(i)
          fdchan(i+nid) = fdd(i) 
       ENDDO
       DO i =1, idatmax
          xxd(i) = xchan(i)
          fdd(i) = fdchan(i)
       ENDDO
    end IF

    !    translation des champs de donnees par rapport
    !    a la nouvelle origine, avec redondance de la
    !       maille a cheval sur les bords

    id0 = 0
    id1 = 0

    DO idat = 1, idatmax
       IF ( xxd( idatmax1- idat ).LT.360.) exit
       id1 = id1 + 1
    ENDDO

    DO idat = 1, idatmax
       IF (xxd(idat).GT.0.) exit
       id0 = id0 + 1
    END DO

    IF( id1 /= 0 ) then
       DO idat = 1, id1
          xxid(idat) = xxd(idatmax - id1 + idat) - 360.
          fxd (idat) = fdd(idatmax - id1 + idat)     
       END DO
       DO idat = 1, idatmax - id1
          xxid(idat + id1) = xxd(idat)
          fxd (idat + id1) = fdd(idat)
       END DO
    end IF

    IF(id0 /= 0) then
       DO idat = 1, idatmax - id0
          xxid(idat) = xxd(idat + id0)
          fxd (idat) = fdd(idat + id0)
       END DO

       DO idat = 1, id0
          xxid (idatmax - id0 + idat) =  xxd(idat) + 360.
          fxd  (idatmax - id0 + idat) =  fdd(idat)   
       END DO
    else 
       DO idat = 1, idatmax
          xxid(idat)  = xxd(idat)
          fxd (idat)  = fdd(idat)
       ENDDO
    end IF
    xxid(idatmax1) = xxid(1) + 360.
    fxd (idatmax1) = fxd(1)

    !   initialisation du champ du modele

    inter_barx(:) = 0.

    ! iteration

    x0   = xim0
    dxm  = 0.
    imod = 1
    idat = 1

    do while (imod <= imodmax)
       do while (xxim(imod).GT.xxid(idat))
          dx   = xxid(idat) - x0
          dxm  = dxm + dx
          inter_barx(imod) = inter_barx(imod) + dx * fxd(idat)
          x0   = xxid(idat)
          idat = idat + 1
       end do
       IF (xxim(imod).LT.xxid(idat)) THEN
          dx   = xxim(imod) - x0
          dxm  = dxm + dx
          inter_barx(imod) = (inter_barx(imod) + dx * fxd(idat)) / dxm
          x0   = xxim(imod)
          dxm  = 0.
          imod = imod + 1
       ELSE
          dx   = xxim(imod) - x0
          dxm  = dxm + dx
          inter_barx(imod) = (inter_barx(imod) + dx * fxd(idat)) / dxm
          x0   = xxim(imod)
          dxm  = 0.
          imod = imod + 1
          idat = idat + 1
       END IF
    end do

  END function inter_barx

  !******************************

  function inter_bary(yjdat, fdat, yjmod)

    ! From dyn3d/inter_bary.F, version 1.1.1.1 2004/05/19 12:53:06
    ! Authors: R. Sadourny, P. Le Van

    ! Interpolation barycentrique basée sur les aires.
    ! Version unidimensionnelle, en latitude.
    ! L'indice 1 correspond à l'interface maille 1 -- maille 2.

    use nr_util, only: assert

    IMPLICIT NONE

    REAL, intent(in):: yjdat(:)
    ! (angles, ordonnées des interfaces des mailles des données, in
    ! degrees, in increasing order)

    REAL, intent(in):: fdat(:) ! champ de données

    REAL, intent(in):: yjmod(:)
    ! (ordonnées des interfaces des mailles du modèle)
    ! (in degrees, in strictly increasing order)

    REAL inter_bary(size(yjmod)) ! champ du modèle

    ! Variables local to the procedure:

    REAL y0, dy, dym 
    INTEGER jdat ! indice du champ de données
    integer jmod ! indice du champ du modèle

    !------------------------------------

    call assert(size(yjdat) == size(fdat), "inter_bary")

    ! Initialisation des variables
    inter_bary(:) = 0.
    y0    = -90.
    dym   = 0.
    jmod  = 1
    jdat  = 1

    do while (jmod <= size(yjmod))
       do while (yjmod(jmod) > yjdat(jdat))
          dy         = yjdat(jdat) - y0
          dym        = dym + dy
          inter_bary(jmod) = inter_bary(jmod) + dy * fdat(jdat)
          y0         = yjdat(jdat)
          jdat       = jdat + 1
       end do
       IF (yjmod(jmod) < yjdat(jdat)) THEN
          dy         = yjmod(jmod) - y0
          dym        = dym + dy
          inter_bary(jmod) = (inter_bary(jmod) + dy * fdat(jdat)) / dym
          y0         = yjmod(jmod)
          dym        = 0.
          jmod       = jmod + 1
       ELSE
          ! {yjmod(jmod) == yjdat(jdat)}
          dy         = yjmod(jmod) - y0
          dym        = dym + dy
          inter_bary(jmod) = (inter_bary(jmod) + dy * fdat(jdat)) / dym
          y0         = yjmod(jmod)
          dym        = 0.
          jmod       = jmod + 1
          jdat       = jdat + 1
       END IF
    end do
    ! Le test de fin suppose que l'interface 0 est commune aux deux
    ! grilles "yjdat" et "yjmod".

  END function inter_bary

  !******************************

  SUBROUTINE ord_coord(xi, xo, decrois)

    ! From dyn3d/ord_coord.F, version 1.1.1.1 2004/05/19 12:53:06
    ! Author : P. Le Van

    ! This procedure receives an array of latitudes.
    ! It converts them to degrees if they are in radians.
    ! If the input latitudes are in decreasing order, the procedure
    ! reverses their order.
    ! Finally, the procedure adds 90° as the last value of the array.

    use nr_util, only: assert_eq, pi

    IMPLICIT NONE

    REAL, intent(in):: xi(:)
    ! (latitude, in degrees or radians, in increasing or decreasing order)
    ! ("xi" should contain latitudes from pole to pole.
    ! "xi" should contain the latitudes of the boundaries of grid
    ! cells, not the centers of grid cells.
    ! So the extreme values should not be 90° and -90°.)

    REAL, intent(out):: xo(:) ! angles in degrees
    LOGICAL, intent(out):: decrois

    ! Variables  local to the procedure:
    INTEGER nmax, i

    !--------------------

    nmax = assert_eq(size(xi), size(xo) - 1, "ord_coord")

    ! Check monotonicity:
    decrois = xi(2) < xi(1)
    DO i = 3, nmax
       IF (decrois .neqv. xi(i) < xi(i-1)) stop &
            '"ord_coord":  latitudes are not monotonic'
    ENDDO

    IF (abs(xi(1)) < pi) then
       ! "xi" contains latitudes in radians
       xo(:nmax) = xi(:) * 180. / pi
    else
       ! "xi" contains latitudes in degrees
       xo(:nmax) = xi(:)
    end IF

    IF (ABS(abs(xo(1)) - 90) < 0.001 .or. ABS(abs(xo(nmax)) - 90) < 0.001) THEN
       print *, "ord_coord"
       PRINT *, '"xi" should contain the latitudes of the boundaries of ' &
            // 'grid cells, not the centers of grid cells.'
       STOP
    ENDIF

    IF (decrois) xo(:nmax) = xo(nmax:1:- 1)
    xo(nmax + 1) = 90.

  END SUBROUTINE ord_coord

  !***********************************

  function ord_coordm(xi)

    ! From dyn3d/ord_coordm.F, version 1.1.1.1 2004/05/19 12:53:06
    ! Author : P. Le Van

    ! This procedure converts to degrees, if necessary, and inverts the
    ! order.

    use nr_util, only: pi

    IMPLICIT NONE

    REAL, intent(in):: xi(:) ! angle, in rad or degrees
    REAL ord_coordm(size(xi)) ! angle, in degrees

    !-----------------------------

    IF (xi(1) < 6.5) THEN
       ! "xi" is in rad
       ord_coordm(:) = xi(size(xi):1:-1) * 180. / pi
    else
       ! "xi" is in degrees
       ord_coordm(:) = xi(size(xi):1:-1)
    ENDIF

  END function ord_coordm

end module inter_barxy_m
