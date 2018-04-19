module inter_barx_m

  implicit none

contains

  function inter_barx(dlonid, fdat, rlonimod)

    ! From dyn3d/inter_barx.F, version 1.1.1.1, 2004/05/19 12:53:06

    ! Authors: Robert Sadourny, P. Le Van

    ! Interpolation barycentrique bas\'ee sur les aires. Version
    ! unidimensionnelle, en longitude.

    use nr_util, only: assert_eq, pi

    REAL, intent(in):: dlonid(:) ! (idatmax)
    ! abscisses des interfaces des mailles donn\'ees
    
    real, intent(in):: fdat(:) ! (idatmax) champ de donn\'ees

    real, intent(in):: rlonimod(:) ! (imodmax)
    ! Abscisses des interfaces des mailles mod\`ele. L'indice 1
    ! correspond a l'interface maille 1 / maille 2. Les abscisses sont
    ! exprim\'ees en degres.

    real inter_barx(size(rlonimod)) ! champ du mod\`ele

    ! Local: 
    INTEGER idatmax, imodmax
    REAL xxid(size(dlonid)+1), xxd(size(dlonid)+1), fdd(size(dlonid)+1)
    REAL fxd(size(dlonid)+1), xchan(size(dlonid)+1), fdchan(size(dlonid)+1)
    REAL xxim(size(rlonimod))
    REAL x0, xim0, dx, dxm
    REAL chmax
    integer idat ! indice du champ de donnees, de 1 a idatmax
    INTEGER imod ! indice du champ du modele, de 1 a imodmax
    integer i, ichang, id0, id1, nid, idatmax1

    !-----------------------------------------------------

    idatmax = assert_eq(size(dlonid), size(fdat), "inter_barx idatmax")
    imodmax = size(rlonimod)

    ! Red\'efinition de l'origine des abscisses \`a l'interface ouest de
    ! la premi\`ere maille du mod\`ele
    DO imod = 1, imodmax
       xxim(imod) = rlonimod(imod)
    ENDDO

    chmax = maxval(xxim)
    IF(chmax < 6.50) THEN
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

    chmax = maxval(xxd(:idatmax))
    IF(chmax < 6.50) THEN
       DO idat = 1, idatmax
          xxd(idat) = xxd(idat) * 180./pi
       ENDDO
    ENDIF

    DO idat = 1, idatmax
       xxd(idat) = AMOD(xxd(idat) - xim0, 360.)
       fdd(idat) = fdat (idat)
    ENDDO

    i = 2
    DO while (xxd(i) >= xxd(i-1) .and. i < idatmax)
       i = i + 1
    ENDDO
    IF (xxd(i) < xxd(i-1)) THEN
       ichang = i
       ! R\'eorganisation des longitudes entre 0 et 360 degr\'es
       nid = idatmax - ichang +1
       DO i = 1, nid
          xchan (i) = xxd(i+ichang -1)
          fdchan(i) = fdd(i+ichang -1)
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

    ! Translation des champs de donn\'ees par rapport \`a la nouvelle
    ! origine, avec redondance de la maille \`a cheval sur les bords

    id0 = 0
    id1 = 0

    DO idat = 1, idatmax
       IF (xxd(idatmax1- idat) < 360.) exit
       id1 = id1 + 1
    ENDDO

    DO idat = 1, idatmax
       IF (xxd(idat) > 0.) exit
       id0 = id0 + 1
    END DO

    IF(id1 /= 0) then
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
          xxid (idatmax - id0 + idat) = xxd(idat) + 360.
          fxd (idatmax - id0 + idat) = fdd(idat)
       END DO
    else
       DO idat = 1, idatmax
          xxid(idat) = xxd(idat)
          fxd (idat) = fdd(idat)
       ENDDO
    end IF
    xxid(idatmax1) = xxid(1) + 360.
    fxd (idatmax1) = fxd(1)

    ! Initialisation du champ du mod\`ele

    inter_barx(:) = 0.

    ! Iteration

    x0 = xim0
    dxm = 0.
    imod = 1
    idat = 1

    do while (imod <= imodmax)
       do while (xxim(imod) > xxid(idat))
          dx = xxid(idat) - x0
          dxm = dxm + dx
          inter_barx(imod) = inter_barx(imod) + dx * fxd(idat)
          x0 = xxid(idat)
          idat = idat + 1
       end do
       IF (xxim(imod) < xxid(idat)) THEN
          dx = xxim(imod) - x0
          dxm = dxm + dx
          inter_barx(imod) = (inter_barx(imod) + dx * fxd(idat)) / dxm
          x0 = xxim(imod)
          dxm = 0.
          imod = imod + 1
       ELSE
          dx = xxim(imod) - x0
          dxm = dxm + dx
          inter_barx(imod) = (inter_barx(imod) + dx * fxd(idat)) / dxm
          x0 = xxim(imod)
          dxm = 0.
          imod = imod + 1
          idat = idat + 1
       END IF
    end do

  END function inter_barx

end module inter_barx_m
