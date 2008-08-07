MODULE dynredem0_m

  IMPLICIT NONE

CONTAINS

  SUBROUTINE dynredem0(fichnom, iday_end, phis)

    ! From dyn3d/dynredem.F, v 1.2 2004/06/22 11:45:30
    ! Ecriture du fichier de redémarrage au format NetCDF (initialisation)

    USE ioipsl, ONLY : ju2ymds, ymds2ju
    USE dimens_m, ONLY : iim, jjm, llm, nqmx
    USE paramet_m, ONLY : iip1, jjp1, llmp1
    USE comconst, ONLY : cpp, daysec, dtvr, g, kappa, omeg, rad
    USE comvert, ONLY : ap, bp, nivsig, nivsigs, pa, preff, presnivs
    USE logic, ONLY : fxyhypb, ysinus
    USE comgeom, ONLY : aire_2d, cu_2d, cv_2d, rlatu, rlatv, rlonu, rlonv
    USE serre, ONLY : clat, clon, dzoomx, dzoomy, grossismx, grossismy, &
         taux, tauy
    USE temps, ONLY : annee_ref, day_ref, itaufin, itau_dyn
    USE ener, ONLY : ang0, etot0, ptot0, stot0, ztot0
    USE iniadvtrac_m, ONLY : tname, ttext
    USE netcdf95, ONLY : nf95_close, nf95_create, nf95_def_dim, &
         nf95_def_var, nf95_enddef, nf95_inq_varid, nf95_put_att, &
         nf95_put_var
    USE netcdf, ONLY : nf90_clobber, nf90_float, nf90_global, &
         nf90_unlimited

    CHARACTER (len=*), INTENT (IN) :: fichnom
    INTEGER, INTENT (IN) :: iday_end
    REAL, INTENT (IN) :: phis(:, :)

    !   Local:

    INTEGER :: iq, l
    INTEGER :: length
    PARAMETER (length=100)
    REAL :: tab_cntrl(length) ! tableau des parametres du run

    !   Variables locales pour NetCDF:

    INTEGER :: dims2(2), dims3(3), dims4(4)
    INTEGER :: idim_index
    INTEGER :: idim_rlonu, idim_rlonv, idim_rlatu, idim_rlatv
    INTEGER :: idim_s, idim_sig
    INTEGER :: idim_tim
    INTEGER :: nid, nvarid

    REAL :: zjulian, hours
    INTEGER :: yyears0, jjour0, mmois0
    CHARACTER (len=30) :: unites

    !-----------------------------------------------------------------------

    PRINT *, 'Call sequence information: dynredem0'

    CALL ymds2ju(annee_ref, 1, iday_end, 0.0, zjulian)
    CALL ju2ymds(zjulian, yyears0, mmois0, jjour0, hours)

    DO l = 1, length
       tab_cntrl(l) = 0.
    END DO
    tab_cntrl(1) = real(iim)
    tab_cntrl(2) = real(jjm)
    tab_cntrl(3) = real(llm)
    tab_cntrl(4) = real(day_ref)
    tab_cntrl(5) = real(annee_ref)
    tab_cntrl(6) = rad
    tab_cntrl(7) = omeg
    tab_cntrl(8) = g
    tab_cntrl(9) = cpp
    tab_cntrl(10) = kappa
    tab_cntrl(11) = daysec
    tab_cntrl(12) = dtvr
    tab_cntrl(13) = etot0
    tab_cntrl(14) = ptot0
    tab_cntrl(15) = ztot0
    tab_cntrl(16) = stot0
    tab_cntrl(17) = ang0
    tab_cntrl(18) = pa
    tab_cntrl(19) = preff

    ! Paramètres  pour le zoom :

    tab_cntrl(20) = clon
    tab_cntrl(21) = clat
    tab_cntrl(22) = grossismx
    tab_cntrl(23) = grossismy

    IF (fxyhypb) THEN
       tab_cntrl(24) = 1.
       tab_cntrl(25) = dzoomx
       tab_cntrl(26) = dzoomy
       tab_cntrl(27) = 0.
       tab_cntrl(28) = taux
       tab_cntrl(29) = tauy
    ELSE
       tab_cntrl(24) = 0.
       tab_cntrl(25) = dzoomx
       tab_cntrl(26) = dzoomy
       tab_cntrl(27) = 0.
       tab_cntrl(28) = 0.
       tab_cntrl(29) = 0.
       IF (ysinus) tab_cntrl(27) = 1.
    END IF

    tab_cntrl(30) = real(iday_end)
    tab_cntrl(31) = real(itau_dyn+itaufin)

    CALL nf95_create(fichnom, nf90_clobber, nid)
    CALL nf95_put_att(nid, nf90_global, 'title', &
         'Fichier de démarrage dynamique')

    ! Definir les dimensions du fichiers:

    CALL nf95_def_dim(nid, 'index', length, idim_index)
    CALL nf95_def_dim(nid, 'rlonu', iip1, idim_rlonu)
    CALL nf95_def_dim(nid, 'rlatu', jjp1, idim_rlatu)
    CALL nf95_def_dim(nid, 'rlonv', iip1, idim_rlonv)
    CALL nf95_def_dim(nid, 'rlatv', jjm, idim_rlatv)
    CALL nf95_def_dim(nid, 'sigs', llm, idim_s)
    CALL nf95_def_dim(nid, 'sig', llmp1, idim_sig)
    CALL nf95_def_dim(nid, 'temps', nf90_unlimited, idim_tim)

    ! Definir et enregistrer certains champs invariants:

    CALL nf95_def_var(nid, 'controle', nf90_float, idim_index, nvarid)
    CALL nf95_put_att(nid, nvarid, 'title', 'Parametres de controle')

    CALL nf95_def_var(nid, 'rlonu', nf90_float, idim_rlonu, nvarid)
    CALL nf95_put_att(nid, nvarid, 'title', 'Longitudes des points U')

    CALL nf95_def_var(nid, 'rlatu', nf90_float, idim_rlatu, nvarid)
    CALL nf95_put_att(nid, nvarid, 'title', 'Latitudes des points U')

    CALL nf95_def_var(nid, 'rlonv', nf90_float, idim_rlonv, nvarid)
    CALL nf95_put_att(nid, nvarid, 'title', 'Longitudes des points V')

    CALL nf95_def_var(nid, 'rlatv', nf90_float, idim_rlatv, nvarid)
    CALL nf95_put_att(nid, nvarid, 'title', 'Latitudes des points V')

    CALL nf95_def_var(nid, 'nivsigs', nf90_float, idim_s, nvarid)
    CALL nf95_put_att(nid, nvarid, 'title', 'Numero naturel des couches s')

    CALL nf95_def_var(nid, 'nivsig', nf90_float, idim_sig, nvarid)
    CALL nf95_put_att(nid, nvarid, 'title', &
         'Numero naturel des couches sigma')

    CALL nf95_def_var(nid, 'ap', nf90_float, idim_sig, nvarid)
    CALL nf95_put_att(nid, nvarid, 'title', 'Coefficient A pour hybride')

    CALL nf95_def_var(nid, 'bp', nf90_float, idim_sig, nvarid)
    CALL nf95_put_att(nid, nvarid, 'title', 'Coefficient B pour hybride')

    CALL nf95_def_var(nid, 'presnivs', nf90_float, idim_s, nvarid)

    ! Coefficients de passage cov. <-> contra. <--> naturel

    dims2(1) = idim_rlonu
    dims2(2) = idim_rlatu
    CALL nf95_def_var(nid, 'cu', nf90_float, dims2, nvarid)
    CALL nf95_put_att(nid, nvarid, 'title', 'Coefficient de passage pour U')

    dims2(1) = idim_rlonv
    dims2(2) = idim_rlatv
    CALL nf95_def_var(nid, 'cv', nf90_float, dims2, nvarid)
    CALL nf95_put_att(nid, nvarid, 'title', 'Coefficient de passage pour V')

    ! Aire de chaque maille:

    dims2(1) = idim_rlonv
    dims2(2) = idim_rlatu
    CALL nf95_def_var(nid, 'aire', nf90_float, dims2, nvarid)
    CALL nf95_put_att(nid, nvarid, 'title', 'Aires de chaque maille')

    ! Geopentiel au sol:

    dims2(1) = idim_rlonv
    dims2(2) = idim_rlatu
    CALL nf95_def_var(nid, 'phisinit', nf90_float, dims2, nvarid)
    CALL nf95_put_att(nid, nvarid, 'title', 'Geopotentiel au sol')

    ! Definir les variables pour pouvoir les enregistrer plus tard:

    CALL nf95_def_var(nid, 'temps', nf90_float, idim_tim, nvarid)
    CALL nf95_put_att(nid, nvarid, 'title', 'Temps de simulation')
    WRITE (unites, 200) yyears0, mmois0, jjour0
200 FORMAT ('days since ', I4, '-', I2.2, '-', I2.2, ' 00:00:00')
    CALL nf95_put_att(nid, nvarid, 'units', unites)


    dims4(1) = idim_rlonu
    dims4(2) = idim_rlatu
    dims4(3) = idim_s
    dims4(4) = idim_tim
    CALL nf95_def_var(nid, 'ucov', nf90_float, dims4, nvarid)
    CALL nf95_put_att(nid, nvarid, 'title', 'Vitesse U')

    dims4(1) = idim_rlonv
    dims4(2) = idim_rlatv
    dims4(3) = idim_s
    dims4(4) = idim_tim
    CALL nf95_def_var(nid, 'vcov', nf90_float, dims4, nvarid)
    CALL nf95_put_att(nid, nvarid, 'title', 'Vitesse V')

    dims4(1) = idim_rlonv
    dims4(2) = idim_rlatu
    dims4(3) = idim_s
    dims4(4) = idim_tim
    CALL nf95_def_var(nid, 'teta', nf90_float, dims4, nvarid)
    CALL nf95_put_att(nid, nvarid, 'title', 'Temperature')

    dims4(1) = idim_rlonv
    dims4(2) = idim_rlatu
    dims4(3) = idim_s
    dims4(4) = idim_tim
    DO iq = 1, nqmx
       CALL nf95_def_var(nid, tname(iq), nf90_float, dims4, nvarid)
       CALL nf95_put_att(nid, nvarid, 'title', ttext(iq))
    END DO

    dims4(1) = idim_rlonv
    dims4(2) = idim_rlatu
    dims4(3) = idim_s
    dims4(4) = idim_tim
    CALL nf95_def_var(nid, 'masse', nf90_float, dims4, nvarid)
    CALL nf95_put_att(nid, nvarid, 'title', 'C est quoi ?')

    dims3(1) = idim_rlonv
    dims3(2) = idim_rlatu
    dims3(3) = idim_tim
    CALL nf95_def_var(nid, 'ps', nf90_float, dims3, nvarid)
    CALL nf95_put_att(nid, nvarid, 'title', 'Pression au sol')

    CALL nf95_enddef(nid)

    CALL nf95_inq_varid(nid, 'controle', nvarid)
    CALL nf95_put_var(nid, nvarid, tab_cntrl)

    CALL nf95_inq_varid(nid, 'rlonu', nvarid)
    CALL nf95_put_var(nid, nvarid, rlonu)

    CALL nf95_inq_varid(nid, 'rlatu', nvarid)
    CALL nf95_put_var(nid, nvarid, rlatu)

    CALL nf95_inq_varid(nid, 'rlonv', nvarid)
    CALL nf95_put_var(nid, nvarid, rlonv)

    CALL nf95_inq_varid(nid, 'rlatv', nvarid)
    CALL nf95_put_var(nid, nvarid, rlatv)

    CALL nf95_inq_varid(nid, 'nivsigs', nvarid)
    CALL nf95_put_var(nid, nvarid, nivsigs)

    CALL nf95_inq_varid(nid, 'nivsig', nvarid)
    CALL nf95_put_var(nid, nvarid, nivsig)

    CALL nf95_inq_varid(nid, 'ap', nvarid)
    CALL nf95_put_var(nid, nvarid, ap)

    CALL nf95_inq_varid(nid, 'bp', nvarid)
    CALL nf95_put_var(nid, nvarid, bp)

    CALL nf95_inq_varid(nid, 'presnivs', nvarid)
    CALL nf95_put_var(nid, nvarid, presnivs)

    CALL nf95_inq_varid(nid, 'cu', nvarid)
    CALL nf95_put_var(nid, nvarid, cu_2d)

    CALL nf95_inq_varid(nid, 'cv', nvarid)
    CALL nf95_put_var(nid, nvarid, cv_2d)

    CALL nf95_inq_varid(nid, 'aire', nvarid)
    CALL nf95_put_var(nid, nvarid, aire_2d)

    CALL nf95_inq_varid(nid, 'phisinit', nvarid)
    CALL nf95_put_var(nid, nvarid, phis)

    CALL nf95_close(nid) ! fermer le fichier

    PRINT *, 'iim, jjm, llm, iday_end', iim, jjm, llm, iday_end
    PRINT *, 'rad, omeg, g, cpp, kappa', rad, omeg, g, cpp, kappa

  END SUBROUTINE dynredem0

END MODULE dynredem0_m
