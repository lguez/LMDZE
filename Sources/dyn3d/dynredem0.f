MODULE dynredem0_m

  IMPLICIT NONE

CONTAINS

  SUBROUTINE dynredem0(fichnom, iday_end, phis)

    ! From dyn3d/dynredem.F, version 1.2 2004/06/22 11:45:30
    ! Ecriture du fichier de redémarrage au format NetCDF (initialisation)

    USE comconst, ONLY: cpp, daysec, dtvr, g, kappa, omeg, rad
    USE comgeom, ONLY: aire_2d, cu_2d, cv_2d, rlatu, rlatv, rlonu, rlonv
    USE dimens_m, ONLY: iim, jjm, llm, nqmx
    USE disvert_m, ONLY: ap, bp, pa, preff, presnivs
    use dynetat0_m, only: day_ref, annee_ref
    USE ener, ONLY: ang0, etot0, ptot0, stot0, ztot0
    USE iniadvtrac_m, ONLY: tname, ttext
    USE ju2ymds_m, ONLY: ju2ymds
    USE netcdf, ONLY: nf90_clobber, nf90_float, nf90_global, nf90_unlimited
    USE netcdf95, ONLY: nf95_close, nf95_create, nf95_def_dim, nf95_def_var, &
         nf95_enddef, nf95_inq_varid, nf95_put_att, nf95_put_var
    USE paramet_m, ONLY: iip1, jjp1, llmp1
    USE serre, ONLY: clat, clon, dzoomx, dzoomy, grossismx, grossismy, taux, &
         tauy
    use ymds2ju_m, only: ymds2ju

    CHARACTER(len=*), INTENT(IN):: fichnom
    INTEGER, INTENT(IN):: iday_end
    REAL, INTENT(IN):: phis(:, :)

    ! Local:

    INTEGER iq, l
    INTEGER, PARAMETER:: length = 100
    REAL tab_cntrl(length) ! tableau des paramètres du run

    ! Pour NetCDF :
    INTEGER idim_index
    INTEGER idim_rlonu, idim_rlonv, idim_rlatu, idim_rlatv
    INTEGER idim_s, idim_sig
    INTEGER dimid_temps
    INTEGER ncid, varid

    REAL zjulian, hours
    INTEGER yyears0, jjour0, mmois0
    CHARACTER(len=30) unites

    !-----------------------------------------------------------------------

    PRINT *, 'Call sequence information: dynredem0'

    CALL ymds2ju(annee_ref, 1, iday_end, 0., zjulian)
    CALL ju2ymds(zjulian, yyears0, mmois0, jjour0, hours)

    DO l = 1, length
       tab_cntrl(l) = 0.
    END DO
    tab_cntrl(1) = iim
    tab_cntrl(2) = jjm
    tab_cntrl(3) = llm
    tab_cntrl(4) = day_ref
    tab_cntrl(5) = annee_ref
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

    ! Paramètres pour le zoom :
    tab_cntrl(20) = clon
    tab_cntrl(21) = clat
    tab_cntrl(22) = grossismx
    tab_cntrl(23) = grossismy
    tab_cntrl(24) = 1.
    tab_cntrl(25) = dzoomx
    tab_cntrl(26) = dzoomy
    tab_cntrl(27) = 0.
    tab_cntrl(28) = taux
    tab_cntrl(29) = tauy

    tab_cntrl(30) = iday_end

    CALL nf95_create(fichnom, nf90_clobber, ncid)
    CALL nf95_put_att(ncid, nf90_global, 'title', &
         'Fichier de démarrage dynamique')

    ! Definir les dimensions du fichiers:

    CALL nf95_def_dim(ncid, 'index', length, idim_index)
    CALL nf95_def_dim(ncid, 'rlonu', iip1, idim_rlonu)
    CALL nf95_def_dim(ncid, 'rlatu', jjp1, idim_rlatu)
    CALL nf95_def_dim(ncid, 'rlonv', iip1, idim_rlonv)
    CALL nf95_def_dim(ncid, 'rlatv', jjm, idim_rlatv)
    CALL nf95_def_dim(ncid, 'sigs', llm, idim_s)
    CALL nf95_def_dim(ncid, 'sig', llmp1, idim_sig)
    CALL nf95_def_dim(ncid, 'temps', nf90_unlimited, dimid_temps)

    ! Definir et enregistrer certains champs invariants:

    CALL nf95_def_var(ncid, 'controle', nf90_float, idim_index, varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Parametres de controle')

    CALL nf95_def_var(ncid, 'rlonu', nf90_float, idim_rlonu, varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Longitudes des points U')

    CALL nf95_def_var(ncid, 'rlatu', nf90_float, idim_rlatu, varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Latitudes des points U')

    CALL nf95_def_var(ncid, 'rlonv', nf90_float, idim_rlonv, varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Longitudes des points V')

    CALL nf95_def_var(ncid, 'rlatv', nf90_float, idim_rlatv, varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Latitudes des points V')

    CALL nf95_def_var(ncid, 'ap', nf90_float, idim_sig, varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Coefficient A pour hybride')

    CALL nf95_def_var(ncid, 'bp', nf90_float, idim_sig, varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Coefficient B pour hybride')

    CALL nf95_def_var(ncid, 'presnivs', nf90_float, idim_s, varid)

    ! Coefficients de passage cov. <-> contra. <--> naturel

    CALL nf95_def_var(ncid, 'cu', nf90_float, (/idim_rlonu, idim_rlatu/), varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Coefficient de passage pour U')

    CALL nf95_def_var(ncid, 'cv', nf90_float, (/idim_rlonv, idim_rlatv/), varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Coefficient de passage pour V')

    ! Aire de chaque maille:

    CALL nf95_def_var(ncid, 'aire', nf90_float, (/idim_rlonv, idim_rlatu/), &
         varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Aires de chaque maille')

    ! Geopentiel au sol:

    CALL nf95_def_var(ncid, 'phisinit', nf90_float, &
         (/idim_rlonv, idim_rlatu/), varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Geopotentiel au sol')

    ! Definir les variables pour pouvoir les enregistrer plus tard:

    CALL nf95_def_var(ncid, 'temps', nf90_float, dimid_temps, varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Temps de simulation')
    WRITE(unites, fmt = 200) yyears0, mmois0, jjour0
200 FORMAT ('days since ', I4, '-', I2.2, '-', I2.2, ' 00:00:00')
    CALL nf95_put_att(ncid, varid, 'units', unites)

    CALL nf95_def_var(ncid, 'ucov', nf90_float, &
         (/idim_rlonu, idim_rlatu, idim_s, dimid_temps/), varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Vitesse U')

    CALL nf95_def_var(ncid, 'vcov', nf90_float, &
         (/idim_rlonv, idim_rlatv, idim_s, dimid_temps/), varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Vitesse V')

    CALL nf95_def_var(ncid, 'teta', nf90_float, &
         (/idim_rlonv, idim_rlatu, idim_s, dimid_temps/), varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Temperature')

    DO iq = 1, nqmx
       CALL nf95_def_var(ncid, tname(iq), nf90_float, &
            (/idim_rlonv, idim_rlatu, idim_s, dimid_temps/), varid)
       CALL nf95_put_att(ncid, varid, 'title', ttext(iq))
    END DO

    CALL nf95_def_var(ncid, 'masse', nf90_float, &
         (/idim_rlonv, idim_rlatu, idim_s, dimid_temps/), varid)
    CALL nf95_put_att(ncid, varid, 'title', 'C est quoi ?')

    CALL nf95_def_var(ncid, 'ps', nf90_float, &
         (/idim_rlonv, idim_rlatu, dimid_temps/), varid)
    CALL nf95_put_att(ncid, varid, 'title', 'Pression au sol')

    CALL nf95_enddef(ncid)

    CALL nf95_inq_varid(ncid, 'controle', varid)
    CALL nf95_put_var(ncid, varid, tab_cntrl)

    CALL nf95_inq_varid(ncid, 'rlonu', varid)
    CALL nf95_put_var(ncid, varid, rlonu)

    CALL nf95_inq_varid(ncid, 'rlatu', varid)
    CALL nf95_put_var(ncid, varid, rlatu)

    CALL nf95_inq_varid(ncid, 'rlonv', varid)
    CALL nf95_put_var(ncid, varid, rlonv)

    CALL nf95_inq_varid(ncid, 'rlatv', varid)
    CALL nf95_put_var(ncid, varid, rlatv)

    CALL nf95_inq_varid(ncid, 'ap', varid)
    CALL nf95_put_var(ncid, varid, ap)

    CALL nf95_inq_varid(ncid, 'bp', varid)
    CALL nf95_put_var(ncid, varid, bp)

    CALL nf95_inq_varid(ncid, 'presnivs', varid)
    CALL nf95_put_var(ncid, varid, presnivs)

    CALL nf95_inq_varid(ncid, 'cu', varid)
    CALL nf95_put_var(ncid, varid, cu_2d)

    CALL nf95_inq_varid(ncid, 'cv', varid)
    CALL nf95_put_var(ncid, varid, cv_2d)

    CALL nf95_inq_varid(ncid, 'aire', varid)
    CALL nf95_put_var(ncid, varid, aire_2d)

    CALL nf95_inq_varid(ncid, 'phisinit', varid)
    CALL nf95_put_var(ncid, varid, phis)

    CALL nf95_close(ncid) ! fermer le fichier

    PRINT *, 'iim, jjm, llm, iday_end', iim, jjm, llm, iday_end
    PRINT *, 'rad, omeg, g, cpp, kappa', rad, omeg, g, cpp, kappa

  END SUBROUTINE dynredem0

END MODULE dynredem0_m
