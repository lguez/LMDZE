module dynredem0_m

  IMPLICIT NONE

contains

  SUBROUTINE dynredem0(fichnom, iday_end, phis)

    ! From dyn3d/dynredem.F, v 1.2 2004/06/22 11:45:30

    ! Ecriture du fichier de redémarrage au format NetCDF (initialisation)

    USE IOIPSL, only: ymds2ju, ju2ymds
    use dimens_m, only: iim, jjm, llm, nqmx
    use paramet_m, only: ip1jmp1, iip1, jjp1, llmp1
    use comconst, only: rad, cpp, daysec, dtvr, kappa, g, omeg
    use comvert, only: pa, bp, ap, nivsigs, preff, presnivs, nivsig
    use logic
    use comgeom
    use serre
    use temps, only: annee_ref, day_ref, itaufin, itau_dyn
    use ener
    use advtrac_m, only: tname, ttext
    use netcdf95, only: nf95_create, nf95_put_att, nf95_def_dim, &
         nf95_def_var, NF95_ENDDEF, NF95_PUT_VAR
    use netcdf, only: NF90_CLOBBER, NF90_GLOBAL, NF90_UNLIMITED, nf90_float

    CHARACTER(len=*), intent(in):: fichnom
    INTEGER, intent(in):: iday_end
    REAL, intent(in):: phis(:, :)

    !   Local:

    include "netcdf.inc"

    INTEGER iq, l
    INTEGER length
    PARAMETER (length = 100)
    REAL tab_cntrl(length) ! tableau des parametres du run
    INTEGER ierr

    !   Variables locales pour NetCDF:

    INTEGER dims2(2), dims3(3), dims4(4)
    INTEGER idim_index
    INTEGER idim_rlonu, idim_rlonv, idim_rlatu, idim_rlatv
    INTEGER idim_s, idim_sig
    INTEGER idim_tim
    INTEGER nid, nvarid

    REAL zjulian, hours
    INTEGER yyears0, jjour0, mmois0
    character(len=30) unites

    !-----------------------------------------------------------------------

    print *, "Call sequence information: dynredem0"

    call ymds2ju(annee_ref, 1, iday_end, 0.0, zjulian)
    call ju2ymds(zjulian, yyears0, mmois0, jjour0, hours)

    DO l=1, length
       tab_cntrl(l) = 0.
    ENDDO
    tab_cntrl(1)  = REAL(iim)
    tab_cntrl(2)  = REAL(jjm)
    tab_cntrl(3)  = REAL(llm)
    tab_cntrl(4)  = REAL(day_ref)
    tab_cntrl(5)  = REAL(annee_ref)
    tab_cntrl(6)  = rad
    tab_cntrl(7)  = omeg
    tab_cntrl(8)  = g
    tab_cntrl(9)  = cpp
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

    tab_cntrl(20)  = clon
    tab_cntrl(21)  = clat
    tab_cntrl(22)  = grossismx
    tab_cntrl(23)  = grossismy

    IF ( fxyhypb )   THEN
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
       IF( ysinus )  tab_cntrl(27) = 1.
    ENDIF

    tab_cntrl(30) = REAL(iday_end)
    tab_cntrl(31) = REAL(itau_dyn + itaufin)

    call nf95_create(fichnom, NF90_CLOBBER, nid)
    call nf95_put_att(nid, NF90_GLOBAL, "title", &
         "Fichier de démarrage dynamique")

    ! Definir les dimensions du fichiers:

    call nf95_def_dim(nid, "index", length, idim_index)
    call NF95_DEF_DIM(nid, "rlonu", iip1, idim_rlonu)
    call NF95_DEF_DIM(nid, "rlatu", jjp1, idim_rlatu)
    call NF95_DEF_DIM(nid, "rlonv", iip1, idim_rlonv)
    call NF95_DEF_DIM(nid, "rlatv", jjm, idim_rlatv)
    call NF95_DEF_DIM(nid, "sigs", llm, idim_s)
    call NF95_DEF_DIM(nid, "sig", llmp1, idim_sig)
    call NF95_DEF_DIM(nid, "temps", NF90_UNLIMITED, idim_tim)

    ! Definir et enregistrer certains champs invariants:

    call nf95_def_var(nid, "controle", NF90_FLOAT, idim_index, nvarid)
    call nf95_put_att(nid, nvarid, "title", "Parametres de controle")
    call NF95_ENDDEF(nid)
    call NF95_PUT_VAR(nid, nvarid, tab_cntrl)

    ierr = NF_REDEF (nid)
    call nf95_def_var(nid, "rlonu", NF90_FLOAT, idim_rlonu, nvarid)
    call nf95_put_att(nid, nvarid, "title", "Longitudes des points U")
    call NF95_ENDDEF(nid)
    call NF95_PUT_VAR(nid, nvarid, rlonu)

    ierr = NF_REDEF (nid)
    call nf95_def_var(nid, "rlatu", NF90_FLOAT, idim_rlatu, nvarid)
    call nf95_put_att(nid, nvarid, "title", "Latitudes des points U")
    call NF95_ENDDEF(nid)
    call NF95_PUT_VAR(nid, nvarid, rlatu)

    ierr = NF_REDEF (nid)
    call nf95_def_var(nid, "rlonv", NF90_FLOAT, idim_rlonv, nvarid)
    call nf95_put_att(nid, nvarid, "title", "Longitudes des points V")
    call NF95_ENDDEF(nid)
    call NF95_PUT_VAR(nid, nvarid, rlonv)

    ierr = NF_REDEF (nid)
    call nf95_def_var(nid, "rlatv", NF90_FLOAT, idim_rlatv, nvarid)
    call nf95_put_att(nid, nvarid, "title", "Latitudes des points V")
    call NF95_ENDDEF(nid)
    call NF95_PUT_VAR(nid, nvarid, rlatv)

    ierr = NF_REDEF (nid)
    call nf95_def_var(nid, "nivsigs", NF90_FLOAT, idim_s, nvarid)
    call nf95_put_att(nid, nvarid, "title", "Numero naturel des couches s")
    call NF95_ENDDEF(nid)
    call NF95_PUT_VAR(nid, nvarid, nivsigs)

    ierr = NF_REDEF (nid)
    call nf95_def_var(nid, "nivsig", NF90_FLOAT, idim_sig, nvarid)
    call nf95_put_att(nid, nvarid, "title", "Numero naturel des couches sigma")
    call NF95_ENDDEF(nid)
    call NF95_PUT_VAR(nid, nvarid, nivsig)

    ierr = NF_REDEF (nid)
    call nf95_def_var(nid, "ap", NF90_FLOAT, idim_sig, nvarid)
    call nf95_put_att(nid, nvarid, "title", "Coefficient A pour hybride")
    call NF95_ENDDEF(nid)
    call NF95_PUT_VAR(nid, nvarid, ap)

    ierr = NF_REDEF (nid)
    call nf95_def_var(nid, "bp", NF90_FLOAT, idim_sig, nvarid)
    call nf95_put_att(nid, nvarid, "title", "Coefficient B pour hybride")
    call NF95_ENDDEF(nid)
    call NF95_PUT_VAR(nid, nvarid, bp)

    ierr = NF_REDEF (nid)
    call nf95_def_var(nid, "presnivs", NF90_FLOAT, idim_s, nvarid)
    call NF95_ENDDEF(nid)
    call NF95_PUT_VAR(nid, nvarid, presnivs)

    ! Coefficients de passage cov. <-> contra. <--> naturel

    ierr = NF_REDEF (nid)
    dims2(1) = idim_rlonu
    dims2(2) = idim_rlatu
    call nf95_def_var(nid, "cu", NF90_FLOAT, dims2, nvarid)
    call nf95_put_att(nid, nvarid, "title", "Coefficient de passage pour U")
    call NF95_ENDDEF(nid)
    call NF95_PUT_VAR(nid, nvarid, cu_2d)

    ierr = NF_REDEF (nid)
    dims2(1) = idim_rlonv
    dims2(2) = idim_rlatv
    call nf95_def_var(nid, "cv", NF90_FLOAT, dims2, nvarid)
    call nf95_put_att(nid, nvarid, "title", "Coefficient de passage pour V")
    call NF95_ENDDEF(nid)
    call NF95_PUT_VAR(nid, nvarid, cv_2d)

    ! Aire de chaque maille:

    ierr = NF_REDEF (nid)
    dims2(1) = idim_rlonv
    dims2(2) = idim_rlatu
    call nf95_def_var(nid, "aire", NF90_FLOAT, dims2, nvarid)
    call nf95_put_att(nid, nvarid, "title", "Aires de chaque maille")
    call NF95_ENDDEF(nid)
    call NF95_PUT_VAR(nid, nvarid, aire_2d)

    ! Geopentiel au sol:

    ierr = NF_REDEF (nid)
    dims2(1) = idim_rlonv
    dims2(2) = idim_rlatu
    call nf95_def_var(nid, "phisinit", NF90_FLOAT, dims2, nvarid)
    call nf95_put_att(nid, nvarid, "title", "Geopotentiel au sol")
    call NF95_ENDDEF(nid)
    call NF95_PUT_VAR(nid, nvarid, phis)

    ! Definir les variables pour pouvoir les enregistrer plus tard:

    ierr = NF_REDEF (nid) ! entrer dans le mode de definition

    call nf95_def_var(nid, "temps", NF90_FLOAT, idim_tim, nvarid)
    call nf95_put_att(nid, nvarid, "title", "Temps de simulation")
    write(unites, 200)yyears0, mmois0, jjour0
200 format('days since ', i4, '-', i2.2, '-', i2.2, ' 00:00:00')
    call nf95_put_att(nid, nvarid, "units", unites)


    dims4(1) = idim_rlonu
    dims4(2) = idim_rlatu
    dims4(3) = idim_s
    dims4(4) = idim_tim
    call nf95_def_var(nid, "ucov", NF90_FLOAT, dims4, nvarid)
    call nf95_put_att(nid, nvarid, "title", "Vitesse U")

    dims4(1) = idim_rlonv
    dims4(2) = idim_rlatv
    dims4(3) = idim_s
    dims4(4) = idim_tim
    call nf95_def_var(nid, "vcov", NF90_FLOAT, dims4, nvarid)
    call nf95_put_att(nid, nvarid, "title", "Vitesse V")

    dims4(1) = idim_rlonv
    dims4(2) = idim_rlatu
    dims4(3) = idim_s
    dims4(4) = idim_tim
    call nf95_def_var(nid, "teta", NF90_FLOAT, dims4, nvarid)
    call nf95_put_att(nid, nvarid, "title", "Temperature")

    dims4(1) = idim_rlonv
    dims4(2) = idim_rlatu
    dims4(3) = idim_s
    dims4(4) = idim_tim
    DO iq=1, nqmx
       call nf95_def_var(nid, tname(iq), NF90_FLOAT, dims4, nvarid)
       call nf95_put_att(nid, nvarid, "title", ttext(iq))
    ENDDO

    dims4(1) = idim_rlonv
    dims4(2) = idim_rlatu
    dims4(3) = idim_s
    dims4(4) = idim_tim
    call nf95_def_var(nid, "masse", NF90_FLOAT, dims4, nvarid)
    call nf95_put_att(nid, nvarid, "title", "C est quoi ?")

    dims3(1) = idim_rlonv
    dims3(2) = idim_rlatu
    dims3(3) = idim_tim
    call nf95_def_var(nid, "ps", NF90_FLOAT, dims3, nvarid)
    call nf95_put_att(nid, nvarid, "title", "Pression au sol")

    ierr = NF_ENDDEF(nid) ! sortir du mode de definition
    ierr = NF_CLOSE(nid) ! fermer le fichier

    PRINT*, 'iim, jjm, llm, iday_end', iim, jjm, llm, iday_end
    PRINT*, 'rad, omeg, g, cpp, kappa', rad, omeg, g, cpp, kappa

  END SUBROUTINE dynredem0

end module dynredem0_m
