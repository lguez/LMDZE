module inithist_m

  ! This module is clean: no C preprocessor directive, no include line

  implicit none

contains

  subroutine inithist(day0, anne0, tstep, nq, fileid, filevid, t_ops, t_wrt)

    ! From inithist.F, version 1.1.1.1 2004/05/19 12:53:05

    ! Routine d'initialisation des écritures des fichiers histoires LMDZ
    ! au format IOIPSL
    ! Appels successifs des routines : histbeg, histhori, histver,
    ! histdef, histend

    ! Entrées :
    ! day0, anne0: date de référence
    ! tstep : durée du pas de temps en secondes
    ! t_ops : fréquence de l'opération pour IOIPSL
    ! t_wrt : fréquence d'écriture sur le fichier
    ! nq : nombre de traceurs

    ! Sorties :
    ! fileid : ID du fichier Netcdf créé
    ! filevid : ID du fichier Netcdf pour la grille v

    ! L. Fairhead, LMD, 03/99

    USE calendar
    use histcom
    use dimens_m
    use paramet_m
    use comconst
    use comvert
    use logic
    use comgeom
    use serre
    use temps
    use ener
    use iniadvtrac_m
    use nr_util, only: pi

    ! Arguments
    integer day0, anne0
    real, intent(in):: tstep, t_ops, t_wrt
    integer fileid, filevid
    integer nq

    ! Variables locales
    real zjulian
    integer iq
    real rlong(iip1, jjp1), rlat(iip1, jjp1)
    integer uhoriid, vhoriid, thoriid, zvertiid
    integer ii, jj
    integer zan, dayref

    !-----------------------------------------------------------------------

    ! Appel a histbeg: creation du fichier netcdf et initialisations diverses

    zan = anne0
    dayref = day0
    CALL ymds2ju(zan, 1, dayref, 0.0, zjulian)

    do jj = 1, jjp1
       do ii = 1, iip1
          rlong(ii, jj) = rlonu(ii) * 180. / pi
          rlat(ii, jj) = rlatu(jj) * 180. / pi
       enddo
    enddo

    call histbeg_totreg("dyn_hist.nc", rlong(:, 1), rlat(1, :), &
         1, iip1, 1, jjp1, &
         itau_dyn, zjulian, tstep, uhoriid, fileid)
    !
    ! Creation du fichier histoire pour la grille en V (oblige pour l'instant,
    ! IOIPSL ne permet pas de grilles avec des nombres de point differents dans
    ! un meme fichier)

    do jj = 1, jjm
       do ii = 1, iip1
          rlong(ii, jj) = rlonv(ii) * 180. / pi
          rlat(ii, jj) = rlatv(jj) * 180. / pi
       enddo
    enddo

    call histbeg_totreg('dyn_histv.nc', rlong(:, 1), rlat(1, :jjm), &
         1, iip1, 1, jjm, &
         itau_dyn, zjulian, tstep, vhoriid, filevid)
    !
    ! Appel a histhori pour rajouter les autres grilles horizontales
    !
    do jj = 1, jjp1
       do ii = 1, iip1
          rlong(ii, jj) = rlonv(ii) * 180. / pi
          rlat(ii, jj) = rlatu(jj) * 180. / pi
       enddo
    enddo

    call histhori_regular(fileid, iip1, rlong, jjp1, rlat, 'scalar', &
         'Grille points scalaires', thoriid)
    !
    ! Appel a histvert pour la grille verticale
    !
    call histvert(fileid, 'sig_s', 'Niveaux sigma', '-', &
         llm, nivsigs, zvertiid)
    ! Pour le fichier V
    call histvert(filevid, 'sig_s', 'Niveaux sigma', '-', &
         llm, nivsigs, zvertiid)
    !
    ! Appels a histdef pour la definition des variables a sauvegarder
    !
    ! Vents U
    !
    call histdef(fileid, 'ucov', 'vents u covariants', 'm/s', &
         iip1, jjp1, uhoriid, llm, 1, llm, zvertiid, &
         'inst(X)', t_ops, t_wrt)
    !
    ! Vents V
    !
    call histdef(filevid, 'vcov', 'vents v covariants', 'm/s', &
         iip1, jjm, vhoriid, llm, 1, llm, zvertiid, &
         'inst(X)', t_ops, t_wrt)

    !
    ! Temperature potentielle
    !
    call histdef(fileid, 'teta', 'temperature potentielle', '-', &
         iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
         'inst(X)', t_ops, t_wrt)
    !
    ! Geopotentiel
    !
    call histdef(fileid, 'phi', 'geopotentiel instantane', '-', &
         iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
         'inst(X)', t_ops, t_wrt)
    !
    ! Traceurs
    !
    DO iq=1, nq
       call histdef(fileid, ttext(iq), ttext(iq), '-', &
            iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
            'inst(X)', t_ops, t_wrt)
    enddo
    !
    ! Masse
    !
    call histdef(fileid, 'masse', 'masse', 'kg', &
         iip1, jjp1, thoriid, 1, 1, 1, -99, &
         'inst(X)', t_ops, t_wrt)
    !
    ! Pression au sol
    !
    call histdef(fileid, 'ps', 'pression naturelle au sol', 'Pa', &
         iip1, jjp1, thoriid, 1, 1, 1, -99, &
         'inst(X)', t_ops, t_wrt)
    !
    ! Pression au sol
    !
    call histdef(fileid, 'phis', 'geopotentiel au sol', '-', &
         iip1, jjp1, thoriid, 1, 1, 1, -99, &
         'inst(X)', t_ops, t_wrt)
    !
    ! Fin
    !
    call histend(fileid)
    call histend(filevid)

  end subroutine inithist

end module inithist_m
