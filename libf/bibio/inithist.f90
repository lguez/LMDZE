module inithist_m

  implicit none

contains

  subroutine inithist(day0, anne0, tstep, nq, t_ops, t_wrt)

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

    ! L. Fairhead, LMD, 03/99

    USE calendar
    use com_io_dyn, only: histid, histvid, histuid
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

    call histbeg_totreg("dyn_histu.nc", rlong(:,1), rlat(1,:), 1, iip1, 1, &
         jjp1, itau_dyn, zjulian, tstep, uhoriid, histuid)

    do jj = 1, jjm
       do ii = 1, iip1
          rlong(ii, jj) = rlonv(ii) * 180. / pi
          rlat(ii, jj) = rlatv(jj) * 180. / pi
       enddo
    enddo

    ! Creation du fichier histoire pour la grille en V (oblige pour l'instant,
    ! IOIPSL ne permet pas de grilles avec des nombres de point differents dans
    ! un meme fichier)

    call histbeg_totreg('dyn_histv.nc', rlong(:, 1), rlat(1, :jjm), &
         1, iip1, 1, jjm, &
         itau_dyn, zjulian, tstep, vhoriid, histvid)
    !
    ! Appel a histhori pour rajouter les autres grilles horizontales
    !
    do jj = 1, jjp1
       do ii = 1, iip1
          rlong(ii, jj) = rlonv(ii) * 180. / pi
          rlat(ii, jj) = rlatu(jj) * 180. / pi
       enddo
    enddo

    call histbeg_totreg("dyn_hist.nc", rlong(:, 1), rlat(1, :), 1, iip1, 1, &
         jjp1, itau_dyn, zjulian, tstep, thoriid, histid)
    !
    ! Appel a histvert pour la grille verticale
    !
      call histvert(histid, 'presnivs', 'Niveaux pression','mb', llm, &
           presnivs/100., zvertiid,'down')
      call histvert(histvid, 'presnivs', 'Niveaux pression','mb', llm, &
           presnivs/100., zvertiid,'down')
      call histvert(histuid, 'presnivs', 'Niveaux pression','mb', llm, &
           presnivs/100., zvertiid,'down')
    !
    ! Appels a histdef pour la definition des variables a sauvegarder
    !
    ! Vents U
    !
    call histdef(histuid, 'u', 'vent u', 'm/s', &
         iip1, jjp1, uhoriid, llm, 1, llm, zvertiid, &
         'inst(X)', t_ops, t_wrt)
    !
    ! Vents V
    !
    call histdef(histvid, 'v', 'vent v', 'm/s', &
         iip1, jjm, vhoriid, llm, 1, llm, zvertiid, &
         'inst(X)', t_ops, t_wrt)

    !
    ! Temperature potentielle
    !
    call histdef(histid, 'teta', 'temperature potentielle', '-', &
         iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
         'inst(X)', t_ops, t_wrt)
    !
    ! Geopotentiel
    !
    call histdef(histid, 'phi', 'geopotentiel', '-', &
         iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
         'inst(X)', t_ops, t_wrt)
    !
    ! Traceurs
    !
    DO iq=1, nq
       call histdef(histid, ttext(iq), ttext(iq), '-', &
            iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
            'inst(X)', t_ops, t_wrt)
    enddo
    !
    ! Masse
    !
    call histdef(histid, 'masse', 'masse', 'kg', &
         iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
         'inst(X)', t_ops, t_wrt)
    !
    ! Pression au sol
    !
    call histdef(histid, 'ps', 'pression naturelle au sol', 'Pa', &
         iip1, jjp1, thoriid, 1, 1, 1, -99, &
         'inst(X)', t_ops, t_wrt)
    !
    ! Geopotentiel au sol
    !
    call histdef(histid, 'phis', 'geopotentiel au sol', '-', &
         iip1, jjp1, thoriid, 1, 1, 1, -99, &
         'inst(X)', t_ops, t_wrt)
    !
    ! Fin
    !
    call histend(histid)
    call histend(histuid)
    call histend(histvid)

  end subroutine inithist

end module inithist_m
