module inithist_m

  implicit none

contains

  subroutine inithist(tstep, nq, t_ops, t_wrt)

    ! From inithist.F, version 1.1.1.1 2004/05/19 12:53:05
    ! L. Fairhead, LMD, 03/99

    ! Routine d'initialisation des écritures des fichiers histoires
    ! LMDZ au format IOIPSL.

    USE com_io_dyn, ONLY: histid, histuid, histvid
    USE dimens_m, ONLY: jjm, llm
    USE disvert_m, ONLY: presnivs
    use dynetat0_m, only: day_ref, annee_ref, rlatu, rlatv, rlonu, rlonv
    USE histbeg_totreg_m, ONLY : histbeg_totreg
    USE histdef_m, ONLY : histdef
    USE histend_m, ONLY : histend
    USE histvert_m, ONLY : histvert
    USE iniadvtrac_m, ONLY: ttext
    USE nr_util, ONLY: pi
    USE paramet_m, ONLY: iip1, jjp1
    USE temps, ONLY: itau_dyn
    USE ymds2ju_m, ONLY: ymds2ju

    real, intent(in):: tstep ! durée du pas de temps en secondes
    integer, intent(in):: nq ! nombre de traceurs
    real, intent(in):: t_ops ! fréquence de l'opération pour IOIPSL
    real, intent(in):: t_wrt ! fréquence d'écriture sur le fichier

    ! Variables locales:
    real zjulian
    integer iq
    real rlong(iip1, jjp1), rlat(iip1, jjp1)
    integer uhoriid, vhoriid, thoriid, zvertiid
    integer ii, jj

    !-----------------------------------------------------------------------

    CALL ymds2ju(annee_ref, 1, day_ref, 0.0, zjulian)

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

    call histbeg_totreg('dyn_histv.nc', rlong(:, 1), rlat(1, :jjm), 1, iip1, &
         1, jjm, itau_dyn, zjulian, tstep, vhoriid, histvid)

    do jj = 1, jjp1
       do ii = 1, iip1
          rlong(ii, jj) = rlonv(ii) * 180. / pi
          rlat(ii, jj) = rlatu(jj) * 180. / pi
       enddo
    enddo

    call histbeg_totreg("dyn_hist.nc", rlong(:, 1), rlat(1, :), 1, iip1, 1, &
         jjp1, itau_dyn, zjulian, tstep, thoriid, histid)

    ! Appel a histvert pour la grille verticale

    call histvert(histid, 'presnivs', 'Niveaux pression','mb', presnivs/100., &
         zvertiid,'down')
    call histvert(histvid, 'presnivs', 'Niveaux pression','mb', &
         presnivs/100., zvertiid,'down')
    call histvert(histuid, 'presnivs', 'Niveaux pression','mb', &
         presnivs/100., zvertiid,'down')

    ! Appels a histdef pour la definition des variables a sauvegarder

    call histdef(histuid, 'u', 'vent u', 'm/s', iip1, jjp1, uhoriid, llm, 1, &
         llm, zvertiid, 'inst(X)', t_ops, t_wrt)
    call histdef(histvid, 'v', 'vent v', 'm/s', iip1, jjm, vhoriid, llm, 1, &
         llm, zvertiid, 'inst(X)', t_ops, t_wrt)
    call histdef(histid, 'teta', 'temperature potentielle', '-', iip1, jjp1, &
         thoriid, llm, 1, llm, zvertiid, 'inst(X)', t_ops, t_wrt)
    call histdef(histid, 'phi', 'geopotentiel', '-', iip1, jjp1, thoriid, &
         llm, 1, llm, zvertiid, 'inst(X)', t_ops, t_wrt)

    ! Traceurs
    DO iq=1, nq
       call histdef(histid, ttext(iq), ttext(iq), '-', iip1, jjp1, thoriid, &
            llm, 1, llm, zvertiid, 'inst(X)', t_ops, t_wrt)
    enddo

    call histdef(histid, 'masse', 'masse', 'kg', iip1, jjp1, thoriid, llm, 1, &
         llm, zvertiid, 'inst(X)', t_ops, t_wrt)
    call histdef(histid, 'ps', 'pression naturelle au sol', 'Pa', iip1, jjp1, &
         thoriid, 1, 1, 1, -99, 'inst(X)', t_ops, t_wrt)
    call histdef(histid, 'phis', 'geopotentiel au sol', '-', iip1, jjp1, &
         thoriid, 1, 1, 1, -99, 'inst(X)', t_ops, t_wrt)

    call histend(histid)
    call histend(histuid)
    call histend(histvid)

  end subroutine inithist

end module inithist_m
