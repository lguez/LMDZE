module inithist_m

  implicit none

  integer histid, histvid, histuid

contains

  subroutine inithist(t_ops, t_wrt)

    ! From inithist.F, version 1.1.1.1, 2004/05/19 12:53:05
    ! L. Fairhead, LMD, 03/99

    ! Routine d'initialisation des écritures des fichiers histoires au
    ! format IOIPSL.

    use comconst, only: dtvr
    USE dimens_m, ONLY: jjm, llm, nqmx
    USE disvert_m, ONLY: presnivs
    use dynetat0_m, only: day_ref, annee_ref, rlatu, rlatv, rlonu, rlonv
    use histbeg_totreg_m, only: histbeg_totreg
    USE histdef_m, ONLY: histdef
    USE histend_m, ONLY: histend
    USE histvert_m, ONLY: histvert
    USE iniadvtrac_m, ONLY: ttext
    USE nr_util, ONLY: pi
    USE paramet_m, ONLY: iip1, jjp1
    USE temps, ONLY: itau_dyn
    use ymds2ju_m, ONLY: ymds2ju

    real, intent(in):: t_ops ! fréquence de l'opération pour IOIPSL
    real, intent(in):: t_wrt ! fréquence d'écriture sur le fichier

    ! Local:
    real julian
    integer iq
    integer uhoriid, vhoriid, thoriid, zvertiid

    !-----------------------------------------------------------------------

    print *, "Call sequence information: inithist"
    CALL ymds2ju(annee_ref, 1, day_ref, 0., julian)

    call histbeg_totreg("dyn_histu.nc", rlonu * 180. / pi, rlatu * 180. / pi, &
         1, iip1, 1, jjp1, itau_dyn, julian, dtvr, uhoriid, histuid)

    ! Creation du fichier histoire pour la grille en V (oblige pour l'instant,
    ! IOIPSL ne permet pas de grilles avec des nombres de point differents dans
    ! un meme fichier)

    call histbeg_totreg('dyn_histv.nc', rlonv * 180. / pi, rlatv * 180. / pi, &
         1, iip1, 1, jjm, itau_dyn, julian, dtvr, vhoriid, histvid)

    call histbeg_totreg("dyn_hist.nc", rlonv * 180. / pi, rlatu * 180. / pi, &
         1, iip1, 1, jjp1, itau_dyn, julian, dtvr, thoriid, histid)

    call histvert(histid, 'presnivs', 'Niveaux pression', 'mb', presnivs/100., &
         zvertiid, 'down')
    call histvert(histvid, 'presnivs', 'Niveaux pression', 'mb', &
         presnivs/100., zvertiid, 'down')
    call histvert(histuid, 'presnivs', 'Niveaux pression', 'mb', &
         presnivs/100., zvertiid, 'down')

    call histdef(histuid, 'u', 'vent u', 'm/s', iip1, jjp1, uhoriid, llm, 1, &
         llm, zvertiid, 'inst(X)', t_ops, t_wrt)
    call histdef(histvid, 'v', 'vent v', 'm/s', iip1, jjm, vhoriid, llm, 1, &
         llm, zvertiid, 'inst(X)', t_ops, t_wrt)
    call histdef(histid, 'temp', 'temperature', 'K', iip1, jjp1, &
         thoriid, llm, 1, llm, zvertiid, 'inst(X)', t_ops, t_wrt)
    call histdef(histid, 'theta', 'temperature potentielle', 'K', iip1, jjp1, &
         thoriid, llm, 1, llm, zvertiid, 'inst(X)', t_ops, t_wrt)
    call histdef(histid, 'phi', 'geopotentiel', '-', iip1, jjp1, thoriid, &
         llm, 1, llm, zvertiid, 'inst(X)', t_ops, t_wrt)

    ! Traceurs
    DO iq = 1, nqmx
       call histdef(histid, ttext(iq), ttext(iq), '-', iip1, jjp1, thoriid, &
            llm, 1, llm, zvertiid, 'inst(X)', t_ops, t_wrt)
    enddo

    call histdef(histid, 'masse', 'masse', 'kg', iip1, jjp1, thoriid, llm, 1, &
         llm, zvertiid, 'inst(X)', t_ops, t_wrt)
    call histdef(histid, 'ps', 'pression naturelle au sol', 'Pa', iip1, jjp1, &
         thoriid, 1, 1, 1, -99, 'inst(X)', t_ops, t_wrt)

    call histend(histid)
    call histend(histuid)
    call histend(histvid)

  end subroutine inithist

end module inithist_m
