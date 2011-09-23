module initdynav_m

  implicit none

contains

  subroutine initdynav(day0, anne0, tstep, nq, fileid, t_ops, t_wrt)

    ! From initdynav.F, version 1.1.1.1, 2004/05/19 12:53:05
    ! L. Fairhead, LMD
    ! Routine d'initialisation des écritures des fichiers histoires LMDZ
    ! au format IOIPSL. Initialisation du fichier histoire moyenne.

    use calendar, ONLY: ymds2ju
    USE comvert, ONLY : nivsigs
    USE comgeom, ONLY : rlatu, rlonv
    USE dimens_m, ONLY : llm
    USE histcom, ONLY: histbeg_totreg, histdef, histend, histvert
    USE iniadvtrac_m, ONLY : ttext
    USE nr_util, ONLY : pi
    USE paramet_m, ONLY : iip1, jjp1
    USE temps, ONLY : itau_dyn

    integer, intent(in):: day0, anne0 ! date de reference
    real, intent(in):: tstep ! frequence d'ecriture
    real, intent(in):: t_ops ! frequence de l'operation pour IOIPSL
    real, intent(in):: t_wrt ! frequence d'ecriture sur le fichier
    integer, intent(out):: fileid ! ID du fichier netcdf cree
    integer, intent(in):: nq ! nombre de traceurs

    ! Variables locales
    integer thoriid, zvertiid
    real zjulian
    integer iq
    real rlong(iip1, jjp1), rlat(iip1, jjp1)
    integer ii, jj
    integer zan, dayref

    !----------------------------------------------------

    print *, "Call sequence information: initdynav"

    zan = anne0
    dayref = day0
    CALL ymds2ju(zan, 1, dayref, 0.0, zjulian)

    do jj = 1, jjp1
       do ii = 1, iip1
          rlong(ii, jj) = rlonv(ii) * 180. / pi
          rlat(ii, jj) = rlatu(jj) * 180. / pi
       enddo
    enddo

    call histbeg_totreg('dyn_hist_ave.nc', rlong(:, 1), rlat(1, :), 1, iip1, &
         1, jjp1, itau_dyn, zjulian, tstep, thoriid, fileid)
    call histvert(fileid, 'sigss', 'Niveaux sigma', 'Pa', llm, nivsigs, &
         zvertiid)

    call histdef(fileid, 'u', 'vents u scalaires moyennes', &
         'm/s', iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
         'ave(X)', t_ops, t_wrt)
    call histdef(fileid, 'v', 'vents v scalaires moyennes', &
         'm/s', iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
         'ave(X)', t_ops, t_wrt)
    call histdef(fileid, 'temp', 'temperature moyennee', 'K', &
         iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
         'ave(X)', t_ops, t_wrt)
    call histdef(fileid, 'theta', 'temperature potentielle', 'K', &
         iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
         'ave(X)', t_ops, t_wrt)
    call histdef(fileid, 'phi', 'geopotentiel moyenne', '-', &
         iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
         'ave(X)', t_ops, t_wrt)

    ! Traceurs
    DO iq = 1, nq
       call histdef(fileid, ttext(iq), ttext(iq), '-', &
            iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
            'ave(X)', t_ops, t_wrt)
    enddo

    call histdef(fileid, 'masse', 'masse', 'kg', &
         iip1, jjp1, thoriid, 1, 1, 1, -99, &
         'ave(X)', t_ops, t_wrt)
    call histdef(fileid, 'ps', 'pression naturelle au sol', 'Pa', &
         iip1, jjp1, thoriid, 1, 1, 1, -99, &
         'ave(X)', t_ops, t_wrt)
    call histdef(fileid, 'phis', 'geopotentiel au sol', '-', &
         iip1, jjp1, thoriid, 1, 1, 1, -99, &
         'ave(X)', t_ops, t_wrt)

    call histend(fileid)

  end subroutine initdynav

end module initdynav_m
