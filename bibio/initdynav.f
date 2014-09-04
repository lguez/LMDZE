module initdynav_m

  implicit none

  integer histaveid

contains

  subroutine initdynav(day0, anne0, tstep, nq, t_ops, t_wrt)

    ! From initdynav.F, version 1.1.1.1, 2004/05/19 12:53:05
    ! L. Fairhead, LMD

    ! Routine d'initialisation des écritures des fichiers histoires au
    ! format IOIPSL. Initialisation du fichier histoire moyenne.

    USE comgeom, ONLY: rlatu, rlonv
    USE dimens_m, ONLY: llm
    USE histbeg_totreg_m, ONLY: histbeg_totreg
    USE histdef_m, ONLY: histdef
    USE histend_m, ONLY: histend
    USE histvert_m, ONLY: histvert
    USE iniadvtrac_m, ONLY: ttext
    USE nr_util, ONLY: pi
    USE paramet_m, ONLY: iip1, jjp1
    USE temps, ONLY: itau_dyn
    use ymds2ju_m, ONLY: ymds2ju

    integer, intent(in):: day0, anne0 ! date de référence
    real, intent(in):: tstep ! fréquence d'écriture
    integer, intent(in):: nq ! nombre de traceurs
    real, intent(in):: t_ops ! fréquence de l'opération pour IOIPSL
    real, intent(in):: t_wrt ! fréquence d'écriture sur le fichier

    ! Variables locales
    integer horiid, zvertiid
    real julian
    integer iq, l

    !----------------------------------------------------

    print *, "Call sequence information: initdynav"

    CALL ymds2ju(anne0, 1, day0, 0., julian)
    call histbeg_totreg('dyn_hist_ave.nc', rlonv * 180. / pi, &
         rlatu * 180. / pi, 1, iip1, 1, jjp1, itau_dyn, julian, tstep, &
         horiid, histaveid)
    call histvert(histaveid, 'sigss', 'Niveaux sigma', '', &
         (/(real(l), l = 1, llm)/), zvertiid)

    call histdef(histaveid, 'u', 'vents u scalaires moyennes', 'm/s', iip1, &
         jjp1, horiid, llm, 1, llm, zvertiid, 'ave(X)', t_ops, t_wrt)
    call histdef(histaveid, 'v', 'vents v scalaires moyennes', 'm/s', iip1, &
         jjp1, horiid, llm, 1, llm, zvertiid, 'ave(X)', t_ops, t_wrt)
    call histdef(histaveid, 'temp', 'temperature moyennee', 'K', iip1, jjp1, &
         horiid, llm, 1, llm, zvertiid, 'ave(X)', t_ops, t_wrt)
    call histdef(histaveid, 'theta', 'temperature potentielle', 'K', iip1, &
         jjp1, horiid, llm, 1, llm, zvertiid, 'ave(X)', t_ops, t_wrt)
    call histdef(histaveid, 'phi', 'geopotentiel moyenne', '-', iip1, jjp1, &
         horiid, llm, 1, llm, zvertiid, 'ave(X)', t_ops, t_wrt)

    ! Traceurs
    DO iq = 1, nq
       call histdef(histaveid, ttext(iq), ttext(iq), '-', iip1, jjp1, &
            horiid, llm, 1, llm, zvertiid, 'ave(X)', t_ops, t_wrt)
    enddo

    call histdef(histaveid, 'masse', 'masse', 'kg', iip1, jjp1, horiid, 1, &
         1, 1, -99, 'ave(X)', t_ops, t_wrt)
    call histdef(histaveid, 'ps', 'pression naturelle au sol', 'Pa', iip1, &
         jjp1, horiid, 1, 1, 1, -99, 'ave(X)', t_ops, t_wrt)
    call histdef(histaveid, 'phis', 'geopotentiel au sol', '-', iip1, jjp1, &
         horiid, 1, 1, 1, -99, 'ave(X)', t_ops, t_wrt)

    call histend(histaveid)

  end subroutine initdynav

end module initdynav_m
