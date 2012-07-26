module writedynav_m

  implicit none

contains

  subroutine writedynav(vcov, ucov, teta, ppk, phi, q, masse, ps, phis, time)

    ! From LMDZ4/libf/bibio/writedynav.F, version 1.1.1.1 2004/05/19 12:53:05
    ! Écriture du fichier histoire au format IOIPSL
    ! L. Fairhead, LMD, 03/99

    ! Appels successifs des routines histwrite

    use covnat_m, only: covnat
    USE histwrite_m, ONLY: histwrite
    USE histsync_m, ONLY: histsync
    USE dimens_m, ONLY: llm
    USE paramet_m, ONLY: iip1, ijp1llm, ip1jm, ip1jmp1, jjp1
    USE comconst, ONLY: cpp
    USE temps, ONLY: itau_dyn
    USE iniadvtrac_m, ONLY: ttext
    use initdynav_m, only: histaveid

    REAL, intent(in):: vcov(ip1jm, llm), ucov(ip1jmp1, llm) ! vents covariants
    REAL, intent(in):: teta(ip1jmp1*llm) ! temperature potentielle
    real, intent(in):: phi(ip1jmp1, llm) ! geopotentiel instantane
    real, intent(in):: ppk(ip1jmp1*llm) 
    REAL, intent(in):: ps(ip1jmp1) ! pression au sol
    real, intent(in):: masse(ip1jmp1, llm) 
    REAL, intent(in):: phis(ip1jmp1) ! geopotentiel au sol
    REAL, intent(in):: q(:, :, :, :) ! (iim + 1, jjm + 1, llm, nqmx) traceurs
    integer, intent(in):: time ! temps de l'ecriture

    ! Variables locales
    integer ndex2d(iip1*jjp1), ndex3d(iip1*jjp1*llm), iq, ii, ll
    real us(ip1jmp1*llm), vs(ip1jmp1*llm)
    real tm(ip1jmp1*llm)
    REAL vnat(ip1jm, llm), unat(ip1jmp1, llm) 
    logical ok_sync
    integer itau_w

    !---------------------------------------------------------------

    ! Initialisations
    ndex3d = 0
    ndex2d = 0
    ok_sync = .TRUE.
    us = 999.999
    vs = 999.999
    tm = 999.999
    vnat = 999.999
    unat = 999.999
    itau_w = itau_dyn + time

    ! Passage aux composantes naturelles du vent
    call covnat(llm, ucov, vcov, unat, vnat)

    ! Appels a histwrite pour l'ecriture des variables a sauvegarder

    ! Vents U scalaire
    call gr_u_scal(llm, unat, us)
    call histwrite(histaveid, 'u', itau_w, us)

    ! Vents V scalaire
    call gr_v_scal(llm, vnat, vs)
    call histwrite(histaveid, 'v', itau_w, vs)

    ! Temperature potentielle moyennee
    call histwrite(histaveid, 'theta', itau_w, teta)

    ! Temperature moyennee
    do ii = 1, ijp1llm
       tm(ii) = teta(ii) * ppk(ii)/cpp
    enddo
    call histwrite(histaveid, 'temp', itau_w, tm)

    ! Geopotentiel
    call histwrite(histaveid, 'phi', itau_w, phi)

    ! Traceurs
    DO iq = 1, size(q, 4)
       call histwrite(histaveid, ttext(iq), itau_w, q(:, :, :, iq))
    enddo

    ! Masse
    call histwrite(histaveid, 'masse', itau_w, masse)

    ! Pression au sol
    call histwrite(histaveid, 'ps', itau_w, ps)

    ! Geopotentiel au sol
    call histwrite(histaveid, 'phis', itau_w, phis)

    if (ok_sync) call histsync(histaveid)

  end subroutine writedynav

end module writedynav_m
