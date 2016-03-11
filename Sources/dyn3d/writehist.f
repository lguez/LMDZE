module writehist_m

  implicit none

contains

  subroutine writehist(time, vcov, ucov, teta, phi, masse, ps)

    ! From writehist.F 1403 2010-07-01 09:02:53Z
    ! Écriture du fichier histoire au format IOIPSL
    ! Appels successifs des routines histwrite
    ! L. Fairhead, LMD, 03/99

    use dimens_m, only: llm
    use com_io_dyn, only: histid, histvid, histuid
    use paramet_m, only: ip1jm, ip1jmp1
    use temps, only: itau_dyn
    use histwrite_m, only: histwrite
    use histsync_m, only: histsync
    use covnat_m, only: covnat

    ! Entree:
    ! time: temps de l'ecriture
    ! vcov: vents v covariants
    ! ucov: vents u covariants
    ! teta: temperature potentielle
    ! phi : geopotentiel instantane
    ! masse: masse
    ! ps :pression au sol

    ! Arguments

    REAL vcov(ip1jm, llm), ucov(ip1jmp1, llm) 
    REAL teta(ip1jmp1, llm), phi(ip1jmp1, llm) 
    REAL, intent(in):: ps(ip1jmp1), masse(ip1jmp1, llm) 
    integer time

    ! This routine needs IOIPSL to work
    ! Variables locales

    logical ok_sync
    integer itau_w
    REAL vnat(ip1jm, llm), unat(ip1jmp1, llm)

    !---------------------------------------------------------------------

    ! Initialisations

    ok_sync =.TRUE.
    itau_w = itau_dyn + time
    ! Passage aux composantes naturelles du vent
    call covnat(llm, ucov, vcov, unat, vnat)

    ! Appels a histwrite pour l'ecriture des variables a sauvegarder

    ! Vents U

    call histwrite(histuid, 'u', itau_w, unat)

    ! Vents V

    call histwrite(histvid, 'v', itau_w, vnat)

    ! Temperature potentielle

    call histwrite(histid, 'teta', itau_w, teta)

    ! Geopotentiel

    call histwrite(histid, 'phi', itau_w, phi)

    ! Masse

    call histwrite(histid, 'masse', itau_w, masse)

    ! Pression au sol
    call histwrite(histid, 'ps', itau_w, ps)

    if (ok_sync) then
       call histsync(histid)
       call histsync(histvid)
       call histsync(histuid)
    endif

  end subroutine writehist

end module writehist_m
