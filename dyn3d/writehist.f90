module writehist_m

  implicit none

contains

  subroutine writehist(vcov, ucov, teta, pk, phi, q, masse, ps, itau_w)

    ! From writehist.F, revision 1403, 2010-07-01 09:02:53
    ! Écriture du fichier histoire au format IOIPSL
    ! L. Fairhead, LMD, 03/99

    use jumble, only: assert

    use covnat_m, only: covnat
    use dimensions, only: llm
    use histsync_m, only: histsync
    use histwrite_m, only: histwrite
    use infotrac_init_m, only: ttext
    use inithist_m, only: histid, histvid, histuid
    use paramet_m, only: ip1jm, ip1jmp1
    use suphec_m, only: rcpd

    ! Vent covariant :
    REAL, intent(in):: vcov(:, :, :) ! (iim + 1, jjm, llm)
    REAL, intent(in):: ucov(:, :, :) ! (iim + 1, jjm + 1, llm)

    REAL, intent(in):: teta(:, :, :) ! (iim + 1, jjm + 1, llm)
    ! temperature potentielle

    real, intent(in):: pk(:, :, :) ! (iim + 1, jjm + 1, llm)

    real, intent(in):: phi(:, :, :) ! (iim + 1, jjm + 1, llm)
    ! geopotential at mid-layer, in m2 s-2

    REAL, intent(in):: q(:, :, :, :) ! (iim + 1, jjm + 1, llm, nqmx) traceurs
    real, intent(in):: masse(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(in):: ps(:, :) ! (iim + 1, jjm + 1) pression au sol
    integer, intent(in):: itau_w ! temps de l'ecriture

    ! Local:
    integer iq
    REAL vnat(ip1jm, llm), unat(ip1jmp1, llm)

    !---------------------------------------------------------------------

    call assert([size(vcov, 1), size(ucov, 1), size(teta, 1), size(phi, 1), &
         size(pk, 1), size(ps, 1), size(masse, 1)] == size(q, 1), &
         "writehist iim")
    call assert([size(vcov, 2) + 1, size(ucov, 2), size(teta, 2), &
         size(phi, 2), size(pk, 2), size(ps, 2), size(masse, 2)] &
         == size(q, 2), "writehist jjm")
    call assert([size(vcov, 3), size(ucov, 3), size(teta, 3), size(phi, 3), &
         size(pk, 3), size(masse, 3), size(q, 3)] == llm, "writehist llm")

    call covnat(llm, ucov, vcov, unat, vnat)

    call histwrite(histuid, 'u', itau_w, unat)
    call histwrite(histvid, 'v', itau_w, vnat)
    call histwrite(histid, 'theta', itau_w, teta)
    call histwrite(histid, 'temp', itau_w, teta * pk / rcpd)
    call histwrite(histid, 'phi', itau_w, phi)

    DO iq = 1, size(q, 4)
       call histwrite(histid, ttext(iq), itau_w, q(:, :, :, iq))
    enddo

    call histwrite(histid, 'masse', itau_w, masse)
    call histwrite(histid, 'ps', itau_w, ps)

    call histsync(histid)
    call histsync(histvid)
    call histsync(histuid)

  end subroutine writehist

end module writehist_m
