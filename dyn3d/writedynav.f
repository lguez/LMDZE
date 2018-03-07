module writedynav_m

  implicit none

contains

  subroutine writedynav(vcov, ucov, teta, pk, phi, q, masse, ps, time)

    ! From LMDZ4/libf/bibio/writedynav.F, version 1.1.1.1 2004/05/19 12:53:05
    ! Écriture du fichier histoire au format IOIPSL
    ! L. Fairhead, LMD, 03/99

    USE comconst, ONLY: cpp
    use covnat_m, only: covnat
    use dimens_m, only: llm, nqmx
    use histsync_m, only: histsync
    use histwrite_m, only: histwrite
    use iniadvtrac_m, only: ttext
    use initdynav_m, only: histaveid
    use nr_util, only: assert
    use paramet_m, only: ip1jm, ip1jmp1
    use temps, only: itau_dyn

    ! Vent covariant :
    REAL, intent(in):: vcov(:, :, :) ! (iim + 1, jjm, llm)
    REAL, intent(in):: ucov(:, :, :) ! (iim + 1, jjm + 1, llm)

    REAL, intent(in):: teta(:, :, :) ! (iim + 1, jjm + 1, llm)
    ! temperature potentielle

    real, intent(in):: pk(:, :, :) ! (iim + 1, jjm + 1, llm)

    real, intent(in):: phi(:, :, :) ! (iim + 1, jjm + 1, llm)
    ! geopotentiel instantann\'e

    REAL, intent(in):: q(:, :, :, :) ! (iim + 1, jjm + 1, llm, nqmx) traceurs
    real, intent(in):: masse(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(in):: ps(:, :) ! (iim + 1, jjm + 1) pression au sol
    integer, intent(in):: time ! temps de l'ecriture

    ! Local:
    integer iq, itau_w
    real us(ip1jmp1*llm), vs(ip1jmp1*llm)
    REAL vnat(ip1jm, llm), unat(ip1jmp1, llm) 

    !---------------------------------------------------------------------

    call assert([size(vcov, 1), size(ucov, 1), size(teta, 1), size(phi, 1), &
         size(pk, 1), size(ps, 1), size(masse, 1)] == size(q, 1), &
         "writedynav iim")
    call assert([size(vcov, 2) + 1, size(ucov, 2), size(teta, 2), &
         size(phi, 2), size(pk, 2), size(ps, 2), size(masse, 2)] &
         == size(q, 2), "writedynav jjm")
    call assert([size(vcov, 3), size(ucov, 3), size(teta, 3), size(phi, 3), &
         size(pk, 3), size(masse, 3), size(q, 3)] == llm, "writedynav llm")
    call assert(size(q, 4) == nqmx, "writedynav nqmx")

    ! Initialisations
    us = 999.999
    vs = 999.999
    
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
    call histwrite(histaveid, 'temp', itau_w, teta * pk / cpp)

    call histwrite(histaveid, 'phi', itau_w, phi)

    ! Traceurs
    DO iq = 1, size(q, 4)
       call histwrite(histaveid, ttext(iq), itau_w, q(:, :, :, iq))
    enddo

    call histwrite(histaveid, 'masse', itau_w, masse)
    call histwrite(histaveid, 'ps', itau_w, ps)

    call histsync(histaveid)

  end subroutine writedynav

end module writedynav_m
