module read_reanalyse_m

  IMPLICIT NONE

contains

  subroutine read_reanalyse(psi, u, v, t, q)

    ! From LMDZ4/libf/dyn3d/read_reanalyse.F, version 1.3, 2005/04/15 12:31:21

    USE conf_guide_m, ONLY: guide_q, guide_t, guide_u, guide_v, ncep
    USE dimens_m, ONLY: iim, jjm, llm
    USE netcdf, ONLY: nf90_nowrite
    USE netcdf95, ONLY: nf95_get_var, nf95_inq_dimid, nf95_inq_varid, &
         nf95_inquire_dimension, nf95_open
    USE paramet_m, ONLY: iip1, jjp1
    use reanalyse2nat_m, only: reanalyse2nat

    real, intent(in):: psi(:, :) ! (iip1, jjp1)
    real, intent(out):: u(:, :, :) ! (iip1, jjp1, llm)
    real, intent(out):: v(:, :, :) ! (iip1, jjm, llm)
    real, intent(out):: t(:, :, :), q(:, :, :) ! (iip1, jjp1, llm)

    ! Local:
    integer, save:: nlevnc
    integer:: timestep = 0
    real pk(iip1, jjp1, llm)
    integer, save:: ncidu, varidu, ncidv, varidv, ncidt, varidt, ncidQ, varidQ
    integer ncid, varid, dimid
    real, allocatable, save:: unc(:, :, :) ! (iip1, jjp1, nlevnc)
    real, allocatable, save:: vnc(:, :, :) ! (iip1, jjm, nlevnc)
    real, allocatable, save:: tnc(:, :, :), Qnc(:, :, :) ! (iip1, jjp1, nlevnc)
    real, allocatable, save:: pl(:) ! (nlevnc)
    logical:: first = .true.
    character(len = 8) name

    ! -----------------------------------------------------------------

    ! Initialisation de la lecture des fichiers

    if (first) then
       print *, 'Intitialisation de read reanalsye'

       ! Vent zonal
       if (guide_u) then
          call nf95_open('u.nc', nf90_nowrite, ncidu)
          call nf95_inq_varid(ncidu, 'UWND', varidu)
       endif

       ! Vent meridien
       if (guide_v) then
          call nf95_open('v.nc', nf90_nowrite, ncidv)
          call nf95_inq_varid(ncidv, 'VWND', varidv)
       endif

       ! Temperature
       if (guide_T) then
          call nf95_open('T.nc', nf90_nowrite, ncidt)
          call nf95_inq_varid(ncidt, 'AIR', varidt)
       endif

       ! Humidite
       if (guide_Q) then
          call nf95_open('hur.nc', nf90_nowrite, ncidQ)
          call nf95_inq_varid(ncidQ, 'RH', varidQ)
       endif

       ! Coordonn\'ee verticale :

       if (guide_u) then
          ncid = ncidu
       else if (guide_v) then
          ncid = ncidv
       else if (guide_T) then
          ncid = ncidt
       else
          ncid = ncidq
       end if

       name = merge('LEVEL   ', 'PRESSURE', ncep)
       call nf95_inq_dimid(ncid, name, dimid)
       call nf95_inquire_dimension(ncid, dimid, nclen = nlevnc)
       call nf95_inq_varid(ncid, name, varid)
       PRINT *, 'nlevnc = ', nlevnc
       allocate(unc(iip1, jjp1, nlevnc), vnc(iip1, jjm, nlevnc))
       allocate(tnc(iip1, jjp1, nlevnc), Qnc(iip1, jjp1, nlevnc), pl(nlevnc))
       call NF95_GET_VAR(ncid, varid, pl)
       pl = 100. * pl ! passage en pascal
       first = .false.
    endif

    ! lecture des champs u, v, T

    timestep = timestep + 1
    unc = 0.
    vnc = 0.
    tnc = 0.
    Qnc = 0.

    ! Vent zonal
    if (guide_u) then
       call NF95_GET_VAR(ncidu, varidu, unc, start = (/1, 1, 1, timestep/))
       call correctbid(iim, jjp1 * nlevnc, unc)
    endif

    ! Temperature
    if (guide_T) then
       call NF95_GET_VAR(ncidt, varidt, tnc, start = (/1, 1, 1, timestep/))
       call correctbid(iim, jjp1 * nlevnc, tnc)
    endif

    ! Humidite
    if (guide_Q) then
       call NF95_GET_VAR(ncidQ, varidQ, Qnc, start = (/1, 1, 1, timestep/))
       call correctbid(iim, jjp1 * nlevnc, Qnc)
    endif

    ! Vent meridien
    if (guide_v) then
       call NF95_GET_VAR(ncidv, varidv, vnc, start = (/1, 1, 1, timestep/))
       call correctbid(iim, jjm * nlevnc, vnc)
    endif

    call reanalyse2nat(nlevnc, psi, unc, vnc, tnc, Qnc, pl, u, v, t, Q, pk)
    call nat2gcm(u, v, t, Q, pk, u, v, t, Q)

  end subroutine read_reanalyse

end module read_reanalyse_m
