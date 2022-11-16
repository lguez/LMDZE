module read_reanalyse_m

  IMPLICIT NONE

contains

  subroutine read_reanalyse(psi, u, v, t, q)

    ! From LMDZ4/libf/dyn3d/read_reanalyse.F, version 1.3, 2005/04/15 12:31:21

    USE conf_guide_m, ONLY: guide_q, guide_t, guide_u, guide_v
    USE dimensions, ONLY: jjm, llm
    use nat2gcm_m, only: nat2gcm
    USE netcdf, ONLY: nf90_nowrite
    USE netcdf95, ONLY: nf95_get_var, nf95_inq_varid, nf95_inquire_dimension, &
         nf95_open, nf95_find_coord
    USE paramet_m, ONLY: iip1, jjp1
    use reanalyse2nat_m, only: reanalyse2nat

    real, intent(in):: psi(:, :) ! (iip1, jjp1)
    real, intent(out):: u(:, :, :) ! (iip1, jjp1, llm)
    real, intent(out):: v(:, :, :) ! (iip1, jjm, llm)
    real, intent(out):: t(:, :, :), q(:, :, :) ! (iip1, jjp1, llm)

    ! Local:
    integer nlevnc
    integer:: timestep = 0
    real pk(iip1, jjp1, llm)
    integer, save:: ncidu, varidu, ncidv, varidv, ncidt, varidt, ncidQ, varidQ
    integer ncid, varid, dimid
    real, allocatable, save:: unc(:, :, :) ! (iip1, jjp1, nlevnc)
    real, allocatable, save:: vnc(:, :, :) ! (iip1, jjm, nlevnc)
    real, allocatable, save:: tnc(:, :, :), Qnc(:, :, :) ! (iip1, jjp1, nlevnc)
    real, allocatable, save:: pl(:) ! (nlevnc)
    real latitude(jjm + 1)
    logical:: first = .true.
    logical, save:: invert_y

    ! -----------------------------------------------------------------

    ! Initialisation de la lecture des fichiers

    if (first) then
       print *, 'Intitialisation de read reanalyse'

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

       call nf95_find_coord(ncid, dimid = dimid, varid = varid, &
            std_name = "plev")
       call nf95_inquire_dimension(ncid, dimid, nclen = nlevnc)
       PRINT *, 'nlevnc = ', nlevnc
       allocate(unc(iip1, jjp1, nlevnc), vnc(iip1, jjm, nlevnc))
       allocate(tnc(iip1, jjp1, nlevnc), Qnc(iip1, jjp1, nlevnc), pl(nlevnc))
       call NF95_GET_VAR(ncid, varid, pl)
       pl = 100. * pl ! passage en pascal

       ! Read latitude values just to know their order:
       call nf95_find_coord(ncid, varid = varid, std_name = "latitude")
       call nf95_get_var(ncid, varid, latitude)
       invert_y = latitude(1) < latitude(2)

       first = .false.
    endif

    ! lecture des champs u, v, T, q

    timestep = timestep + 1

    ! Vent zonal
    if (guide_u) then
       call NF95_GET_VAR(ncidu, varidu, unc, start = (/1, 1, 1, timestep/))
    else
       unc = 0.
    end if

    ! Temperature
    if (guide_T) then
       call NF95_GET_VAR(ncidt, varidt, tnc, start = (/1, 1, 1, timestep/))
    else
       tnc = 0.
    end if

    ! Humidite
    if (guide_Q) then
       call NF95_GET_VAR(ncidQ, varidQ, Qnc, start = (/1, 1, 1, timestep/))
    else
       Qnc = 0.
    end if

    ! Vent meridien
    if (guide_v) then
       call NF95_GET_VAR(ncidv, varidv, vnc, start = (/1, 1, 1, timestep/))
    else
       vnc = 0.
    end if

    call reanalyse2nat(invert_y, psi, unc, vnc, tnc, Qnc, pl, u, v, t, Q, pk)
    call nat2gcm(pk, u, v, t)

  end subroutine read_reanalyse

end module read_reanalyse_m
