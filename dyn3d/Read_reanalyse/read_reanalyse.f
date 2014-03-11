module read_reanalyse_m

  IMPLICIT NONE

contains

  subroutine read_reanalyse(timestep, psi, u, v, t, q, masse, nlevnc)

    ! From LMDZ4/libf/dyn3d/read_reanalyse.F, version 1.3, 2005/04/15 12:31:21

    USE conf_guide_m, ONLY: guide_q, guide_t, guide_u, guide_v, ncep
    USE dimens_m, ONLY: iim, jjm, llm
    USE netcdf, ONLY: nf90_get_var, nf90_inq_varid, nf90_nowrite, nf90_open
    USE paramet_m, ONLY: iip1, jjp1
    use reanalyse2nat_m, only: reanalyse2nat

    integer timestep
    real, intent(in):: psi(iip1, jjp1)
    real u(iip1, jjp1, llm), v(iip1, jjm, llm)
    real t(iip1, jjp1, llm), q(iip1, jjp1, llm)
    real masse(iip1, jjp1, llm)
    integer nlevnc

    ! Local:

    integer l
    real pk(iip1, jjp1, llm)
    integer, save:: ncidu, varidu, ncidv, varidv, ncidt, varidt
    integer, save:: ncidpl
    integer, save:: varidpl, ncidQ, varidQ
    real unc(iip1, jjp1, nlevnc), vnc(iip1, jjm, nlevnc)
    real tnc(iip1, jjp1, nlevnc)
    real Qnc(iip1, jjp1, nlevnc)
    real pl(nlevnc)
    integer start(4), count(4), status
    real rcode
    logical:: first = .true.

    ! -----------------------------------------------------------------

    !   Initialisation de la lecture des fichiers

    if (first) then
       ncidpl=-99
       print *, 'Intitialisation de read reanalsye'

       ! Vent zonal
       if (guide_u) then
          rcode=nf90_open('u.nc', nf90_nowrite, ncidu)
          rcode = nf90_inq_varid(ncidu, 'UWND', varidu)
          print *, 'ncidu, varidu', ncidu, varidu
          if (ncidpl.eq.-99) ncidpl=ncidu
       endif

       ! Vent meridien
       if (guide_v) then
          rcode=nf90_open('v.nc', nf90_nowrite, ncidv)
          rcode = nf90_inq_varid(ncidv, 'VWND', varidv)
          print *, 'ncidv, varidv', ncidv, varidv
          if (ncidpl.eq.-99) ncidpl=ncidv
       endif

       ! Temperature
       if (guide_T) then
          rcode=nf90_open('T.nc', nf90_nowrite, ncidt)
          rcode = nf90_inq_varid(ncidt, 'AIR', varidt)
          print *, 'ncidt, varidt', ncidt, varidt
          if (ncidpl.eq.-99) ncidpl=ncidt
       endif

       ! Humidite
       if (guide_Q) then
          rcode=nf90_open('hur.nc', nf90_nowrite, ncidQ)
          rcode = nf90_inq_varid(ncidQ, 'RH', varidQ)
          print *, 'ncidQ, varidQ', ncidQ, varidQ
          if (ncidpl.eq.-99) ncidpl=ncidQ
       endif

       ! Coordonnee verticale
       if (ncep) then
          print *, 'Vous etes entrain de lire des donnees NCEP'
          rcode = nf90_inq_varid(ncidpl, 'LEVEL', varidpl)
       else
          print *, 'Vous etes entrain de lire des donnees ECMWF'
          rcode = nf90_inq_varid(ncidpl, 'PRESSURE', varidpl)
       endif
       print *, 'ncidu, varidpl', ncidu, varidpl
    endif
    print *, 'ok1'

    ! Niveaux de pression
    print *, 'WARNING!!! Il n y a pas de test de coherence'
    print *, 'sur le nombre de niveaux verticaux dans le fichier nc'
    status=NF90_GET_VAR(ncidpl, varidpl, pl)
    !  passage en pascal
    pl(:)=100.*pl(:)
    if (first) then
       do l=1, nlevnc
          print *, 'PL(', l, ')=', pl(l)
       enddo
    endif

    !   lecture des champs u, v, T

    !  dimensions pour les champs scalaires et le vent zonal

    start(1)=1
    start(2)=1
    start(3)=1
    start(4)=timestep

    count(1)=iip1
    count(2)=jjp1
    count(3)=nlevnc
    count(4)=1

    ! mise a zero des tableaux

    unc(:, :, :)=0.
    vnc(:, :, :)=0.
    tnc(:, :, :)=0.
    Qnc(:, :, :)=0.

    !  Vent zonal

    if (guide_u) then
       print *, 'avant la lecture de UNCEP nd de niv:', nlevnc
       status=NF90_GET_VAR(ncidu, varidu, unc, start, count)
       print *, 'WARNING!!! Correction bidon pour palier a un '
       print *, 'probleme dans la creation des fichiers nc'
       call correctbid(iim, jjp1*nlevnc, unc)
       call dump2d(iip1, jjp1, unc, 'UNC COUCHE 1 ')
    endif

    !  Temperature

    print *, 'ncidt=', ncidt, 'varidt=', varidt, 'start=', start
    print *, 'count=', count
    if (guide_T) then
       status=NF90_GET_VAR(ncidt, varidt, tnc, start, count)
       call dump2d(iip1, jjp1, tnc, 'TNC COUCHE 1 AAA ')
       call correctbid(iim, jjp1*nlevnc, tnc)
       call dump2d(iip1, jjp1, tnc, 'TNC COUCHE 1 BBB ')
    endif

    !  Humidite

    if (guide_Q) then
       status=NF90_GET_VAR(ncidQ, varidQ, Qnc, start, count)
       call correctbid(iim, jjp1*nlevnc, Qnc)
       call dump2d(iip1, jjp1, Qnc, 'QNC COUCHE 1 ')
    endif

    count(2)=jjm
    !  Vent meridien

    if (guide_v) then
       status=NF90_GET_VAR(ncidv, varidv, vnc, start, count)
       call correctbid(iim, jjm*nlevnc, vnc)
       call dump2d(iip1, jjm, vnc, 'VNC COUCHE 1 ')
    endif

    start(3)=timestep
    start(4)=0
    count(2)=jjp1
    count(3)=1
    count(4)=0

    !  Interpolation verticale sur les niveaux modele

    call reanalyse2nat(nlevnc, psi, unc, vnc, tnc, Qnc, pl, u, v, t, Q, &
         masse, pk)

    call dump2d(iip1, jjm, v, 'V COUCHE APRES ')

    !  Passage aux variables du modele (vents covariants, temperature
    !  potentielle et humidite specifique)

    call nat2gcm(u, v, t, Q, pk, u, v, t, Q)
    print *, 'TIMESTEP ', timestep
    first=.false.

  end subroutine read_reanalyse

end module read_reanalyse_m
