module interfoce_lim_m

  implicit none

contains

  SUBROUTINE interfoce_lim(itime, dtime, jour,  &
       klon, nisurf, knon, knindex,  &
       debut,  &
       lmt_sst, pctsrf_new)

    ! Cette routine sert d'interface entre le modele atmospherique et
    ! un fichier de conditions aux limites

    ! L. Fairhead 02/2000

    use abort_gcm_m, only: abort_gcm
    use indicesol

    integer, intent(IN) :: itime ! numero du pas de temps courant
    real , intent(IN) :: dtime ! pas de temps de la physique (en s)
    integer, intent(IN) :: jour ! jour a lire dans l'annee
    integer, intent(IN) :: nisurf ! index de la surface a traiter (1 = sol continental)
    integer, intent(IN) :: knon ! nombre de points dans le domaine a traiter
    integer, intent(IN) :: klon ! taille de la grille
    integer, dimension(klon), intent(in) :: knindex ! index des points de la surface a traiter
    logical, intent(IN) :: debut ! logical: 1er appel a la physique (initialisation)

    ! Parametres de sortie
    ! output:
    ! lmt_sst SST lues dans le fichier de CL
    ! pctsrf_new sous-maille fractionnelle
    real, intent(out), dimension(klon) :: lmt_sst
    real, intent(out), dimension(klon, nbsrf) :: pctsrf_new

    ! Variables locales
    integer :: ii
    INTEGER, save :: lmt_pas ! frequence de lecture des conditions limites
    ! (en pas de physique)
    logical, save :: deja_lu ! pour indiquer que le jour a lire a deja
    ! lu pour une surface precedente
    integer, save :: jour_lu
    integer :: ierr
    character (len = 20) :: modname = 'interfoce_lim'
    character (len = 80) :: abort_message
    logical, save :: newlmt = .TRUE.
    logical, save :: check = .FALSE.
    ! Champs lus dans le fichier de CL
    real, allocatable , save, dimension(:) :: sst_lu, rug_lu, nat_lu
    real, allocatable , save, dimension(:, :) :: pct_tmp

    ! quelques variables pour netcdf

    include "netcdf.inc"
    integer :: nid, nvarid
    integer, dimension(2) :: start, epais

    ! --------------------------------------------------

    if (debut .and. .not. allocated(sst_lu)) then
       lmt_pas = nint(86400./dtime * 1.0) ! pour une lecture une fois par jour
       jour_lu = jour - 1
       allocate(sst_lu(klon))
       allocate(nat_lu(klon))
       allocate(pct_tmp(klon, nbsrf))
    endif

    if ((jour - jour_lu) /= 0) deja_lu = .false.

    if (check) write(*, *)modname, ' :: jour, jour_lu, deja_lu', jour, jour_lu, &
         deja_lu
    if (check) write(*, *)modname, ' :: itime, lmt_pas ', itime, lmt_pas, dtime

    ! Tester d'abord si c'est le moment de lire le fichier
    if (mod(itime-1, lmt_pas) == 0 .and. .not. deja_lu) then

       ! Ouverture du fichier

       ierr = NF_OPEN ('limit.nc', NF_NOWRITE, nid)
       if (ierr.NE.NF_NOERR) then
          abort_message &
               = 'Pb d''ouverture du fichier de conditions aux limites'
          call abort_gcm(modname, abort_message, 1)
       endif

       ! La tranche de donnees a lire:

       start(1) = 1
       start(2) = jour
       epais(1) = klon
       epais(2) = 1

       if (newlmt) then

          ! Fraction "ocean"

          ierr = NF_INQ_VARID(nid, 'FOCE', nvarid)
          if (ierr /= NF_NOERR) then
             abort_message = 'Le champ <FOCE> est absent'
             call abort_gcm(modname, abort_message, 1)
          endif
          ierr = NF_GET_VARA_REAL(nid, nvarid, start, epais, pct_tmp(1, is_oce))
          if (ierr /= NF_NOERR) then
             abort_message = 'Lecture echouee pour <FOCE>'
             call abort_gcm(modname, abort_message, 1)
          endif

          ! Fraction "glace de mer"

          ierr = NF_INQ_VARID(nid, 'FSIC', nvarid)
          if (ierr /= NF_NOERR) then
             abort_message = 'Le champ <FSIC> est absent'
             call abort_gcm(modname, abort_message, 1)
          endif
          ierr = NF_GET_VARA_REAL(nid, nvarid, start, epais, pct_tmp(1, is_sic))
          if (ierr /= NF_NOERR) then
             abort_message = 'Lecture echouee pour <FSIC>'
             call abort_gcm(modname, abort_message, 1)
          endif

          ! Fraction "terre"

          ierr = NF_INQ_VARID(nid, 'FTER', nvarid)
          if (ierr /= NF_NOERR) then
             abort_message = 'Le champ <FTER> est absent'
             call abort_gcm(modname, abort_message, 1)
          endif
          ierr = NF_GET_VARA_REAL(nid, nvarid, start, epais, pct_tmp(1, is_ter))
          if (ierr /= NF_NOERR) then
             abort_message = 'Lecture echouee pour <FTER>'
             call abort_gcm(modname, abort_message, 1)
          endif

          ! Fraction "glacier terre"

          ierr = NF_INQ_VARID(nid, 'FLIC', nvarid)
          if (ierr /= NF_NOERR) then
             abort_message = 'Le champ <FLIC> est absent'
             call abort_gcm(modname, abort_message, 1)
          endif
          ierr = NF_GET_VARA_REAL(nid, nvarid, start, epais, pct_tmp(1, is_lic))
          if (ierr /= NF_NOERR) then
             abort_message = 'Lecture echouee pour <FLIC>'
             call abort_gcm(modname, abort_message, 1)
          endif

       else ! on en est toujours a rnatur

          ierr = NF_INQ_VARID(nid, 'NAT', nvarid)
          if (ierr /= NF_NOERR) then
             abort_message = 'Le champ <NAT> est absent'
             call abort_gcm(modname, abort_message, 1)
          endif
          ierr = NF_GET_VARA_REAL(nid, nvarid, start, epais, nat_lu)
          if (ierr /= NF_NOERR) then
             abort_message = 'Lecture echouee pour <NAT>'
             call abort_gcm(modname, abort_message, 1)
          endif

          ! Remplissage des fractions de surface
          ! nat = 0, 1, 2, 3 pour ocean, terre, glacier, seaice

          pct_tmp = 0.0
          do ii = 1, klon
             pct_tmp(ii, nint(nat_lu(ii)) + 1) = 1.
          enddo


          ! On se retrouve avec ocean en 1 et terre en 2 alors qu'on veut le contraire

          pctsrf_new = pct_tmp
          pctsrf_new (:, 2)= pct_tmp (:, 1)
          pctsrf_new (:, 1)= pct_tmp (:, 2)
          pct_tmp = pctsrf_new
       endif ! fin test sur newlmt

       ! Lecture SST

       ierr = NF_INQ_VARID(nid, 'SST', nvarid)
       if (ierr /= NF_NOERR) then
          abort_message = 'Le champ <SST> est absent'
          call abort_gcm(modname, abort_message, 1)
       endif
       ierr = NF_GET_VARA_REAL(nid, nvarid, start, epais, sst_lu)
       if (ierr /= NF_NOERR) then
          abort_message = 'Lecture echouee pour <SST>'
          call abort_gcm(modname, abort_message, 1)
       endif


       ! Fin de lecture

       ierr = NF_CLOSE(nid)
       deja_lu = .true.
       jour_lu = jour
    endif

    ! Recopie des variables dans les champs de sortie

    lmt_sst = 999999999.
    do ii = 1, knon
       lmt_sst(ii) = sst_lu(knindex(ii))
    enddo

    pctsrf_new(:, is_oce) = pct_tmp(:, is_oce)
    pctsrf_new(:, is_sic) = pct_tmp(:, is_sic)

  END SUBROUTINE interfoce_lim

end module interfoce_lim_m
