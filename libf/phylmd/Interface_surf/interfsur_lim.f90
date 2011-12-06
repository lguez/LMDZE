module interfsur_lim_m

  implicit none

contains

  SUBROUTINE interfsur_lim(itime, dtime, jour,  &
       klon, nisurf, knon, knindex,  &
       debut,  &
       lmt_alb, lmt_rug)

    ! Cette routine sert d'interface entre le modèle atmosphérique et
    ! un fichier de conditions aux limites.

    ! L. Fairhead 02/2000

    use abort_gcm_m, only: abort_gcm

    ! Parametres d'entree
    ! input:
    ! itime numero du pas de temps courant
    ! dtime pas de temps de la physique (en s)
    ! jour jour a lire dans l'annee
    ! nisurf index de la surface a traiter (1 = sol continental)
    ! knon nombre de points dans le domaine a traiter
    ! knindex index des points de la surface a traiter
    ! klon taille de la grille
    ! debut logical: 1er appel a la physique (initialisation)
    integer, intent(IN) :: itime
    real , intent(IN) :: dtime
    integer, intent(IN) :: jour
    integer, intent(IN) :: nisurf
    integer, intent(IN) :: knon
    integer, intent(IN) :: klon
    integer, dimension(klon), intent(in) :: knindex
    logical, intent(IN) :: debut

    ! Parametres de sortie
    ! output:
    ! lmt_sst SST lues dans le fichier de CL
    ! lmt_alb Albedo lu
    ! lmt_rug longueur de rugosité lue
    ! pctsrf_new sous-maille fractionnelle
    real, intent(out), dimension(klon) :: lmt_alb
    real, intent(out), dimension(klon) :: lmt_rug

    ! Variables locales
    integer :: ii
    integer, save :: lmt_pas ! frequence de lecture des conditions limites
    ! (en pas de physique)
    logical, save :: deja_lu_sur! pour indiquer que le jour a lire a deja
    ! lu pour une surface precedente
    integer, save :: jour_lu_sur
    integer :: ierr
    character (len = 20) :: modname = 'interfsur_lim'
    character (len = 80) :: abort_message
    logical, save :: newlmt = .false.
    logical, save :: check = .false.
    ! Champs lus dans le fichier de CL
    real, allocatable , save, dimension(:) :: alb_lu, rug_lu

    ! quelques variables pour netcdf

    include "netcdf.inc"
    integer , save :: nid, nvarid
    integer, dimension(2), save :: start, epais

    !------------------------------------------------------------

    if (debut) then
       lmt_pas = nint(86400./dtime * 1.0) ! pour une lecture une fois par jour
       jour_lu_sur = jour - 1
       allocate(alb_lu(klon))
       allocate(rug_lu(klon))
    endif

    if ((jour - jour_lu_sur) /= 0) deja_lu_sur = .false.

    if (check) write(*, *)modname, ':: jour_lu_sur, deja_lu_sur', jour_lu_sur, &
         deja_lu_sur
    if (check) write(*, *)modname, ':: itime, lmt_pas', itime, lmt_pas

    ! Tester d'abord si c'est le moment de lire le fichier
    if (mod(itime-1, lmt_pas) == 0 .and. .not. deja_lu_sur) then

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

       ! Lecture Albedo

       ierr = NF_INQ_VARID(nid, 'ALB', nvarid)
       if (ierr /= NF_NOERR) then
          abort_message = 'Le champ <ALB> est absent'
          call abort_gcm(modname, abort_message, 1)
       endif
       ierr = NF_GET_VARA_REAL(nid, nvarid, start, epais, alb_lu)
       if (ierr /= NF_NOERR) then
          abort_message = 'Lecture echouee pour <ALB>'
          call abort_gcm(modname, abort_message, 1)
       endif

       ! Lecture rugosité

       ierr = NF_INQ_VARID(nid, 'RUG', nvarid)
       if (ierr /= NF_NOERR) then
          abort_message = 'Le champ <RUG> est absent'
          call abort_gcm(modname, abort_message, 1)
       endif
       ierr = NF_GET_VARA_REAL(nid, nvarid, start, epais, rug_lu)
       if (ierr /= NF_NOERR) then
          abort_message = 'Lecture echouee pour <RUG>'
          call abort_gcm(modname, abort_message, 1)
       endif


       ! Fin de lecture

       ierr = NF_CLOSE(nid)
       deja_lu_sur = .true.
       jour_lu_sur = jour
    endif

    ! Recopie des variables dans les champs de sortie

!!$ lmt_alb = 0.0
!!$ lmt_rug = 0.0
    lmt_alb = 999999.
    lmt_rug = 999999.
    DO ii = 1, knon
       lmt_alb(ii) = alb_lu(knindex(ii))
       lmt_rug(ii) = rug_lu(knindex(ii))
    enddo

  END SUBROUTINE interfsur_lim

end module interfsur_lim_m
