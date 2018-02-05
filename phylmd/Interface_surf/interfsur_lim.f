module interfsur_lim_m

  implicit none

contains

  SUBROUTINE interfsur_lim(dtime, jour, knindex, debut, albedo, z0_new)

    ! Cette routine sert d'interface entre le mod\`ele atmosph\'erique et
    ! un fichier de conditions aux limites.

    ! Laurent FAIRHEAD, February 2000

    USE dimphy, ONLY: klon
    use netcdf, only: NF90_NOWRITE
    use netcdf95, only: NF95_close, NF95_GET_VAR, NF95_INQ_VARID, NF95_OPEN
    use time_phylmdz, only: itap

    real, intent(IN):: dtime ! pas de temps de la physique (en s)
    integer, intent(IN):: jour ! jour a lire dans l'annee

    integer, intent(in):: knindex(:) ! (knon)
    ! index des points de la surface \`a traiter

    logical, intent(IN):: debut ! premier appel \`a la physique (initialisation)
    real, intent(out):: albedo(:) ! (knon) albedo lu
    real, intent(out):: z0_new(:) ! (klon) longueur de rugosit\'e lue

    ! Local:

    integer knon ! nombre de points dans le domaine a traiter

    integer, save:: lmt_pas ! frequence de lecture des conditions limites
    ! (en pas de physique)

    logical, save:: deja_lu_sur
    ! jour \`a lire d\'ej\`a lu pour une surface pr\'ec\'edente

    integer, save:: jour_lu_sur

    ! Champs lus dans le fichier de conditions aux limites :
    real, allocatable, save:: alb_lu(:), rug_lu(:)

    integer ncid, varid

    !------------------------------------------------------------

    knon = size(knindex)

    if (debut) then
       lmt_pas = nint(86400. / dtime) ! pour une lecture une fois par jour
       jour_lu_sur = jour - 1
       allocate(alb_lu(klon))
       allocate(rug_lu(klon))
    endif

    if (jour - jour_lu_sur /= 0) deja_lu_sur = .false.

    ! Tester d'abord si c'est le moment de lire le fichier
    if (mod(itap - 1, lmt_pas) == 0 .and. .not. deja_lu_sur) then
       call NF95_OPEN('limit.nc', NF90_NOWRITE, ncid)

       ! Lecture Albedo
       call NF95_INQ_VARID(ncid, 'ALB', varid)
       call NF95_GET_VAR(ncid, varid, alb_lu, start=(/1, jour/))

       ! Lecture rugosit\'e
       call NF95_INQ_VARID(ncid, 'RUG', varid)
       call NF95_GET_VAR(ncid, varid, rug_lu, start=(/1, jour/))

       call NF95_CLOSE(ncid)
       deja_lu_sur = .true.
       jour_lu_sur = jour
    endif

    ! Recopie des variables dans les champs de sortie
    albedo = alb_lu(knindex)
    z0_new(:knon) = rug_lu(knindex)
    z0_new(knon + 1:) = 999999.

  END SUBROUTINE interfsur_lim

end module interfsur_lim_m
