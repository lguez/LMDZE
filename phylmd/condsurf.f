! $Header: /home/cvsroot/LMDZ4/libf/phylmd/condsurf.F,v 1.2 2005/12/01
! 11:27:29 fairhead Exp $

SUBROUTINE condsurf(jour, jourvrai, lmt_bils)
  USE dimens_m
  USE indicesol
  USE dimphy
  USE temps
  USE clesphys2, ONLY: ok_limitvrai
  USE netcdf
  IMPLICIT NONE

  ! I. Musat 05.2005

  ! Lire chaque jour le bilan de chaleur au sol issu
  ! d'un run atmospherique afin de l'utiliser dans
  ! dans un run "slab" ocean
  ! -----------------------------------------
  ! jour     : input  , numero du jour a lire
  ! jourvrai : input  , vrai jour de la simulation

  ! lmt_bils: bilan chaleur au sol (a utiliser pour "slab-ocean")

  INTEGER nid, nvarid
  INTEGER debut(2)
  INTEGER epais(2)


  INTEGER nannemax
  PARAMETER (nannemax=60)

  INTEGER jour
  INTEGER jourvrai
  REAL lmt_bils(klon) !bilan chaleur au sol

  ! Variables locales:
  INTEGER ig, i, kt, ierr
  LOGICAL ok
  INTEGER anneelim, anneemax
  CHARACTER *20 fich
  ! c
  ! c   .....................................................................
  ! c
  ! c    Pour lire le fichier limit correspondant vraiment  a l'annee de la
  ! c     simulation en cours , il suffit de mettre  ok_limitvrai = .TRUE.
  ! c
  ! c
  ! ......................................................................


  IF (jour<0 .OR. jour>(360-1)) THEN
    PRINT *, 'Le jour demande n est pas correct: ', jour
    STOP 1
  END IF

  anneelim = annee_ref
  anneemax = annee_ref + nannemax


  IF (ok_limitvrai) THEN
    DO kt = 1, nannemax
      IF (jourvrai<=(kt-1)*360+359) THEN
        WRITE (fich, '("limit",i4,".nc")') anneelim
        ! PRINT *,' Fichier  Limite ',fich
        GO TO 100
      END IF
      anneelim = anneelim + 1
    END DO

    PRINT *, ' PBS ! Le jour a lire sur le fichier limit ne se '
    PRINT *, ' trouve pas sur les ', nannemax, ' annees a partir de '
    PRINT *, ' l annee de debut', annee_ref
    STOP 1

100 CONTINUE

  ELSE

    WRITE (fich, '("limitNEW.nc")')
    ! PRINT *,' Fichier  Limite ',fich
  END IF

  ! Ouvrir le fichier en format NetCDF:

  ierr = nf90_open(fich, nf90_nowrite, nid)
  IF (ierr/=nf90_noerr) THEN
    WRITE (6, *) ' Pb d''ouverture du fichier ', fich
    WRITE (6, *) ' Le fichier limit ', fich, ' (avec 4 chiffres , pour'
    WRITE (6, *) '       l an 2000 )  ,  n existe  pas !  '
    WRITE (6, *) ' ierr = ', ierr
    STOP 1
  END IF
  ! DO k = 1, jour
  ! La tranche de donnees a lire:

  debut(1) = 1
  debut(2) = jour
  epais(1) = klon
  epais(2) = 1

  ! Bilan flux de chaleur au sol:

  ierr = nf90_inq_varid(nid, 'BILS', nvarid)
  IF (ierr/=nf90_noerr) THEN
    PRINT *, 'condsurf: Le champ <BILS> est absent'
    STOP 1
  END IF
  ! PRINT*,'debut,epais',debut,epais
  ierr = nf90_get_var(nid, nvarid, lmt_bils, debut, epais)
  IF (ierr/=nf90_noerr) THEN
    PRINT *, 'condsurf: Lecture echouee pour <BILS>'
    STOP 1
  END IF
  ! ENDDO !k = 1, jour

  ! Fermer le fichier:

  ierr = nf90_close(nid)


  ! PRINT*, 'lmt_bils est lu pour jour: ', jour

  RETURN
END SUBROUTINE condsurf
