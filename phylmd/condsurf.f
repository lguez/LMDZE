c $Header: /home/cvsroot/LMDZ4/libf/phylmd/condsurf.F,v 1.2 2005/12/01 11:27:29 fairhead Exp $
c
      SUBROUTINE condsurf( jour, jourvrai, lmt_bils )
      use dimens_m
      use indicesol
      use dimphy
      use temps
      use clesphys2, only: ok_limitvrai
            use netcdf
      IMPLICIT none
c
c I. Musat 05.2005
c
c Lire chaque jour le bilan de chaleur au sol issu 
c d'un run atmospherique afin de l'utiliser dans
c dans un run "slab" ocean 
c -----------------------------------------
c jour     : input  , numero du jour a lire
c jourvrai : input  , vrai jour de la simulation  
c
c lmt_bils: bilan chaleur au sol (a utiliser pour "slab-ocean")
c
      INTEGER nid, nvarid
      INTEGER debut(2)
      INTEGER epais(2)
c
c
      INTEGER     nannemax
      PARAMETER ( nannemax = 60 )
c
      INTEGER jour
      INTEGER jourvrai
      REAL lmt_bils(klon) !bilan chaleur au sol
c
c Variables locales:
      INTEGER ig, i, kt, ierr
      LOGICAL ok
      INTEGER anneelim,anneemax
      CHARACTER*20 fich
cc
cc   .....................................................................
cc
cc    Pour lire le fichier limit correspondant vraiment  a l'annee de la
cc     simulation en cours , il suffit de mettre  ok_limitvrai = .TRUE.
cc
cc   ......................................................................
c
c
      IF (jour.LT.0 .OR. jour.GT.(360-1)) THEN
         PRINT*,'Le jour demande n est pas correct: ', jour
         STOP 1
      ENDIF
c
       anneelim  = annee_ref
       anneemax  = annee_ref + nannemax
c
c
       IF( ok_limitvrai )       THEN
          DO  kt = 1, nannemax
           IF(jourvrai.LE. (kt-1)*360 + 359  )  THEN
              WRITE(fich,'("limit",i4,".nc")') anneelim
c             PRINT *,' Fichier  Limite ',fich
              GO TO 100
             ENDIF
           anneelim = anneelim + 1
          ENDDO

         PRINT *,' PBS ! Le jour a lire sur le fichier limit ne se '
         PRINT *,' trouve pas sur les ',nannemax,' annees a partir de '
         PRINT *,' l annee de debut', annee_ref
         stop 1
c
100     CONTINUE
c
       ELSE
     
            WRITE(fich,'("limitNEW.nc")') 
c           PRINT *,' Fichier  Limite ',fich
       ENDIF
c
c Ouvrir le fichier en format NetCDF:
c
      ierr = NF90_OPEN (fich, NF90_NOWRITE,nid)
      IF (ierr.NE.NF90_NOERR) THEN
        WRITE(6,*)' Pb d''ouverture du fichier ', fich
        WRITE(6,*)' Le fichier limit ',fich,' (avec 4 chiffres , pour' 
        WRITE(6,*)'       l an 2000 )  ,  n existe  pas !  ' 
        WRITE(6,*)' ierr = ', ierr
        stop 1
      ENDIF
c     DO k = 1, jour
c La tranche de donnees a lire:
c
      debut(1) = 1
      debut(2) = jour
      epais(1) = klon
      epais(2) = 1
c
c Bilan flux de chaleur au sol:
c
      ierr = NF90_INQ_VARID (nid, "BILS", nvarid)
      IF (ierr .NE. NF90_NOERR) THEN
         PRINT*, "condsurf: Le champ <BILS> est absent"
         stop 1
      ENDIF
c     PRINT*,'debut,epais',debut,epais
      ierr = NF90_GET_VAR(nid, nvarid,lmt_bils,debut,epais)
      IF (ierr .NE. NF90_NOERR) THEN
         PRINT*, "condsurf: Lecture echouee pour <BILS>"
         stop 1
      ENDIF
c     ENDDO !k = 1, jour
c
c Fermer le fichier:
c
      ierr = NF90_CLOSE(nid)
c
c
c     PRINT*, 'lmt_bils est lu pour jour: ', jour
c
      RETURN
      END
