
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/initrrnpb.F,v 1.2 2004/06/22
! 11:45:33 lmdzadmin Exp $

SUBROUTINE initrrnpb(ftsol, pctsrf, masktr, fshtr, hsoltr, tautr, vdeptr, &
    scavtr)
  USE dimens_m
  USE indicesol
  USE dimphy
  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): AA + CG (LGGE/CNRS) Date 24-06-94
  ! Objet: initialisation des constantes des traceurs
  ! AA Revison pour le controle avec la temperature du sol
  ! AA
  ! AA   it = 1 radon ss controle de ts
  ! AA   it = 2 plomb ss controle de ts
  ! ======================================================================
  ! Arguments:
  ! nbtr------input-I- nombre de vrais traceurs (sans l'eau)
  ! ftsol-------input-R- Temperature du sol (Kelvin)
  ! pctsrf-----input-R-  Nature de sol (pourcentage de sol)
  ! masktr---output-R- Masque reservoir de sol traceur (1 = reservoir)
  ! fshtr----output-R- Flux surfacique de production dans le sol
  ! hsoltr---output-R- Epaisseur du reservoir de sol
  ! tautr----output-R- Constante de decroissance du traceur
  ! vdeptr---output-R- Vitesse de depot sec dans la couche Brownienne
  ! scavtr---output-R- Coefficient de lessivage
  ! ======================================================================
  ! ======================================================================

  INTEGER i, it
  REAL pctsrf(klon, nbsrf) !Pourcentage de sol (f(nature du sol))
  REAL ftsol(klon, nbsrf) ! Temperature du sol pour le controle Rn
  ! ! le cas echeant
  REAL masktr(klon, nbtr) ! Masque de l'echange avec la surface
  ! (possible => 1 )
  REAL fshtr(klon, nbtr) ! Flux surfacique dans le reservoir de sol
  REAL hsoltr(nbtr) ! Epaisseur equivalente du reservoir de sol
  REAL tautr(nbtr) ! Constante de decroissance radioactive
  REAL vdeptr(nbtr) ! Vitesse de depot sec dans la couche Brownienne
  REAL scavtr(nbtr) ! Coefficient de lessivage
  REAL s

  PRINT *, 'Call sequence information: initrrnpb'
  PRINT *, 'nbtr= ', nbtr
  PRINT *, 'nbsrf= ', nbsrf
  PRINT *, 'klon= ', klon

  ! Puis les initialisation specifiques a chaque traceur (pour le moment,
  ! Rn222)


  ! Radon it = 1

  IF (nbtr<=0) STOP 'initrrnpb pas glop pas glop'
  it = 1
  s = 1.E4 !  Source: atome par m2
  hsoltr(it) = 0.1 ! Hauteur equivalente du reservoir :
  ! 1 m * porosite 0.1
  tautr(it) = 4.765E5 ! Decroissance du radon, secondes
  ! AA
  ! tautr(it) = 4.765E55  ! Decroissance du radon,infinie
  ! AA
  vdeptr(it) = 0. ! Pas de depot sec pour le radon
  scavtr(it) = 0. ! Pas de lessivage pour le radon

  PRINT *, '-------------- SOURCE DU RADON ------------------------ '
  PRINT *, 'it = ', it
  PRINT *, 'Source : ', s
  PRINT *, 'Hauteur equivalente du reservoir de sol: ', hsoltr(it)
  PRINT *, 'Decroissance (s): ', tautr(it)
  PRINT *, 'Vitesse de depot sec: ', vdeptr(it)
  PRINT *, 'Facteur de lessivage: ', scavtr(it)

  DO i = 1, klon
    masktr(i, it) = 0.
    IF (nint(pctsrf(i,1))==1) masktr(i, it) = 1.
    fshtr(i, it) = s*masktr(i, it)

    ! AA
    ! AA POur l'instant le pctsrf(i,3) = 1.0
    ! AA lorsqu'il ya de la terre mias ne prend aucune autre valeur
    ! AA il n'est donc pas necessaire de multiplier fshtr par pctsrf
    ! AA

  END DO

  ! 210Pb it = 2

  IF (nbtr<=1) STOP 'initrrnpb pas glop pas glop'
  it = 2
  s = 0. !  Pas de source !!!
  hsoltr(it) = 10. ! Hauteur equivalente du reservoir
  ! a partir duquel le
  ! depot Brownien a lieu
  tautr(it) = 1.028E9 ! Decroissance du Pb210, secondes
  vdeptr(it) = 1.E-3 ! 1 mm/s pour le 210Pb
  scavtr(it) = .5 ! Lessivage du Pb210
  DO i = 1, klon
    masktr(i, it) = 1. ! Le depot sec peut avoir lieu partout
    fshtr(i, it) = s*masktr(i, it)
  END DO
  PRINT *, '-------------- SOURCE DU PLOMB ------------------------ '
  PRINT *, 'it = ', it
  PRINT *, 'Source : ', s
  PRINT *, 'Hauteur equivalente du reservoir : ', hsoltr(it)
  PRINT *, 'Decroissance (s): ', tautr(it)
  PRINT *, 'Vitesse de depot sec: ', vdeptr(it)
  PRINT *, 'Facteur de lessivage: ', scavtr(it)

  WRITE (*, *) 'initialisation rnpb ok'

  RETURN
END SUBROUTINE initrrnpb
