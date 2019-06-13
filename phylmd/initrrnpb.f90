module initrrnpb_m

  IMPLICIT NONE

contains

  SUBROUTINE initrrnpb(pctsrf, masktr, fshtr, hsoltr, tautr, vdeptr, scavtr)

    ! From LMDZ4/libf/phylmd/initrrnpb.F, version 1.2, 2004/06/22 11:45:33

    ! Authors: AA + Christophe Genthon (LGGE/CNRS)
    ! Date: June 1994
    ! Objet : initialisation des constantes des traceurs

    USE dimphy, ONLY: klon

    REAL, intent(in):: pctsrf(:, :) ! (klon, nbsrf) 
    ! pourcentage de sol (f(nature du sol))
    ! Nature de sol (pourcentage de sol)

    REAL masktr(:, :) ! (klon, nqmx - 2)
    ! Masque de l'echange avec la surface
    ! (possible => 1)
    ! masktr---output-R- Masque reservoir de sol traceur (1 = reservoir)

    REAL fshtr(:, :) ! (klon, nqmx - 2)
    ! Flux surfacique dans le reservoir de sol
    ! fshtr----output-R- Flux surfacique de production dans le sol
    
    REAL hsoltr(:) ! (nqmx - 2)
    ! Epaisseur equivalente du reservoir de sol
    ! hsoltr---output-R- Epaisseur du reservoir de sol

    REAL tautr(:) ! (nqmx - 2)
    ! Constante de decroissance radioactive
    ! tautr----output-R- Constante de decroissance du traceur

    REAL vdeptr(:) ! (nqmx - 2)
    ! Vitesse de depot sec dans la couche Brownienne
    ! vdeptr---output-R- Vitesse de depot sec dans la couche Brownienne

    REAL scavtr(:) ! (nqmx - 2)
    ! Coefficient de lessivage
    ! scavtr---output-R- Coefficient de lessivage

    ! Local:
    INTEGER i
    REAL s

    !-------------------------------------------------------------

    PRINT *, 'Call sequence information: initrrnpb'

    ! Initialisations spécifiques à chaque traceur (pour le moment, Rn222)

    ! Radon

    s = 1.E4 ! Source: atome par m2

    hsoltr(1) = 0.1
    ! Hauteur equivalente du reservoir :
    ! 1 m * porosite 0.1

    tautr(1) = 4.765E5 ! Decroissance du radon, secondes
    vdeptr(1) = 0. ! Pas de depot sec pour le radon
    scavtr(1) = 0. ! Pas de lessivage pour le radon

    PRINT *, 'Source du radon :'
    PRINT *, 'Source : ', s
    PRINT *, 'Hauteur equivalente du reservoir de sol: ', hsoltr(1)
    PRINT *, 'Decroissance (s): ', tautr(1)
    PRINT *, 'Vitesse de depot sec: ', vdeptr(1)
    PRINT *, 'Facteur de lessivage: ', scavtr(1)

    DO i = 1, klon
       masktr(i, 1) = 0.
       IF (nint(pctsrf(i, 1))==1) masktr(i, 1) = 1.
       fshtr(i, 1) = s*masktr(i, 1)

       ! Pour l'instant le pctsrf(i, 3) = 1.
       ! lorsqu'il y a de la terre mais ne prend aucune autre valeur
       ! il n'est donc pas nécessaire de multiplier fshtr par pctsrf
    END DO

    ! Pb210

    s = 0. ! pas de source

    hsoltr(2) = 10.
    ! hauteur équivalente du réservoir à partir duquel le dépôt
    ! Brownien a lieu

    tautr(2) = 1.028E9 ! Decroissance du Pb210, secondes
    vdeptr(2) = 1.E-3 ! 1 mm/s pour le 210Pb
    scavtr(2) = .5 ! Lessivage du Pb210
    DO i = 1, klon
       masktr(i, 2) = 1. ! Le depot sec peut avoir lieu partout
       fshtr(i, 2) = s*masktr(i, 2)
    END DO
    PRINT *, 'Source du plomb :'
    PRINT *, 'Source : ', s
    PRINT *, 'Hauteur equivalente du reservoir : ', hsoltr(2)
    PRINT *, 'Decroissance (s): ', tautr(2)
    PRINT *, 'Vitesse de depot sec: ', vdeptr(2)
    PRINT *, 'Facteur de lessivage: ', scavtr(2)

  END SUBROUTINE initrrnpb

end module initrrnpb_m
