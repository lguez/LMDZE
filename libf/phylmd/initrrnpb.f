!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/initrrnpb.F,v 1.2 2004/06/22 11:45:33 lmdzadmin Exp $
!
      SUBROUTINE  initrrnpb(ftsol,pctsrf,masktr,fshtr,hsoltr,tautr
     .                   ,vdeptr,scavtr)
      use dimens_m
      use indicesol
      use dimphy
      IMPLICIT none
c======================================================================
c Auteur(s): AA + CG (LGGE/CNRS) Date 24-06-94
c Objet: initialisation des constantes des traceurs
CAA Revison pour le controle avec la temperature du sol
cAA
CAA   it = 1 radon ss controle de ts
cAA   it = 2 plomb ss controle de ts  
c======================================================================
c Arguments:
c nbtr------input-I- nombre de vrais traceurs (sans l'eau)
c ftsol-------input-R- Temperature du sol (Kelvin)
c pctsrf-----input-R-  Nature de sol (pourcentage de sol)
c masktr---output-R- Masque reservoir de sol traceur (1 = reservoir)
c fshtr----output-R- Flux surfacique de production dans le sol
c hsoltr---output-R- Epaisseur du reservoir de sol
c tautr----output-R- Constante de decroissance du traceur
c vdeptr---output-R- Vitesse de depot sec dans la couche Brownienne
c scavtr---output-R- Coefficient de lessivage
c======================================================================
c======================================================================
C
      INTEGER i, it
      REAL pctsrf(klon,nbsrf) !Pourcentage de sol (f(nature du sol))
      REAL ftsol(klon,nbsrf)  ! Temperature du sol pour le controle Rn
c                             ! le cas echeant
      REAL masktr(klon,nbtr)  ! Masque de l'echange avec la surface
c                                 (possible => 1 )
      REAL fshtr(klon,nbtr)  ! Flux surfacique dans le reservoir de sol
      REAL hsoltr(nbtr)      ! Epaisseur equivalente du reservoir de sol
      REAL tautr(nbtr)       ! Constante de decroissance radioactive
      REAL vdeptr(nbtr)      ! Vitesse de depot sec dans la couche Brownienne
      REAL scavtr(nbtr)      ! Coefficient de lessivage
      REAL s
C
      print *, "Call sequence information: initrrnpb"
      print*,'nbtr= ',nbtr
      print*,'nbsrf= ',nbsrf
      print*,'klon= ',klon
C
C Puis les initialisation specifiques a chaque traceur (pour le moment, Rn222)
C
C
C Radon it = 1
c
      IF ( nbtr .LE. 0 ) STOP 'initrrnpb pas glop pas glop' 
      it = 1
      s = 1.E4  !  Source: atome par m2
      hsoltr(it) = 0.1      ! Hauteur equivalente du reservoir : 
c                              1 m * porosite 0.1
      tautr(it) = 4.765E5  ! Decroissance du radon, secondes
cAA
c      tautr(it) = 4.765E55  ! Decroissance du radon,infinie
cAA
      vdeptr(it) = 0. ! Pas de depot sec pour le radon
      scavtr(it) = 0. ! Pas de lessivage pour le radon

      print*, '-------------- SOURCE DU RADON ------------------------ '
      print*,'it = ',it
      print*,'Source : ', s
      print*,'Hauteur equivalente du reservoir de sol: ',hsoltr(it) 
      print*,'Decroissance (s): ', tautr(it)
      print*,'Vitesse de depot sec: ',vdeptr(it) 
      print*,'Facteur de lessivage: ',scavtr(it)

      DO i = 1,klon
        masktr(i,it) = 0.
        IF ( NINT(pctsrf(i,1)) .EQ. 1 ) masktr(i,it) = 1.
        fshtr(i,it) = s * masktr(i,it)

cAA
cAA POur l'instant le pctsrf(i,3) = 1.0 
cAA lorsqu'il ya de la terre mias ne prend aucune autre valeur
cAA il n'est donc pas necessaire de multiplier fshtr par pctsrf
cAA 

      END DO
C
C 210Pb it = 2
C
      IF ( nbtr .LE. 1 ) STOP 'initrrnpb pas glop pas glop' 
      it = 2
      s = 0. !  Pas de source !!!
      hsoltr(it) = 10.     ! Hauteur equivalente du reservoir 
c                              a partir duquel le
c                              depot Brownien a lieu
      tautr(it) = 1.028E9 ! Decroissance du Pb210, secondes
      vdeptr(it) = 1.E-3 ! 1 mm/s pour le 210Pb
      scavtr(it) =  .5   ! Lessivage du Pb210
      DO i = 1,klon
        masktr(i,it) = 1. ! Le depot sec peut avoir lieu partout
        fshtr(i,it) = s * masktr(i,it)
      END DO
      print*, '-------------- SOURCE DU PLOMB ------------------------ '
      print*,'it = ',it
      print*,'Source : ', s
      print*,'Hauteur equivalente du reservoir : ',hsoltr(it) 
      print*,'Decroissance (s): ', tautr(it)
      print*,'Vitesse de depot sec: ',vdeptr(it) 
      print*,'Facteur de lessivage: ',scavtr(it)
c
      WRITE(*,*) 'initialisation rnpb ok'
c
      RETURN
      END
