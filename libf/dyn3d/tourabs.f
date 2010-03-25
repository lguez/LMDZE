      SUBROUTINE tourabs ( ntetaSTD,vcov, ucov, vorabs )

c=======================================================================
c
c   Modif:  I. Musat (28/10/04)
c   -------
c   adaptation du code tourpot.F pour le calcul de la vorticite absolue
c   cf. P. Le Van
c
c   Objet: 
c   ------
c
c    *******************************************************************
c    .............  calcul de la vorticite absolue     .................
c    *******************************************************************
c
c     ntetaSTD, vcov,ucov      sont des argum. d'entree pour le s-pg .
c             vorabs            est  un argum.de sortie pour le s-pg .
c
c=======================================================================

      use dimens_m
      use paramet_m
      use comconst
      use logic
      use comgeom
      use filtreg_m, only: filtreg

      IMPLICIT NONE
c
      INTEGER ntetaSTD
      REAL vcov( ip1jm,ntetaSTD ), ucov( ip1jmp1,ntetaSTD )
      REAL vorabs( ip1jm,ntetaSTD )
c
c variables locales
      INTEGER l, ij, i, j
      REAL  rot( ip1jm,ntetaSTD )



c  ... vorabs = ( Filtre( d(vcov)/dx - d(ucov)/dy ) + fext ) ..



c    ........  Calcul du rotationnel du vent V  puis filtrage  ........

      DO 5 l = 1,ntetaSTD

      DO 2 i = 1, iip1
      DO 2 j = 1, jjm
c
       ij=i+(j-1)*iip1
       IF(ij.LE.ip1jm - 1) THEN
c
        IF(cv(ij).EQ.0..OR.cv(ij+1).EQ.0..OR.
     $     cu(ij).EQ.0..OR.cu(ij+iip1).EQ.0.) THEN
         rot( ij,l ) = 0.
         continue
        ELSE
         rot( ij,l ) = (vcov(ij+1,l)/cv(ij+1)-vcov(ij,l)/cv(ij))/
     $                 (2.*pi*RAD*cos(rlatv(j)))*float(iim)
     $                +(ucov(ij+iip1,l)/cu(ij+iip1)-ucov(ij,l)/cu(ij))/
     $                 (pi*RAD)*(float(jjm)-1.)
c
        ENDIF
       ENDIF !(ij.LE.ip1jm - 1) THEN
   2  CONTINUE

c    ....  correction pour  rot( iip1,j,l )  .....
c    ....     rot(iip1,j,l) = rot(1,j,l)    .....

CDIR$ IVDEP

      DO 3 ij = iip1, ip1jm, iip1
      rot( ij,l ) = rot( ij -iim, l )
   3  CONTINUE

   5  CONTINUE


      CALL  filtreg( rot, jjm, ntetaSTD, 2, 1, .FALSE., 1 )


      DO 10 l = 1, ntetaSTD

      DO 6 ij = 1, ip1jm - 1
      vorabs( ij,l ) = ( rot(ij,l) + fext(ij)*unsairez(ij) )
   6  CONTINUE

c    ..... correction pour  vorabs( iip1,j,l)  .....
c    ....   vorabs(iip1,j,l)= vorabs(1,j,l) ....
CDIR$ IVDEP
      DO 8 ij = iip1, ip1jm, iip1
      vorabs( ij,l ) = vorabs( ij -iim,l )
   8  CONTINUE

  10  CONTINUE

      RETURN
      END
