!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/massbar.F,v 1.1.1.1 2004/05/19 12:53:05 lmdzadmin Exp $
!
      SUBROUTINE massbar(  masse, massebx, masseby )
c
c **********************************************************************
c
c  Calcule les moyennes en x et  y de la masse d'air dans chaque maille.
c **********************************************************************
c    Auteurs : P. Le Van , Fr. Hourdin  .
c   ..........
c
c  ..  masse                 est  un argum. d'entree  pour le s-pg ...
c  ..  massebx,masseby      sont des argum. de sortie pour le s-pg ...
c     
c
c     IMPLICIT NONE
c
      use dimens_m
      use paramet_m
      use comconst
      use comgeom
c
      REAL    masse( ip1jmp1,llm ), massebx( ip1jmp1,llm )  ,
     *      masseby(   ip1jm,llm )
c
c
c   Methode pour calculer massebx et masseby .
c   ----------------------------------------
c
c    A chaque point scalaire P (i,j) est affecte 4 coefficients d'aires
c       alpha1(i,j)  calcule  au point ( i+1/4,j-1/4 )
c       alpha2(i,j)  calcule  au point ( i+1/4,j+1/4 )
c       alpha3(i,j)  calcule  au point ( i-1/4,j+1/4 )
c       alpha4(i,j)  calcule  au point ( i-1/4,j-1/4 )
c
c    Avec  alpha1(i,j) = aire(i+1/4,j-1/4)/ aire(i,j)        
c
c    N.B .  Pour plus de details, voir s-pg  ...  iniconst ...
c
c
c
c   alpha4 .         . alpha1    . alpha4
c    (i,j)             (i,j)       (i+1,j)
c
c             P .        U .          . P
c           (i,j)       (i,j)         (i+1,j)
c
c   alpha3 .         . alpha2    .alpha3 
c    (i,j)              (i,j)     (i+1,j)
c
c             V .        Z .          . V
c           (i,j)
c
c   alpha4 .         . alpha1    .alpha4
c   (i,j+1)            (i,j+1)   (i+1,j+1) 
c
c             P .        U .          . P
c          (i,j+1)                    (i+1,j+1)
c
c
c
c                       On  a :
c
c    massebx(i,j) = masse(i  ,j) * ( alpha1(i  ,j) + alpha2(i,j))   +
c                   masse(i+1,j) * ( alpha3(i+1,j) + alpha4(i+1,j) )
c     localise  au point  ... U (i,j) ...
c
c    masseby(i,j) = masse(i,j  ) * ( alpha2(i,j  ) + alpha3(i,j  )  +
c                   masse(i,j+1) * ( alpha1(i,j+1) + alpha4(i,j+1)  
c     localise  au point  ... V (i,j) ...
c
c
c=======================================================================

      DO   100    l = 1 , llm
c
        DO  ij = 1, ip1jmp1 - 1
         massebx(ij,l) =  masse( ij, l) * alpha1p2( ij  )     + 
     *                   masse(ij+1, l) * alpha3p4(ij+1 )
        ENDDO

c    .... correction pour massebx( iip1,j) .....
c    ...    massebx(iip1,j)= massebx(1,j) ...
c
CDIR$ IVDEP
        DO  ij = iip1, ip1jmp1, iip1
         massebx( ij,l ) = massebx( ij - iim,l )
        ENDDO


         DO  ij = 1,ip1jm
         masseby( ij,l ) = masse(  ij   , l ) * alpha2p3(   ij    )  +
     *                     masse(ij+iip1, l ) * alpha1p4( ij+iip1 )
         ENDDO

100   CONTINUE
c
      RETURN
      END
