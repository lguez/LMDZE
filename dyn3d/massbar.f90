SUBROUTINE massbar(  masse, massebx, masseby )

  ! From LMDZ4/libf/dyn3d/massbar.F,v 1.1.1.1 2004/05/19 12:53:05

  !
  ! **********************************************************************
  !
  !  Calcule les moyennes en x et  y de la masse d'air dans chaque maille.
  ! **********************************************************************
  !    Auteurs : P. Le Van , Fr. Hourdin  .
  !   ..........
  !
  !  ..  masse                 est  un argum. d'entree  pour le s-pg ...
  !  ..  massebx,masseby      sont des argum. de sortie pour le s-pg ...
  !     
  !
  !     IMPLICIT NONE
  !
  use dimens_m
  use paramet_m
  use comconst
  use comgeom
  !
  REAL    masse( ip1jmp1,llm ), massebx( ip1jmp1,llm )  , &
       masseby(   ip1jm,llm )
  !
  !
  ! Méthode pour calculer massebx et masseby. A chaque point scalaire
  ! P(i, j) sont affectés quatre coefficients d'aire.

  ! alpha1(i, j) calculé au point (i + 1/4, j - 1/4)
  ! alpha2(i, j) calculé au point (i + 1/4, j + 1/4)
  ! alpha3(i, j) calculé au point (i - 1/4, j + 1/4)
  ! alpha4(i, j) calculé au point (i - 1/4, j - 1/4)

  ! Avec alpha1(i, j) = aire(i + 1/4, j - 1/4)/ aire(i, j) 

  ! Pour plus de détails, voir sous-programme "iniconst" et
  ! "massbar.txt".


  DO       l = 1 , llm
     !
     DO  ij = 1, ip1jmp1 - 1
        massebx(ij,l) =  masse( ij, l) * alpha1p2( ij  )     +  &
             masse(ij+1, l) * alpha3p4(ij+1 )
     ENDDO

     !    .... correction pour massebx( iip1,j) .....
     !    ...    massebx(iip1,j)= massebx(1,j) ...
     !
     !DIR$ IVDEP
     DO  ij = iip1, ip1jmp1, iip1
        massebx( ij,l ) = massebx( ij - iim,l )
     ENDDO


     DO  ij = 1,ip1jm
        masseby( ij,l ) = masse(  ij   , l ) * alpha2p3(   ij    )  + &
             masse(ij+iip1, l ) * alpha1p4( ij+iip1 )
     ENDDO

  end DO

END SUBROUTINE massbar
