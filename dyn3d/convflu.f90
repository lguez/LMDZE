SUBROUTINE convflu( xflu,yflu,nbniv,convfl )

  ! From LMDZ4/libf/dyn3d/convflu.F, version 1.1.1.1 2004/05/19 12:53:05

  !  P. Le Van

  ! Calcule la (convergence horiz. * aire locale) du flux ayant pour
  ! composantes xflu et yflu ,variables extensives .

  ! nbniv   est le nombre de niveaux vert. de  xflu et de yflu.

  use dimensions
  use paramet_m
  use comgeom

  IMPLICIT NONE

  integer, intent(in):: nbniv
  REAL, intent(in):: xflu( ip1jmp1,nbniv ), yflu( ip1jm,nbniv )
  real, intent(out):: convfl( ip1jmp1,nbniv )

  real convpn,convps
  INTEGER    l,ij
  REAL       SSUM

  !------------------------------------------------------------------

  DO  l = 1,nbniv
     DO  ij = iip2, ip1jm - 1
        convfl( ij + 1,l ) =  xflu(   ij,l ) - xflu( ij +  1,l )   + &
             yflu(ij +1,l ) - yflu( ij -iim,l )
     end DO

     ! correction pour  convfl( 1,j,l)
     ! convfl(1,j,l)= convfl(iip1,j,l)

     DO  ij = iip2,ip1jm,iip1
        convfl( ij,l ) = convfl( ij + iim,l )
     end DO

     ! calcul aux pôles

     convpn =   SSUM( iim, yflu(     1    ,l ),  1 )
     convps = - SSUM( iim, yflu( ip1jm-iim,l ),  1 )
     DO  ij = 1,iip1
        convfl(     ij   ,l ) = convpn * aire(   ij     ) / apoln
        convfl( ij+ ip1jm,l ) = convps * aire( ij+ ip1jm) / apols
     end DO
  end DO

END SUBROUTINE convflu
