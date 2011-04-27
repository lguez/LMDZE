SUBROUTINE covnat(klevel,ucov, vcov, unat, vnat )

  ! From LMDZ4/libf/dyn3d/covnat.F,v 1.1.1.1 2004/05/19 12:53:07

  use dimens_m
  use paramet_m
  use comgeom

  IMPLICIT NONE

  !   Auteur:  F Hourdin Phu LeVan
  !   Objet:
  !    calcul des compos. naturelles a partir des comp.covariantes

  INTEGER klevel
  REAL ucov( ip1jmp1,klevel ),  vcov( ip1jm,klevel )
  REAL unat( ip1jmp1,klevel ), vnat( ip1jm,klevel )
  INTEGER   l,ij

  !------------------------------------------------------------------

  DO l = 1,klevel
     DO ij = 1, iip1
        unat (ij,l) =0.
     END DO

     DO ij = iip2, ip1jm
        unat( ij,l ) = ucov( ij,l ) / cu(ij)
     ENDDO

     DO ij = ip1jm+1, ip1jmp1  
        unat (ij,l) =0.
     END DO

     DO ij = 1,ip1jm
        vnat( ij,l ) = vcov( ij,l ) / cv(ij)
     ENDDO
  ENDDO

END SUBROUTINE covnat
