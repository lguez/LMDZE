module covnat_m

  IMPLICIT NONE

contains

  SUBROUTINE covnat(klevel, ucov, vcov, unat, vnat)

    ! From LMDZ4/libf/dyn3d/covnat.F, version 1.1.1.1 2004/05/19 12:53:07

    USE comgeom, ONLY: cu, cv
    USE paramet_m, ONLY: iip1, iip2, ip1jm, ip1jmp1

    ! Authors: F. Hourdin, Phu Le Van

    ! Objet : calcul des composantes naturelles du vent  à partir des
    ! composantes covariantes.

    INTEGER, intent(in):: klevel
    REAL, intent(in):: ucov(ip1jmp1, klevel), vcov(ip1jm, klevel)
    REAL, intent(out):: unat(ip1jmp1, klevel), vnat(ip1jm, klevel)

    ! Local:
    INTEGER l, ij

    !------------------------------------------------------------------

    DO l = 1, klevel
       DO ij = 1, iip1
          unat(ij, l) =0.
       END DO

       DO ij = iip2, ip1jm
          unat(ij, l) = ucov(ij, l) / cu(ij)
       ENDDO

       DO ij = ip1jm+1, ip1jmp1 
          unat(ij, l) =0.
       END DO

       DO ij = 1, ip1jm
          vnat(ij, l) = vcov(ij, l) / cv(ij)
       ENDDO
    ENDDO

  END SUBROUTINE covnat

end module covnat_m
