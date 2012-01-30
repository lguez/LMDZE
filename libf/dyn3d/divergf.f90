module divergf_m

  IMPLICIT NONE

contains

  SUBROUTINE divergf(klevel, x, y, div)

    ! From libf/dyn3d/divergf.F, v 1.1.1.1 2004/05/19 12:53:05

    ! P. Le Van
    ! Calcule la divergence à tous les niveaux d'un vecteur de
    ! composantes x et y. x et y sont des composantes covariantes.

    USE dimens_m, ONLY: iim
    USE paramet_m, ONLY: iip1, iip2, ip1jm, ip1jmi1, ip1jmp1, jjp1
    USE comgeom, ONLY: apoln, apols, cuvsurcv, cvusurcu, unsaire
    USE filtreg_m, ONLY: filtreg

    INTEGER, intent(in):: klevel
    REAL, intent(in):: x(ip1jmp1, klevel), y(ip1jm, klevel)
    real, intent(out):: div(ip1jmp1, klevel) ! in (unit of x, y) m-2

    ! Variables locales :

    INTEGER l, ij
    REAL aiy1(iim) , aiy2(iim)
    REAL sumypn, sumyps

    !------------------------------------------------------------

    DO l = 1, klevel
       DO ij = iip2, ip1jm - 1
          div(ij + 1, l) = cvusurcu(ij+1) * x(ij+1, l) &
               - cvusurcu(ij) * x(ij , l) + cuvsurcv(ij-iim) * y(ij-iim, l) &
               - cuvsurcv(ij+1) * y(ij+1, l) 
       ENDDO

       DO ij = iip2, ip1jm, iip1
          div(ij, l) = div(ij + iim, l)
       ENDDO

       ! Calcul aux pôles 

       DO ij = 1, iim
          aiy1(ij) = cuvsurcv(ij) * y(ij , l)
          aiy2(ij) = cuvsurcv(ij+ ip1jmi1) * y(ij+ ip1jmi1, l)
       ENDDO
       sumypn = SUM(aiy1) / apoln
       sumyps = SUM(aiy2) / apols

       DO ij = 1, iip1
          div(ij , l) = - sumypn
          div(ij + ip1jm, l) = sumyps
       ENDDO
    end DO

    CALL filtreg(div, jjp1, klevel, 2, 2, .TRUE., 1)

    DO l = 1, klevel
       DO ij = iip2, ip1jm
          div(ij, l) = div(ij, l) * unsaire(ij) 
       ENDDO
    ENDDO

  END SUBROUTINE divergf

end module divergf_m
