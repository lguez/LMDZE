module gradiv2_m

  IMPLICIT NONE

contains

  SUBROUTINE gradiv2(klevel, xcov, ycov, ld, gdx, gdy, cdivu)

    ! From LMDZ4/libf/dyn3d/gradiv2.F, version 1.1.1.1 2004/05/19 12:53:07

    ! P. Le Van
    ! calcul de grad div du vecteur v
    ! xcov et ycov etant les composantes covariantes de v
    ! xcont, ycont et ld sont des arguments d'entree pour le sous-programme
    ! gdx et gdy sont des arguments de sortie pour le sous-programme

    use dimens_m
    use paramet_m
    use comgeom
    use filtreg_m, only: filtreg

    ! variables en arguments 

    INTEGER klevel
    REAL xcov( ip1jmp1,klevel), ycov( ip1jm,klevel)
    integer, intent(in):: ld
    REAL gdx( ip1jmp1,klevel), gdy( ip1jm,klevel)
    real, intent(in):: cdivu

    ! variables locales

    REAL div(ip1jmp1,llm)
    REAL nugrads
    INTEGER l,ij,iter

    !--------------------------------------------------------------

    CALL SCOPY( ip1jmp1 * klevel, xcov, 1, gdx, 1)
    CALL SCOPY( ip1jm * klevel, ycov, 1, gdy, 1)

    CALL divergf( klevel, gdx, gdy, div)

    IF( ld.GT.1) THEN
       CALL laplacien ( klevel, div, div)

       ! Iteration de l'operateur laplacien_gam
       DO iter = 1, ld -2
          CALL laplacien_gam ( klevel,cuvscvgam1,cvuscugam1,unsair_gam1, &
               unsapolnga1, unsapolsga1, div, div)
       ENDDO
    ENDIF

    CALL filtreg( div, jjp1, klevel, 2, 1, .TRUE., 1)
    CALL grad ( klevel, div, gdx, gdy)
    nugrads = (-1.)**ld * cdivu

    DO l = 1, klevel
       DO ij = 1, ip1jmp1
          gdx( ij,l) = gdx( ij,l) * nugrads
       ENDDO
       DO ij = 1, ip1jm
          gdy( ij,l) = gdy( ij,l) * nugrads
       ENDDO
    ENDDO

  END SUBROUTINE gradiv2

end module gradiv2_m
