
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/laplacien_rotgam.F,v 1.1.1.1
! 2004/05/19 12:53:05 lmdzadmin Exp $

SUBROUTINE laplacien_rotgam(klevel, rotin, rotout)

  ! P. Le Van

  ! ************************************************************
  ! ... calcul de  (rotat x nxgrad)_gam  du rotationnel rotin ..
  ! ************************************************************
  ! klevel et teta  sont des arguments  d'entree pour le s-prog
  ! divgra     est  un argument  de sortie pour le s-prog

  USE dimensions
  USE paramet_m
  USE comgeom
  IMPLICIT NONE



  ! .............   variables  en  arguments    ...........

  INTEGER, INTENT (IN) :: klevel
  REAL rotin(ip1jm, klevel), rotout(ip1jm, klevel)

  ! ............     variables   locales     ...............

  INTEGER l, ij
  REAL ghy(ip1jm, llm), ghx(ip1jmp1, llm)
  ! ........................................................



  CALL nxgrad_gam(klevel, rotin, ghx, ghy)
  CALL rotat_nfil(klevel, ghx, ghy, rotout)

  DO l = 1, klevel
    DO ij = 1, ip1jm
      rotout(ij, l) = rotout(ij, l)*unsairz_gam(ij)
    END DO
  END DO

  RETURN
END SUBROUTINE laplacien_rotgam
