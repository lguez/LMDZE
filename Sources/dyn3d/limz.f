
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/limz.F,v 1.1.1.1 2004/05/19
! 12:53:07 lmdzadmin Exp $

SUBROUTINE limz(s0, sz, sm, pente_max)

  ! Auteurs:   P.Le Van, F.Hourdin, F.Forget

  ! ********************************************************************
  ! Shema  d'advection " pseudo amont " .
  ! ********************************************************************
  ! nq,iq,q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....


  ! --------------------------------------------------------------------
  USE dimens_m
  USE paramet_m
  USE comconst
  USE disvert_m
  USE conf_gcm_m
  USE comgeom
  IMPLICIT NONE



  ! Arguments:
  ! ----------
  REAL pente_max
  REAL s0(ip1jmp1, llm), sm(ip1jmp1, llm)
  REAL sz(ip1jmp1, llm)

  ! Local
  ! ---------

  INTEGER ij, l

  REAL q(ip1jmp1, llm)
  REAL dzq(ip1jmp1, llm)

  REAL dzqw(ip1jmp1)
  REAL adzqw(ip1jmp1), dzqmax(ip1jmp1)

  LOGICAL first
  SAVE first

  REAL ssum
  INTEGER ismax, ismin
  EXTERNAL ssum, convflu, ismin, ismax

  DATA first/.TRUE./


  DO l = 1, llm
    DO ij = 1, ip1jmp1
      q(ij, l) = s0(ij, l)/sm(ij, l)
      dzq(ij, l) = sz(ij, l)/sm(ij, l)
    END DO
  END DO

  ! calcul de la pente en haut et en bas de la maille
  DO ij = 1, ip1jmp1
    DO l = 1, llm - 1
      dzqw(l) = q(ij, l+1) - q(ij, l)
    END DO
    dzqw(llm) = 0.

    DO l = 1, llm
      adzqw(l) = abs(dzqw(l))
    END DO

    ! calcul de la pente maximum dans la maille en valeur absolue

    DO l = 2, llm - 1
      dzqmax(l) = pente_max*min(adzqw(l-1), adzqw(l))
    END DO

    ! calcul de la pente avec limitation

    DO l = 2, llm - 1
      IF (dzqw(l-1)*dzqw(l)>0. .AND. dzq(ij,l)*dzqw(l)>0.) THEN
        dzq(ij, l) = sign(min(abs(dzq(ij,l)),dzqmax(l)), dzq(ij,l))
      ELSE
        ! extremum local
        dzq(ij, l) = 0.
      END IF
    END DO

    DO l = 1, llm
      sz(ij, l) = dzq(ij, l)*sm(ij, l)
    END DO

  END DO

  RETURN
END SUBROUTINE limz
