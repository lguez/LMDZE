
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/limx.F,v 1.1.1.1 2004/05/19
! 12:53:06 lmdzadmin Exp $

SUBROUTINE limx(s0, sx, sm, pente_max)

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
  REAL sx(ip1jmp1, llm)

  ! Local
  ! ---------

  INTEGER ij, l

  REAL q(ip1jmp1, llm)
  REAL dxq(ip1jmp1, llm)

  REAL dxqu(ip1jmp1)
  REAL adxqu(ip1jmp1), dxqmax(ip1jmp1)

  EXTERNAL convflu

  DO l = 1, llm
    DO ij = 1, ip1jmp1
      q(ij, l) = s0(ij, l)/sm(ij, l)
      dxq(ij, l) = sx(ij, l)/sm(ij, l)
    END DO
  END DO

  ! calcul de la pente a droite et a gauche de la maille

  DO l = 1, llm
    DO ij = iip2, ip1jm - 1
      dxqu(ij) = q(ij+1, l) - q(ij, l)
    END DO
    DO ij = iip1 + iip1, ip1jm, iip1
      dxqu(ij) = dxqu(ij-iim)
    END DO

    DO ij = iip2, ip1jm
      adxqu(ij) = abs(dxqu(ij))
    END DO

    ! calcul de la pente maximum dans la maille en valeur absolue

    DO ij = iip2 + 1, ip1jm
      dxqmax(ij) = pente_max*min(adxqu(ij-1), adxqu(ij))
    END DO

    DO ij = iip1 + iip1, ip1jm, iip1
      dxqmax(ij-iim) = dxqmax(ij)
    END DO

    ! calcul de la pente avec limitation

    DO ij = iip2 + 1, ip1jm
      IF (dxqu(ij-1)*dxqu(ij)>0. .AND. dxq(ij,l)*dxqu(ij)>0.) THEN
        dxq(ij, l) = sign(min(abs(dxq(ij,l)),dxqmax(ij)), dxq(ij,l))
      ELSE
        ! extremum local
        dxq(ij, l) = 0.
      END IF
    END DO
    DO ij = iip1 + iip1, ip1jm, iip1
      dxq(ij-iim, l) = dxq(ij, l)
    END DO

    DO ij = 1, ip1jmp1
      sx(ij, l) = dxq(ij, l)*sm(ij, l)
    END DO

  END DO

  RETURN
END SUBROUTINE limx
