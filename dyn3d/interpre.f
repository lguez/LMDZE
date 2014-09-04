
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/interpre.F,v 1.1.1.1 2004/05/19
! 12:53:07 lmdzadmin Exp $

SUBROUTINE interpre(q, qppm, w, fluxwppm, masse, apppm, bpppm, massebx, &
    masseby, pbaru, pbarv, unatppm, vnatppm, psppm)

  USE dimens_m
  USE paramet_m
  USE comconst
  USE disvert_m
  USE conf_gcm_m
  USE conf_gcm_m
  USE comgeom
  USE temps
  IMPLICIT NONE

  ! ---------------------------------------------------
  ! Arguments
  REAL apppm(llm+1), bpppm(llm+1)
  REAL q(iip1, jjp1, llm), qppm(iim, jjp1, llm)
  ! ---------------------------------------------------
  REAL masse(iip1, jjp1, llm)
  REAL massebx(iip1, jjp1, llm), masseby(iip1, jjm, llm)
  REAL w(iip1, jjp1, llm+1)
  REAL fluxwppm(iim, jjp1, llm)
  REAL, INTENT (IN) :: pbaru(iip1, jjp1, llm)
  REAL, INTENT (IN) :: pbarv(iip1, jjm, llm)
  REAL unatppm(iim, jjp1, llm)
  REAL vnatppm(iim, jjp1, llm)
  REAL psppm(iim, jjp1)
  ! ---------------------------------------------------
  ! Local
  REAL vnat(iip1, jjp1, llm)
  REAL unat(iip1, jjp1, llm)
  REAL fluxw(iip1, jjp1, llm)
  REAL smass(iip1, jjp1)
  ! ----------------------------------------------------
  INTEGER l, i, j

  ! CALCUL DE LA PRESSION DE SURFACE
  ! Les coefficients ap et bp sont passés en common
  ! Calcul de la pression au sol en mb optimisée pour
  ! la vectorialisation

  DO j = 1, jjp1
    DO i = 1, iip1
      smass(i, j) = 0.
    END DO
  END DO

  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iip1
        smass(i, j) = smass(i, j) + masse(i, j, l)
      END DO
    END DO
  END DO

  DO j = 1, jjp1
    DO i = 1, iim
      psppm(i, j) = smass(i, j)/aire_2d(i, j)*g*0.01
    END DO
  END DO

  ! RECONSTRUCTION DES CHAMPS CONTRAVARIANTS
  ! Le programme ppm3d travaille avec les composantes
  ! de vitesse et pas les flux, on doit donc passer de l'un à l'autre
  ! Dans le même temps, on fait le changement d'orientation du vent en v
  DO l = 1, llm
    DO j = 1, jjm
      DO i = 1, iip1
        vnat(i, j, l) = -pbarv(i, j, l)/masseby(i, j, l)*cv_2d(i, j)
      END DO
    END DO
    DO i = 1, iim
      vnat(i, jjp1, l) = 0.
    END DO
    DO j = 1, jjp1
      DO i = 1, iip1
        unat(i, j, l) = pbaru(i, j, l)/massebx(i, j, l)*cu_2d(i, j)
      END DO
    END DO
  END DO

  ! CALCUL DU FLUX MASSIQUE VERTICAL
  ! Flux en l=1 (sol) nul
  fluxw = 0.
  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iip1
        fluxw(i, j, l) = w(i, j, l)*g*0.01/aire_2d(i, j)
      END DO
    END DO
  END DO

  ! INVERSION DES NIVEAUX
  ! le programme ppm3d travaille avec une 3ème coordonnée inversée par
  ! rapport
  ! de celle du LMDZ: z=1<=>niveau max, z=llm+1<=>surface
  ! On passe donc des niveaux du LMDZ à ceux de Lin

  DO l = 1, llm + 1
    apppm(l) = ap(llm+2-l)
    bpppm(l) = bp(llm+2-l)
  END DO

  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iim
        unatppm(i, j, l) = unat(i, j, llm-l+1)
        vnatppm(i, j, l) = vnat(i, j, llm-l+1)
        fluxwppm(i, j, l) = fluxw(i, j, llm-l+1)
        qppm(i, j, l) = q(i, j, llm-l+1)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE interpre






