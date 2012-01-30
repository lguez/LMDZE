SUBROUTINE dteta1(teta, pbaru, pbarv, dteta)

  ! From LMDZ4/libf/dyn3d/dteta1.F, version 1.1.1.1 2004/05/19 12:53:06
  ! Auteurs : P. Le Van, F. Forget

  ! Calcul du terme de convergence horizontale du flux d'enthalpie
  ! potentielle.

  ! dteta est un argument de sortie pour le s-pg

  use dimens_m
  use paramet_m
  use conf_gcm_m
  use filtreg_m, only: filtreg

  IMPLICIT NONE

  REAL, intent(in):: teta(ip1jmp1, llm), pbaru(ip1jmp1, llm), pbarv(ip1jm, llm)
  REAL dteta(ip1jmp1, llm)
  INTEGER l, ij

  REAL hbyv(ip1jm, llm), hbxu(ip1jmp1, llm)

  !----------------------------------------------------------------

  DO l = 1, llm
     DO ij = iip2, ip1jm - 1
        hbxu(ij, l) = pbaru(ij, l) * 0.5 * (teta(ij, l) + teta(ij + 1, l))
     end DO

     DO ij = iip1+ iip1, ip1jm, iip1
        hbxu(ij, l) = hbxu(ij - iim, l)
     end DO

     DO ij = 1, ip1jm
        hbyv(ij, l)= pbarv(ij, l) * 0.5 * (teta(ij, l) + teta(ij + iip1, l))
     end DO
  end DO

  CALL convflu(hbxu, hbyv, llm, dteta)

  ! stockage dans dh de la convergence horizont. filtree' du flux
  ! d'enthalpie potentielle
  CALL filtreg(dteta, jjp1, llm, 2, 2, .true., 1)

END SUBROUTINE dteta1
