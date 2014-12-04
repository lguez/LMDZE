SUBROUTINE alboc_cd(rmu0, albedo)
  ! From LMDZ4/libf/phylmd/albedo.F,v 1.2 2005/02/07 15:00:52
  USE dimens_m
  USE dimphy
  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS)
  ! date: 19940624
  ! Calculer l'albedo sur l'ocean en fonction de l'angle zenithal moyen
  ! Formule due a Larson and Barkstrom (1977) Proc. of the symposium
  ! on radiation in the atmosphere, 19-28 August 1976, science Press,
  ! 1977 pp 451-453, ou These de 3eme cycle de Sylvie Joussaume.

  ! Arguments
  ! rmu0    (in): cosinus de l'angle solaire zenithal
  ! albedo (out): albedo de surface de l'ocean
  ! ======================================================================
  REAL rmu0(klon), albedo(klon)

  REAL fmagic ! un facteur magique pour regler l'albedo
  ! cc      PARAMETER (fmagic=0.7)
  ! ccIM => a remplacer
  ! PARAMETER (fmagic=1.32)
  PARAMETER (fmagic=1.0)
  ! PARAMETER (fmagic=0.7)

  REAL fauxo
  INTEGER i
  ! ccIM
  LOGICAL ancien_albedo
  PARAMETER (ancien_albedo=.FALSE.)
  ! SAVE albedo

  IF (ancien_albedo) THEN

    DO i = 1, klon

      rmu0(i) = max(rmu0(i), 0.0)

      fauxo = (1.47-acos(rmu0(i)))/0.15
      albedo(i) = fmagic*(.03+.630/(1.+fauxo*fauxo))
      albedo(i) = max(min(albedo(i),0.60), 0.04)
    END DO

    ! nouvel albedo

  ELSE

    DO i = 1, klon
      rmu0(i) = max(rmu0(i), 0.0)
      ! IM:orig albedo(i) = 0.058/(rmu0(i) + 0.30)
      albedo(i) = fmagic*0.058/(rmu0(i)+0.30)
      albedo(i) = max(min(albedo(i),0.60), 0.04)
    END DO

  END IF

  RETURN
END SUBROUTINE alboc_cd
