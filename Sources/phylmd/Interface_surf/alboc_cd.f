module alboc_cd_m

  IMPLICIT NONE

contains

  SUBROUTINE alboc_cd(rmu0, albedo)

    ! From LMDZ4/libf/phylmd/albedo.F, version 1.2 2005/02/07 15:00:52

    ! Author: Z. X. Li (LMD/CNRS)
    ! Date: 1994/06/24

    ! Calculer l'alb\'edo sur l'oc\'ean en fonction de l'angle
    ! z\'enithal moyen. Formule due \`a Larson and Barkstrom,
    ! Proceedings of the symposium on radiation in the atmosphere,
    ! 19-28 August 1976, science Press, 1977, pages 451-453, ou
    ! th\`ese de 3\`eme cycle de Sylvie Joussaume.

    REAL, intent(in):: rmu0(:) ! cosinus de l'angle solaire z\'enithal
    real, intent(out):: albedo(:) ! alb\'edo de surface de l'oc\'ean

    ! Local:
    REAL, PARAMETER:: fmagic = 1. ! facteur magique pour r\'egler l'alb\'edo

    !----------------------------------------------------------

    albedo = max(min(fmagic * 0.058 / (max(rmu0, 0.) + 0.3), 0.6), 0.04)

  END SUBROUTINE alboc_cd

end module alboc_cd_m
