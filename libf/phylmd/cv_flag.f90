
!==================================================================
      SUBROUTINE cv_flag
            use cvflag
      implicit none


! -- si .TRUE., on rend la gravite plus explicite et eventuellement
! differente de 10.0 dans convect3:
      cvflag_grav = .TRUE.

      return
      end
