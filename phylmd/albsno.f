module albsno_m

  ! From phylmd/interface_surf.F90,v 1.8 2005/05/25 13:10:09

  IMPLICIT none

contains

  SUBROUTINE albsno(klon, knon,dtime,agesno,alb_neig_grid, precip_snow)

    INTEGER :: klon, knon
    INTEGER, PARAMETER :: nvm = 8
    REAL   :: dtime
    REAL, dimension(klon,nvm) :: veget
    REAL, DIMENSION(klon) :: alb_neig_grid, agesno, precip_snow

    INTEGER :: i, nv

    REAL, DIMENSION(nvm),SAVE :: init, decay
    REAL :: as
    DATA init /0.55, 0.14, 0.18, 0.29, 0.15, 0.15, 0.14, 0./
    DATA decay/0.30, 0.67, 0.63, 0.45, 0.40, 0.14, 0.06, 1./

    veget = 0.
    veget(:,1) = 1.     ! desert partout
    DO i = 1, knon
       alb_neig_grid(i) = 0.0
    ENDDO
    DO nv = 1, nvm
       DO i = 1, knon
          as = init(nv)+decay(nv)*EXP(-agesno(i)/5.)
          alb_neig_grid(i) = alb_neig_grid(i) + veget(i,nv)*as
       ENDDO
    ENDDO
    !
    !! modilation en fonction de l'age de la neige
    !
    DO i = 1, knon
       agesno(i)  = (agesno(i) + (1.-agesno(i)/50.)*dtime/86400.)&
            &             * EXP(-1.*MAX(0.0,precip_snow(i))*dtime/0.3)
       agesno(i) =  MAX(agesno(i),0.0)
    ENDDO

  END SUBROUTINE albsno

end module albsno_m
