module ustarhb_m

  IMPLICIT none

contains

  pure function ustarhb(u, v, cd_m)

    ! From LMDZ4/libf/phylmd/ustarhb.F, version 1.1, 2004/06/22 11:45:35

    ! Laurent Li (LMD/CNRS), le 30 septembre 1998
    ! Couche limite non-locale. Adaptation du code du CCM3.
    ! Code non test\'e, donc \`a ne pas utiliser.

    ! Nonlocal scheme that determines eddy diffusivities based on a
    ! diagnosed boundary layer height and a turbulent velocity scale.
    ! Also countergradient effects for heat and moisture are included.

    ! For more information, see Holtslag, A.A.M. and B.A. Boville, 1993:
    ! Local versus nonlocal boundary-layer diffusion in a global climate
    ! model. J. of Climate, vol. 6, 1825-1842.

    REAL, intent(in):: u(:), v(:) ! wind in first layer (m/s)
    REAL, intent(in):: cd_m(:) ! coefficient de friction au sol pour vitesse
    REAL ustarhb(size(u))

    ! Local:
    INTEGER i
    REAL zxmod, taux, tauy

    !---------------------------------------------------------------

    DO i = 1, size(u)
       zxmod = 1. + SQRT(u(i)**2 + v(i)**2)
       taux = u(i) * zxmod * cd_m(i)
       tauy = v(i) * zxmod * cd_m(i)
       ustarhb(i) = SQRT(taux**2 + tauy**2)
    ENDDO

  end function ustarhb

end module ustarhb_m
