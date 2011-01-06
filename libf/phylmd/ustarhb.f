!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/ustarhb.F,v 1.1 2004/06/22 11:45:35 lmdzadmin Exp $
!
      SUBROUTINE ustarhb(knon,u,v,cd_m, ustar)
      use dimens_m
      use dimphy
      use SUPHEC_M
      use yoethf_m
      use fcttre
      IMPLICIT none
c======================================================================
c Laurent Li (LMD/CNRS), le 30 septembre 1998
c Couche limite non-locale. Adaptation du code du CCM3.
c Code non teste, donc a ne pas utiliser.
c======================================================================
c Nonlocal scheme that determines eddy diffusivities based on a
c diagnosed boundary layer height and a turbulent velocity scale.
c Also countergradient effects for heat and moisture are included.
c
c For more information, see Holtslag, A.A.M., and B.A. Boville, 1993:
c Local versus nonlocal boundary-layer diffusion in a global climate
c model. J. of Climate, vol. 6, 1825-1842.
c======================================================================
c
c Arguments:
c
      INTEGER knon ! nombre de points a calculer
      REAL u(klon,klev) ! vitesse U (m/s)
      REAL v(klon,klev) ! vitesse V (m/s)
      REAL cd_m(klon) ! coefficient de friction au sol pour vitesse
      REAL ustar(klon)
c
      INTEGER i, k
      REAL zxt, zxq, zxu, zxv, zxmod, taux, tauy
      REAL zx_alf1, zx_alf2 ! parametres pour extrapolation
      LOGICAL unssrf(klon)  ! unstb pbl w/lvls within srf pbl lyr
      LOGICAL unsout(klon)  ! unstb pbl w/lvls in outer pbl lyr
      LOGICAL check(klon)   ! True=>chk if Richardson no.>critcal
c
      DO i = 1, knon
        zx_alf1 = 1.0
        zx_alf2 = 1.0 - zx_alf1
        zxu = u(i,1)*zx_alf1+u(i,2)*zx_alf2
        zxv = v(i,1)*zx_alf1+v(i,2)*zx_alf2
        zxmod = 1.0+SQRT(zxu**2+zxv**2)
        taux = zxu *zxmod*cd_m(i)
        tauy = zxv *zxmod*cd_m(i)
        ustar(i) = SQRT(taux**2+tauy**2)
c       print*,'Ust ',zxu,zxmod,taux,ustar(i)
      ENDDO
c
      return
      end
