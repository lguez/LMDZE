!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/atm2geo.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
C
      SUBROUTINE atm2geo ( im, jm, pte, ptn, plon, plat, pxx, pyy, pzz )
cc
cc Change wind local atmospheric coordinates to
cc geocentric
cc
c$$$      INCLUDE 'param.h'
c
      INTEGER, INTENT (in)              :: im, jm
      REAL, DIMENSION (im,jm), INTENT (in) :: pte, ptn
      REAL, DIMENSION (im,jm), INTENT (in) :: plon, plat
      REAL, DIMENSION (im,jm), INTENT(out) :: pxx, pyy, pzz
c
      REAL, PARAMETER :: rpi = 3.141592653E0
      REAL, PARAMETER :: rad = rpi / 180.0E0
c
      REAL, DIMENSION (im,jm) :: zsinlon, zcoslon
      REAL, DIMENSION (im,jm) :: zsinlat, zcoslat
c
      LOGICAL, SAVE :: linit = .FALSE.
c
c$$$      IF ( .NOT. linit ) THEN 
          zsinlon = SIN (rad * plon)
          zcoslon = COS (rad * plon)
          zsinlat = SIN (rad * plat)
          zcoslat = COS (rad * plat)
          linit = .TRUE.
c$$$      ENDIF 
c
      pxx = - zsinlon * pte - zsinlat * zcoslon * ptn
      pyy =   zcoslon * pte - zsinlat * zsinlon * ptn
      pzz =   zcoslat * ptn
c
c Value at North Pole
      pxx ( :,  1) = - ptn ( 1, 1)
      pyy ( :,  1) = - pte ( 1, 1)
      pzz ( :,  1) = 0.0
c Value at South Pole
      pxx ( :, jm) = + ptn ( 1, jm)
      pyy ( :, jm) = + pte ( 1, jm)
      pzz ( :, jm) = 0.0
c
      RETURN 
      END SUBROUTINE atm2geo
