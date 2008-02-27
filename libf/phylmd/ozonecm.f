!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/ozonecm.F,v 1.3 2005/06/06 13:16:33 fairhead Exp $
!
      SUBROUTINE ozonecm(rjour, rlat, paprs, o3)
      use dimens_m
      use dimphy
      use clesphys
      use YOMCST
      IMPLICIT none
C
C The ozone climatology is based on an analytic formula which fits the
C Krueger and Mintzner (1976) profile, as well as the variations with
C altitude and latitude of the maximum ozone concentrations and the total
C column ozone concentration of Keating and Young (1986). The analytic
C formula have been established by J-F Royer (CRNM, Meteo France), who
C also provided us the code.
C
C A. J. Krueger and R. A. Minzner, A Mid-Latitude Ozone Model for the
C 1976 U.S. Standard Atmosphere, J. Geophys. Res., 81, 4477, (1976).
C
C Keating, G. M. and D. F. Young, 1985: Interim reference models for the
C middle atmosphere, Handbook for MAP, vol. 16, 205-229.
C

      real, intent(in):: rjour
      REAL, intent(in):: rlat(klon)
      real, intent(in):: paprs(klon,klev+1)
      REAL, intent(out):: o3(klon,klev)   ! ozone concentration in kg/kg

      REAL tozon
      real pi, pl
      INTEGER i, k
C----------------------------------------------------------
      REAL field(klon,klev+1)
      REAL ps
      PARAMETER (ps=101325.0)
      REAL an, unit, zo3q3
      SAVE an, unit, zo3q3
      REAL mu,gms, zslat, zsint, zcost, z, ppm, qpm, a
      REAL asec, bsec, aprim, zo3a3
C----------------------------------------------------------
c         data an /365.25/   (meteo)
      DATA an /360.00/
      DATA unit /2.1415e-05/
      DATA zo3q3 /4.0e-08/

      pi = 4.0 * ATAN(1.0)
      DO k = 1, klev
      DO i = 1, klon
      zslat = SIN(pi/180.*rlat(i))
      zsint = SIN(2.*pi*(rjour+15.)/an)
      zcost = COS(2.*pi*(rjour+15.)/an)
      z = 0.0531+zsint*(-0.001595+0.009443*zslat) +
     .  zcost*(-0.001344-0.00346*zslat) +
     .  zslat**2*(.056222+zslat**2*(-.037609
     . +.012248*zsint+.00521*zcost+.008890*zslat))
      zo3a3 = zo3q3/ps/2.
      z = z-zo3q3*ps
      gms = z
      mu = ABS(sin(pi/180.*rlat(i)))
      ppm = 800.-(500.*zslat+150.*zcost)*zslat
      qpm = 1.74e-5-(7.5e-6*zslat+1.7e-6*zcost)*zslat
      bsec = 2650.+5000.*zslat**2
      a = 4.0*(bsec)**(3./2.)*(ppm)**(3./2.)*(1.0+(bsec/ps)**(3./2.))
      a = a/(bsec**(3./2.)+ppm**(3./2.))**2
      aprim = (2.666666*qpm*ppm-a*gms)/(1.0-a)
      aprim = amax1(0.,aprim)
      asec = (gms-aprim)*(1.0+(bsec/ps)**(3./2.))
      asec = AMAX1(0.0,asec)
      aprim = gms-asec/(1.+(bsec/ps)**(3./2.))
      pl = paprs(i,k)
      tozon = aprim/(1.+3.*(ppm/pl)**2)+asec/(1.+(bsec/pl)**(3./2.))
     .  + zo3a3*pl*pl
      tozon = tozon / 9.81  ! en kg/m**2
      tozon = tozon / unit  ! en kg/m**2  > u dobson (10e-5 m)
      tozon = tozon / 1000. ! en cm
      field(i,k) = tozon
      ENDDO
      ENDDO
c
      DO i = 1, klon
         field(i,klev+1) = 0.0
      ENDDO
      DO k = 1, klev
        DO i = 1, klon
          o3(i,k) = field(i,k) - field(i,k+1)
          IF (.not. bug_ozone) then
c         convert o3 into kg/kg         
            o3(i,k)=MAX(o3(i,k),1.0e-12)*RG/46.6968
     .               /(paprs(i,k)-paprs(i,k+1))
          ENDIF
        ENDDO
      ENDDO
c
      RETURN
      END
