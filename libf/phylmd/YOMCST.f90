module YOMCST

  implicit none

  ! From phylmd/YOMCST.h, version 1.2 2005/06/06 13:16:33

  ! A1.0 Fundamental constants
  REAL RPI,RCLUM,RHPLA,RKBOL,RNAVO
  ! A1.1 Astronomical constants
  REAL RDAY,REA,REPSM,RSIYEA,RSIDAY,ROMEGA

  ! A1.1.bis Constantes concernant l'orbite de la Terre:
  REAL R_ecc ! excentricité
  real R_peri ! équinoxe
  real R_incl ! inclinaison

  ! A1.2 Geoide
  REAL RA,RG,R1SA
  ! A1.3 Radiation
  !     REAL RSIGMA,RI0
  REAL RSIGMA
  ! A1.4 Thermodynamic gas phase
  REAL R,RMD,RMO3,RMV,RD,RV,RCPD,RCPV,RCVD,RCVV
  REAL RKAPPA,RETV
  ! A1.5,6 Thermodynamic liquid,solid phases
  REAL RCW,RCS
  ! A1.7 Thermodynamic transition of phase
  REAL RLVTT,RLSTT,RLMLT,RTT,RATM
  ! A1.8 Curve of saturation
  REAL RESTT,RALPW,RBETW,RGAMW,RALPS,RBETS,RGAMS
  REAL RALPD,RBETD,RGAMD

  save

end module YOMCST
