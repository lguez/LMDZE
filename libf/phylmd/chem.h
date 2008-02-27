!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/chem.h,v 1.1.1.1 2004/05/19 12:53:08 lmdzadmin Exp $
!
      INTEGER idms, iso2, iso4, ih2s, idmso, imsa, ih2o2
      PARAMETER (idms=1, iso2=2, iso4=3)
      PARAMETER (ih2s=4, idmso=5, imsa=6, ih2o2=7)

      REAL n_avogadro, masse_s, masse_so4, rho_water, rho_ice
      PARAMETER (n_avogadro=6.02E23)                  !--molec mol-1
      PARAMETER (masse_s=32.0)                        !--g mol-1
      PARAMETER (masse_so4=96.0)                      !--g mol-1
      PARAMETER (rho_water=1000.0)                    !--kg m-3
      PARAMETER (rho_ice=500.0)                       !--kg m-3

