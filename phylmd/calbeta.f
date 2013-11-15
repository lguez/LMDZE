      SUBROUTINE calbeta(dtime,indice,knon,snow,qsol,
     .                    vbeta,vcal,vdif)
      use dimens_m
      use indicesol
      use dimphy
      use conf_gcm_m
      use SUPHEC_M
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) (adaptation du GCM du LMD)
c date: 19940414
c======================================================================
c
c Calculer quelques parametres pour appliquer la couche limite
c ------------------------------------------------------------
      REAL tau_gl ! temps de relaxation pour la glace de mer
ccc      PARAMETER (tau_gl=86400.0*30.0)
      PARAMETER (tau_gl=86400.0*5.0)
      REAL mx_eau_sol
      PARAMETER (mx_eau_sol=150.0)
c
      REAL calsol, calsno, calice ! epaisseur du sol: 0.15 m
      PARAMETER (calsol=1.0/(2.5578E+06*0.15))
      PARAMETER (calsno=1.0/(2.3867E+06*0.15))
      PARAMETER (calice=1.0/(5.1444E+06*0.15))
C
      INTEGER i
c
      REAL dtime
      REAL snow(klon), qsol(klon)
      INTEGER indice, knon
C
      REAL vbeta(klon)
      REAL vcal(klon)
      REAL vdif(klon)
C

      IF (indice.EQ.is_oce) THEN
      DO i = 1, knon
          vcal(i)   = 0.0
          vbeta(i)  = 1.0
          vdif(i) = 0.0
      ENDDO
      ENDIF
c
      IF (indice.EQ.is_sic) THEN
      DO i = 1, knon
          vcal(i) = calice
          IF (snow(i) .GT. 0.0) vcal(i) = calsno
          vbeta(i)  = 1.0
          vdif(i) = 1.0/tau_gl
ccc          vdif(i) = calice/tau_gl ! c'etait une erreur
      ENDDO
      ENDIF
c
      IF (indice.EQ.is_ter) THEN
      DO i = 1, knon
          vcal(i) = calsol
          IF (snow(i) .GT. 0.0) vcal(i) = calsno
          vbeta(i)  = MIN(2.0*qsol(i)/mx_eau_sol, 1.0)
          vdif(i) = 0.0
      ENDDO
      ENDIF
c
      IF (indice.EQ.is_lic) THEN
      DO i = 1, knon
          vcal(i) = calice
          IF (snow(i) .GT. 0.0) vcal(i) = calsno
          vbeta(i)  = 1.0
          vdif(i) = 0.0
      ENDDO
      ENDIF
c
      RETURN
      END
