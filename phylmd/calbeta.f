SUBROUTINE calbeta(dtime, indice, knon, snow, qsol, vbeta, vcal, vdif)
  USE dimens_m
  USE indicesol
  USE dimphy
  USE conf_gcm_m
  USE suphec_m
  IMPLICIT NONE
  ! ======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS) (adaptation du GCM du LMD)
  ! date: 19940414
  ! ======================================================================

  ! Calculer quelques parametres pour appliquer la couche limite
  ! ------------------------------------------------------------
  REAL tau_gl ! temps de relaxation pour la glace de mer
  ! cc      PARAMETER (tau_gl=86400.0*30.0)
  PARAMETER (tau_gl=86400.0*5.0)
  REAL mx_eau_sol
  PARAMETER (mx_eau_sol=150.0)

  REAL calsol, calsno, calice ! epaisseur du sol: 0.15 m
  PARAMETER (calsol=1.0/(2.5578E+06*0.15))
  PARAMETER (calsno=1.0/(2.3867E+06*0.15))
  PARAMETER (calice=1.0/(5.1444E+06*0.15))

  INTEGER i

  REAL dtime
  REAL snow(klon), qsol(klon)
  INTEGER indice, knon

  REAL vbeta(klon)
  REAL vcal(klon)
  REAL vdif(klon)


  IF (indice==is_oce) THEN
    DO i = 1, knon
      vcal(i) = 0.0
      vbeta(i) = 1.0
      vdif(i) = 0.0
    END DO
  END IF

  IF (indice==is_sic) THEN
    DO i = 1, knon
      vcal(i) = calice
      IF (snow(i)>0.0) vcal(i) = calsno
      vbeta(i) = 1.0
      vdif(i) = 1.0/tau_gl
      ! cc          vdif(i) = calice/tau_gl ! c'etait une erreur
    END DO
  END IF

  IF (indice==is_ter) THEN
    DO i = 1, knon
      vcal(i) = calsol
      IF (snow(i)>0.0) vcal(i) = calsno
      vbeta(i) = min(2.0*qsol(i)/mx_eau_sol, 1.0)
      vdif(i) = 0.0
    END DO
  END IF

  IF (indice==is_lic) THEN
    DO i = 1, knon
      vcal(i) = calice
      IF (snow(i)>0.0) vcal(i) = calsno
      vbeta(i) = 1.0
      vdif(i) = 0.0
    END DO
  END IF

  RETURN
END SUBROUTINE calbeta
