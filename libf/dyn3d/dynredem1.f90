SUBROUTINE dynredem1(fichnom, time, vcov, ucov, teta, q, masse, ps)

  ! From dyn3d/dynredem.F, v 1.2 2004/06/22 11:45:30

  ! Ecriture du fichier de redémarrage au format NetCDF

  USE dimens_m, ONLY : llm, nqmx
  USE paramet_m, ONLY : ip1jm, ip1jmp1
  USE temps, ONLY : itaufin, itau_dyn
  USE abort_gcm_m, ONLY : abort_gcm
  USE iniadvtrac_m, ONLY : tname
  use netcdf, only: nf90_open, nf90_write, nf90_noerr, nf90_inq_varid

  IMPLICIT NONE

  REAL, INTENT (IN) :: vcov(ip1jm, llm), ucov(ip1jmp1, llm), teta(ip1jmp1, llm)
  REAL, INTENT (IN) :: ps(ip1jmp1), masse(ip1jmp1, llm)
  REAL, INTENT (IN) :: q(ip1jmp1, llm, nqmx)
  CHARACTER(len=*), INTENT (IN) :: fichnom

  INCLUDE 'netcdf.inc'

  REAL :: time
  INTEGER :: nid, nvarid
  INTEGER :: ierr
  INTEGER :: iq
  INTEGER :: length
  PARAMETER (length=100)
  REAL :: tab_cntrl(length) ! tableau des parametres du run
  CHARACTER (len=20) :: modname
  CHARACTER (len=80) :: abort_message
  INTEGER :: nb = 0

  !---------------------------------------------------------

  PRINT *, 'Call sequence information: dynredem1'

  modname = 'dynredem1'
  ierr = nf90_open(fichnom, nf90_write, nid)
  IF (ierr/=nf90_noerr) THEN
     PRINT *, 'Pb. d ouverture ' // fichnom
     STOP 1
  END IF

  !  Ecriture/extension de la coordonnee temps

  nb = nb + 1
  ierr = nf90_inq_varid(nid, 'temps', nvarid)
  IF (ierr/=nf90_noerr) THEN
     PRINT *, nf_strerror(ierr)
     abort_message = 'Variable temps n est pas definie'
     CALL abort_gcm(modname, abort_message, ierr)
  END IF
  ierr = nf_put_var1_real(nid, nvarid, nb, time)
  PRINT *, 'Enregistrement pour ', nb, time


  !  Re-ecriture du tableau de controle, itaufin n'est plus defini quand
  !  on passe dans dynredem0
  ierr = nf90_inq_varid(nid, 'controle', nvarid)
  IF (ierr/=nf90_noerr) THEN
     abort_message = 'dynredem1: Le champ <controle> est absent'
     ierr = 1
     CALL abort_gcm(modname, abort_message, ierr)
  END IF
  ierr = nf_get_var_real(nid, nvarid, tab_cntrl)
  tab_cntrl(31) = real(itau_dyn+itaufin)
  ierr = nf_put_var_real(nid, nvarid, tab_cntrl)

  !  Ecriture des champs

  ierr = nf90_inq_varid(nid, 'ucov', nvarid)
  IF (ierr/=nf90_noerr) THEN
     PRINT *, 'Variable ucov n est pas definie'
     STOP 1
  END IF
  ierr = nf_put_var_real(nid, nvarid, ucov)

  ierr = nf90_inq_varid(nid, 'vcov', nvarid)
  IF (ierr/=nf90_noerr) THEN
     PRINT *, 'Variable vcov n est pas definie'
     STOP 1
  END IF
  ierr = nf_put_var_real(nid, nvarid, vcov)

  ierr = nf90_inq_varid(nid, 'teta', nvarid)
  IF (ierr/=nf90_noerr) THEN
     PRINT *, 'Variable teta n est pas definie'
     STOP 1
  END IF
  ierr = nf_put_var_real(nid, nvarid, teta)

  DO iq = 1, nqmx
     ierr = nf90_inq_varid(nid, tname(iq), nvarid)
     IF (ierr/=nf90_noerr) THEN
        PRINT *, 'Variable  ', tname(iq), 'n''est pas définie'
        STOP 1
     END IF
     ierr = nf_put_var_real(nid, nvarid, q(1, 1, iq))
  END DO

  ierr = nf90_inq_varid(nid, 'masse', nvarid)
  IF (ierr/=nf90_noerr) THEN
     PRINT *, 'Variable masse n est pas definie'
     STOP 1
  END IF
  ierr = nf_put_var_real(nid, nvarid, masse)

  ierr = nf90_inq_varid(nid, 'ps', nvarid)
  IF (ierr/=nf90_noerr) THEN
     PRINT *, 'Variable ps n est pas definie'
     STOP 1
  END IF
  ierr = nf_put_var_real(nid, nvarid, ps)

  ierr = nf_close(nid)

END SUBROUTINE dynredem1
