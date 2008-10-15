SUBROUTINE dynredem1(fichnom, time, vcov, ucov, teta, q, masse, ps)

  ! From dyn3d/dynredem.F, v 1.2 2004/06/22 11:45:30

  ! Ecriture du fichier de redémarrage au format NetCDF

  USE dimens_m, ONLY : iim, jjm, llm, nqmx
  USE temps, ONLY : itaufin, itau_dyn
  USE iniadvtrac_m, ONLY : tname
  use netcdf, only: nf90_open, nf90_write, nf90_noerr, nf90_put_var, &
       nf90_get_var, nf90_close
  use netcdf95, only: nf95_inq_varid

  IMPLICIT NONE

  REAL, INTENT (IN) :: vcov(iim + 1, jjm, llm), ucov(iim+1, jjm+1, llm)
  REAL, INTENT (IN) :: teta(iim+1, jjm+1, llm)
  REAL, INTENT (IN) :: ps(iim+1, jjm+1), masse(iim+1, jjm+1, llm)
  REAL, INTENT (IN) :: q(iim+1, jjm+1, llm, nqmx)
  CHARACTER(len=*), INTENT (IN) :: fichnom

  REAL :: time
  INTEGER :: nid, nvarid
  INTEGER :: ierr
  INTEGER :: iq
  INTEGER :: length
  PARAMETER (length=100)
  REAL :: tab_cntrl(length) ! tableau des parametres du run
  INTEGER :: nb = 0

  !---------------------------------------------------------

  PRINT *, 'Call sequence information: dynredem1'

  ierr = nf90_open(fichnom, nf90_write, nid)
  IF (ierr/=nf90_noerr) THEN
     PRINT *, 'Pb. d ouverture ' // fichnom
     STOP 1
  END IF

  !  Ecriture/extension de la coordonnee temps

  nb = nb + 1
  call nf95_inq_varid(nid, 'temps', nvarid)
  ierr = nf90_put_var(nid, nvarid, time, (/nb/))
  PRINT *, 'Enregistrement pour ', nb, time


  !  Re-ecriture du tableau de controle, itaufin n'est plus defini quand
  !  on passe dans dynredem0
  call nf95_inq_varid(nid, 'controle', nvarid)
  ierr = nf90_get_var(nid, nvarid, tab_cntrl)
  tab_cntrl(31) = real(itau_dyn+itaufin)
  ierr = nf90_put_var(nid, nvarid, tab_cntrl)

  !  Ecriture des champs

  call nf95_inq_varid(nid, 'ucov', nvarid)
  ierr = nf90_put_var(nid, nvarid, ucov)

  call nf95_inq_varid(nid, 'vcov', nvarid)
  ierr = nf90_put_var(nid, nvarid, vcov)

  call nf95_inq_varid(nid, 'teta', nvarid)
  ierr = nf90_put_var(nid, nvarid, teta)

  DO iq = 1, nqmx
     call nf95_inq_varid(nid, tname(iq), nvarid)
     ierr = nf90_put_var(nid, nvarid, q(:, :, :, iq))
  END DO

  call nf95_inq_varid(nid, 'masse', nvarid)
  ierr = nf90_put_var(nid, nvarid, masse)

  call nf95_inq_varid(nid, 'ps', nvarid)
  ierr = nf90_put_var(nid, nvarid, ps)

  ierr = nf90_close(nid)

END SUBROUTINE dynredem1
