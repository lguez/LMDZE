SUBROUTINE dynredem1(fichnom, time, vcov, ucov, teta, q, masse, ps)

  ! From dyn3d/dynredem.F, v 1.2 2004/06/22 11:45:30

  ! Ecriture du fichier de red�marrage au format NetCDF

  USE dimens_m, ONLY : iim, jjm, llm, nqmx
  USE temps, ONLY : itaufin, itau_dyn
  USE iniadvtrac_m, ONLY : tname
  use netcdf, only: nf90_write
  use netcdf95, only: nf95_close, nf95_gw_var, nf95_inq_varid, nf95_open, &
       nf95_put_var

  IMPLICIT NONE

  CHARACTER(len=*), INTENT (IN) :: fichnom
  REAL, INTENT (IN):: time
  REAL, INTENT (IN) :: vcov(iim + 1, jjm, llm), ucov(iim+1, jjm+1, llm)
  REAL, INTENT (IN) :: teta(iim+1, jjm+1, llm)
  REAL, INTENT (IN) :: q(iim+1, jjm+1, llm, nqmx)
  REAL, INTENT (IN) :: ps(iim+1, jjm+1), masse(iim+1, jjm+1, llm)

  ! Variables local to the procedure:
  INTEGER nid, nvarid, ierr
  INTEGER iq
  REAL, pointer:: tab_cntrl(:) ! tableau des param�tres du run
  INTEGER :: nb = 0

  !---------------------------------------------------------

  PRINT *, 'Call sequence information: dynredem1'

  call nf95_open(fichnom, nf90_write, nid)

  ! �criture/extension de la coordonn�e temps
  nb = nb + 1
  call nf95_inq_varid(nid, 'temps', nvarid)
  call nf95_put_var(nid, nvarid, time, start=(/nb/))
  PRINT *, 'Enregistrement pour :'
  print *, "nb = ", nb
  print *, "time = ", time

  ! R�criture du tableau de contr�le, "itaufin" n'est plus d�fini quand
  ! on passe dans "dynredem0"
  call nf95_inq_varid(nid, 'controle', nvarid)
  call nf95_gw_var(nid, nvarid, tab_cntrl)
  tab_cntrl(31) = real(itau_dyn + itaufin)
  call nf95_put_var(nid, nvarid, tab_cntrl)
  deallocate(tab_cntrl) ! pointer

  ! �criture des champs

  call nf95_inq_varid(nid, 'ucov', nvarid)
  call nf95_put_var(nid, nvarid, ucov)

  call nf95_inq_varid(nid, 'vcov', nvarid)
  call nf95_put_var(nid, nvarid, vcov)

  call nf95_inq_varid(nid, 'teta', nvarid)
  call nf95_put_var(nid, nvarid, teta)

  DO iq = 1, nqmx
     call nf95_inq_varid(nid, tname(iq), nvarid)
     call nf95_put_var(nid, nvarid, q(:, :, :, iq))
  END DO

  call nf95_inq_varid(nid, 'masse', nvarid)
  call nf95_put_var(nid, nvarid, masse)

  call nf95_inq_varid(nid, 'ps', nvarid)
  call nf95_put_var(nid, nvarid, ps)

  call nf95_close(nid)

END SUBROUTINE dynredem1
