SUBROUTINE dynredem1(fichnom,time, vcov,ucov,teta,q,nq,masse,ps)

  ! From dyn3d/dynredem.F,v 1.2 2004/06/22 11:45:30

  !  Ecriture du fichier de redémarrage au format NetCDF

  use dimens_m
  use paramet_m
  use comgeom
  use temps
  use abort_gcm_m, only: abort_gcm
  use advtrac_m, only: tname

  IMPLICIT NONE
  include "netcdf.inc"

  INTEGER nq
  REAL, intent(in):: vcov(ip1jm,llm),ucov(ip1jmp1,llm) 
  REAL teta(ip1jmp1,llm)                   
  REAL, intent(in):: ps(ip1jmp1)
  real masse(ip1jmp1,llm)                   
  REAL, intent(in):: q(ip1jmp1,llm,nq)
  CHARACTER(len=*) fichnom

  REAL time
  INTEGER nid, nvarid
  INTEGER ierr
  INTEGER iq
  INTEGER length
  PARAMETER (length = 100)
  REAL tab_cntrl(length) ! tableau des parametres du run
  character(len=20) modname
  character(len=80) abort_message

  INTEGER nb
  SAVE nb
  DATA nb / 0 /

  !---------------------------------------------------------

  print *, "Call sequence information: dynredem1"

  modname = 'dynredem1'
  ierr = NF_OPEN(fichnom, NF_WRITE, nid)
  IF (ierr .NE. NF_NOERR) THEN
     PRINT*, "Pb. d ouverture "//fichnom
     stop 1
  ENDIF

  !  Ecriture/extension de la coordonnee temps

  nb = nb + 1
  ierr = NF_INQ_VARID(nid, "temps", nvarid)
  IF (ierr .NE. NF_NOERR) THEN
     print *, NF_STRERROR(ierr)
     abort_message='Variable temps n est pas definie'
     CALL abort_gcm(modname,abort_message,ierr)
  ENDIF
  ierr = NF_PUT_VAR1_REAL (nid,nvarid,nb,time)
  PRINT*, "Enregistrement pour ", nb, time

  !
  !  Re-ecriture du tableau de controle, itaufin n'est plus defini quand
  !  on passe dans dynredem0
  ierr = NF_INQ_VARID (nid, "controle", nvarid)
  IF (ierr .NE. NF_NOERR) THEN
     abort_message="dynredem1: Le champ <controle> est absent"
     ierr = 1
     CALL abort_gcm(modname,abort_message,ierr)
  ENDIF
  ierr = NF_GET_VAR_REAL(nid, nvarid, tab_cntrl)
  tab_cntrl(31) = REAL(itau_dyn + itaufin)
  ierr = NF_PUT_VAR_REAL (nid,nvarid,tab_cntrl)

  !  Ecriture des champs
  !
  ierr = NF_INQ_VARID(nid, "ucov", nvarid)
  IF (ierr .NE. NF_NOERR) THEN
     PRINT*, "Variable ucov n est pas definie"
     stop 1
  ENDIF
  ierr = NF_PUT_VAR_REAL (nid,nvarid,ucov)

  ierr = NF_INQ_VARID(nid, "vcov", nvarid)
  IF (ierr .NE. NF_NOERR) THEN
     PRINT*, "Variable vcov n est pas definie"
     stop 1
  ENDIF
  ierr = NF_PUT_VAR_REAL (nid,nvarid,vcov)

  ierr = NF_INQ_VARID(nid, "teta", nvarid)
  IF (ierr .NE. NF_NOERR) THEN
     PRINT*, "Variable teta n est pas definie"
     stop 1
  ENDIF
  ierr = NF_PUT_VAR_REAL (nid,nvarid,teta)

  IF(nq.GE.1) THEN
     do iq=1,nq   
        ierr = NF_INQ_VARID(nid, tname(iq), nvarid)
        IF (ierr .NE. NF_NOERR) THEN
           PRINT*, "Variable  ", tname(iq), "n'est pas définie"
           stop 1
        ENDIF
        ierr = NF_PUT_VAR_REAL (nid,nvarid,q(1,1,iq))
     ENDDO
  ENDIF
  !
  ierr = NF_INQ_VARID(nid, "masse", nvarid)
  IF (ierr .NE. NF_NOERR) THEN
     PRINT*, "Variable masse n est pas definie"
     stop 1
  ENDIF
  ierr = NF_PUT_VAR_REAL (nid,nvarid,masse)
  !
  ierr = NF_INQ_VARID(nid, "ps", nvarid)
  IF (ierr .NE. NF_NOERR) THEN
     PRINT*, "Variable ps n est pas definie"
     stop 1
  ENDIF
  ierr = NF_PUT_VAR_REAL (nid,nvarid,ps)

  ierr = NF_CLOSE(nid)

END SUBROUTINE dynredem1
