!$Header: /home/ioipsl/CVSROOT/IOIPSL/src/Attic/ioipslmpp.f90,v 2.0 2004/04/05 14:50:16 adm Exp $
!-
MODULE ioipslmpp
!---------------------------------------------------------------------
  USE errioipsl, ONLY : histerr
!-
  IMPLICIT NONE
!-
  PRIVATE
  PUBLIC :: ioipslmpp_file, ioipslmpp_addatt
!-
  LOGICAL,SAVE :: ison_mpp=.FALSE., lock=.FALSE.
!-
! Number of distributed dimension for mpp
!-
  INTEGER,PARAMETER :: jpp=4
!-
  INTEGER,SAVE :: pe_number, pe_total_number
  INTEGER,SAVE,DIMENSION(jpp) :: &
 & domain_global_size, domain_local_size, domain_abs_first, &
 & domain_abs_last, domain_halo_start_size, domain_halo_end_size
!-
CONTAINS
!-
!===
!-
SUBROUTINE ioipslmpp_file (filename)
!---------------------------------------------------------------------
!- Update the netCDF file to include the PE number on which this
!- copy of IOIPSL runs.
!- This routine is called by histbeg and not by user anyway
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*),INTENT(inout) :: filename
!-
  INTEGER :: il
  CHARACTER(LEN=3) :: str
!---------------------------------------------------------------------
  IF (ison_mpp) THEN
    WRITE(str,'(I3.3)') pe_number
!-- Tester la taille de la chaine
    il = INDEX(filename,'.nc')
    filename = filename(1:il-1)//'_'//str//'.nc'
  ENDIF
!-
! This as to be done after ioipslmpp
!-
  lock = .TRUE.
!---------------------------------------------------------------------
END SUBROUTINE ioipslmpp_file
!-
!===
!-
SUBROUTINE ioipslmpp_addatt (fid)
!---------------------------------------------------------------------
!- Adds the attributed to the netCDF file.
!- This routine is called by histend and not by user anyway
!---------------------------------------------------------------------
  USE netcdf
!-
  IMPLICIT NONE
!-
  INTEGER,INTENT(in) :: fid
!-
  INTEGER :: iret
!---------------------------------------------------------------------
  IF (ison_mpp) THEN
    iret = NF90_PUT_ATT (fid,NF90_GLOBAL, &
 &      'PE_number',pe_number)
    iret = NF90_PUT_ATT (fid,NF90_GLOBAL, &
 &      'PE_total_number',pe_total_number)
    iret = NF90_PUT_ATT (fid,NF90_GLOBAL, &
 &      'DOMAIN_global_size',domain_global_size(1:jpp))
    iret = NF90_PUT_ATT (fid,NF90_GLOBAL, &
 &      'DOMAIN_local_size',domain_local_size(1:jpp))
    iret = NF90_PUT_ATT (fid,NF90_GLOBAL, &
 &      'DOMAIN_absolute_first_point_number',domain_abs_first(1:jpp))
    iret = NF90_PUT_ATT (fid,NF90_GLOBAL, &
 &      'DOMAIN_absolute_last_point_number',domain_abs_last(1:jpp))
    iret = NF90_PUT_ATT (fid,NF90_GLOBAL, &
 &      'DOMAIN_start_halo_size',domain_halo_start_size(1:jpp))
    iret = NF90_PUT_ATT (fid,NF90_GLOBAL, &
 &      'DOMAIN_end_halo_size',domain_halo_end_size(1:jpp))
  ENDIF
!-
! This as to be done after ioipslmpp
!-
  lock = .TRUE.
!------------------------------
END SUBROUTINE ioipslmpp_addatt
!-
!===
!-
END MODULE ioipslmpp
