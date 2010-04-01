!$Header: /home/ioipsl/CVSROOT/IOIPSL/src/Attic/ioipslmpp.f90,v 2.0 2004/04/05 14:50:16 adm Exp $
!-
MODULE ioipslmpp
!---------------------------------------------------------------------
  USE errioipsl, ONLY : histerr
!-
  IMPLICIT NONE
!-
  PRIVATE
  PUBLIC :: ioipsl_inimpp, ioipslmpp_file, ioipslmpp_addatt
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
SUBROUTINE ioipsl_inimpp &
 & (petotnb, penb, pglo, ploc, pabsf, pabsl, phals, phale)
!---------------------------------------------------------------------
!- This routine sets up the MPP activity of IOIPSL.
!- It will store all the PE information and allow it to be stored
!- in the netCDF file and change the file names.
!-
!- INPUT
!-
!- penb    : process number
!- petotnb : total number of process
!- pglo(1) : total number of points in first direction
!- pglo(2) : total number of points in second direction
!- ploc(1) : local number of points in first direction
!- ploc(2) : local number of points in second direction
!- pabsf(1) : absolute position of first local point for
!-            first dimension
!- pabsf(2) : absolute position of first local point for
!-            second dimension
!- pabsl(1) : absolute position of last local point for
!-            first dimension
!- pabsl(2) : absolute position of last local point for
!-            second dimension
!- phals(1) : start halo size in first direction
!- phals(2) : start halo size in second direction
!- phale(1) : end halo size in first direction
!- phale(2) : end halo size in second direction
!- phale(2) : end halo size in second direction
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(in) :: penb, petotnb
  INTEGER,DIMENSION(:),INTENT(in) :: &
 & pglo, ploc, pabsf, pabsl, phals, phale
!---------------------------------------------------------------------
  IF (lock) THEN
    CALL histerr (3,'ioipslmpp','ioipslmpp called to late', &
 &                'please call ioipslmpp before first histbeg','')
  ELSE
!-- Take note of the fact that we are on an MPP.
    ison_mpp=.TRUE.
!--
    pe_number = penb
    pe_total_number = petotnb
    domain_global_size(:) = pglo(:)
    domain_local_size(:) = ploc(:)
    domain_abs_first(:) = pabsf(:)
    domain_abs_last(:) = pabsl(:)
    domain_halo_start_size(:) = phals(:)
    domain_halo_end_size(:) = phale(:)
!-- Lock this information into the module
!-- so that it does not get changed
    lock = .TRUE.
  ENDIF
!---------------------------
END SUBROUTINE ioipsl_inimpp
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
