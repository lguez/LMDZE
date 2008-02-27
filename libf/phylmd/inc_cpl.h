!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/inc_cpl.h,v 1.2 2005/04/29 14:03:26 fairhead Exp $
!
!
!
! -- inc_cpl.h   1998-04
!    **********
!@
!@  Contents : variables describing pipe and field names
!@  --------
!@
!@ -- cl_write  : for fields to write
!@
!@ -- cl_read  : for fields to read
!@
!     -------------------------------------------------------------------
!
      INTEGER jpread, jpwrit
      PARAMETER (jpread=0, jpwrit=1)
      CHARACTER(len=8) cl_writ(jpmaxfld), cl_read(jpmaxfld)
      CHARACTER(len=8) cl_f_writ(jpmaxfld), cl_f_read(jpmaxfld)
      COMMON / comcpl / cl_writ, cl_read, cl_f_writ, cl_f_read
!     -------------------------------------------------------------------
