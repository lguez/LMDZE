!!
!! @author Jacques Bellier, Marie-Alice Foujols, Jan Polcher
!! @Version : $Revision: 2.0 $, $Date: 2004/04/05 14:47:47 $
!! 
!! $Header: /home/ioipsl/CVSROOT/IOIPSL/src/def.prec,v 2.0 2004/04/05 14:47:47 adm Exp $
!! 
MODULE defprec 
!---------------------------------------------------------------------
!! This module should be used by every modules
!! to keep the right precision for every variable
!!
!! Set default precision for computation
!---------------------------------------------------------------------
  INTEGER,PARAMETER :: i_4=SELECTED_INT_KIND(9) 
  INTEGER,PARAMETER :: i_8=SELECTED_INT_KIND(13)
  INTEGER,PARAMETER :: r_4=SELECTED_REAL_KIND(6,37) 
  INTEGER,PARAMETER :: r_8=SELECTED_REAL_KIND(15,307)
  INTEGER,PARAMETER :: i_std=i_4, r_std=r_4
!-----------------
END MODULE defprec
