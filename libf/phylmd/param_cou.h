!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/param_cou.h,v 1.2 2005/02/07 16:41:55 fairhead Exp $
!
! $Id: param_cou.h,v 1.2 2005/02/07 16:41:55 fairhead Exp $
!
! -- param_cou.h
!
        INTEGER jpmaxfld
       PARAMETER(jpmaxfld = 40)        ! Maximum number of fields exchanged
                                        ! between ocean and atmosphere
! -- LOOP
       INTEGER jpflda2o1
       PARAMETER(jpflda2o1 = 13)        ! Number of fields exchanged from
                                         ! atmosphere to ocean via flx.F
! -- LOOP
       INTEGER jpflda2o2
       PARAMETER(jpflda2o2 = 6)         ! Number of fields exchanged from
                                         ! atmosphere to ocean via tau.F
!
       INTEGER jpfldo2a
       PARAMETER(jpfldo2a = 4)          ! Number of fields exchanged from
                                         ! ocean to atmosphere
!
