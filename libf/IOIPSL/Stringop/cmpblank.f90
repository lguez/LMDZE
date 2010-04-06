module cmpblank_m
contains
   SUBROUTINE cmpblank (str)
!---------------------------------------------------------------------
!-
!---------------------------------------------------------------------
   CHARACTER(LEN=*),INTENT(inout) :: str
!-
   INTEGER :: lcc,ipb
!---------------------------------------------------------------------
   lcc = LEN_TRIM(str)
   ipb = 1
   DO
     IF (ipb >= lcc)   EXIT
     IF (str(ipb:ipb+1) == '  ') THEN
       str(ipb+1:) = str(ipb+2:lcc)
       lcc = lcc-1
     ELSE
       ipb = ipb+1
     ENDIF
   ENDDO
!-------------------------
   END SUBROUTINE cmpblank
 end module cmpblank_m
