module strlowercase_m

implicit none

contains

!=
   SUBROUTINE strlowercase (str)
!---------------------------------------------------------------------
!- Converts a string into lowercase
!---------------------------------------------------------------------
   IMPLICIT NONE
!-
   CHARACTER(LEN=*) :: str
!-
   INTEGER :: i,ic
!---------------------------------------------------------------------
   DO i=1,LEN_TRIM(str)
     ic = IACHAR(str(i:i))
     IF ( (ic >= 65) .AND. (ic <= 90) )   str(i:i) = ACHAR(ic+32)
   ENDDO
!-----------------------------
   END SUBROUTINE strlowercase

 end module strlowercase_m
