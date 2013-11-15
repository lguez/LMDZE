module nocomma_m

implicit none

contains

!=
   SUBROUTINE nocomma (str)
!---------------------------------------------------------------------
!-
!---------------------------------------------------------------------
   IMPLICIT NONE
!-
   CHARACTER(LEN=*) :: str
!-
   INTEGER :: i
!---------------------------------------------------------------------
   DO i=1,LEN_TRIM(str)
     IF (str(i:i) == ',')   str(i:i) = ' '
   ENDDO
!------------------------
   END SUBROUTINE nocomma

 end module nocomma_m
