module gensig_m

implicit none

contains

!=
   SUBROUTINE gensig (str, sig)
!---------------------------------------------------------------------
!- Generate a signature from the first 30 characters of the string
!- This signature is not unique and thus when one looks for the
!- one needs to also verify the string.
!---------------------------------------------------------------------
   IMPLICIT NONE
!-
   CHARACTER(LEN=*) :: str
   INTEGER          :: sig
!-
   INTEGER :: i
   INTEGER, DIMENSION(30) :: prime=(/1,2,3,5,7,11,13,17,19,23,29,31,37,41,43, &
        47,53,59,61,67,71,73,79,83,89,97,101,103,107,109/)
!---------------------------------------------------------------------
   sig = 0
   DO i=1,MIN(len_trim(str),30)
      sig = sig  + prime(i)*IACHAR(str(i:i))
   ENDDO
!-----------------------------
 END SUBROUTINE gensig

end module gensig_m
