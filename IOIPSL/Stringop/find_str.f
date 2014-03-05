module find_str_m

implicit none

contains

!=
   SUBROUTINE find_str (nb_str,str_tab,str_len_tab,str,pos)
!---------------------------------------------------------------------
!- This subroutine looks for a string in a table
!---------------------------------------------------------------------
!- INPUT
!-   nb_str      : length of table
!-   str_tab     : Table  of strings
!-   str_len_tab : Table  of string-length
!-   str         : Target we are looking for
!- OUTPUT
!-   pos         : -1 if str not found, else value in the table
!---------------------------------------------------------------------
   IMPLICIT NONE
!-
   INTEGER :: nb_str
   CHARACTER(LEN=*),DIMENSION(nb_str) :: str_tab
   INTEGER,DIMENSION(nb_str) :: str_len_tab
   CHARACTER(LEN=*) :: str
   INTEGER :: pos
!-
   INTEGER :: i,il
!---------------------------------------------------------------------
   pos = -1
   il = LEN_TRIM(str)
   IF ( nb_str > 0 ) THEN
      DO i=1,nb_str
         IF (     (INDEX(str_tab(i),str(1:il)) > 0) &
              .AND.(str_len_tab(i) == il) ) THEN
            pos = i
            EXIT
         ENDIF
      ENDDO
   ENDIF
!-------------------------
   END SUBROUTINE find_str

 end module find_str_m
