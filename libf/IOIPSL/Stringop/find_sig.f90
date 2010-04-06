module find_sig_m

implicit none

contains

!=
!------------------
   SUBROUTINE find_sig (nb_sig, str_tab, str, sig_tab, sig, pos)
!---------------------------------------------------------------------
!- Find the string signature in a list of signatures
!---------------------------------------------------------------------
!- INPUT
!-   nb_sig      : length of table of signatures
!-   str_tab     : Table of strings
!-   str         : Target string we are looking for
!-   sig_tab     : Table of signatures
!-   sig         : Target signature we are looking for
!- OUTPUT
!-   pos         : -1 if str not found, else value in the table
!---------------------------------------------------------------------
   IMPLICIT NONE
!-
   INTEGER :: nb_sig
   CHARACTER(LEN=*),DIMENSION(nb_sig) :: str_tab
   CHARACTER(LEN=*) :: str
   INTEGER, DIMENSION(nb_sig) :: sig_tab
   INTEGER :: sig
!-
   INTEGER :: pos
   INTEGER, DIMENSION(nb_sig) :: loczeros
!-
   INTEGER :: il, my_len
   INTEGER, DIMENSION(1) :: minpos
!---------------------------------------------------------------------
!-
   pos = -1
   il = LEN_TRIM(str)
!-
   IF ( nb_sig > 0 ) THEN
      !
      loczeros = ABS(sig_tab(1:nb_sig)-sig)
      !
      IF ( COUNT(loczeros < 1) == 1 ) THEN
         !
         minpos = MINLOC(loczeros)
         my_len = LEN_TRIM(str_tab(minpos(1)))
         IF ( (INDEX(str_tab(minpos(1)),str(1:il)) > 0) &
                 .AND.(my_len == il) ) THEN
            pos = minpos(1)
         ENDIF
         !
      ELSE IF ( COUNT(loczeros < 1) > 1 ) THEN
         !
         DO WHILE (COUNT(loczeros < 1) >= 1 .AND. pos < 0 )
            minpos = MINLOC(loczeros)
            my_len = LEN_TRIM(str_tab(minpos(1)))
            IF ( (INDEX(str_tab(minpos(1)),str(1:il)) > 0) &
                 .AND.(my_len == il) ) THEN
               pos = minpos(1)
            ELSE
               loczeros(minpos(1)) = 99999
            ENDIF
         ENDDO
         !
      ENDIF
      !
   ENDIF
!-
 END SUBROUTINE find_sig

end module find_sig_m
