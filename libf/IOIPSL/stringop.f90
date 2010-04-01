MODULE stringop

!$Header: /home/ioipsl/CVSROOT/IOIPSL/src/stringop.f90,v 2.0 2004/04/05 14:47:51 adm Exp $
!
!---------------------------------------------------------------------
!-
  INTEGER, DIMENSION(30) :: &
       & prime=(/1,2,3,5,7,11,13,17,19,23,29,31,37,41,43, &
       & 47,53,59,61,67,71,73,79,83,89,97,101,103,107,109/)
!-
!---------------------------------------------------------------------
CONTAINS
!=
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
!=
   INTEGER FUNCTION cntpos (c_c,l_c,c_r,l_r)
!---------------------------------------------------------------------
!- Finds number of occurences of c_r in c_c
!---------------------------------------------------------------------
   IMPLICIT NONE
!-
   CHARACTER(LEN=*),INTENT(in) :: c_c
   INTEGER,INTENT(IN) :: l_c
   CHARACTER(LEN=*),INTENT(in) :: c_r
   INTEGER,INTENT(IN) :: l_r
!-
   INTEGER :: ipos,indx
!---------------------------------------------------------------------
   cntpos = 0
   ipos   = 1
   DO
     indx = INDEX(c_c(ipos:l_c),c_r(1:l_r))
     IF (indx > 0) THEN
       cntpos = cntpos+1
       ipos   = ipos+indx+l_r-1
     ELSE
       EXIT
     ENDIF
   ENDDO
!---------------------
   END FUNCTION cntpos
!=
   INTEGER FUNCTION findpos (c_c,l_c,c_r,l_r)
!---------------------------------------------------------------------
!- Finds position of c_r in c_c
!---------------------------------------------------------------------
   IMPLICIT NONE
!-
   CHARACTER(LEN=*),INTENT(in) :: c_c
   INTEGER,INTENT(IN) :: l_c
   CHARACTER(LEN=*),INTENT(in) :: c_r
   INTEGER,INTENT(IN) :: l_r
!---------------------------------------------------------------------
    findpos = INDEX(c_c(1:l_c),c_r(1:l_r))
    IF (findpos == 0)   findpos=-1
!----------------------
   END FUNCTION findpos
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
!=
   SUBROUTINE struppercase (str)
!---------------------------------------------------------------------
!- Converts a string into uppercase
!---------------------------------------------------------------------
   IMPLICIT NONE
!-
   CHARACTER(LEN=*) :: str
!-
   INTEGER :: i,ic
!---------------------------------------------------------------------
   DO i=1,LEN_TRIM(str)
     ic = IACHAR(str(i:i))
     IF ( (ic >= 97) .AND. (ic <= 122) )   str(i:i) = ACHAR(ic-32)
   ENDDO
!-----------------------------
   END SUBROUTINE struppercase
!=
!------------------
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
!---------------------------------------------------------------------
   sig = 0
   DO i=1,MIN(len_trim(str),30)
      sig = sig  + prime(i)*IACHAR(str(i:i))
   ENDDO
!-----------------------------
 END SUBROUTINE gensig
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
!=
!------------------
END MODULE stringop
