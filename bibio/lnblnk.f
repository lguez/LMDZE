
! $Header: /home/cvsroot/LMDZ4/libf/bibio/lnblnk.F,v 1.1.1.1 2004/05/19
! 12:53:05 lmdzadmin Exp $

INTEGER FUNCTION lnblnk(letter)

  ! --------------------------------------------------------
  ! Fonction qui determine la longeur d'un string sans les
  ! blancs qui suivent. Le critere pour determiner la fin du
  ! string est, trois blancs de suite
  ! ---------------------------------------------------------
  ! ARGUMENTS
  ! +++++++++
  ! letter: CHARACTER*xxx (xxx < imax)
  ! le string dont on determine la longuer
  ! lnblnk: INTEGER
  ! le nombre de characteres

  ! PARAMETER
  ! +++++++++
  ! imax : INTEGER
  ! le nombre maximale de character que peut contenir le string
  ! a traiter

  IMPLICIT NONE
  INTEGER i, imax
  PARAMETER (imax=256)
  CHARACTER *256 letter

  i = 0

10 i = i + 1
  IF (letter(i:i+3)=='   ') GO TO 20
  GO TO 10

20 lnblnk = i - 1

  RETURN
END FUNCTION lnblnk

