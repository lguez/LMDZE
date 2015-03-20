MODULE errioipsl

  ! From errioipsl.f90, version 2.0, 2004/04/05 14:47:47

  IMPLICIT NONE

CONTAINS

  SUBROUTINE histerr(plev, pcname, pstr1, pstr2, pstr3)

    INTEGER, intent(in):: plev
    !- plev   : Category of message to be reported to the user
    !-          1 = Note to the user
    !-          2 = Warning to the user
    !-          3 = Fatal error

    CHARACTER(LEN=*), intent(in):: pcname 
    ! name of subroutine which has called histerr

    CHARACTER(LEN=*), intent(in):: pstr1, pstr2, pstr3
    ! strings containing the explanations to the user

    ! Local:
    CHARACTER(LEN=30), DIMENSION(3) :: pemsg = &
         (/ "NOTE TO THE USER FROM ROUTINE ", &
         "WARNING FROM ROUTINE          ", &
         "FATAL ERROR FROM ROUTINE      " /)

    !---------------------------------------------------------------------

    IF ((plev >= 1).AND.(plev <= 3)) THEN
       print '(A, " ", A)', TRIM(pemsg(plev)), TRIM(pcname)
       print '(" --> ", a)', pstr1
       print '(" --> ", a)', pstr2
       print '(" --> ", a)', pstr3
    ENDIF
    IF (plev == 3) THEN
       print *, 'Fatal error from IOIPSL'
       STOP 1
    ENDIF

  END SUBROUTINE histerr

END MODULE errioipsl
