MODULE getparam

  ! From dyn3d/getparam.F90, version 1.1.1.1 2004/05/19 12:53:07

   USE getincom

   INTERFACE getpar
     MODULE PROCEDURE getparamr, getparami, getparaml
   END INTERFACE

   private getparamr, getparami, getparaml
   INTEGER, PARAMETER :: out_eff=99

CONTAINS

  SUBROUTINE ini_getparam(fichier)
    !
    IMPLICIT NONE
    !
    CHARACTER*(*) :: fichier
    open(out_eff, file=fichier, status='unknown', form='formatted') 
  END SUBROUTINE ini_getparam

  !**********************************************************

  SUBROUTINE fin_getparam
    !
    IMPLICIT NONE
    !
    close(out_eff)

  END SUBROUTINE fin_getparam

  !**********************************************************

  SUBROUTINE getparamr(TARGET, def_val, ret_val, comment)
    !
    IMPLICIT NONE
    !
    !   Get a real scalar. We first check if we find it
    !   in the database and if not we get it from the run.def
    !
    !   getinr1d and getinr2d are written on the same pattern
    !
    CHARACTER*(*) :: TARGET
    REAL :: def_val
    REAL :: ret_val
    CHARACTER*(*) :: comment

    ret_val=def_val
    call getin(TARGET, ret_val)

    write(out_eff, *) '******'
    write(out_eff, *) comment
    write(out_eff, *) TARGET, '=', ret_val

  END SUBROUTINE getparamr

  !**********************************************************

  SUBROUTINE getparami(TARGET, def_val, ret_val, comment)
    !
    IMPLICIT NONE
    !
    !   Get a real scalar. We first check if we find it
    !   in the database and if not we get it from the run.def
    !
    !   getinr1d and getinr2d are written on the same pattern
    !
    CHARACTER*(*) :: TARGET
    INTEGER :: def_val
    INTEGER :: ret_val
    CHARACTER*(*) :: comment

    ret_val=def_val
    call getin(TARGET, ret_val)

    write(out_eff, *) '***'
    write(out_eff, *) '*** ', comment, ' ***'
    write(out_eff, *) comment
    write(out_eff, *) TARGET, '=', ret_val

  END SUBROUTINE getparami

  !**********************************************************

  SUBROUTINE getparaml(TARGET, def_val, ret_val, comment)
    !
    IMPLICIT NONE
    !
    !   Get a real scalar. We first check if we find it
    !   in the database and if not we get it from the run.def
    !
    !   getinr1d and getinr2d are written on the same pattern
    !
    CHARACTER*(*) :: TARGET
    LOGICAL :: def_val
    LOGICAL :: ret_val
    CHARACTER*(*) :: comment

    ret_val=def_val
    call getin(TARGET, ret_val)

    write(out_eff, *) '***'
    write(out_eff, *) '*** ', comment, ' ***'
    write(out_eff, *) TARGET, '=', ret_val

  END SUBROUTINE getparaml

END MODULE getparam
