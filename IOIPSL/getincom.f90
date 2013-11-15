MODULE getincom

  ! From getincom.f90, version 2.0 2004/04/05 14:47:48

  use gensig_m, only: gensig
  use find_sig_m, only: find_sig
  use getincom2, only: nb_keys, keysig, keystr, getfill, getdbwl, getdbrl, &
       getfili, getdbwi, getdbri, getfilr, getdbwr, getdbrr

  IMPLICIT NONE

  PRIVATE
  PUBLIC getin

  INTERFACE getin
     MODULE PROCEDURE getinrs, getinis, getinls
  END INTERFACE

CONTAINS

  SUBROUTINE getinrs(MY_TARGET, ret_val)

    ! Get a real scalar. We first check whether we find it in the
    ! database and if not we get it from "run.def". "getinr1d" and
    ! "getinr2d" are written on the same pattern.

    CHARACTER(LEN=*) MY_TARGET
    REAL ret_val

    ! Local:
    REAL, DIMENSION(1):: tmp_ret_val
    INTEGER:: target_sig, pos, status = 0, fileorig

    !--------------------------------------------------------------------

    ! Compute the signature of the target
    CALL gensig(MY_TARGET, target_sig)

    ! Do we have this my_target in our database ?

    ! "find_sig" should not be called if "keystr" and "keysig" are not
    ! allocated.
    ! Avoid this problem with a test on "nb_keys":
    if (nb_keys > 0) then
       CALL find_sig(nb_keys, keystr, my_target, keysig, target_sig, pos)
    else
       pos = -1
    end if

    tmp_ret_val(1) = ret_val

    IF (pos < 0) THEN
       ! Get the information out of the file
       CALL getfilr(MY_TARGET, status, fileorig, tmp_ret_val)
       ! Put the data into the database
       CALL getdbwr(MY_TARGET, target_sig, status, fileorig, 1, tmp_ret_val)
    ELSE
       ! Get the value out of the database
       CALL getdbrr (pos, 1, MY_TARGET, tmp_ret_val)
    ENDIF
    ret_val = tmp_ret_val(1)

  END SUBROUTINE getinrs

  !****************************

  SUBROUTINE getinis(MY_TARGET, ret_val)

    ! Get a interer scalar. We first check if we find it
    ! in the database and if not we get it from the run.def

    ! getini1d and getini2d are written on the same pattern


    CHARACTER(LEN=*) :: MY_TARGET
    INTEGER :: ret_val

    INTEGER, DIMENSION(1) :: tmp_ret_val
    INTEGER :: target_sig, pos, status=0, fileorig


    ! Compute the signature of the target

    CALL gensig(MY_TARGET, target_sig)

    ! Do we have this target in our database ?

    CALL find_sig (nb_keys, keystr, my_target, keysig, target_sig, pos)

    tmp_ret_val(1) = ret_val

    IF (pos < 0) THEN
       ! Ge the information out of the file
       CALL getfili(MY_TARGET, status, fileorig, tmp_ret_val)
       ! Put the data into the database
       CALL getdbwi(MY_TARGET, target_sig, status, fileorig, 1, tmp_ret_val)
    ELSE
       ! Get the value out of the database
       CALL getdbri (pos, 1, MY_TARGET, tmp_ret_val)
    ENDIF
    ret_val = tmp_ret_val(1)

  END SUBROUTINE getinis

  !****************************

  !=== LOGICAL INTERFACES

  SUBROUTINE getinls(MY_TARGET, ret_val)

    ! Get a logical scalar. We first check if we find it
    ! in the database and if not we get it from the run.def

    ! getinl1d and getinl2d are written on the same pattern


    CHARACTER(LEN=*) :: MY_TARGET
    LOGICAL :: ret_val

    LOGICAL, DIMENSION(1) :: tmp_ret_val
    INTEGER :: target_sig, pos, status=0, fileorig


    ! Compute the signature of the target

    CALL gensig(MY_TARGET, target_sig)

    ! Do we have this target in our database ?

    if (nb_keys > 0) then
       CALL find_sig(nb_keys, keystr, my_target, keysig, target_sig, pos)
    else
       pos = -1
    end if

    tmp_ret_val(1) = ret_val

    IF (pos < 0) THEN
       ! Ge the information out of the file
       CALL getfill(MY_TARGET, status, fileorig, tmp_ret_val)
       ! Put the data into the database
       CALL getdbwl(MY_TARGET, target_sig, status, fileorig, 1, tmp_ret_val)
    ELSE
       ! Get the value out of the database
       CALL getdbrl (pos, 1, MY_TARGET, tmp_ret_val)
    ENDIF
    ret_val = tmp_ret_val(1)

  END SUBROUTINE getinls

END MODULE getincom
