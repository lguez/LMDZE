MODULE getincom

  ! From getincom.f90, version 2.0 2004/04/05 14:47:48

  use gensig_m, only: gensig
  use find_sig_m, only: find_sig
  use getincom2, only: nb_keys, keysig, keystr, getfill, getdbwl, getdbrl, &
       getfilc, getdbwc, getdbrc, getfili, getdbwi, getdbri, getfilr, &
       getdbwr, getdbrr

  IMPLICIT NONE

  PRIVATE
  PUBLIC getin

  INTERFACE getin
     MODULE PROCEDURE getinrs, getinr1d, getinr2d, getinis, getini1d, &
          getini2d, getincs, getinc1d, getinc2d, getinls, getinl1d, getinl2d
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

  SUBROUTINE getinr1d(MY_TARGET, ret_val)

    ! See getinrs for details. It is the same thing but for a vector


    CHARACTER(LEN=*) :: MY_TARGET
    REAL, DIMENSION(:) :: ret_val

    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: tmp_ret_val
    INTEGER, SAVE :: tmp_ret_size = 0
    INTEGER :: target_sig, pos, size_of_in, status=0, fileorig


    ! Compute the signature of the target

    CALL gensig(MY_TARGET, target_sig)

    ! Do we have this target in our database ?

    CALL find_sig (nb_keys, keystr, my_target, keysig, target_sig, pos)

    size_of_in = SIZE(ret_val)
    IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
       ALLOCATE (tmp_ret_val(size_of_in))
    ELSE IF (size_of_in > tmp_ret_size) THEN
       DEALLOCATE (tmp_ret_val)
       ALLOCATE (tmp_ret_val(size_of_in))
       tmp_ret_size = size_of_in
    ENDIF
    tmp_ret_val(1:size_of_in) = ret_val(1:size_of_in)

    IF (pos < 0) THEN
       ! Ge the information out of the file
       CALL getfilr(MY_TARGET, status, fileorig, tmp_ret_val)
       ! Put the data into the database
       CALL getdbwr &
            & (MY_TARGET, target_sig, status, fileorig, size_of_in, tmp_ret_val)
    ELSE
       ! Get the value out of the database
       CALL getdbrr (pos, size_of_in, MY_TARGET, tmp_ret_val)
    ENDIF
    ret_val(1:size_of_in) = tmp_ret_val(1:size_of_in)

  END SUBROUTINE getinr1d

  !****************************

  SUBROUTINE getinr2d(MY_TARGET, ret_val)

    ! See getinrs for details. It is the same thing but for a matrix


    CHARACTER(LEN=*) :: MY_TARGET
    REAL, DIMENSION(:, :) :: ret_val

    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: tmp_ret_val
    INTEGER, SAVE :: tmp_ret_size = 0
    INTEGER :: target_sig, pos, size_of_in, size_1, size_2, status=0, fileorig
    INTEGER :: jl, jj, ji


    ! Compute the signature of the target

    CALL gensig(MY_TARGET, target_sig)

    ! Do we have this target in our database ?

    CALL find_sig (nb_keys, keystr, my_target, keysig, target_sig, pos)

    size_of_in = SIZE(ret_val)
    size_1 = SIZE(ret_val, 1)
    size_2 = SIZE(ret_val, 2)
    IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
       ALLOCATE (tmp_ret_val(size_of_in))
    ELSE IF (size_of_in > tmp_ret_size) THEN
       DEALLOCATE (tmp_ret_val)
       ALLOCATE (tmp_ret_val(size_of_in))
       tmp_ret_size = size_of_in
    ENDIF

    jl=0
    DO jj=1, size_2
       DO ji=1, size_1
          jl=jl+1
          tmp_ret_val(jl) = ret_val(ji, jj)
       ENDDO
    ENDDO

    IF (pos < 0) THEN
       ! Ge the information out of the file
       CALL getfilr(MY_TARGET, status, fileorig, tmp_ret_val)
       ! Put the data into the database
       CALL getdbwr &
            & (MY_TARGET, target_sig, status, fileorig, size_of_in, tmp_ret_val)
    ELSE
       ! Get the value out of the database
       CALL getdbrr (pos, size_of_in, MY_TARGET, tmp_ret_val)
    ENDIF

    jl=0
    DO jj=1, size_2
       DO ji=1, size_1
          jl=jl+1
          ret_val(ji, jj) = tmp_ret_val(jl)
       ENDDO
    ENDDO

  END SUBROUTINE getinr2d

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

  SUBROUTINE getini1d(MY_TARGET, ret_val)

    ! See getinis for details. It is the same thing but for a vector


    CHARACTER(LEN=*) :: MY_TARGET
    INTEGER, DIMENSION(:) :: ret_val

    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: tmp_ret_val
    INTEGER, SAVE :: tmp_ret_size = 0
    INTEGER :: target_sig, pos, size_of_in, status=0, fileorig


    ! Compute the signature of the target

    CALL gensig(MY_TARGET, target_sig)

    ! Do we have this target in our database ?

    CALL find_sig (nb_keys, keystr, my_target, keysig, target_sig, pos)

    size_of_in = SIZE(ret_val)
    IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
       ALLOCATE (tmp_ret_val(size_of_in))
    ELSE IF (size_of_in > tmp_ret_size) THEN
       DEALLOCATE (tmp_ret_val)
       ALLOCATE (tmp_ret_val(size_of_in))
       tmp_ret_size = size_of_in
    ENDIF
    tmp_ret_val(1:size_of_in) = ret_val(1:size_of_in)

    IF (pos < 0) THEN
       ! Ge the information out of the file
       CALL getfili(MY_TARGET, status, fileorig, tmp_ret_val)
       ! Put the data into the database
       CALL getdbwi &
            & (MY_TARGET, target_sig, status, fileorig, size_of_in, tmp_ret_val)
    ELSE
       ! Get the value out of the database
       CALL getdbri (pos, size_of_in, MY_TARGET, tmp_ret_val)
    ENDIF
    ret_val(1:size_of_in) = tmp_ret_val(1:size_of_in)

  END SUBROUTINE getini1d

  !****************************

  SUBROUTINE getini2d(MY_TARGET, ret_val)

    ! See getinis for details. It is the same thing but for a matrix


    CHARACTER(LEN=*) :: MY_TARGET
    INTEGER, DIMENSION(:, :) :: ret_val

    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: tmp_ret_val
    INTEGER, SAVE :: tmp_ret_size = 0
    INTEGER :: target_sig, pos, size_of_in, size_1, size_2, status=0, fileorig
    INTEGER :: jl, jj, ji


    ! Compute the signature of the target

    CALL gensig(MY_TARGET, target_sig)

    ! Do we have this target in our database ?

    CALL find_sig (nb_keys, keystr, my_target, keysig, target_sig, pos)

    size_of_in = SIZE(ret_val)
    size_1 = SIZE(ret_val, 1)
    size_2 = SIZE(ret_val, 2)
    IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
       ALLOCATE (tmp_ret_val(size_of_in))
    ELSE IF (size_of_in > tmp_ret_size) THEN
       DEALLOCATE (tmp_ret_val)
       ALLOCATE (tmp_ret_val(size_of_in))
       tmp_ret_size = size_of_in
    ENDIF

    jl=0
    DO jj=1, size_2
       DO ji=1, size_1
          jl=jl+1
          tmp_ret_val(jl) = ret_val(ji, jj)
       ENDDO
    ENDDO

    IF (pos < 0) THEN
       ! Ge the information out of the file
       CALL getfili(MY_TARGET, status, fileorig, tmp_ret_val)
       ! Put the data into the database
       CALL getdbwi &
            & (MY_TARGET, target_sig, status, fileorig, size_of_in, tmp_ret_val)
    ELSE
       ! Get the value out of the database
       CALL getdbri (pos, size_of_in, MY_TARGET, tmp_ret_val)
    ENDIF

    jl=0
    DO jj=1, size_2
       DO ji=1, size_1
          jl=jl+1
          ret_val(ji, jj) = tmp_ret_val(jl)
       ENDDO
    ENDDO

  END SUBROUTINE getini2d

  !****************************

  !=== CHARACTER INTERFACES

  SUBROUTINE getincs(MY_TARGET, ret_val)

    ! Get a CHARACTER scalar. We first check if we find it
    ! in the database and if not we get it from the run.def

    ! getinc1d and getinc2d are written on the same pattern


    CHARACTER(LEN=*) :: MY_TARGET
    CHARACTER(LEN=*) :: ret_val

    CHARACTER(LEN=100), DIMENSION(1) :: tmp_ret_val
    INTEGER :: target_sig, pos, status=0, fileorig


    ! Compute the signature of the target

    CALL gensig(MY_TARGET, target_sig)

    ! Do we have this target in our database ?

    CALL find_sig (nb_keys, keystr, my_target, keysig, target_sig, pos)

    tmp_ret_val(1) = ret_val

    IF (pos < 0) THEN
       ! Ge the information out of the file
       CALL getfilc(MY_TARGET, status, fileorig, tmp_ret_val)
       ! Put the data into the database
       CALL getdbwc(MY_TARGET, target_sig, status, fileorig, 1, tmp_ret_val)
    ELSE
       ! Get the value out of the database
       CALL getdbrc (pos, 1, MY_TARGET, tmp_ret_val)
    ENDIF
    ret_val = tmp_ret_val(1)

  END SUBROUTINE getincs

  !****************************

  SUBROUTINE getinc1d(MY_TARGET, ret_val)

    ! See getincs for details. It is the same thing but for a vector


    CHARACTER(LEN=*) :: MY_TARGET
    CHARACTER(LEN=*), DIMENSION(:) :: ret_val

    CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE, SAVE :: tmp_ret_val
    INTEGER, SAVE :: tmp_ret_size = 0
    INTEGER :: target_sig, pos, size_of_in, status=0, fileorig


    ! Compute the signature of the target

    CALL gensig(MY_TARGET, target_sig)

    ! Do we have this target in our database ?

    CALL find_sig (nb_keys, keystr, my_target, keysig, target_sig, pos)

    size_of_in = SIZE(ret_val)
    IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
       ALLOCATE (tmp_ret_val(size_of_in))
    ELSE IF (size_of_in > tmp_ret_size) THEN
       DEALLOCATE (tmp_ret_val)
       ALLOCATE (tmp_ret_val(size_of_in))
       tmp_ret_size = size_of_in
    ENDIF
    tmp_ret_val(1:size_of_in) = ret_val(1:size_of_in)

    IF (pos < 0) THEN
       ! Ge the information out of the file
       CALL getfilc(MY_TARGET, status, fileorig, tmp_ret_val)
       ! Put the data into the database
       CALL getdbwc &
            & (MY_TARGET, target_sig, status, fileorig, size_of_in, tmp_ret_val)
    ELSE
       ! Get the value out of the database
       CALL getdbrc (pos, size_of_in, MY_TARGET, tmp_ret_val)
    ENDIF
    ret_val(1:size_of_in) = tmp_ret_val(1:size_of_in)

  END SUBROUTINE getinc1d

  !****************************

  SUBROUTINE getinc2d(MY_TARGET, ret_val)

    ! See getincs for details. It is the same thing but for a matrix


    CHARACTER(LEN=*) :: MY_TARGET
    CHARACTER(LEN=*), DIMENSION(:, :) :: ret_val

    CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE, SAVE :: tmp_ret_val
    INTEGER, SAVE :: tmp_ret_size = 0
    INTEGER :: target_sig, pos, size_of_in, size_1, size_2, status=0, fileorig
    INTEGER :: jl, jj, ji


    ! Compute the signature of the target

    CALL gensig(MY_TARGET, target_sig)

    ! Do we have this target in our database ?

    CALL find_sig (nb_keys, keystr, my_target, keysig, target_sig, pos)

    size_of_in = SIZE(ret_val)
    size_1 = SIZE(ret_val, 1)
    size_2 = SIZE(ret_val, 2)
    IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
       ALLOCATE (tmp_ret_val(size_of_in))
    ELSE IF (size_of_in > tmp_ret_size) THEN
       DEALLOCATE (tmp_ret_val)
       ALLOCATE (tmp_ret_val(size_of_in))
       tmp_ret_size = size_of_in
    ENDIF

    jl=0
    DO jj=1, size_2
       DO ji=1, size_1
          jl=jl+1
          tmp_ret_val(jl) = ret_val(ji, jj)
       ENDDO
    ENDDO

    IF (pos < 0) THEN
       ! Ge the information out of the file
       CALL getfilc(MY_TARGET, status, fileorig, tmp_ret_val)
       ! Put the data into the database
       CALL getdbwc &
            & (MY_TARGET, target_sig, status, fileorig, size_of_in, tmp_ret_val)
    ELSE
       ! Get the value out of the database
       CALL getdbrc (pos, size_of_in, MY_TARGET, tmp_ret_val)
    ENDIF

    jl=0
    DO jj=1, size_2
       DO ji=1, size_1
          jl=jl+1
          ret_val(ji, jj) = tmp_ret_val(jl)
       ENDDO
    ENDDO

  END SUBROUTINE getinc2d

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

  !****************************

  SUBROUTINE getinl1d(MY_TARGET, ret_val)

    ! See getinls for details. It is the same thing but for a vector


    CHARACTER(LEN=*) :: MY_TARGET
    LOGICAL, DIMENSION(:) :: ret_val

    LOGICAL, DIMENSION(:), ALLOCATABLE, SAVE :: tmp_ret_val
    INTEGER, SAVE :: tmp_ret_size = 0
    INTEGER :: target_sig, pos, size_of_in, status=0, fileorig


    ! Compute the signature of the target

    CALL gensig(MY_TARGET, target_sig)

    ! Do we have this target in our database ?

    CALL find_sig (nb_keys, keystr, my_target, keysig, target_sig, pos)

    size_of_in = SIZE(ret_val)
    IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
       ALLOCATE (tmp_ret_val(size_of_in))
    ELSE IF (size_of_in > tmp_ret_size) THEN
       DEALLOCATE (tmp_ret_val)
       ALLOCATE (tmp_ret_val(size_of_in))
       tmp_ret_size = size_of_in
    ENDIF
    tmp_ret_val(1:size_of_in) = ret_val(1:size_of_in)

    IF (pos < 0) THEN
       ! Ge the information out of the file
       CALL getfill(MY_TARGET, status, fileorig, tmp_ret_val)
       ! Put the data into the database
       CALL getdbwl &
            & (MY_TARGET, target_sig, status, fileorig, size_of_in, tmp_ret_val)
    ELSE
       ! Get the value out of the database
       CALL getdbrl (pos, size_of_in, MY_TARGET, tmp_ret_val)
    ENDIF
    ret_val(1:size_of_in) = tmp_ret_val(1:size_of_in)

  END SUBROUTINE getinl1d

  !****************************

  SUBROUTINE getinl2d(MY_TARGET, ret_val)

    ! See getinls for details. It is the same thing but for a matrix


    CHARACTER(LEN=*) :: MY_TARGET
    LOGICAL, DIMENSION(:, :) :: ret_val

    LOGICAL, DIMENSION(:), ALLOCATABLE, SAVE :: tmp_ret_val
    INTEGER, SAVE :: tmp_ret_size = 0
    INTEGER :: target_sig, pos, size_of_in, size_1, size_2, status=0, fileorig
    INTEGER :: jl, jj, ji


    ! Compute the signature of the target

    CALL gensig(MY_TARGET, target_sig)

    ! Do we have this target in our database ?

    CALL find_sig (nb_keys, keystr, my_target, keysig, target_sig, pos)

    size_of_in = SIZE(ret_val)
    size_1 = SIZE(ret_val, 1)
    size_2 = SIZE(ret_val, 2)
    IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
       ALLOCATE (tmp_ret_val(size_of_in))
    ELSE IF (size_of_in > tmp_ret_size) THEN
       DEALLOCATE (tmp_ret_val)
       ALLOCATE (tmp_ret_val(size_of_in))
       tmp_ret_size = size_of_in
    ENDIF

    jl=0
    DO jj=1, size_2
       DO ji=1, size_1
          jl=jl+1
          tmp_ret_val(jl) = ret_val(ji, jj)
       ENDDO
    ENDDO

    IF (pos < 0) THEN
       ! Ge the information out of the file
       CALL getfill(MY_TARGET, status, fileorig, tmp_ret_val)
       ! Put the data into the database
       CALL getdbwl &
            & (MY_TARGET, target_sig, status, fileorig, size_of_in, tmp_ret_val)
    ELSE
       ! Get the value out of the database
       CALL getdbrl (pos, size_of_in, MY_TARGET, tmp_ret_val)
    ENDIF

    jl=0
    DO jj=1, size_2
       DO ji=1, size_1
          jl=jl+1
          ret_val(ji, jj) = tmp_ret_val(jl)
       ENDDO
    ENDDO

  END SUBROUTINE getinl2d

END MODULE getincom
