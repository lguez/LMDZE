MODULE getincom

  ! From getincom.f90,v 2.0 2004/04/05 14:47:48

  USE stringop, ONLY : findpos,nocomma,cmpblank,strlowercase,gensig,find_sig

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: getin, getin_dump

  INTERFACE getin
     MODULE PROCEDURE getinrs, getinr1d, getinr2d, &
          &                   getinis, getini1d, getini2d, &
          &                   getincs, getinc1d, getinc2d, &
          &                   getinls, getinl1d, getinl2d
  END INTERFACE

  INTEGER,PARAMETER :: max_files=100
  CHARACTER(LEN=100),DIMENSION(max_files),SAVE :: filelist
  INTEGER,SAVE      :: nbfiles

  INTEGER,PARAMETER :: max_lines=500
  INTEGER,SAVE :: nb_lines
  CHARACTER(LEN=100),DIMENSION(max_lines),SAVE :: fichier
  INTEGER,DIMENSION(max_lines),SAVE :: targetsiglist,fromfile,compline
  CHARACTER(LEN=30),DIMENSION(max_lines),SAVE  :: targetlist

  ! The data base of parameters

  INTEGER,PARAMETER :: memslabs=200
  INTEGER,PARAMETER :: compress_lim = 20

  INTEGER,SAVE :: nb_keys=0
  INTEGER,SAVE :: keymemsize=0
  INTEGER,SAVE,ALLOCATABLE :: keysig(:)
  CHARACTER(LEN=30),SAVE,ALLOCATABLE :: keystr(:)

  ! keystatus definition
  ! keystatus = 1 : Value comes from run.def
  ! keystatus = 2 : Default value is used
  ! keystatus = 3 : Some vector elements were taken from default

  INTEGER,SAVE,ALLOCATABLE :: keystatus(:)

  ! keytype definition
  ! keytype = 1 : Interger
  ! keytype = 2 : Real
  ! keytype = 3 : Character
  ! keytype = 4 : Logical

  INTEGER,SAVE,ALLOCATABLE :: keytype(:)

  ! Allow compression for keys (only for integer and real)
  ! keycompress < 0 : not compresses
  ! keycompress > 0 : number of repeat of the value

  INTEGER,SAVE,ALLOCATABLE :: keycompress(:)
  INTEGER,SAVE,ALLOCATABLE :: keyfromfile(:)

  INTEGER,SAVE,ALLOCATABLE :: keymemstart(:)
  INTEGER,SAVE,ALLOCATABLE :: keymemlen(:)

  INTEGER,SAVE,ALLOCATABLE :: intmem(:)
  INTEGER,SAVE             :: intmemsize=0, intmempos=0
  REAL,SAVE,ALLOCATABLE :: realmem(:)
  INTEGER,SAVE          :: realmemsize=0, realmempos=0
  CHARACTER(LEN=100),SAVE,ALLOCATABLE :: charmem(:)
  INTEGER,SAVE             :: charmemsize=0, charmempos=0
  LOGICAL,SAVE,ALLOCATABLE :: logicmem(:)
  INTEGER,SAVE             :: logicmemsize=0, logicmempos=0

CONTAINS

  !=== REAL INTERFACES

  SUBROUTINE getinrs (TARGET,ret_val)
    !---------------------------------------------------------------------
    !-  Get a real scalar. We first check if we find it
    !-  in the database and if not we get it from the run.def
  
    !-  getinr1d and getinr2d are written on the same pattern
    !---------------------------------------------------------------------
  
    CHARACTER(LEN=*) :: TARGET
    REAL :: ret_val
  
    REAL,DIMENSION(1) :: tmp_ret_val
    INTEGER :: target_sig, pos, status=0, fileorig
    !---------------------------------------------------------------------
  
    ! Compute the signature of the target
  
    CALL gensig (TARGET,target_sig)
  
    ! Do we have this target in our database ?

    ! Modification by Lionel GUEZ, April 4th, 2007
    ! "find_sig" should not be called if "keystr" and "keysig" are not
    ! allocated.
    ! Avoid this problem with a test on "nb_keys":
    if (nb_keys > 0) then
       CALL find_sig(nb_keys,keystr,target,keysig,target_sig,pos)
    else
       pos = -1
    end if
  
    tmp_ret_val(1) = ret_val
  
    IF (pos < 0) THEN
       !-- Get the information out of the file
       CALL getfilr (TARGET,status,fileorig,tmp_ret_val)
       !-- Put the data into the database
       CALL getdbwr (TARGET,target_sig,status,fileorig,1,tmp_ret_val)
    ELSE
       !-- Get the value out of the database
       CALL getdbrr (pos,1,TARGET,tmp_ret_val)
    ENDIF
    ret_val = tmp_ret_val(1)
    !---------------------
  END SUBROUTINE getinrs

  !****************************

  SUBROUTINE getinr1d (TARGET,ret_val)
    !---------------------------------------------------------------------
    !- See getinrs for details. It is the same thing but for a vector
    !---------------------------------------------------------------------
  
    CHARACTER(LEN=*) :: TARGET
    REAL,DIMENSION(:) :: ret_val
  
    REAL,DIMENSION(:),ALLOCATABLE,SAVE :: tmp_ret_val
    INTEGER,SAVE :: tmp_ret_size = 0
    INTEGER :: target_sig, pos, size_of_in, status=0, fileorig
    !---------------------------------------------------------------------
  
    ! Compute the signature of the target
  
    CALL gensig (TARGET,target_sig)
  
    ! Do we have this target in our database ?
  
    CALL find_sig (nb_keys,keystr,target,keysig,target_sig,pos)
  
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
       !-- Ge the information out of the file
       CALL getfilr (TARGET,status,fileorig,tmp_ret_val)
       !-- Put the data into the database
       CALL getdbwr &
            &   (TARGET,target_sig,status,fileorig,size_of_in,tmp_ret_val)
    ELSE
       !-- Get the value out of the database
       CALL getdbrr (pos,size_of_in,TARGET,tmp_ret_val)
    ENDIF
    ret_val(1:size_of_in) = tmp_ret_val(1:size_of_in)
    !----------------------
  END SUBROUTINE getinr1d

  !****************************

  SUBROUTINE getinr2d (TARGET,ret_val)
    !---------------------------------------------------------------------
    !- See getinrs for details. It is the same thing but for a matrix
    !---------------------------------------------------------------------
  
    CHARACTER(LEN=*) :: TARGET
    REAL,DIMENSION(:,:) :: ret_val
  
    REAL,DIMENSION(:),ALLOCATABLE,SAVE :: tmp_ret_val
    INTEGER,SAVE :: tmp_ret_size = 0
    INTEGER :: target_sig,pos,size_of_in,size_1,size_2,status=0,fileorig
    INTEGER :: jl, jj, ji
    !---------------------------------------------------------------------
  
    ! Compute the signature of the target
  
    CALL gensig (TARGET,target_sig)
  
    ! Do we have this target in our database ?
  
    CALL find_sig (nb_keys,keystr,target,keysig,target_sig,pos)
  
    size_of_in = SIZE(ret_val)
    size_1 = SIZE(ret_val,1)
    size_2 = SIZE(ret_val,2)
    IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
       ALLOCATE (tmp_ret_val(size_of_in))
    ELSE IF (size_of_in > tmp_ret_size) THEN
       DEALLOCATE (tmp_ret_val)
       ALLOCATE (tmp_ret_val(size_of_in))
       tmp_ret_size = size_of_in
    ENDIF
  
    jl=0
    DO jj=1,size_2
       DO ji=1,size_1
          jl=jl+1
          tmp_ret_val(jl) = ret_val(ji,jj)
       ENDDO
    ENDDO
  
    IF (pos < 0) THEN
       !-- Ge the information out of the file
       CALL getfilr (TARGET,status,fileorig,tmp_ret_val)
       !-- Put the data into the database
       CALL getdbwr &
            &   (TARGET,target_sig,status,fileorig,size_of_in,tmp_ret_val)
    ELSE
       !-- Get the value out of the database
       CALL getdbrr (pos,size_of_in,TARGET,tmp_ret_val)
    ENDIF
  
    jl=0
    DO jj=1,size_2
       DO ji=1,size_1
          jl=jl+1
          ret_val(ji,jj) = tmp_ret_val(jl)
       ENDDO
    ENDDO
    !----------------------
  END SUBROUTINE getinr2d

  !****************************

  SUBROUTINE getfilr (TARGET,status,fileorig,ret_val)
    !---------------------------------------------------------------------
    !- Subroutine that will extract from the file the values
    !- attributed to the keyword target
  
    !- REALS
    !- -----
  
    !- target   : in  : CHARACTER(LEN=*)  target for which we will
    !-                                    look in the file
    !- status   : out : INTEGER tells us from where we obtained the data
    !- fileorig : out : The index of the file from which the key comes
    !- ret_val  : out : REAL(nb_to_ret) values read
    !---------------------------------------------------------------------
  
    CHARACTER(LEN=*) :: TARGET
    INTEGER :: status, fileorig
    REAL,DIMENSION(:) :: ret_val
  
    INTEGER :: nb_to_ret
    INTEGER :: it, pos, len_str, epos, ppos, int_tmp, status_cnt
    CHARACTER(LEN=3)  :: cnt, tl, dl
    CHARACTER(LEN=10) :: fmt
    CHARACTER(LEN=30) :: full_target
    CHARACTER(LEN=80) :: str_READ, str_READ_lower, str_tmp
    INTEGER :: full_target_sig
    REAL :: compvalue
  
    INTEGER,SAVE :: max_len = 0
    LOGICAL,SAVE,ALLOCATABLE :: found(:)
    LOGICAL :: def_beha
    LOGICAL :: compressed = .FALSE.
    !---------------------------------------------------------------------
    nb_to_ret = SIZE(ret_val)
    CALL getin_read
  
    ! Get the variables and memory we need
  
    IF (max_len == 0) THEN
       ALLOCATE(found(nb_to_ret))
       max_len = nb_to_ret
    ENDIF
    IF (max_len < nb_to_ret) THEN
       DEALLOCATE(found)
       ALLOCATE(found(nb_to_ret))
       max_len = nb_to_ret
    ENDIF
    found(:) = .FALSE.
  
    ! See what we find in the files read
  
    DO it=1,nb_to_ret
       !---
     
       !-- First try the target as it is
       !---
       full_target = TARGET(1:len_TRIM(target))
       CALL gensig (full_target,full_target_sig)
       CALL find_sig (nb_lines,targetlist,full_target, &
            &                 targetsiglist,full_target_sig,pos)
       !---
       !-- Another try
       !---
       IF (pos < 0) THEN
          WRITE(cnt,'(I3.3)') it
          full_target = TARGET(1:len_TRIM(target))//'__'//cnt
          CALL gensig (full_target,full_target_sig)
          CALL find_sig (nb_lines,targetlist,full_target, &
               &                   targetsiglist,full_target_sig,pos)
       ENDIF
       !---
       !-- A priori we dont know from which file the target could come.
       !-- Thus by default we attribute it to the first file :
       !---
       fileorig = 1
       !--
       IF (pos > 0) THEN
          !----
          found(it) = .TRUE.
          fileorig = fromfile(pos)
          !-----
          !---- DECODE
          !-----
          str_READ = TRIM(ADJUSTL(fichier(pos)))
          str_READ_lower = str_READ
          CALL strlowercase (str_READ_lower)
          !----
          IF (    (     (INDEX(str_READ_lower,'def') == 1)     &
               &             .AND.(LEN_TRIM(str_READ_lower) == 3)   )    &
               &        .OR.(     (INDEX(str_READ_lower,'default') == 1) &
               &             .AND.(LEN_TRIM(str_READ_lower) == 7)   )   ) THEN
             def_beha = .TRUE.
          ELSE
             def_beha = .FALSE.
             len_str = LEN_TRIM(str_READ)
             epos = INDEX(str_READ,'e')
             ppos = INDEX(str_READ,'.')
             !------
             IF (epos > 0) THEN
                WRITE(tl,'(I3.3)') len_str
                WRITE(dl,'(I3.3)') epos-ppos-1
                fmt='(e'//tl//'.'//dl//')'
                READ(str_READ,fmt) ret_val(it)
             ELSE IF (ppos > 0) THEN
                WRITE(tl,'(I3.3)') len_str
                WRITE(dl,'(I3.3)') len_str-ppos
                fmt='(f'//tl//'.'//dl//')'
                READ(str_READ,fmt) ret_val(it)
             ELSE
                WRITE(tl,'(I3.3)') len_str
                fmt = '(I'//tl//')'
                READ(str_READ,fmt) int_tmp
                ret_val(it) = REAL(int_tmp)
             ENDIF
          ENDIF
          !----
          targetsiglist(pos) = -1
          !-----
          !---- Is this the value of a compressed field ?
          !-----
          IF (compline(pos) > 0) THEN
             IF (compline(pos) == nb_to_ret) THEN
                compressed = .TRUE.
                compvalue = ret_val(it)
             ELSE
                WRITE(*,*) 'WARNING from getfilr'
                WRITE(*,*) 'For key ',TRIM(TARGET), &
                     & ' we have a compressed field but which does not have the right size.'
                WRITE(*,*) 'We will try to fix that '
                compressed = .TRUE.
                compvalue = ret_val(it)
             ENDIF
          ENDIF
       ELSE
          found(it) = .FALSE.
       ENDIF
    ENDDO
    !--
    ! If this is a compressed field then we will uncompress it
    !--
    IF (compressed) THEN
       DO it=1,nb_to_ret
          IF (.NOT. found(it)) THEN
             ret_val(it) = compvalue
             found(it) = .TRUE.
          ENDIF
       ENDDO
    ENDIF
  
    ! Now we get the status for what we found
  
    IF (def_beha) THEN
       status = 2
       WRITE(*,*) 'USING DEFAULT BEHAVIOUR FOR ',TRIM(TARGET)
    ELSE
       status_cnt = 0
       DO it=1,nb_to_ret
          IF (.NOT. found(it)) THEN
             status_cnt = status_cnt+1
             IF (nb_to_ret > 1) THEN
                WRITE(str_tmp,'(a,"__",I3.3)') TRIM(TARGET),it
             ELSE
                str_tmp = TRIM(TARGET)
             ENDIF
             WRITE(*,*) 'USING DEFAULTS : ',TRIM(str_tmp),'=',ret_val(it)
          ENDIF
       ENDDO
       !---
       IF (status_cnt == 0) THEN
          status = 1
       ELSE IF (status_cnt == nb_to_ret) THEN
          status = 2
       ELSE
          status = 3
       ENDIF
    ENDIF
    !---------------------
  END SUBROUTINE getfilr

  !=== INTEGER INTERFACES

  SUBROUTINE getinis (TARGET,ret_val)
    !---------------------------------------------------------------------
    !- Get a interer scalar. We first check if we find it
    !- in the database and if not we get it from the run.def
  
    !- getini1d and getini2d are written on the same pattern
    !---------------------------------------------------------------------
  
    CHARACTER(LEN=*) :: TARGET
    INTEGER :: ret_val
  
    INTEGER,DIMENSION(1) :: tmp_ret_val
    INTEGER :: target_sig, pos, status=0, fileorig
    !---------------------------------------------------------------------
  
    ! Compute the signature of the target
  
    CALL gensig (TARGET,target_sig)
  
    ! Do we have this target in our database ?
  
    CALL find_sig (nb_keys,keystr,target,keysig,target_sig,pos)
  
    tmp_ret_val(1) = ret_val
  
    IF (pos < 0) THEN
       !-- Ge the information out of the file
       CALL getfili (TARGET,status,fileorig,tmp_ret_val)
       !-- Put the data into the database
       CALL getdbwi (TARGET,target_sig,status,fileorig,1,tmp_ret_val)
    ELSE
       !-- Get the value out of the database
       CALL getdbri (pos,1,TARGET,tmp_ret_val)
    ENDIF
    ret_val = tmp_ret_val(1)
    !---------------------
  END SUBROUTINE getinis

  !****************************

  SUBROUTINE getini1d (TARGET,ret_val)
    !---------------------------------------------------------------------
    !- See getinis for details. It is the same thing but for a vector
    !---------------------------------------------------------------------
  
    CHARACTER(LEN=*) :: TARGET
    INTEGER,DIMENSION(:) :: ret_val
  
    INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: tmp_ret_val
    INTEGER,SAVE :: tmp_ret_size = 0
    INTEGER :: target_sig, pos, size_of_in, status=0, fileorig
    !---------------------------------------------------------------------
  
    ! Compute the signature of the target
  
    CALL gensig (TARGET,target_sig)
  
    ! Do we have this target in our database ?
  
    CALL find_sig (nb_keys,keystr,target,keysig,target_sig,pos)
  
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
       !-- Ge the information out of the file
       CALL getfili (TARGET,status,fileorig,tmp_ret_val)
       !-- Put the data into the database
       CALL getdbwi &
            &   (TARGET,target_sig,status,fileorig,size_of_in,tmp_ret_val)
    ELSE
       !-- Get the value out of the database
       CALL getdbri (pos,size_of_in,TARGET,tmp_ret_val)
    ENDIF
    ret_val(1:size_of_in) = tmp_ret_val(1:size_of_in)
    !----------------------
  END SUBROUTINE getini1d

  !****************************

  SUBROUTINE getini2d (TARGET,ret_val)
    !---------------------------------------------------------------------
    !- See getinis for details. It is the same thing but for a matrix
    !---------------------------------------------------------------------
  
    CHARACTER(LEN=*) :: TARGET
    INTEGER,DIMENSION(:,:) :: ret_val
  
    INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: tmp_ret_val
    INTEGER,SAVE :: tmp_ret_size = 0
    INTEGER :: target_sig,pos,size_of_in,size_1,size_2,status=0,fileorig
    INTEGER :: jl, jj, ji
    !---------------------------------------------------------------------
  
    ! Compute the signature of the target
  
    CALL gensig (TARGET,target_sig)
  
    ! Do we have this target in our database ?
  
    CALL find_sig (nb_keys,keystr,target,keysig,target_sig,pos)
  
    size_of_in = SIZE(ret_val)
    size_1 = SIZE(ret_val,1)
    size_2 = SIZE(ret_val,2)
    IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
       ALLOCATE (tmp_ret_val(size_of_in))
    ELSE IF (size_of_in > tmp_ret_size) THEN
       DEALLOCATE (tmp_ret_val)
       ALLOCATE (tmp_ret_val(size_of_in))
       tmp_ret_size = size_of_in
    ENDIF
  
    jl=0
    DO jj=1,size_2
       DO ji=1,size_1
          jl=jl+1
          tmp_ret_val(jl) = ret_val(ji,jj)
       ENDDO
    ENDDO
  
    IF (pos < 0) THEN
       !-- Ge the information out of the file
       CALL getfili (TARGET,status,fileorig,tmp_ret_val)
       !-- Put the data into the database
       CALL getdbwi &
            &   (TARGET,target_sig,status,fileorig,size_of_in,tmp_ret_val)
    ELSE
       !-- Get the value out of the database
       CALL getdbri (pos,size_of_in,TARGET,tmp_ret_val)
    ENDIF
  
    jl=0
    DO jj=1,size_2
       DO ji=1,size_1
          jl=jl+1
          ret_val(ji,jj) = tmp_ret_val(jl)
       ENDDO
    ENDDO
    !----------------------
  END SUBROUTINE getini2d

  !****************************

  SUBROUTINE getfili (TARGET,status,fileorig,ret_val)
    !---------------------------------------------------------------------
    !- Subroutine that will extract from the file the values
    !- attributed to the keyword target
  
    !- INTEGER
    !- -------
  
    !- target   : in  : CHARACTER(LEN=*)  target for which we will
    !-                                    look in the file
    !- status   : out : INTEGER tells us from where we obtained the data
    !- fileorig : out : The index of the file from which the key comes
    !- ret_val  : out : INTEGER(nb_to_ret) values read
    !---------------------------------------------------------------------
  
    CHARACTER(LEN=*) :: TARGET
    INTEGER :: status, fileorig
    INTEGER :: ret_val(:)
  
    INTEGER :: nb_to_ret
    INTEGER :: it, pos, len_str, status_cnt
    CHARACTER(LEN=3)  :: cnt, chlen
    CHARACTER(LEN=10) ::  fmt
    CHARACTER(LEN=30) :: full_target
    CHARACTER(LEN=80) :: str_READ, str_READ_lower, str_tmp
    INTEGER :: full_target_sig
    INTEGER :: compvalue
  
    INTEGER,SAVE :: max_len = 0
    LOGICAL,SAVE,ALLOCATABLE :: found(:)
    LOGICAL :: def_beha
    LOGICAL :: compressed = .FALSE.
    !---------------------------------------------------------------------
    nb_to_ret = SIZE(ret_val)
    CALL getin_read
  
    ! Get the variables and memory we need
  
    IF (max_len == 0) THEN
       ALLOCATE(found(nb_to_ret))
       max_len = nb_to_ret
    ENDIF
    IF (max_len < nb_to_ret) THEN
       DEALLOCATE(found)
       ALLOCATE(found(nb_to_ret))
       max_len = nb_to_ret
    ENDIF
    found(:) = .FALSE.
  
    ! See what we find in the files read
  
    DO it=1,nb_to_ret
       !---
       !-- First try the target as it is
       !---
       full_target = TARGET(1:len_TRIM(target))
       CALL gensig (full_target,full_target_sig)
       CALL find_sig (nb_lines,targetlist,full_target, &
            &                 targetsiglist,full_target_sig,pos)
       !---
       !-- Another try
       !---
       IF (pos < 0) THEN
          WRITE(cnt,'(I3.3)') it
          full_target = TARGET(1:len_TRIM(target))//'__'//cnt
          CALL gensig (full_target,full_target_sig)
          CALL find_sig (nb_lines,targetlist,full_target, &
               &                   targetsiglist,full_target_sig,pos)
       ENDIF
       !---
       !-- A priori we dont know from which file the target could come.
       !-- Thus by default we attribute it to the first file :
       !---
       fileorig = 1
     
       IF (pos > 0) THEN
          !-----
          found(it) = .TRUE.
          fileorig = fromfile(pos)
          !-----
          !---- DECODE
          !----
          str_READ = TRIM(ADJUSTL(fichier(pos)))
          str_READ_lower = str_READ
          CALL strlowercase (str_READ_lower)
          !-----
          IF (    (     (INDEX(str_READ_lower,'def') == 1)     &
               &             .AND.(LEN_TRIM(str_READ_lower) == 3)   )    &
               &        .OR.(     (INDEX(str_READ_lower,'default') == 1) &
               &             .AND.(LEN_TRIM(str_READ_lower) == 7)   )   ) THEN
             def_beha = .TRUE.
          ELSE
             def_beha = .FALSE.
             len_str = LEN_TRIM(str_READ)
             WRITE(chlen,'(I3.3)') len_str
             fmt = '(I'//chlen//')'
             READ(str_READ,fmt) ret_val(it)
          ENDIF
          !-----
          targetsiglist(pos) = -1
          !-----
          !---- Is this the value of a compressed field ?
          !-----
          IF (compline(pos) > 0) THEN
             IF (compline(pos) == nb_to_ret) THEN
                compressed = .TRUE.
                compvalue = ret_val(it)
             ELSE
                WRITE(*,*) 'WARNING from getfilr'
                WRITE(*,*) 'For key ',TRIM(TARGET), &
                     & ' we have a compressed field but which does not have the right size.'
                WRITE(*,*) 'We will try to fix that '
                compressed = .TRUE.
                compvalue = ret_val(it)
             ENDIF
          ENDIF
       ELSE
          found(it) = .FALSE.
       ENDIF
    ENDDO
  
    ! If this is a compressed field then we will uncompress it
  
    IF (compressed) THEN
       DO it=1,nb_to_ret
          IF (.NOT. found(it)) THEN
             ret_val(it) = compvalue
             found(it) = .TRUE.
          ENDIF
       ENDDO
    ENDIF
  
    ! Now we get the status for what we found
  
    IF (def_beha) THEN
       status = 2
       WRITE(*,*) 'USING DEFAULT BEHAVIOUR FOR ',TRIM(TARGET)
    ELSE
       status_cnt = 0
       DO it=1,nb_to_ret
          IF (.NOT. found(it)) THEN
             status_cnt = status_cnt+1
             IF (nb_to_ret > 1) THEN
                WRITE(str_tmp,'(a,"__",I3.3)') TRIM(TARGET),it
             ELSE
                str_tmp = TRIM(TARGET)
             ENDIF
             WRITE(*,*) 'USING DEFAULTS : ',TRIM(str_tmp),'=',ret_val(it)
          ENDIF
       ENDDO
       !---
       IF (status_cnt == 0) THEN
          status = 1
       ELSE IF (status_cnt == nb_to_ret) THEN
          status = 2
       ELSE
          status = 3
       ENDIF
    ENDIF
    !---------------------
  END SUBROUTINE getfili

  !=== CHARACTER INTERFACES

  SUBROUTINE getincs (TARGET,ret_val)
    !---------------------------------------------------------------------
    !- Get a CHARACTER scalar. We first check if we find it
    !- in the database and if not we get it from the run.def
  
    !- getinc1d and getinc2d are written on the same pattern
    !---------------------------------------------------------------------
  
    CHARACTER(LEN=*) :: TARGET
    CHARACTER(LEN=*) :: ret_val
  
    CHARACTER(LEN=100),DIMENSION(1) :: tmp_ret_val
    INTEGER :: target_sig, pos, status=0, fileorig
    !---------------------------------------------------------------------
  
    ! Compute the signature of the target
  
    CALL gensig (TARGET,target_sig)
  
    ! Do we have this target in our database ?
  
    CALL find_sig (nb_keys,keystr,target,keysig,target_sig,pos)
  
    tmp_ret_val(1) = ret_val
  
    IF (pos < 0) THEN
       !-- Ge the information out of the file
       CALL getfilc (TARGET,status,fileorig,tmp_ret_val)
       !-- Put the data into the database
       CALL getdbwc (TARGET,target_sig,status,fileorig,1,tmp_ret_val)
    ELSE
       !-- Get the value out of the database
       CALL getdbrc (pos,1,TARGET,tmp_ret_val)
    ENDIF
    ret_val = tmp_ret_val(1)
    !---------------------
  END SUBROUTINE getincs

  !****************************

  SUBROUTINE getinc1d (TARGET,ret_val)
    !---------------------------------------------------------------------
    !- See getincs for details. It is the same thing but for a vector
    !---------------------------------------------------------------------
  
    CHARACTER(LEN=*) :: TARGET
    CHARACTER(LEN=*),DIMENSION(:) :: ret_val
  
    CHARACTER(LEN=100),DIMENSION(:),ALLOCATABLE,SAVE :: tmp_ret_val
    INTEGER,SAVE :: tmp_ret_size = 0
    INTEGER :: target_sig, pos, size_of_in, status=0, fileorig
    !---------------------------------------------------------------------
  
    ! Compute the signature of the target
  
    CALL gensig (TARGET,target_sig)
  
    ! Do we have this target in our database ?
  
    CALL find_sig (nb_keys,keystr,target,keysig,target_sig,pos)
  
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
       !-- Ge the information out of the file
       CALL getfilc (TARGET,status,fileorig,tmp_ret_val)
       !-- Put the data into the database
       CALL getdbwc &
            &   (TARGET,target_sig,status,fileorig,size_of_in,tmp_ret_val)
    ELSE
       !-- Get the value out of the database
       CALL getdbrc (pos,size_of_in,TARGET,tmp_ret_val)
    ENDIF
    ret_val(1:size_of_in) = tmp_ret_val(1:size_of_in)
    !----------------------
  END SUBROUTINE getinc1d

  !****************************

  SUBROUTINE getinc2d (TARGET,ret_val)
    !---------------------------------------------------------------------
    !- See getincs for details. It is the same thing but for a matrix
    !---------------------------------------------------------------------
  
    CHARACTER(LEN=*) :: TARGET
    CHARACTER(LEN=*),DIMENSION(:,:) :: ret_val
  
    CHARACTER(LEN=100),DIMENSION(:),ALLOCATABLE,SAVE :: tmp_ret_val
    INTEGER,SAVE :: tmp_ret_size = 0
    INTEGER :: target_sig,pos,size_of_in,size_1,size_2,status=0,fileorig
    INTEGER :: jl,jj,ji
    !---------------------------------------------------------------------
  
    ! Compute the signature of the target
  
    CALL gensig (TARGET,target_sig)
  
    ! Do we have this target in our database ?
  
    CALL find_sig (nb_keys,keystr,target,keysig,target_sig,pos)
  
    size_of_in = SIZE(ret_val)
    size_1 = SIZE(ret_val,1)
    size_2 = SIZE(ret_val,2)
    IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
       ALLOCATE (tmp_ret_val(size_of_in))
    ELSE IF (size_of_in > tmp_ret_size) THEN
       DEALLOCATE (tmp_ret_val)
       ALLOCATE (tmp_ret_val(size_of_in))
       tmp_ret_size = size_of_in
    ENDIF
  
    jl=0
    DO jj=1,size_2
       DO ji=1,size_1
          jl=jl+1
          tmp_ret_val(jl) = ret_val(ji,jj)
       ENDDO
    ENDDO
  
    IF (pos < 0) THEN
       !-- Ge the information out of the file
       CALL getfilc (TARGET,status,fileorig,tmp_ret_val)
       !-- Put the data into the database
       CALL getdbwc &
            &   (TARGET,target_sig,status,fileorig,size_of_in,tmp_ret_val)
    ELSE
       !-- Get the value out of the database
       CALL getdbrc (pos,size_of_in,TARGET,tmp_ret_val)
    ENDIF
  
    jl=0
    DO jj=1,size_2
       DO ji=1,size_1
          jl=jl+1
          ret_val(ji,jj) = tmp_ret_val(jl)
       ENDDO
    ENDDO
    !----------------------
  END SUBROUTINE getinc2d

  !****************************

  SUBROUTINE getfilc (TARGET,status,fileorig,ret_val)
    !---------------------------------------------------------------------
    !- Subroutine that will extract from the file the values
    !- attributed to the keyword target
  
    !- CHARACTER
    !- ---------
  
    !- target   : in  : CHARACTER(LEN=*)  target for which we will
    !-                                    look in the file
    !- status   : out : INTEGER tells us from where we obtained the data
    !- fileorig : out : The index of the file from which the key comes
    !- ret_val  : out : CHARACTER(nb_to_ret) values read
    !---------------------------------------------------------------------
  
  
    CHARACTER(LEN=*) :: TARGET
    INTEGER :: status, fileorig
    CHARACTER(LEN=*),DIMENSION(:) :: ret_val
  
    INTEGER :: nb_to_ret
    INTEGER :: it, pos, len_str, status_cnt
    CHARACTER(LEN=3)  :: cnt
    CHARACTER(LEN=30) :: full_target
    CHARACTER(LEN=80) :: str_READ, str_READ_lower, str_tmp
    INTEGER :: full_target_sig
  
    INTEGER,SAVE :: max_len = 0
    LOGICAL,DIMENSION(:),SAVE,ALLOCATABLE :: found
    LOGICAL :: def_beha
    !---------------------------------------------------------------------
    nb_to_ret = SIZE(ret_val)
    CALL getin_read
  
    ! Get the variables and memory we need
  
    IF (max_len == 0) THEN
       ALLOCATE(found(nb_to_ret))
       max_len = nb_to_ret
    ENDIF
    IF (max_len < nb_to_ret) THEN
       DEALLOCATE(found)
       ALLOCATE(found(nb_to_ret))
       max_len = nb_to_ret
    ENDIF
    found(:) = .FALSE.
  
    ! See what we find in the files read
  
    DO it=1,nb_to_ret
       !---
       !-- First try the target as it is
       full_target = TARGET(1:len_TRIM(target))
       CALL gensig (full_target,full_target_sig)
       CALL find_sig (nb_lines,targetlist,full_target, &
            &                 targetsiglist,full_target_sig,pos)
       !---
       !-- Another try
       !---
       IF (pos < 0) THEN
          WRITE(cnt,'(I3.3)') it
          full_target = TARGET(1:len_TRIM(target))//'__'//cnt
          CALL gensig (full_target,full_target_sig)
          CALL find_sig (nb_lines,targetlist,full_target, &
               &                   targetsiglist,full_target_sig,pos)
       ENDIF
       !---
       !-- A priori we dont know from which file the target could come.
       !-- Thus by default we attribute it to the first file :
       !---
       fileorig = 1
       !---
       IF (pos > 0) THEN
          !-----
          found(it) = .TRUE.
          fileorig = fromfile(pos)
          !-----
          !---- DECODE
          !-----
          str_READ = TRIM(ADJUSTL(fichier(pos)))
          str_READ_lower = str_READ
          CALL strlowercase (str_READ_lower)
          !-----
          IF (    (     (INDEX(str_READ_lower,'def') == 1)     &
               &             .AND.(LEN_TRIM(str_READ_lower) == 3)   )    &
               &        .OR.(     (INDEX(str_READ_lower,'default') == 1) &
               &             .AND.(LEN_TRIM(str_READ_lower) == 7)   )   ) THEN
             def_beha = .TRUE.
          ELSE
             def_beha = .FALSE.
             len_str = LEN_TRIM(str_READ)
             ret_val(it) = str_READ(1:len_str)
          ENDIF
          !-----
          targetsiglist(pos) = -1
          !-----
       ELSE
          found(it) = .FALSE.
       ENDIF
    ENDDO
  
    ! Now we get the status for what we found
  
    IF (def_beha) THEN
       status = 2
       WRITE(*,*) 'USING DEFAULT BEHAVIOUR FOR ',TRIM(TARGET)
    ELSE
       status_cnt = 0
       DO it=1,nb_to_ret
          IF (.NOT. found(it)) THEN
             status_cnt = status_cnt+1
             IF (nb_to_ret > 1) THEN
                WRITE(str_tmp,'(a,"__",I3.3)') TRIM(TARGET),it
             ELSE
                str_tmp = TARGET(1:len_TRIM(target))
             ENDIF
             WRITE(*,*) 'USING DEFAULTS : ',TRIM(str_tmp),'=',ret_val(it)
          ENDIF
       ENDDO
     
       IF (status_cnt == 0) THEN
          status = 1
       ELSE IF (status_cnt == nb_to_ret) THEN
          status = 2
       ELSE
          status = 3
       ENDIF
    ENDIF
    !---------------------
  END SUBROUTINE getfilc

  !=== LOGICAL INTERFACES

  SUBROUTINE getinls (TARGET,ret_val)
    !---------------------------------------------------------------------
    !- Get a logical scalar. We first check if we find it
    !- in the database and if not we get it from the run.def
  
    !- getinl1d and getinl2d are written on the same pattern
    !---------------------------------------------------------------------
  
    CHARACTER(LEN=*) :: TARGET
    LOGICAL :: ret_val
  
    LOGICAL,DIMENSION(1) :: tmp_ret_val
    INTEGER :: target_sig, pos, status=0, fileorig
    !---------------------------------------------------------------------
  
    ! Compute the signature of the target
  
    CALL gensig (TARGET,target_sig)
  
    ! Do we have this target in our database ?
  
    if (nb_keys > 0) then
       CALL find_sig(nb_keys,keystr,target,keysig,target_sig,pos)
    else
       pos = -1
    end if
  
    tmp_ret_val(1) = ret_val
  
    IF (pos < 0) THEN
       !-- Ge the information out of the file
       CALL getfill (TARGET,status,fileorig,tmp_ret_val)
       !-- Put the data into the database
       CALL getdbwl (TARGET,target_sig,status,fileorig,1,tmp_ret_val)
    ELSE
       !-- Get the value out of the database
       CALL getdbrl (pos,1,TARGET,tmp_ret_val)
    ENDIF
    ret_val = tmp_ret_val(1)
    !---------------------
  END SUBROUTINE getinls

  !****************************

  SUBROUTINE getinl1d (TARGET,ret_val)
    !---------------------------------------------------------------------
    !- See getinls for details. It is the same thing but for a vector
    !---------------------------------------------------------------------
  
    CHARACTER(LEN=*) :: TARGET
    LOGICAL,DIMENSION(:) :: ret_val
  
    LOGICAL,DIMENSION(:),ALLOCATABLE,SAVE :: tmp_ret_val
    INTEGER,SAVE :: tmp_ret_size = 0
    INTEGER :: target_sig, pos, size_of_in, status=0, fileorig
    !---------------------------------------------------------------------
  
    ! Compute the signature of the target
  
    CALL gensig (TARGET,target_sig)
  
    ! Do we have this target in our database ?
  
    CALL find_sig (nb_keys,keystr,target,keysig,target_sig,pos)
  
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
       !-- Ge the information out of the file
       CALL getfill (TARGET,status,fileorig,tmp_ret_val)
       !-- Put the data into the database
       CALL getdbwl &
            &   (TARGET,target_sig,status,fileorig,size_of_in,tmp_ret_val)
    ELSE
       !-- Get the value out of the database
       CALL getdbrl (pos,size_of_in,TARGET,tmp_ret_val)
    ENDIF
    ret_val(1:size_of_in) = tmp_ret_val(1:size_of_in)
    !----------------------
  END SUBROUTINE getinl1d

  !****************************

  SUBROUTINE getinl2d (TARGET,ret_val)
    !---------------------------------------------------------------------
    !- See getinls for details. It is the same thing but for a matrix
    !---------------------------------------------------------------------
  
    CHARACTER(LEN=*) :: TARGET
    LOGICAL,DIMENSION(:,:) :: ret_val
  
    LOGICAL,DIMENSION(:),ALLOCATABLE,SAVE :: tmp_ret_val
    INTEGER,SAVE :: tmp_ret_size = 0
    INTEGER :: target_sig,pos,size_of_in,size_1,size_2,status=0,fileorig
    INTEGER :: jl,jj,ji
    !---------------------------------------------------------------------
  
    ! Compute the signature of the target
  
    CALL gensig (TARGET,target_sig)
  
    ! Do we have this target in our database ?
  
    CALL find_sig (nb_keys,keystr,target,keysig,target_sig,pos)
  
    size_of_in = SIZE(ret_val)
    size_1 = SIZE(ret_val,1)
    size_2 = SIZE(ret_val,2)
    IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
       ALLOCATE (tmp_ret_val(size_of_in))
    ELSE IF (size_of_in > tmp_ret_size) THEN
       DEALLOCATE (tmp_ret_val)
       ALLOCATE (tmp_ret_val(size_of_in))
       tmp_ret_size = size_of_in
    ENDIF
  
    jl=0
    DO jj=1,size_2
       DO ji=1,size_1
          jl=jl+1
          tmp_ret_val(jl) = ret_val(ji,jj)
       ENDDO
    ENDDO
  
    IF (pos < 0) THEN
       !-- Ge the information out of the file
       CALL getfill (TARGET,status,fileorig,tmp_ret_val)
       !-- Put the data into the database
       CALL getdbwl &
            &   (TARGET,target_sig,status,fileorig,size_of_in,tmp_ret_val)
    ELSE
       !-- Get the value out of the database
       CALL getdbrl (pos,size_of_in,TARGET,tmp_ret_val)
    ENDIF
  
    jl=0
    DO jj=1,size_2
       DO ji=1,size_1
          jl=jl+1
          ret_val(ji,jj) = tmp_ret_val(jl)
       ENDDO
    ENDDO
    !----------------------
  END SUBROUTINE getinl2d

  !****************************

  SUBROUTINE getfill (TARGET,status,fileorig,ret_val)
    !---------------------------------------------------------------------
    !- Subroutine that will extract from the file the values
    !- attributed to the keyword target
  
    !- LOGICAL
    !- -------
  
    !- target   : in  : CHARACTER(LEN=*)  target for which we will
    !-                                    look in the file
    !- status   : out : INTEGER tells us from where we obtained the data
    !- fileorig : out : The index of the file from which the key comes
    !- ret_val  : out : LOGICAL(nb_to_ret) values read
    !---------------------------------------------------------------------
  
    CHARACTER(LEN=*) :: TARGET
    INTEGER :: status, fileorig
    LOGICAL,DIMENSION(:) :: ret_val
  
    INTEGER :: nb_to_ret
    INTEGER :: it, pos, len_str, ipos_tr, ipos_fl, status_cnt
    CHARACTER(LEN=3)  :: cnt
    CHARACTER(LEN=30) :: full_target
    CHARACTER(LEN=80) :: str_READ, str_READ_lower, str_tmp
    INTEGER :: full_target_sig
  
    INTEGER,SAVE :: max_len = 0
    LOGICAL,DIMENSION(:),SAVE,ALLOCATABLE :: found
    LOGICAL :: def_beha
    !---------------------------------------------------------------------
    nb_to_ret = SIZE(ret_val)
    CALL getin_read
  
    ! Get the variables and memory we need
  
    IF (max_len == 0) THEN
       ALLOCATE(found(nb_to_ret))
       max_len = nb_to_ret
    ENDIF
    IF (max_len < nb_to_ret) THEN
       DEALLOCATE(found)
       ALLOCATE(found(nb_to_ret))
       max_len = nb_to_ret
    ENDIF
    found(:) = .FALSE.
  
    ! See what we find in the files read
  
    DO it=1,nb_to_ret
       !---
       !-- First try the target as it is
       !---
       full_target = TARGET(1:len_TRIM(target))
       CALL gensig (full_target,full_target_sig)
       CALL find_sig (nb_lines,targetlist,full_target, &
            &                 targetsiglist,full_target_sig,pos)
       !---
       !-- Another try
       !---
       IF (pos < 0) THEN
          WRITE(cnt,'(I3.3)') it
          full_target = TARGET(1:len_TRIM(target))//'__'//cnt
          CALL gensig (full_target,full_target_sig)
          CALL find_sig (nb_lines,targetlist,full_target, &
               &                   targetsiglist,full_target_sig,pos)
       ENDIF
       !---
       !-- A priori we dont know from which file the target could come.
       !-- Thus by default we attribute it to the first file :
       !---
       fileorig = 1
       !---
       IF (pos > 0) THEN
          !-----
          found(it) = .TRUE.
          fileorig = fromfile(pos)
          !-----
          !---- DECODE
          !-----
          str_READ = TRIM(ADJUSTL(fichier(pos)))
          str_READ_lower = str_READ
          CALL strlowercase (str_READ_lower)
          !-----
          IF (    (     (INDEX(str_READ_lower,'def') == 1)     &
               &             .AND.(LEN_TRIM(str_READ_lower) == 3)   )    &
               &        .OR.(     (INDEX(str_READ_lower,'default') == 1) &
               &             .AND.(LEN_TRIM(str_READ_lower) == 7)   )   ) THEN
             def_beha = .TRUE.
          ELSE
             def_beha = .FALSE.
             len_str = LEN_TRIM(str_READ)
             ipos_tr = -1
             ipos_fl = -1
             !-------
             ipos_tr = MAX(INDEX(str_READ,'tru'),INDEX(str_READ,'TRU'), &
                  &                    INDEX(str_READ,'y'),INDEX(str_READ,'Y'))
             ipos_fl = MAX(INDEX(str_READ,'fal'),INDEX(str_READ,'FAL'), &
                  &                    INDEX(str_READ,'n'),INDEX(str_READ,'N'))
             !-------
             IF (ipos_tr > 0) THEN
                ret_val(it) = .TRUE.
             ELSE IF (ipos_fl > 0) THEN
                ret_val(it) = .FALSE.
             ELSE
                WRITE(*,*) "ERROR : getfill : TARGET ", &
                     &                   TRIM(TARGET)," is not of logical value"
                STOP 'getinl'
             ENDIF
          ENDIF
          !-----
          targetsiglist(pos) = -1
          !-----
       ELSE
        
          found(it) = .FALSE.
        
       ENDIF
     
    ENDDO
  
    ! Now we get the status for what we found
  
    IF (def_beha) THEN
       status = 2
       WRITE(*,*) 'USING DEFAULT BEHAVIOUR FOR ',TRIM(TARGET)
    ELSE
       status_cnt = 0
       DO it=1,nb_to_ret
          IF (.NOT. found(it)) THEN
             status_cnt = status_cnt+1
             IF (nb_to_ret > 1) THEN
                WRITE(str_tmp,'(a,"__",I3.3)') TRIM(TARGET),it
             ELSE
                str_tmp = TRIM(TARGET)
             ENDIF
             WRITE(*,*) 'USING DEFAULTS : ',TRIM(str_tmp),'=',ret_val(it)
          ENDIF
       ENDDO
       !---
       IF (status_cnt == 0) THEN
          status = 1
       ELSE IF (status_cnt == nb_to_ret) THEN
          status = 2
       ELSE
          status = 3
       ENDIF
    ENDIF
    !---------------------
  END SUBROUTINE getfill

  !****************************

  SUBROUTINE getin_read
    !---------------------------------------------------------------------
  
    INTEGER,SAVE :: allread=0
    INTEGER,SAVE :: current
    !---------------------------------------------------------------------
    IF (allread == 0) THEN
       !-- Allocate a first set of memory.
       CALL getin_allockeys
       CALL getin_allocmem (1,0)
       CALL getin_allocmem (2,0)
       CALL getin_allocmem (3,0)
       CALL getin_allocmem (4,0)
       !-- Start with reading the files
       nbfiles = 1
       filelist(1) = 'run.def'
       current = 1
       nb_lines = 0
       !--
       DO WHILE (current <= nbfiles)
          CALL getin_readdef (current)
          current = current+1
       ENDDO
       allread = 1
       CALL getin_checkcohe ()
    ENDIF
    !------------------------
  END SUBROUTINE getin_read

  !****************************

  SUBROUTINE getin_readdef(current)
    !---------------------------------------------------------------------
    !- This subroutine will read the files and only keep the
    !- the relevant information. The information is kept as it
    !- found in the file. The data will be analysed later.
    !---------------------------------------------------------------------
  
    INTEGER :: current
  
    CHARACTER(LEN=100) :: READ_str, NEW_str, new_key, last_key, key_str
    CHARACTER(LEN=3) :: cnt
    INTEGER :: nb_lastkey
  
    INTEGER :: eof, ptn, len_str, i, it, iund
    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    eof = 0
    ptn = 1
    nb_lastkey = 0
  
    IF (check) THEN
       WRITE(*,*) 'getin_readdef : Open file ',TRIM(filelist(current))
    ENDIF
  
    OPEN (22,file=filelist(current),ERR=9997,STATUS="OLD")
  
    DO WHILE (eof /= 1)
       !---
       CALL getin_skipafew (22,READ_str,eof,nb_lastkey)
       len_str = LEN_TRIM(READ_str)
       ptn = INDEX(READ_str,'=')
       !---
       IF (ptn > 0) THEN
          !---- Get the target
          key_str = TRIM(ADJUSTL(READ_str(1:ptn-1)))
          !---- Make sure that if a vector keyword has the right length
          iund =  INDEX(key_str,'__')
          IF (iund > 0) THEN
             SELECTCASE( len_trim(key_str)-iund )
             CASE(2)
                READ(key_str(iund+2:len_trim(key_str)),'(I1)') it
             CASE(3)
                READ(key_str(iund+2:len_trim(key_str)),'(I2)') it
             CASE(4)
                READ(key_str(iund+2:len_trim(key_str)),'(I3)') it
             CASE DEFAULT
                it = -1
             END SELECT
             IF (it > 0) THEN
                WRITE(cnt,'(I3.3)') it
                key_str = key_str(1:iund+1)//cnt
             ELSE
                WRITE(*,*) &
                     &          'getin_readdef : A very strange key has just been found'
                WRITE(*,*) 'getin_readdef : ',key_str(1:len_TRIM(key_str))
                STOP 'getin_readdef'
             ENDIF
          ENDIF
          !---- Prepare the content
          NEW_str = TRIM(ADJUSTL(READ_str(ptn+1:len_str)))
          CALL nocomma (NEW_str)
          CALL cmpblank (NEW_str)
          NEW_str  = TRIM(ADJUSTL(NEW_str))
          IF (check) THEN
             WRITE(*,*) &
                  &        '--> getin_readdef : ',TRIM(key_str),' :: ',TRIM(NEW_str)
          ENDIF
          !---- Decypher the content of NEW_str
        
          !---- This has to be a new key word, thus :
          nb_lastkey = 0
          !----
          CALL getin_decrypt (current,key_str,NEW_str,last_key,nb_lastkey)
          !----
       ELSE IF (len_str > 0) THEN
          !---- Prepare the key if we have an old one to which
          !---- we will add the line just read
          IF (nb_lastkey > 0) THEN
             iund =  INDEX(last_key,'__')
             IF (iund > 0) THEN
                !-------- We only continue a keyword, thus it is easy
                key_str = last_key(1:iund-1)
             ELSE
                IF (nb_lastkey /= 1) THEN
                   WRITE(*,*) &
                        &   'getin_readdef : An error has occured. We can not have a scalar'
                   WRITE(*,*) 'getin_readdef : keywod and a vector content'
                   STOP 'getin_readdef'
                ENDIF
                !-------- The last keyword needs to be transformed into a vector.
                targetlist(nb_lines) = &
                     &          last_key(1:MIN(len_trim(last_key),30))//'__001'
                CALL gensig (targetlist(nb_lines),targetsiglist(nb_lines))
                key_str = last_key(1:len_TRIM(last_key))
             ENDIF
          ENDIF
          !---- Prepare the content
          NEW_str = TRIM(ADJUSTL(READ_str(1:len_str)))
          CALL getin_decrypt (current,key_str,NEW_str,last_key,nb_lastkey)
       ELSE
          !---- If we have an empty line the the keyword finishes
          nb_lastkey = 0
          IF (check) THEN
             WRITE(*,*) 'getin_readdef : Have found an emtpy line '
          ENDIF
       ENDIF
    ENDDO
  
    CLOSE(22)
  
    IF (check) THEN
       OPEN (22,file='run.def.test')
       DO i=1,nb_lines
          WRITE(22,*) targetlist(i)," : ",fichier(i)
       ENDDO
       CLOSE(22)
    ENDIF
  
    RETURN
  
9997 WRITE(*,*) "getin_readdef : Could not open file ", &
         & TRIM(filelist(current))
    !---------------------------
  END SUBROUTINE getin_readdef

  !****************************

  SUBROUTINE getin_decrypt(current,key_str,NEW_str,last_key,nb_lastkey)
    !---------------------------------------------------------------------
    !- This subroutine is going to decypher the line.
    !- It essentialy checks how many items are included and
    !- it they can be attached to a key.
    !---------------------------------------------------------------------
  
    ! ARGUMENTS
  
    INTEGER :: current, nb_lastkey
    CHARACTER(LEN=*) :: key_str, NEW_str, last_key
  
    ! LOCAL
  
    INTEGER :: len_str, blk, nbve, starpos
    CHARACTER(LEN=100) :: tmp_str, new_key, mult
    CHARACTER(LEN=3)   :: cnt, chlen
    CHARACTER(LEN=10)  :: fmt
    !---------------------------------------------------------------------
    len_str = LEN_TRIM(NEW_str)
    blk = INDEX(NEW_str(1:len_str),' ')
    tmp_str = NEW_str(1:len_str)
  
    ! If the key is a new file then we take it up. Else
    ! we save the line and go on.
  
    IF (INDEX(key_str,'INCLUDEDEF') > 0) THEN
       DO WHILE (blk > 0)
          IF (nbfiles+1 > max_files) THEN
             WRITE(*,*) 'FATAL ERROR : Too many files to include'
             STOP 'getin_readdef'
          ENDIF
          !-----
          nbfiles = nbfiles+1
          filelist(nbfiles) = tmp_str(1:blk)
          !-----
          tmp_str = TRIM(ADJUSTL(tmp_str(blk+1:LEN_TRIM(tmp_str))))
          blk = INDEX(tmp_str(1:LEN_TRIM(tmp_str)),' ')
       ENDDO
       !---
       IF (nbfiles+1 > max_files) THEN
          WRITE(*,*) 'FATAL ERROR : Too many files to include'
          STOP 'getin_readdef'
       ENDIF
       !---
       nbfiles =  nbfiles+1
       filelist(nbfiles) = TRIM(ADJUSTL(tmp_str))
       !---
       last_key = 'INCLUDEDEF'
       nb_lastkey = 1
    ELSE
     
       !-- We are working on a new line of input
     
       nb_lines = nb_lines+1
       IF (nb_lines > max_lines) THEN
          WRITE(*,*) &
               &      'Too many line in the run.def files. You need to increase'
          WRITE(*,*) 'the parameter max_lines in the module getincom.'
          STOP 'getin_decrypt'
       ENDIF
     
       !-- First we solve the issue of conpressed information. Once
       !-- this is done all line can be handled in the same way.
     
       starpos = INDEX(NEW_str(1:len_str),'*')
       IF ( (starpos > 0).AND.(tmp_str(1:1) /= '"') &
            &                    .AND.(tmp_str(1:1) /= "'") ) THEN
          !-----
          IF (INDEX(key_str(1:len_TRIM(key_str)),'__') > 0) THEN
             WRITE(*,*) 'ERROR : getin_decrypt'
             WRITE(*,*) &
                  &         'We can not have a compressed field of values for in a'
             WRITE(*,*) &
                  &         'vector notation. If a target is of the type TARGET__1'
             WRITE(*,*) 'then only a scalar value is allowed'
             WRITE(*,*) 'The key at fault : ',key_str(1:len_TRIM(key_str))
             STOP 'getin_decrypt'
          ENDIF
        
          !---- Read the multiplied
        
          mult = TRIM(ADJUSTL(NEW_str(1:starpos-1)))
          !---- Construct the new string and its parameters
          NEW_str = TRIM(ADJUSTL(NEW_str(starpos+1:len_str)))
          len_str = LEN_TRIM(NEW_str)
          blk = INDEX(NEW_str(1:len_str),' ')
          IF (blk > 1) THEN
             WRITE(*,*) &
                  &       'This is a strange behavior of getin_decrypt you could report'
          ENDIF
          WRITE(chlen,'(I3.3)') LEN_TRIM(mult)
          fmt = '(I'//chlen//')'
          READ(mult,fmt) compline(nb_lines)
          !---
       ELSE
          compline(nb_lines) = -1
       ENDIF
     
       !-- If there is no space wthin the line then the target is a scalar
       !-- or the element of a properly written vector.
       !-- (ie of the type TARGET__1)
     
       IF (    (blk <= 1) &
            &      .OR.(tmp_str(1:1) == '"') &
            &      .OR.(tmp_str(1:1) == "'") ) THEN
        
          IF (nb_lastkey == 0) THEN
             !------ Save info of current keyword as a scalar
             !------ if it is not a continuation
             targetlist(nb_lines) = key_str(1:MIN(len_trim(key_str),30))
             last_key = key_str(1:MIN(len_trim(key_str),30))
             nb_lastkey = 1
          ELSE
             !------ We are continuing a vector so the keyword needs
             !------ to get the underscores
             WRITE(cnt,'(I3.3)') nb_lastkey+1
             targetlist(nb_lines) = &
                  &        key_str(1:MIN(len_trim(key_str),25))//'__'//cnt
             last_key = key_str(1:MIN(len_trim(key_str),25))//'__'//cnt
             nb_lastkey = nb_lastkey+1
          ENDIF
          !-----
          CALL gensig (targetlist(nb_lines),targetsiglist(nb_lines))
          fichier(nb_lines) = NEW_str(1:len_str)
          fromfile(nb_lines) = current
       ELSE
        
          !---- If there are blanks whithin the line then we are dealing
          !---- with a vector and we need to split it in many entries
          !---- with the TRAGET__1 notation.
          !----
          !---- Test if the targer is not already a vector target !
        
          IF (INDEX(TRIM(key_str),'__') > 0) THEN
             WRITE(*,*) 'ERROR : getin_decrypt'
             WRITE(*,*) 'We have found a mixed vector notation'
             WRITE(*,*) 'If a target is of the type TARGET__1'
             WRITE(*,*) 'then only a scalar value is allowed'
             WRITE(*,*) 'The key at fault : ',key_str(1:len_TRIM(key_str))
             STOP 'getin_decrypt'
          ENDIF
        
          nbve = nb_lastkey
          nbve = nbve+1
          WRITE(cnt,'(I3.3)') nbve
        
          DO WHILE (blk > 0)
           
             !------ Save the content of target__nbve
           
             fichier(nb_lines) = tmp_str(1:blk)
             new_key =  key_str(1:MIN(len_trim(key_str),25))//'__'//cnt
             targetlist(nb_lines) = new_key(1:MIN(len_trim(new_key),30))
             CALL gensig (targetlist(nb_lines),targetsiglist(nb_lines))
             fromfile(nb_lines) = current
           
             tmp_str = TRIM(ADJUSTL(tmp_str(blk+1:LEN_TRIM(tmp_str))))
             blk = INDEX(TRIM(tmp_str),' ')
           
             nb_lines = nb_lines+1
             IF (nb_lines > max_lines) THEN
                WRITE(*,*) &
                     &          'Too many line in the run.def files. You need to increase'
                WRITE(*,*) 'the parameter max_lines in the module getincom.'
                STOP 'getin_decrypt'
             ENDIF
             nbve = nbve+1
             WRITE(cnt,'(I3.3)') nbve
           
          ENDDO
        
          !---- Save the content of the last target
        
          fichier(nb_lines) = tmp_str(1:LEN_TRIM(tmp_str))
          new_key = key_str(1:MIN(len_trim(key_str),25))//'__'//cnt
          targetlist(nb_lines) = new_key(1:MIN(len_trim(new_key),30))
          CALL gensig (targetlist(nb_lines),targetsiglist(nb_lines))
          fromfile(nb_lines) = current
        
          last_key = key_str(1:MIN(len_trim(key_str),25))//'__'//cnt
          nb_lastkey = nbve
        
       ENDIF
     
    ENDIF
    !---------------------------
  END SUBROUTINE getin_decrypt

  !****************************

  SUBROUTINE getin_checkcohe ()
    !---------------------------------------------------------------------
    !- This subroutine checks for redundancies.
    !---------------------------------------------------------------------
  
    ! Arguments
  
  
    ! LOCAL
  
    INTEGER :: line,i,sig
    INTEGER :: found
    CHARACTER(LEN=30) :: str
    !---------------------------------------------------------------------
    DO line=1,nb_lines-1
     
       CALL find_sig &
            &    (nb_lines-line,targetlist(line+1:nb_lines),targetlist(line), &
            &     targetsiglist(line+1:nb_lines),targetsiglist(line),found)
       !---
       !-- IF we have found it we have a problem to solve.
       !---
       IF (found > 0) THEN
          WRITE(*,*) 'COUNT : ', &
               &  COUNT(ABS(targetsiglist(line+1:nb_lines)-targetsiglist(line)) < 1)
          !-----
          WRITE(*,*) &
               & 'getin_checkcohe : Found a problem on key ',targetlist(line)
          WRITE(*,*) &
               & 'getin_checkcohe : The following values were encoutered :'
          WRITE(*,*) &
               & '                ',TRIM(targetlist(line)), &
               &               targetsiglist(line),' == ',fichier(line)
          WRITE(*,*) &
               & '                ',TRIM(targetlist(line+found)), &
               &               targetsiglist(line+found),' == ',fichier(line+found)
          WRITE(*,*) &
               & 'getin_checkcohe : We will keep only the last value'
          !-----
          targetsiglist(line) = 1
       ENDIF
    ENDDO
  
  END SUBROUTINE getin_checkcohe

  !****************************

  SUBROUTINE getin_skipafew (unit,out_string,eof,nb_lastkey)
    !---------------------------------------------------------------------
  
    INTEGER :: unit, eof, nb_lastkey
    CHARACTER(LEN=100) :: dummy
    CHARACTER(LEN=100) :: out_string
    CHARACTER(LEN=1) :: first
    !---------------------------------------------------------------------
    first="#"
    eof = 0
    out_string = "    "
  
    DO WHILE (first == "#")
       READ (unit,'(a100)',ERR=9998,END=7778) dummy
       dummy = TRIM(ADJUSTL(dummy))
       first=dummy(1:1)
       IF (first == "#") THEN
          nb_lastkey = 0
       ENDIF
    ENDDO
    out_string=dummy
  
    RETURN
  
9998 WRITE(*,*) " GETIN_SKIPAFEW : Error while reading file "
    STOP 'getin_skipafew'
  
7778 eof = 1
    !----------------------------
  END SUBROUTINE getin_skipafew

  !=== INTEGER database INTERFACE

  SUBROUTINE getdbwi &
       &  (target,target_sig,status,fileorig,size_of_in,tmp_ret_val)
    !---------------------------------------------------------------------
    !- Write the INTEGER data into the data base
    !---------------------------------------------------------------------
  
    CHARACTER(LEN=*) :: target
    INTEGER :: target_sig, status, fileorig, size_of_in
    INTEGER,DIMENSION(:) :: tmp_ret_val
    !---------------------------------------------------------------------
  
    ! First check if we have sufficiant space for the new key
  
    IF (nb_keys+1 > keymemsize) THEN
       CALL getin_allockeys ()
    ENDIF
  
    ! Fill out the items of the data base
  
    nb_keys = nb_keys+1
    keysig(nb_keys) = target_sig
    keystr(nb_keys) = target(1:MIN(len_trim(target),30))
    keystatus(nb_keys) = status
    keytype(nb_keys) = 1
    keyfromfile(nb_keys) = fileorig
  
    ! Can we compress the data base entry ?
  
    IF (     (MINVAL(tmp_ret_val) == MAXVAL(tmp_ret_val)) &
         &    .AND.(size_of_in > compress_lim)) THEN
       keymemstart(nb_keys) = intmempos+1
       keycompress(nb_keys) = size_of_in
       keymemlen(nb_keys) = 1
    ELSE
       keymemstart(nb_keys) = intmempos+1
       keycompress(nb_keys) = -1
       keymemlen(nb_keys) = size_of_in
    ENDIF
  
    ! Before writing the actual size lets see if we have the space
  
    IF (keymemstart(nb_keys)+keymemlen(nb_keys) > intmemsize) THEN
       CALL getin_allocmem (1,keymemlen(nb_keys))
    ENDIF
  
    intmem(keymemstart(nb_keys): &
         &       keymemstart(nb_keys)+keymemlen(nb_keys)-1) = &
         &  tmp_ret_val(1:keymemlen(nb_keys))
    intmempos = keymemstart(nb_keys)+keymemlen(nb_keys)-1
    !---------------------
  END SUBROUTINE getdbwi

  !****************************

  SUBROUTINE getdbri (pos,size_of_in,target,tmp_ret_val)
    !---------------------------------------------------------------------
    !- Read the required variables in the database for INTEGERS
    !---------------------------------------------------------------------
  
    INTEGER :: pos, size_of_in
    CHARACTER(LEN=*) :: target
    INTEGER,DIMENSION(:) :: tmp_ret_val
    !---------------------------------------------------------------------
    IF (keytype(pos) /= 1) THEN
       WRITE(*,*) 'FATAL ERROR : Wrong data type for keyword ',target
       STOP 'getdbri'
    ENDIF
  
    IF (keycompress(pos) > 0) THEN
       IF ( keycompress(pos) /= size_of_in .OR. keymemlen(pos) /= 1 ) THEN
          WRITE(*,*) &
               &      'FATAL ERROR : Wrong compression length for keyword ',target
          STOP 'getdbri'
       ELSE
          tmp_ret_val(1:size_of_in) = intmem(keymemstart(pos))
       ENDIF
    ELSE
       IF (keymemlen(pos) /= size_of_in) THEN
          WRITE(*,*) 'FATAL ERROR : Wrong array length for keyword ',target
          STOP 'getdbri'
       ELSE
          tmp_ret_val(1:size_of_in) = &
               &      intmem(keymemstart(pos):keymemstart(pos)+keymemlen(pos)-1)
       ENDIF
    ENDIF
    !---------------------
  END SUBROUTINE getdbri

  !=== REAL database INTERFACE

  SUBROUTINE getdbwr &
       &  (target,target_sig,status,fileorig,size_of_in,tmp_ret_val)
    !---------------------------------------------------------------------
    !- Write the REAL data into the data base
    !---------------------------------------------------------------------
  
    CHARACTER(LEN=*) :: target
    INTEGER :: target_sig, status, fileorig, size_of_in
    REAL,DIMENSION(:) :: tmp_ret_val
    !---------------------------------------------------------------------
  
    ! First check if we have sufficiant space for the new key
  
    IF (nb_keys+1 > keymemsize) THEN
       CALL getin_allockeys ()
    ENDIF
  
    ! Fill out the items of the data base
  
    nb_keys = nb_keys+1
    keysig(nb_keys) = target_sig
    keystr(nb_keys) = target(1:MIN(len_trim(target),30))
    keystatus(nb_keys) = status
    keytype(nb_keys) = 2
    keyfromfile(nb_keys) = fileorig
  
    ! Can we compress the data base entry ?
  
    IF (     (MINVAL(tmp_ret_val) == MAXVAL(tmp_ret_val)) &
         &     .AND.(size_of_in > compress_lim)) THEN
       keymemstart(nb_keys) = realmempos+1
       keycompress(nb_keys) = size_of_in
       keymemlen(nb_keys) = 1
    ELSE
       keymemstart(nb_keys) = realmempos+1
       keycompress(nb_keys) = -1
       keymemlen(nb_keys) = size_of_in
    ENDIF
  
    ! Before writing the actual size lets see if we have the space
  
    IF (keymemstart(nb_keys)+keymemlen(nb_keys) > realmemsize) THEN
       CALL getin_allocmem (2,keymemlen(nb_keys))
    ENDIF
  
    realmem(keymemstart(nb_keys): &
         &        keymemstart(nb_keys)+keymemlen(nb_keys)-1) = &
         &  tmp_ret_val(1:keymemlen(nb_keys))
    realmempos = keymemstart(nb_keys)+keymemlen(nb_keys)-1
    !---------------------
  END SUBROUTINE getdbwr

  !****************************

  SUBROUTINE getdbrr (pos,size_of_in,target,tmp_ret_val)
    !---------------------------------------------------------------------
    !- Read the required variables in the database for REALS
    !---------------------------------------------------------------------
  
    INTEGER :: pos, size_of_in
    CHARACTER(LEN=*) :: target
    REAL,DIMENSION(:) :: tmp_ret_val
    !---------------------------------------------------------------------
    IF (keytype(pos) /= 2) THEN
       WRITE(*,*) 'FATAL ERROR : Wrong data type for keyword ',target
       STOP 'getdbrr'
    ENDIF
  
    IF (keycompress(pos) > 0) THEN
       IF (    (keycompress(pos) /= size_of_in) &
            &      .OR.(keymemlen(pos) /= 1) ) THEN
          WRITE(*,*) &
               &      'FATAL ERROR : Wrong compression length for keyword ',target
          STOP 'getdbrr'
       ELSE
          tmp_ret_val(1:size_of_in) = realmem(keymemstart(pos))
       ENDIF
    ELSE
       IF (keymemlen(pos) /= size_of_in) THEN
          WRITE(*,*) 'FATAL ERROR : Wrong array length for keyword ',target
          STOP 'getdbrr'
       ELSE
          tmp_ret_val(1:size_of_in) = &
               &      realmem(keymemstart(pos):keymemstart(pos)+keymemlen(pos)-1)
       ENDIF
    ENDIF
    !---------------------
  END SUBROUTINE getdbrr

  !=== CHARACTER database INTERFACE

  SUBROUTINE getdbwc &
       &  (target,target_sig,status,fileorig,size_of_in,tmp_ret_val)
    !---------------------------------------------------------------------
    !- Write the CHARACTER data into the data base
    !---------------------------------------------------------------------
  
    CHARACTER(LEN=*) :: target
    INTEGER :: target_sig,status,fileorig,size_of_in
    CHARACTER(LEN=*),DIMENSION(:) :: tmp_ret_val
    !---------------------------------------------------------------------
  
    ! First check if we have sufficiant space for the new key
  
    IF (nb_keys+1 > keymemsize) THEN
       CALL getin_allockeys ()
    ENDIF
  
    ! Fill out the items of the data base
  
    nb_keys = nb_keys+1
    keysig(nb_keys) = target_sig
    keystr(nb_keys) = target(1:MIN(len_trim(target),30))
    keystatus(nb_keys) = status
    keytype(nb_keys) = 3
    keyfromfile(nb_keys) = fileorig
    keymemstart(nb_keys) = charmempos+1
    keymemlen(nb_keys) = size_of_in
  
    ! Before writing the actual size lets see if we have the space
  
    IF (keymemstart(nb_keys)+keymemlen(nb_keys) > realmemsize) THEN
       CALL getin_allocmem (3,keymemlen(nb_keys))
    ENDIF
  
    charmem(keymemstart(nb_keys): &
         &        keymemstart(nb_keys)+keymemlen(nb_keys)-1) = &
         &  tmp_ret_val(1:keymemlen(nb_keys))
    charmempos = keymemstart(nb_keys)+keymemlen(nb_keys)-1
    !---------------------
  END SUBROUTINE getdbwc

  !****************************

  SUBROUTINE getdbrc(pos,size_of_in,target,tmp_ret_val)
    !---------------------------------------------------------------------
    !- Read the required variables in the database for CHARACTER
    !---------------------------------------------------------------------
  
    INTEGER :: pos, size_of_in
    CHARACTER(LEN=*) :: target
    CHARACTER(LEN=*),DIMENSION(:) :: tmp_ret_val
    !---------------------------------------------------------------------
    IF (keytype(pos) /= 3) THEN
       WRITE(*,*) 'FATAL ERROR : Wrong data type for keyword ',target
       STOP 'getdbrc'
    ENDIF
  
    IF (keymemlen(pos) /= size_of_in) THEN
       WRITE(*,*) 'FATAL ERROR : Wrong array length for keyword ',target
       STOP 'getdbrc'
    ELSE
       tmp_ret_val(1:size_of_in) = &
            &    charmem(keymemstart(pos):keymemstart(pos)+keymemlen(pos)-1)
    ENDIF
    !---------------------
  END SUBROUTINE getdbrc

  !=== LOGICAL database INTERFACE

  SUBROUTINE getdbwl &
       &  (target,target_sig,status,fileorig,size_of_in,tmp_ret_val)
    !---------------------------------------------------------------------
    !- Write the LOGICAL data into the data base
    !---------------------------------------------------------------------
  
    CHARACTER(LEN=*) :: target
    INTEGER :: target_sig, status, fileorig, size_of_in
    LOGICAL,DIMENSION(:) :: tmp_ret_val
    !---------------------------------------------------------------------
  
    ! First check if we have sufficiant space for the new key
  
    IF (nb_keys+1 > keymemsize) THEN
       CALL getin_allockeys ()
    ENDIF
  
    ! Fill out the items of the data base
  
    nb_keys = nb_keys+1
    keysig(nb_keys) = target_sig
    keystr(nb_keys) = target(1:MIN(len_trim(target),30))
    keystatus(nb_keys) = status
    keytype(nb_keys) = 4
    keyfromfile(nb_keys) = fileorig
    keymemstart(nb_keys) = logicmempos+1
    keymemlen(nb_keys) = size_of_in
  
    ! Before writing the actual size lets see if we have the space
  
    IF (keymemstart(nb_keys)+keymemlen(nb_keys) > logicmemsize) THEN
       CALL getin_allocmem (4,keymemlen(nb_keys))
    ENDIF
  
    logicmem(keymemstart(nb_keys): &
         &         keymemstart(nb_keys)+keymemlen(nb_keys)-1) = &
         &  tmp_ret_val(1:keymemlen(nb_keys))
    logicmempos = keymemstart(nb_keys)+keymemlen(nb_keys)-1
    !---------------------
  END SUBROUTINE getdbwl

  !****************************

  SUBROUTINE getdbrl(pos,size_of_in,target,tmp_ret_val)
    !---------------------------------------------------------------------
    !- Read the required variables in the database for LOGICALS
    !---------------------------------------------------------------------
  
    INTEGER :: pos, size_of_in
    CHARACTER(LEN=*) :: target
    LOGICAL,DIMENSION(:) :: tmp_ret_val
    !---------------------------------------------------------------------
    IF (keytype(pos) /= 4) THEN
       WRITE(*,*) 'FATAL ERROR : Wrong data type for keyword ',target
       STOP 'getdbrl'
    ENDIF
  
    IF (keymemlen(pos) /= size_of_in) THEN
       WRITE(*,*) 'FATAL ERROR : Wrong array length for keyword ',target
       STOP 'getdbrl'
    ELSE
       tmp_ret_val(1:size_of_in) = &
            &    logicmem(keymemstart(pos):keymemstart(pos)+keymemlen(pos)-1)
    ENDIF
    !---------------------
  END SUBROUTINE getdbrl

  !****************************

  SUBROUTINE getin_allockeys

    INTEGER,ALLOCATABLE :: tmp_int(:)
    CHARACTER(LEN=100),ALLOCATABLE :: tmp_str(:)

    !---------------------------------------------------------------------

    !!print *, "Call sequence information: getin_allockeys"
    ! Either nothing exists in these arrays and it is easy to do

    IF (keymemsize == 0) THEN
       ALLOCATE(keysig(memslabs))
       ALLOCATE(keystr(memslabs))
       ALLOCATE(keystatus(memslabs))
       ALLOCATE(keytype(memslabs))
       ALLOCATE(keycompress(memslabs))
       ALLOCATE(keyfromfile(memslabs))
       ALLOCATE(keymemstart(memslabs))
       ALLOCATE(keymemlen(memslabs))
       nb_keys = 0
       keymemsize = memslabs
       keycompress(:) = -1
    ELSE
       !-- There is something already in the memory,
       !-- we need to transfer and reallocate.
       ALLOCATE(tmp_str(keymemsize))

       ALLOCATE(tmp_int(keymemsize))
       tmp_int(1:keymemsize) = keysig(1:keymemsize)

       DEALLOCATE(keysig)
       ALLOCATE(keysig(keymemsize+memslabs))
       keysig(1:keymemsize) = tmp_int(1:keymemsize)

       tmp_str(1:keymemsize) = keystr(1:keymemsize)
       DEALLOCATE(keystr)
       ALLOCATE(keystr(keymemsize+memslabs))
       keystr(1:keymemsize) = tmp_str(1:keymemsize)

       tmp_int(1:keymemsize) = keystatus(1:keymemsize)
       DEALLOCATE(keystatus)
       ALLOCATE(keystatus(keymemsize+memslabs))
       keystatus(1:keymemsize) = tmp_int(1:keymemsize)

       tmp_int(1:keymemsize) = keytype(1:keymemsize)
       DEALLOCATE(keytype)
       ALLOCATE(keytype(keymemsize+memslabs))
       keytype(1:keymemsize) = tmp_int(1:keymemsize)

       tmp_int(1:keymemsize) = keycompress(1:keymemsize)
       DEALLOCATE(keycompress)
       ALLOCATE(keycompress(keymemsize+memslabs))
       keycompress(:) = -1
       keycompress(1:keymemsize) = tmp_int(1:keymemsize)

       tmp_int(1:keymemsize) = keyfromfile(1:keymemsize)
       DEALLOCATE(keyfromfile)
       ALLOCATE(keyfromfile(keymemsize+memslabs))
       keyfromfile(1:keymemsize) = tmp_int(1:keymemsize)

       tmp_int(1:keymemsize) = keymemstart(1:keymemsize)
       DEALLOCATE(keymemstart)
       ALLOCATE(keymemstart(keymemsize+memslabs))
       keymemstart(1:keymemsize) = tmp_int(1:keymemsize)

       tmp_int(1:keymemsize) = keymemlen(1:keymemsize)
       DEALLOCATE(keymemlen)
       ALLOCATE(keymemlen(keymemsize+memslabs))
       keymemlen(1:keymemsize) = tmp_int(1:keymemsize)

       keymemsize = keymemsize+memslabs

       DEALLOCATE(tmp_int)
       DEALLOCATE(tmp_str)
    ENDIF

  END SUBROUTINE getin_allockeys

  !****************************

  SUBROUTINE getin_allocmem (type,len_wanted)
    !---------------------------------------------------------------------
    !- Allocate the memory of the data base for all 4 types of memory
  
    !- 1 = INTEGER
    !- 2 = REAL
    !- 3 = CHAR
    !- 4 = LOGICAL
    !---------------------------------------------------------------------
  
    INTEGER :: type, len_wanted
  
    INTEGER,ALLOCATABLE :: tmp_int(:)
    CHARACTER(LEN=100),ALLOCATABLE :: tmp_char(:)
    REAL,ALLOCATABLE :: tmp_real(:)
    LOGICAL,ALLOCATABLE :: tmp_logic(:)
    INTEGER :: ier
    !---------------------------------------------------------------------
    SELECT CASE (type)
    CASE(1)
       IF (intmemsize == 0) THEN
          ALLOCATE(intmem(memslabs),stat=ier)
          IF (ier /= 0) THEN
             WRITE(*,*) &
                  &    'getin_allocmem : Unable to allocate db-memory intmem to ', &
                  &    memslabs
             STOP
          ENDIF
          intmemsize=memslabs
       ELSE
          ALLOCATE(tmp_int(intmemsize),stat=ier)
          IF (ier /= 0) THEN
             WRITE(*,*) &
                  &    'getin_allocmem : Unable to allocate tmp_int to ', &
                  &    intmemsize
             STOP
          ENDIF
          tmp_int(1:intmemsize) = intmem(1:intmemsize)
          DEALLOCATE(intmem)
          ALLOCATE(intmem(intmemsize+MAX(memslabs,len_wanted)),stat=ier)
          IF (ier /= 0) THEN
             WRITE(*,*) &
                  &    'getin_allocmem : Unable to re-allocate db-memory intmem to ', &
                  &    intmemsize+MAX(memslabs,len_wanted)
             STOP
          ENDIF
          intmem(1:intmemsize) = tmp_int(1:intmemsize)
          intmemsize = intmemsize+MAX(memslabs,len_wanted)
          DEALLOCATE(tmp_int)
       ENDIF
    CASE(2)
       IF (realmemsize == 0) THEN
          ALLOCATE(realmem(memslabs),stat=ier)
          IF (ier /= 0) THEN
             WRITE(*,*) &
                  &    'getin_allocmem : Unable to allocate db-memory realmem to ', &
                  &    memslabs
             STOP
          ENDIF
          realmemsize =  memslabs
       ELSE
          ALLOCATE(tmp_real(realmemsize),stat=ier)
          IF (ier /= 0) THEN
             WRITE(*,*) &
                  &    'getin_allocmem : Unable to allocate tmp_real to ', &
                  &    realmemsize
             STOP
          ENDIF
          tmp_real(1:realmemsize) = realmem(1:realmemsize)
          DEALLOCATE(realmem)
          ALLOCATE(realmem(realmemsize+MAX(memslabs,len_wanted)),stat=ier)
          IF (ier /= 0) THEN
             WRITE(*,*) &
                  &    'getin_allocmem : Unable to re-allocate db-memory realmem to ', &
                  &    realmemsize+MAX(memslabs,len_wanted)
             STOP
          ENDIF
          realmem(1:realmemsize) = tmp_real(1:realmemsize)
          realmemsize = realmemsize+MAX(memslabs,len_wanted)
          DEALLOCATE(tmp_real)
       ENDIF
    CASE(3)
       IF (charmemsize == 0) THEN
          ALLOCATE(charmem(memslabs),stat=ier)
          IF (ier /= 0) THEN
             WRITE(*,*) &
                  &    'getin_allocmem : Unable to allocate db-memory charmem to ', &
                  &    memslabs
             STOP
          ENDIF
          charmemsize = memslabs
       ELSE
          ALLOCATE(tmp_char(charmemsize),stat=ier)
          IF (ier /= 0) THEN
             WRITE(*,*) &
                  &    'getin_allocmem : Unable to allocate tmp_char to ', &
                  &    charmemsize
             STOP
          ENDIF
          tmp_char(1:charmemsize) = charmem(1:charmemsize)
          DEALLOCATE(charmem)
          ALLOCATE(charmem(charmemsize+MAX(memslabs,len_wanted)),stat=ier)
          IF (ier /= 0) THEN
             WRITE(*,*) &
                  &    'getin_allocmem : Unable to re-allocate db-memory charmem to ', &
                  &    charmemsize+MAX(memslabs,len_wanted)
             STOP
          ENDIF
          charmem(1:charmemsize) = tmp_char(1:charmemsize)
          charmemsize = charmemsize+MAX(memslabs,len_wanted)
          DEALLOCATE(tmp_char)
       ENDIF
    CASE(4)
       IF (logicmemsize == 0) THEN
          ALLOCATE(logicmem(memslabs),stat=ier)
          IF (ier /= 0) THEN
             WRITE(*,*) &
                  &    'getin_allocmem : Unable to allocate db-memory logicmem to ', &
                  &    memslabs
             STOP
          ENDIF
          logicmemsize = memslabs
       ELSE
          ALLOCATE(tmp_logic(logicmemsize),stat=ier)
          IF (ier /= 0) THEN
             WRITE(*,*) &
                  &    'getin_allocmem : Unable to allocate tmp_logic to ', &
                  &    logicmemsize
             STOP
          ENDIF
          tmp_logic(1:logicmemsize) = logicmem(1:logicmemsize)
          DEALLOCATE(logicmem)
          ALLOCATE(logicmem(logicmemsize+MAX(memslabs,len_wanted)),stat=ier)
          IF (ier /= 0) THEN
             WRITE(*,*) &
                  &    'getin_allocmem : Unable to re-allocate db-memory logicmem to ', &
                  &    logicmemsize+MAX(memslabs,len_wanted)
             STOP
          ENDIF
          logicmem(1:logicmemsize) = tmp_logic(1:logicmemsize)
          logicmemsize = logicmemsize+MAX(memslabs,len_wanted)
          DEALLOCATE(tmp_logic)
       ENDIF
    CASE DEFAULT
       WRITE(*,*) 'getin_allocmem : Unknown type : ',type
       STOP
    END SELECT
    !----------------------------
  END SUBROUTINE getin_allocmem

  !****************************

  SUBROUTINE getin_dump (fileprefix)
    !---------------------------------------------------------------------
    !- This subroutine will dump the content of the database into  file
    !- which has the same format as the run.def. The idea is that the user
    !- can see which parameters were used and re-use the file for another
    !- run.
  
    !- The argument file allows the user to change the name of the file
    !- in which the data will be archived
    !---------------------------------------------------------------------
  
    CHARACTER(*),OPTIONAL :: fileprefix
  
    CHARACTER(LEN=80) :: usedfileprefix = "used"
    INTEGER :: ikey,if,iff,iv
    CHARACTER(LEN=3) :: tmp3
    CHARACTER(LEN=100) :: tmp_str, used_filename
    LOGICAL :: check = .FALSE.
    !---------------------------------------------------------------------
    IF (PRESENT(fileprefix)) THEN
       usedfileprefix = fileprefix(1:MIN(len_trim(fileprefix),80))
    ENDIF
  
    DO if=1,nbfiles
       !---
       used_filename = TRIM(usedfileprefix)//'_'//TRIM(filelist(if))
       IF (check) THEN
          WRITE(*,*) &
               &      'GETIN_DUMP : opens file : ',TRIM(used_filename),' if = ',if
          WRITE(*,*) 'GETIN_DUMP : NUMBER OF KEYS : ',nb_keys
       ENDIF
       OPEN(unit=76,file=used_filename)
     
       !-- If this is the first file we need to add the list
       !-- of file which belong to it
     
       IF ( (if == 1) .AND. (nbfiles > 1) ) THEN
          WRITE(76,*) '# '
          WRITE(76,*) '# This file is linked to the following files :'
          WRITE(76,*) '# '
          DO iff=2,nbfiles
             WRITE(76,*) 'INCLUDEDEF = ',TRIM(filelist(iff))
          ENDDO
          WRITE(76,*) '# '
       ENDIF
       !---
       DO ikey=1,nb_keys
        
          !---- Is this key form this file ?
        
          IF (keyfromfile(ikey) == if) THEN
           
             !---- Write some comments
           
             WRITE(76,*) '#'
             SELECT CASE (keystatus(ikey))
             CASE(1)
                WRITE(76,*) '# Values of ', &
                     &          TRIM(keystr(ikey)),' comes from the run.def.'
             CASE(2)
                WRITE(76,*) '# Values of ', &
                     &          TRIM(keystr(ikey)),' are all defaults.'
             CASE(3)
                WRITE(76,*) '# Values of ', &
                     &          TRIM(keystr(ikey)),' are a mix of run.def and defaults.'
             CASE DEFAULT
                WRITE(76,*) '# Dont know from where the value of ', &
                     &          TRIM(keystr(ikey)),' comes.'
             END SELECT
             WRITE(76,*) '#'
           
             !---- Write the values
           
             SELECT CASE (keytype(ikey))
              
             CASE(1)
                IF (keymemlen(ikey) == 1) THEN
                   IF (keycompress(ikey) < 0) THEN
                      WRITE(76,*) &
                           &              TRIM(keystr(ikey)),' = ',intmem(keymemstart(ikey))
                   ELSE
                      WRITE(76,*) &
                           &              TRIM(keystr(ikey)),' = ',keycompress(ikey), &
                           &              ' * ',intmem(keymemstart(ikey))
                   ENDIF
                ELSE
                   DO iv=0,keymemlen(ikey)-1
                      WRITE(tmp3,'(I3.3)') iv+1
                      WRITE(76,*) &
                           &              TRIM(keystr(ikey)),'__',tmp3, &
                           &              ' = ',intmem(keymemstart(ikey)+iv)
                   ENDDO
                ENDIF
              
             CASE(2)
                IF (keymemlen(ikey) == 1) THEN
                   IF (keycompress(ikey) < 0) THEN
                      WRITE(76,*) &
                           &              TRIM(keystr(ikey)),' = ',realmem(keymemstart(ikey))
                   ELSE
                      WRITE(76,*) &
                           &              TRIM(keystr(ikey)),' = ',keycompress(ikey),&
                           & ' * ',realmem(keymemstart(ikey))
                   ENDIF
                ELSE
                   DO iv=0,keymemlen(ikey)-1
                      WRITE(tmp3,'(I3.3)') iv+1
                      WRITE(76,*) &
                           &              TRIM(keystr(ikey)),'__',tmp3, &
                           &              ' = ',realmem(keymemstart(ikey)+iv)
                   ENDDO
                ENDIF
             CASE(3)
                IF (keymemlen(ikey) == 1) THEN
                   tmp_str = charmem(keymemstart(ikey))
                   WRITE(76,*) TRIM(keystr(ikey)),' = ',TRIM(tmp_str)
                ELSE
                   DO iv=0,keymemlen(ikey)-1
                      WRITE(tmp3,'(I3.3)') iv+1
                      tmp_str = charmem(keymemstart(ikey)+iv)
                      WRITE(76,*) &
                           &              TRIM(keystr(ikey)),'__',tmp3,' = ',TRIM(tmp_str)
                   ENDDO
                ENDIF
             CASE(4)
                IF (keymemlen(ikey) == 1) THEN
                   IF (logicmem(keymemstart(ikey))) THEN
                      WRITE(76,*) TRIM(keystr(ikey)),' = TRUE '
                   ELSE
                      WRITE(76,*) TRIM(keystr(ikey)),' = FALSE '
                   ENDIF
                ELSE
                   DO iv=0,keymemlen(ikey)-1
                      WRITE(tmp3,'(I3.3)') iv+1
                      IF (logicmem(keymemstart(ikey)+iv)) THEN
                         WRITE(76,*) TRIM(keystr(ikey)),'__',tmp3,' = TRUE '
                      ELSE
                         WRITE(76,*) TRIM(keystr(ikey)),'__',tmp3,' = FALSE '
                      ENDIF
                   ENDDO
                ENDIF
              
             CASE DEFAULT
                WRITE(*,*) &
                     &          'FATAL ERROR : Unknown type for variable ', &
                     &          TRIM(keystr(ikey))
                STOP 'getin_dump'
             END SELECT
          ENDIF
       ENDDO
     
       CLOSE(unit=76)
     
    ENDDO
    !------------------------
  END SUBROUTINE getin_dump

END MODULE getincom
