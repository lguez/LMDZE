MODULE histcom_var

  implicit none

  ! Fixed parameter
  INTEGER, PARAMETER:: nb_files_max=20, nb_var_max=400, &
       nb_hax_max=5, nb_zax_max=10, nbopp_max=10
  REAL, PARAMETER:: missing_val = 1e20

  INTEGER:: bufftmp_max(nb_files_max) = 1

  ! Time variables
  INTEGER, SAVE:: itau0(nb_files_max)=0
  REAL, DIMENSION(nb_files_max), SAVE::date0, deltat

  ! Counter of elements
  INTEGER, SAVE:: nb_files=0
  INTEGER, DIMENSION(nb_files_max), SAVE:: nb_var=0, nb_tax=0

  ! NETCDF IDs for files and axes
  INTEGER, DIMENSION(nb_files_max), SAVE:: ncdf_ids, xid, yid, tid
  CHARACTER(LEN=500):: assc_file = ''

  ! General definitions in the NETCDF file
  INTEGER, DIMENSION(nb_files_max, 2), SAVE:: full_size=0, slab_ori, slab_sz

  ! The horizontal axes
  INTEGER, SAVE:: nb_hax(nb_files_max)=0
  CHARACTER(LEN=25), SAVE:: hax_name(nb_files_max, nb_hax_max, 2)

  ! The vertical axes
  INTEGER, SAVE:: nb_zax(nb_files_max)=0
  INTEGER, DIMENSION(nb_files_max, nb_zax_max), SAVE:: &
       zax_size, zax_ids, zax_name_length
  CHARACTER(LEN=20), SAVE:: zax_name(nb_files_max, nb_zax_max)

  ! Informations on each variable
  INTEGER, DIMENSION(nb_files_max, nb_var_max), SAVE:: &
       name_length, nbopp
  CHARACTER(LEN=20), DIMENSION(nb_files_max, nb_var_max), SAVE:: &
       name, unit_name
  CHARACTER(LEN=80), DIMENSION(nb_files_max, nb_var_max), SAVE:: &
       title, fullop
  CHARACTER(LEN=7), SAVE:: topp(nb_files_max, nb_var_max)
  CHARACTER(LEN=7), SAVE:: sopps(nb_files_max, nb_var_max, nbopp_max)
  REAL, SAVE:: scal(nb_files_max, nb_var_max, nbopp_max)
  ! Sizes of the associated grid and zommed area
  INTEGER, DIMENSION(nb_files_max, nb_var_max, 3), SAVE:: &
       scsize, zorig, zsize
  ! Sizes for the data as it goes through the various math operations
  INTEGER, SAVE:: datasz_in(nb_files_max, nb_var_max, 3) = -1
  INTEGER, SAVE:: datasz_max(nb_files_max, nb_var_max) = -1

  INTEGER, DIMENSION(nb_files_max, nb_var_max), SAVE:: var_haxid, var_zaxid, &
       var_axid, ncvar_ids

  REAL, SAVE:: minmax(nb_files_max, nb_var_max, 2)

  REAL, DIMENSION(nb_files_max, nb_var_max), SAVE:: &
       freq_opp, freq_wrt
  INTEGER, DIMENSION(nb_files_max, nb_var_max), SAVE:: &
       last_opp, last_wrt, last_opp_chk, last_wrt_chk, nb_opp, nb_wrt, point

  ! Book keeping for the buffers
  INTEGER:: buff_pos = 0
  REAL, ALLOCATABLE, SAVE:: buffer(:)
  LOGICAL:: zoom(nb_files_max) = .FALSE., regular(nb_files_max) = .TRUE.

  ! Book keeping of the axes

  INTEGER, DIMENSION(nb_files_max, nb_var_max), SAVE:: tdimid, tax_last, &
       tax_name_length
  CHARACTER(LEN=40), DIMENSION(nb_files_max, nb_var_max), SAVE:: tax_name

  ! A list of functions which require special action
  ! (Needs to be updated when functions are added
  !  but they are well located here)

  CHARACTER(LEN=120):: indchfun = 'scatter, fill, gather, coll', &
       fuchnbout = 'scatter, fill'
  ! Some configurable variables with locks
  CHARACTER(LEN=80):: model_name = 'An IPSL model'
  LOGICAL:: lock_modname = .FALSE.

END MODULE histcom_var
