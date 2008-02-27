module in_out

  implicit none

  private prep_file, go_column, s_pr_mat, d_pr_mat
  interface pr_matrix
     module procedure s_pr_mat, d_pr_mat
  end interface

contains

  !***********************************************************

  integer function new_unit()

    logical opened, exist

    !------------------------------------------------------

    new_unit = 0
    do
       inquire(unit=new_unit, opened=opened, exist=exist)
       if (exist .and. .not. opened) exit
       new_unit = new_unit + 1
    end do

  end function new_unit

  !***********************************************************

  function csvread(file, first_r, first_c, last_r, last_c)

    ! Reads comma-separated numeric values in a file. The
    ! last column and/or last row parameters may be 0. This is
    ! interpreted as "last in the file".

    real, pointer:: csvread(:,:)

    character(len=*), intent(in):: file
    integer, intent(in), optional:: first_r ! (first row to read)
    integer, intent(in), optional:: first_c ! (first column to read)
    integer, intent(in), optional:: last_r ! (last row to read)
    integer, intent(in), optional:: last_c ! (last column to read)

    ! Variables local to the subprogram:
    integer i, unit
    integer f_r_loc ! (first row to read, local variable)
    integer f_c_loc ! (first column to read, local variable)
    integer l_r_loc ! (last row to read, local variable)
    integer l_c_loc ! (last column to read, local variable)

    !------------------------------------------------------

    print *, 'Reading data from file "' // file // '"'
    unit = new_unit()
    open(unit, file=file, status='old', action='read', position='rewind')

    call prep_file(unit, first_r, first_c, last_r, last_c, f_r_loc, &
         f_c_loc, l_r_loc, l_c_loc)

    allocate(csvread(l_r_loc - f_r_loc + 1, l_c_loc - f_c_loc + 1))

    do i = 1, l_r_loc - f_r_loc + 1
       call go_column(unit, f_c_loc)
       read(unit, fmt=*) csvread(i, :)
    end do
    ! (no implicit loop in read to allow partial reading of a line)

    close(unit)

  end function csvread

  !***********************************************************

  function csvread_dp(file, first_r, first_c, last_r, last_c)

    ! Reads comma-separated numeric values from a file, into a
    ! double precision array. The last column and/or last row parameters may be
    ! 0. This is interpreted as "last in the file".

    double precision, pointer:: csvread_dp(:,:)

    character(len=*), intent(in):: file
    integer, intent(in), optional:: first_r ! (first row to read)
    integer, intent(in), optional:: first_c ! (first column to read)
    integer, intent(in), optional:: last_r ! (last row to read)
    integer, intent(in), optional:: last_c ! (last column to read)

    ! Variables local to the subprogram:
    integer i, unit
    integer f_r_loc ! (first row to read, local variable)
    integer f_c_loc ! (first column to read, local variable)
    integer l_r_loc ! (last row to read, local variable)
    integer l_c_loc ! (last column to read, local variable)

    !------------------------------------------------------

    print *, 'Reading data from file "' // file // '"'
    unit = new_unit()
    open(unit, file=file, status='old', action='read', position='rewind')

    call prep_file(unit, first_r, first_c, last_r, last_c, f_r_loc, &
         f_c_loc, l_r_loc, l_c_loc)

    allocate(csvread_dp(l_r_loc - f_r_loc + 1, l_c_loc - f_c_loc + 1))

    do i = 1, l_r_loc - f_r_loc + 1
       call go_column(unit, f_c_loc)
       read(unit, fmt=*) csvread_dp(i, :)
    end do
    ! (no implicit loop in read to allow partial reading of a line)

    close(unit)

  end function csvread_dp

  !***********************************************************

  subroutine go_column(unit, column)

    ! This subroutine is used by the various versions of "csvread". On
    ! the current line of file, it advances to the input column. Columns are
    ! assumend to be separated by commas.

    integer, intent(in):: unit ! logical unit for input file
    integer, intent(in):: column

    ! Variables local to the subprogram:
    integer j
    character c

    !------------------------------------------------------

    ! Skip columns before "column" :
    j = 1 ! column index
    do while (j <= column - 1)
       read(unit, fmt='(a)', advance='no') c
       if (c == ',') j = j + 1
    end do

  end subroutine go_column

  !***********************************************************

  subroutine prep_file(unit, first_r, first_c, last_r, last_c, f_r_not_opt, &
       f_c_not_opt, l_r_not_opt, l_c_not_opt)

    ! This subroutine is used by the various versions of "csvread". It
    ! fills non-optional arguments: first and last row, first and last
    ! column which will actually be read, taking information from the
    ! file itself if necessary. It also positions the input file on the
    ! first row to read.

    integer, intent(in):: unit ! logical unit for input file
    integer, intent(in), optional:: first_r ! (first row to read)
    integer, intent(in), optional:: first_c ! (first column to read)
    integer, intent(in), optional:: last_r ! (last row to read)
    integer, intent(in), optional:: last_c ! (last column to read)
    integer, intent(out):: f_r_not_opt ! (first row to read, not optional)
    integer, intent(out):: f_c_not_opt ! (first column to read, not optional)
    integer, intent(out):: l_r_not_opt ! (last row to read, not optional)
    integer, intent(out):: l_c_not_opt ! (last column to read, not optional)

    ! Variables local to the subprogram:
    integer iostat, i
    character c

    !------------------------------------------------------

    f_r_not_opt = opt_merge(first_r, 1)
    f_c_not_opt = opt_merge(first_c, 1)
    l_r_not_opt = opt_merge(last_r, 0)
    l_c_not_opt = opt_merge(last_c, 0)

    if (l_r_not_opt == 0) then
       ! Count the number of lines in the file:
       i = 0
       do
          read(unit, fmt=*, iostat=iostat)
          if (iostat /= 0) exit
          i = i + 1
       end do
       l_r_not_opt = i
       if (l_r_not_opt == 0) stop 'Empty file.'

       rewind(unit)
    end if

    ! Go to first row to read:
    do i = 1, f_r_not_opt - 1
       read(unit, fmt=*)
    end do

    if (l_c_not_opt == 0) then
       ! Count the number of values per line:
       i = 0
       do
          read(unit, fmt='(a)', advance='no', iostat=iostat) c
          if (iostat /= 0) exit
          if (c == ',') i = i + 1
       end do
       l_c_not_opt = i + 1

       backspace(unit)
    end if

    print *, 'Reading column(s) ', f_c_not_opt, ':', l_c_not_opt, &
         ', row(s) ', f_r_not_opt, ':', l_r_not_opt

  end subroutine prep_file

  !***********************************************************
  
  integer function opt_merge(param, default)

    ! Analogous to the intrinsic procedure "merge" : merges an
    ! optional parameter and a default value depending on the
    ! presence of the optional parameter.

    integer, intent(in), optional:: param
    integer, intent(in):: default

    !--------------

    if (present(param)) then
       opt_merge = param
    else
       opt_merge = default
    end if

  end function opt_merge

  !***********************************************************

  subroutine s_pr_mat(name, a)

    ! This subroutine prints a rank 2 real matrix.

    character(len=*), intent(in):: name
    real, intent(in):: a(:,:)

    character(len=20) fmt
    integer n_lines, n_col, i

    !-----------------

    n_lines = size(a, 1)
    n_col = size(a, 2)
    if (n_lines <= 10 .and. n_col <= 5) then
       print *, name, ":"
       write(unit=fmt, fmt='("(1p, ", i0, "(g10.3: 1X))")') n_col
       do i = 1, n_lines
          print fmt, a(i, :)
       end do
    else
       print *, '"', name, '" is too big to print.'
    end if

  end subroutine s_pr_mat

  !***********************************************************

  subroutine d_pr_mat(name, a)

    ! This subroutine prints a rank 2 double precision matrix.

    character(len=*), intent(in):: name
    double precision, intent(in):: a(:,:)

    character(len=20) fmt
    integer n_lines, n_col, i

    !-----------------

    n_lines = size(a, 1)
    n_col = size(a, 2)
    if (n_lines <= 10 .and. n_col <= 5) then
       print *, name, ":"
       write(unit=fmt, fmt='("(1p, ", i0, "(g8.1: 1X))")') n_col
       do i = 1, n_lines
          print fmt, a(i, :)
       end do
    else
       print *, '"', name, '" is too big to print.'
    end if

  end subroutine d_pr_mat

end module in_out
