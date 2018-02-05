module histdef_m

  USE histcom_var, ONLY: nb_files_max, nb_var_max

  implicit none

  INTEGER:: buff_pos = 0
  INTEGER, SAVE:: point(nb_files_max, nb_var_max)
  private nb_files_max, nb_var_max

contains

  SUBROUTINE histdef(fileid, varname, ptitle, unit, xsize, ysize, &
       horiid, pzsize, oriz, szz, zid, opp, pfreq_opp, pfreq_wrt)

    ! With this subroutine each variable to be archived on the history
    ! tape should be declared. It gives the user the choice of
    ! operation to be performed on the variable, the frequency of
    ! this operation and the frequency of the archiving.

    USE buildop_m, ONLY: buildop
    USE errioipsl, ONLY: histerr
    USE find_str_m, ONLY: find_str
    use histbeg_totreg_m, only: deltat
    USE histcom_var, ONLY: freq_opp, freq_wrt, fullop, full_size, itau0, &
         last_opp, last_opp_chk, last_wrt, last_wrt_chk, missing_val, name, &
         name_length, nbopp, nbopp_max, nb_hax, nb_opp, nb_tax, nb_var, &
         nb_wrt, nb_zax, scal, scsize, slab_ori, slab_sz, sopps, &
         tax_last, tax_name, tax_name_length, title, topp, unit_name, &
         var_axid, var_haxid, var_zaxid, zax_name, zax_size, zorig, zsize
    USE ioget_calendar_m, ONLY: ioget_calendar_real

    INTEGER, INTENT(IN):: fileid
    ! (ID of the file the variable should be archived in)

    CHARACTER(len=*), INTENT(IN):: varname
    ! (name of the variable, short and easy to remember)

    CHARACTER(len=*), INTENT(IN):: ptitle ! full name of the variable
    CHARACTER(len=*), INTENT(IN):: unit ! units of the variable

    ! The next 3 arguments give the size of that data
    ! that will be passed to histwrite. The zoom will be
    ! done there with the horizontal information obtained
    ! in "histbeg" and the vertical information to follow.
    INTEGER, INTENT(IN):: xsize, ysize ! Sizes in X and Y directions
    INTEGER, INTENT(IN):: horiid ! ID of the horizontal axis

    ! The next two arguments give the vertical zoom to use.

    INTEGER, INTENT(IN):: pzsize
    ! (Size in Z direction (If 1 then no axis is declared for this
    ! variable and zid is not used)

    INTEGER, INTENT(IN):: oriz ! Off set of the zoom
    INTEGER, INTENT(IN):: szz ! Size of the zoom

    INTEGER, INTENT(IN):: zid
    ! (ID of the vertical axis to use. It has to have the size of the zoom.)

    CHARACTER(len=*), INTENT(IN):: opp
    ! Operation to be performed. The following options exist today:
    ! inst: keeps instantaneous values for writting
    ! ave: Computes the average from call between writes

    REAL, INTENT(IN):: pfreq_opp ! Frequency of this operation (in seconds)

    REAL, INTENT(IN):: pfreq_wrt
    ! (Frequency at which the variable should be written, in seconds)

    ! Local:

    INTEGER:: iv, i, nb
    CHARACTER(len=70):: str70, str71, str72
    CHARACTER(len=20):: tmp_name
    CHARACTER(len=20):: str20, tab_str20(nb_var_max)
    INTEGER:: tab_str20_length(nb_var_max)
    CHARACTER(len=40):: str40, tab_str40(nb_var_max)
    INTEGER:: tab_str40_length(nb_var_max)
    CHARACTER(len=10):: str10
    CHARACTER(len=80):: tmp_str80
    CHARACTER(len=7):: tmp_topp, tmp_sopp(nbopp_max)
    CHARACTER(len=120):: ex_topps
    REAL:: tmp_scal(nbopp_max), un_an, un_jour, test_fopp, test_fwrt
    INTEGER:: pos, buff_sz

    !---------------------------------------------------------------------

    ex_topps = 'ave, inst, t_min, t_max, t_sum, once, never, l_max, l_min'

    nb_var(fileid) = nb_var(fileid) + 1
    iv = nb_var(fileid)

    IF (iv>nb_var_max) THEN
       CALL histerr(3, 'histdef', &
            'Table of variables too small. You should increase nb_var_max', &
            'in M_HISTCOM.f90 in order to accomodate all these variables', ' ')
    END IF

    ! 1.0 Transfer informations on the variable to the common
    !     and verify that it does not already exist

    IF (iv>1) THEN
       str20 = varname
       nb = iv - 1
       tab_str20(1:nb) = name(fileid, 1:nb)
       tab_str20_length(1:nb) = name_length(fileid, 1:nb)
       CALL find_str(nb, tab_str20, tab_str20_length, str20, pos)
    ELSE
       pos = 0
    END IF

    IF (pos>0) THEN
       str70 = 'Variable already exists'
       WRITE (str71, '("Check variable  ", a, " in file", I3)') str20, &
            fileid
       str72 = 'Can also be a wrong file ID in another declaration'
       CALL histerr(3, 'histdef', str70, str71, str72)
    END IF

    name(fileid, iv) = varname
    name_length(fileid, iv) = len_trim(name(fileid, iv))
    title(fileid, iv) = ptitle
    unit_name(fileid, iv) = unit
    tmp_name = name(fileid, iv)

    ! 1.1 decode the operations

    fullop(fileid, iv) = opp
    tmp_str80 = opp
    CALL buildop(tmp_str80, ex_topps, tmp_topp, nbopp_max, missing_val, &
         tmp_sopp, tmp_scal, nbopp(fileid, iv))

    topp(fileid, iv) = tmp_topp
    DO i = 1, nbopp(fileid, iv)
       sopps(fileid, iv, i) = tmp_sopp(i)
       scal(fileid, iv, i) = tmp_scal(i)
    END DO

    ! 1.2 If we have an even number of operations
    !     then we need to add identity

    IF (2*int(nbopp(fileid, iv)/2.0)==nbopp(fileid, iv)) THEN
       nbopp(fileid, iv) = nbopp(fileid, iv) + 1
       sopps(fileid, iv, nbopp(fileid, iv)) = 'ident'
       scal(fileid, iv, nbopp(fileid, iv)) = missing_val
    END IF

    ! 2.0 Put the size of the variable in the common and check

    scsize(fileid, iv, :) = (/ xsize, ysize, pzsize/)

    zorig(fileid, iv, 1:3) = (/ slab_ori(fileid, 1), slab_ori(fileid, 2), &
         oriz/)

    zsize(fileid, iv, 1:3) = (/ slab_sz(fileid, 1), slab_sz(fileid, 2), &
         szz/)

    ! Is the size of the full array the same as that of the coordinates  ?

    IF ((xsize>full_size(fileid, 1)) .OR. (ysize>full_size(fileid, &
         2))) THEN

       str70 = 'The size of the variable is different ' // &
            'from the one of the coordinates'
       WRITE (str71, '("Size of coordinates:", 2I4)') full_size(fileid, 1), &
            full_size(fileid, 2)
       WRITE (str72, '("Size declared for variable ", a, ":", 2I4)') &
            trim(tmp_name), xsize, ysize
       CALL histerr(3, 'histdef', str70, str71, str72)
    END IF

    ! Is the size of the zoom smaler than the coordinates ?

    IF ((full_size(fileid, 1)<slab_sz(fileid, 1)) .OR. (full_size(fileid, &
         2)<slab_sz(fileid, 2))) THEN
       str70 = 'Size of variable should be greater or equal &
            &to those of the zoom'
       WRITE (str71, '("Size of XY zoom:", 2I4)') slab_sz(fileid, 1), &
            slab_sz(fileid, 1)
       WRITE (str72, '("Size declared for variable ", a, ":", 2I4)') &
            trim(tmp_name), xsize, ysize
       CALL histerr(3, 'histdef', str70, str71, str72)
    END IF

    ! 2.1 We store the horizontal grid information with minimal
    !     and a fall back onto the default grid

    IF (horiid>0 .AND. horiid<=nb_hax(fileid)) THEN
       var_haxid(fileid, iv) = horiid
    ELSE
       var_haxid(fileid, iv) = 1
       CALL histerr(2, 'histdef', &
            'We use the default grid for variable as an invalide', &
            'ID was provided for variable: ', varname)
    END IF

    ! 2.2 Check the vertical coordinates if needed

    IF (szz>1) THEN

       ! Does the vertical coordinate exist ?

       IF (zid>nb_zax(fileid)) THEN
          WRITE (str70, '("The vertical coordinate chosen for variable ", a)' &
               ) trim(tmp_name)
          str71 = ' Does not exist.'
          CALL histerr(3, 'histdef', str70, str71, ' ')
       END IF

       ! Is the vertical size of the variable equal to that of the axis ?

       IF (szz/=zax_size(fileid, zid)) THEN
          str20 = zax_name(fileid, zid)
          WRITE (str71, '("Size of zoom in z:", I4)') szz
          WRITE (str72, '("Size declared for axis ", a, ":", I4)') &
               trim(str20), zax_size(fileid, zid)
          CALL histerr(3, 'histdef', 'The size of the zoom does not ' &
               // 'correspond to the size of the chosen vertical axis', &
               str71, str72)
       END IF

       ! Is the zoom smaler that the total size of the variable ?

       IF (pzsize<szz) THEN
          str20 = zax_name(fileid, zid)
          str70 = 'The vertical size of variable ' // &
               'is smaller than that of the zoom.'
          WRITE (str71, '("Declared vertical size of data:", I5)') pzsize
          WRITE (str72, '("Size of zoom for variable ", a, " = ", I5)') &
               trim(tmp_name), szz
          CALL histerr(3, 'histdef', str70, str71, str72)
       END IF
       var_zaxid(fileid, iv) = zid
    ELSE
       var_zaxid(fileid, iv) = -99
    END IF

    ! 3.0 Determine the position of the variable in the buffer
    !     If it is instantaneous output then we do not use the buffer

    ! 3.1 We get the size of the arrays histwrite will get and check
    !     that they fit into the tmp_buffer

    buff_sz = zsize(fileid, iv, 1)*zsize(fileid, iv, 2)*zsize(fileid, iv, 3)

    ! 3.2 move the pointer of the buffer array for operation
    !     which need bufferisation

    IF ((trim(tmp_topp)/='inst') .AND. (trim(tmp_topp)/='once') .AND. ( &
         trim(tmp_topp)/='never')) THEN
       point(fileid, iv) = buff_pos + 1
       buff_pos = buff_pos + buff_sz
    END IF

    ! 4.0 Transfer the frequency of the operations and check
    !     for validity. We have to pay attention to negative values
    !     of the frequency which indicate monthly time-steps.
    !     The strategy is to bring it back to seconds for the tests

    freq_opp(fileid, iv) = pfreq_opp
    freq_wrt(fileid, iv) = pfreq_wrt

    CALL ioget_calendar_real(un_an, un_jour)
    IF (pfreq_opp<0) THEN
       CALL ioget_calendar_real(un_an)
       test_fopp = pfreq_opp*(-1.)*un_an/12.*un_jour
    ELSE
       test_fopp = pfreq_opp
    END IF
    IF (pfreq_wrt<0) THEN
       CALL ioget_calendar_real(un_an)
       test_fwrt = pfreq_wrt*(-1.)*un_an/12.*un_jour
    ELSE
       test_fwrt = pfreq_wrt
    END IF

    ! 4.1 Frequency of operations and output should be larger than deltat !

    IF (test_fopp<deltat(fileid)) THEN
       str70 = 'Frequency of operations should be larger than deltat'
       WRITE (str71, '("It is not the case for variable ", a, ":", F10.4)') &
            trim(tmp_name), pfreq_opp
       str72 = 'PATCH: frequency set to deltat'

       CALL histerr(2, 'histdef', str70, str71, str72)

       freq_opp(fileid, iv) = deltat(fileid)
    END IF

    IF (test_fwrt<deltat(fileid)) THEN
       str70 = 'Frequency of output should be larger than deltat'
       WRITE (str71, '("It is not the case for variable ", a, ":", F10.4)') &
            trim(tmp_name), pfreq_wrt
       str72 = 'PATCH: frequency set to deltat'

       CALL histerr(2, 'histdef', str70, str71, str72)

       freq_wrt(fileid, iv) = deltat(fileid)
    END IF

    ! 4.2 First the existence of the operation is tested and then
    !     its compatibility with the choice of frequencies

    IF (trim(tmp_topp)=='inst') THEN
       IF (test_fopp/=test_fwrt) THEN
          str70 = 'For instantaneous output the frequency ' // &
               'of operations and output'
          WRITE (str71, &
               '("should be the same, this was not case for variable ", a)') &
               trim(tmp_name)
          str72 = 'PATCH: The smalest frequency of both is used'
          CALL histerr(2, 'histdef', str70, str71, str72)
          IF (test_fopp<test_fwrt) THEN
             freq_opp(fileid, iv) = pfreq_opp
             freq_wrt(fileid, iv) = pfreq_opp
          ELSE
             freq_opp(fileid, iv) = pfreq_wrt
             freq_wrt(fileid, iv) = pfreq_wrt
          END IF
       END IF
    ELSE IF (index(ex_topps, trim(tmp_topp))>0) THEN
       IF (test_fopp>test_fwrt) THEN
          str70 = 'For averages the frequency of operations ' // &
               'should be smaller or equal'
          WRITE (str71, &
               '("to that of output. It is not the case for variable ", a)') &
               trim(tmp_name)
          str72 = 'PATCH: The output frequency is used for both'
          CALL histerr(2, 'histdef', str70, str71, str72)
          freq_opp(fileid, iv) = pfreq_wrt
       END IF
    ELSE
       WRITE (str70, '("Operation on variable ", a, " is unknown")') &
            trim(tmp_name)
       WRITE (str71, '("operation requested is:", a)') tmp_topp
       WRITE (str72, '("File ID:", I3)') fileid
       CALL histerr(3, 'histdef', str70, str71, str72)
    END IF

    ! 5.0 Initialize other variables of the common

    last_opp(fileid, iv) = itau0(fileid)
    ! - freq_opp(fileid, iv)/2./deltat(fileid)
    last_wrt(fileid, iv) = itau0(fileid)
    ! - freq_wrt(fileid, iv)/2./deltat(fileid)
    last_opp_chk(fileid, iv) = itau0(fileid)
    ! - freq_opp(fileid, iv)/2./deltat(fileid)
    last_wrt_chk(fileid, iv) = itau0(fileid)
    ! - freq_wrt(fileid, iv)/2./deltat(fileid)
    nb_opp(fileid, iv) = 0
    nb_wrt(fileid, iv) = 0

    ! 6.0 Get the time axis for this variable

    IF (freq_wrt(fileid, iv)>0) THEN
       WRITE (str10, '(I8.8)') int(freq_wrt(fileid, iv))
       str40 = trim(tmp_topp) // '_' // trim(str10)
    ELSE
       WRITE (str10, '(I2.2, "month")') abs(int(freq_wrt(fileid, iv)))
       str40 = trim(tmp_topp) // '_' // trim(str10)
    END IF

    DO i = 1, nb_tax(fileid)
       tab_str40(i) = tax_name(fileid, i)
       tab_str40_length(i) = tax_name_length(fileid, i)
    END DO

    CALL find_str(nb_tax(fileid), tab_str40, tab_str40_length, str40, pos)

    ! No time axis for once, l_max, l_min or never operation

    IF ((trim(tmp_topp)/='once') .AND. (trim(tmp_topp)/='never') .AND. ( &
         trim(tmp_topp)/='l_max') .AND. (trim(tmp_topp)/='l_min')) THEN
       IF (pos<0) THEN
          nb_tax(fileid) = nb_tax(fileid) + 1
          tax_name(fileid, nb_tax(fileid)) = str40
          tax_name_length(fileid, nb_tax(fileid)) = len_trim(str40)
          tax_last(fileid, nb_tax(fileid)) = 0
          var_axid(fileid, iv) = nb_tax(fileid)
       ELSE
          var_axid(fileid, iv) = pos
       END IF
    ELSE
       var_axid(fileid, iv) = -99
    END IF

    ! 7.0 prepare frequence of writing and operation
    !     for never or once operation

    IF ((trim(tmp_topp)=='once') .OR. (trim(tmp_topp)=='never')) THEN
       freq_opp(fileid, iv) = 0.
       freq_wrt(fileid, iv) = 0.
    END IF

  END SUBROUTINE histdef

end module histdef_m
