module unit_nml_m

  implicit none

  integer, protected:: unit_nml
  ! logical unit number for file containing used namelists

contains

  subroutine set_unit_nml

    use jumble, only: new_unit

    !--------------------------------------------------------------

    call new_unit(unit_nml)
    
  end subroutine set_unit_nml
  
end module unit_nml_m
