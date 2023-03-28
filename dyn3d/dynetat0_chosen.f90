module dynetat0_chosen_m

  IMPLICIT NONE

  integer, protected:: day_ref  = 1 ! jour de l'ann\'ee de l'\'etat initial
  ! (= 350 si 20 d\'ecembre par exemple)

  integer, protected:: annee_ref = 1998
  ! Ann\'ee de l'\'etat initial (avec 4 chiffres)

  REAL, protected, save:: clon ! longitude of the center of the zoom, in rad
  real, protected, save:: clat ! latitude of the center of the zoom, in rad

  real, protected, save:: grossismx, grossismy
  ! facteurs de grossissement du zoom, selon la longitude et la latitude
  ! = 2 si 2 fois, = 3 si 3 fois, etc.

  real, protected, save:: dzoomx, dzoomy
  ! extensions en longitude et latitude de la zone du zoom (fractions
  ! de la zone totale)

  real, protected, save:: taux, tauy
  ! raideur de la transition de l'int\'erieur \`a l'ext\'erieur du zoom
  
contains

  SUBROUTINE dynetat0_chosen(ncid_start)

    ! This procedure reads the initial state of the atmosphere. Values
    ! that were chosen in ce0l.

    ! Libraries:
    use jumble, only: assert
    use netcdf95, only: nf95_inq_varid, NF95_Gw_VAR

    use conf_gcm_m, only: raz_date
    use dimensions, only: iim, jjm, llm
    use unit_nml_m, only: unit_nml

    integer, intent(in):: ncid_start

    ! Local: 
    REAL, allocatable:: tab_cntrl(:) ! tableau des param\`etres du run
    INTEGER varid

    namelist /dynetat0_nml/ day_ref, annee_ref

    !-----------------------------------------------------------------------

    print *, "Call sequence information: dynetat0_chosen"

    call nf95_inq_varid(ncid_start, "controle", varid)
    call NF95_Gw_VAR(ncid_start, varid, tab_cntrl)

    call assert(int(tab_cntrl(1)) == iim, "dynetat0_chosen tab_cntrl iim") 
    call assert(int(tab_cntrl(2)) == jjm, "dynetat0_chosen tab_cntrl jjm") 
    call assert(int(tab_cntrl(3)) == llm, "dynetat0_chosen tab_cntrl llm") 

    clon = tab_cntrl(20)
    clat = tab_cntrl(21)
    grossismx = tab_cntrl(22)
    grossismy = tab_cntrl(23)
    dzoomx = tab_cntrl(25)
    dzoomy = tab_cntrl(26)
    taux = tab_cntrl(28)
    tauy = tab_cntrl(29)

    if (raz_date) then
       print *, "Enter namelist 'dynetat0_nml'."
       read(unit = *, nml = dynetat0_nml)
       write(unit_nml, nml = dynetat0_nml)
    else
       day_ref = tab_cntrl(4)
       annee_ref = tab_cntrl(5)
    end if

  END SUBROUTINE dynetat0_chosen

  !********************************************************************

  subroutine read_serre

    use jumble, only: assert, pi

    use unit_nml_m, only: unit_nml

    REAL:: clon_deg = 0. ! longitude of the center of the zoom, in degrees
    real:: clat_deg = 0. ! latitude of the center of the zoom, in degrees

    namelist /serre_nml/ clon_deg, clat_deg, grossismx, grossismy, dzoomx, &
         dzoomy, taux, tauy
    namelist /dynetat0_nml/ day_ref, annee_ref

    !-------------------------------------------------

    ! Default values:
    grossismx = 1.
    grossismy = 1.
    dzoomx = 0.2
    dzoomy = 0.2
    taux = 3.
    tauy = 3. 

    print *, "Enter namelist 'serre_nml'."
    read(unit=*, nml=serre_nml)
    write(unit_nml, nml=serre_nml)

    call assert(grossismx >= 1. .and. grossismy >= 1., "read_serre grossism")
    call assert(dzoomx > 0., dzoomx < 1., dzoomy < 1., &
         "read_serre dzoomx dzoomy")
    clon = clon_deg / 180. * pi
    clat = clat_deg / 180. * pi

    print *, "Enter namelist 'dynetat0_nml'."
    read(unit = *, nml = dynetat0_nml)
    write(unit_nml, nml = dynetat0_nml)

  end subroutine read_serre

end module dynetat0_chosen_m
