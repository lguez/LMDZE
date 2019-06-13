module dynetat0_chosen_m

  IMPLICIT NONE

  integer, protected, save:: day_ref ! jour de l'ann\'ee de l'\'etat initial
  ! (= 350 si 20 d\'ecembre par exemple)

  integer, protected, save:: annee_ref
  ! Annee de l'etat initial (avec 4 chiffres)

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
  
  real, protected, save:: pa ! in Pa

contains

  SUBROUTINE dynetat0_chosen

    ! This procedure reads the initial state of the atmosphere. Values
    ! that were chosen in ce0l.

    ! Libraries:
    use netcdf, only: NF90_NOWRITE
    use netcdf95, only: nf95_open, nf95_inq_varid, NF95_CLOSE, NF95_Gw_VAR
    use nr_util, only: assert

    use comconst, only: dtvr
    use conf_gcm_m, only: raz_date
    use dimensions, only: iim, jjm, llm
    use unit_nml_m, only: unit_nml

    ! Local: 
    REAL, allocatable:: tab_cntrl(:) ! tableau des param\`etres du run
    INTEGER ncid, varid

    namelist /dynetat0_nml/ day_ref, annee_ref

    !-----------------------------------------------------------------------

    print *, "Call sequence information: dynetat0_chosen"

    ! Fichier \'etat initial :
    call nf95_open("start.nc", NF90_NOWRITE, ncid)

    call nf95_inq_varid(ncid, "controle", varid)
    call NF95_Gw_VAR(ncid, varid, tab_cntrl)

    call assert(int(tab_cntrl(1)) == iim, "dynetat0_chosen tab_cntrl iim") 
    call assert(int(tab_cntrl(2)) == jjm, "dynetat0_chosen tab_cntrl jjm") 
    call assert(int(tab_cntrl(3)) == llm, "dynetat0_chosen tab_cntrl llm") 

    IF (dtvr /= tab_cntrl(12)) THEN
       print *, 'Warning: the time steps from day_step and "start.nc" ' // &
            'are different.'
       print *, 'dtvr from day_step: ', dtvr
       print *, 'dtvr from "start.nc": ', tab_cntrl(12)
       print *, 'Using the value from day_step.'
    ENDIF

    pa = tab_cntrl(18)

    clon = tab_cntrl(20)
    clat = tab_cntrl(21)
    grossismx = tab_cntrl(22)
    grossismy = tab_cntrl(23)
    dzoomx = tab_cntrl(25)
    dzoomy = tab_cntrl(26)
    taux = tab_cntrl(28)
    tauy = tab_cntrl(29)

    if (raz_date) then
       ! Default values:
       day_ref = 1
       annee_ref = 1998
       
       print *, "Enter namelist 'dynetat0_nml'."
       read(unit = *, nml = dynetat0_nml)
       write(unit_nml, nml = dynetat0_nml)
    else
       day_ref = tab_cntrl(4)
       annee_ref = tab_cntrl(5)
    end if

    call NF95_CLOSE(ncid)

  END SUBROUTINE dynetat0_chosen

  !********************************************************************

  subroutine read_serre

    use unit_nml_m, only: unit_nml
    use nr_util, only: assert, pi

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

    ! Default values:
    day_ref = 1
    annee_ref = 1998

    print *, "Enter namelist 'dynetat0_nml'."
    read(unit = *, nml = dynetat0_nml)
    write(unit_nml, nml = dynetat0_nml)

    pa = 5e4

  end subroutine read_serre

end module dynetat0_chosen_m
