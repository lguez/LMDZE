module logic

  implicit none

  LOGICAL:: purmats= .FALSE.
  ! Help = Choix du schema d'integration temporel.
  ! y = pure Matsuno sinon c'est du Matsuno-leapfrog

  logical:: fxyhypb = .TRUE.
  ! (fonction f(y) � d�riv�e tangente hyperbolique, sinon � d�riv�e
  ! sinuso�dale)

  logical:: ysinus = .TRUE.
  ! (Fonction f(y) avec y = Sin(latit.) si = .true. sinon y = latit.)

  logical:: ok_guide= .FALSE.
  ! Help = Guidage

  INTEGER:: iflag_phys = 1
  ! contr�le l'appel � la physique :
  ! 0 : pas de physique
  ! 1 : physique normale (appel � phylmd, phymars...) (default)
  ! 2 : rappel Newtonien pour la temp�rature + friction au sol

contains

  subroutine read_logic

    namelist /logic_nml/ purmats, fxyhypb, ysinus, ok_guide, iflag_phys

    !------------------------------------

    print *, "Enter namelist 'logic_nml'."
    read(unit=*, nml=logic_nml)
    write(unit=*, nml=logic_nml)

  end subroutine read_logic

end module logic
