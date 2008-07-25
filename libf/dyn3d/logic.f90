module logic

  implicit none

  LOGICAL:: purmats= .FALSE.
  ! Help = Choix du schema d'integration temporel.
  ! y = pure Matsuno sinon c'est du Matsuno-leapfrog

  logical:: fxyhypb = .TRUE.
  ! (fonction f(y) à dérivée tangente hyperbolique, sinon à dérivée
  ! sinusoïdale)

  logical:: ysinus = .TRUE.
  ! (Fonction f(y) avec y = Sin(latit.) si = .true. sinon y = latit.)

  logical:: ok_guide= .FALSE.
  ! Help = Guidage

  INTEGER:: iflag_phys = 1
  ! contrôle l'appel à la physique :
  ! 0 : pas de physique
  ! 1 : physique normale (appel à phylmd, phymars...) (default)
  ! 2 : rappel Newtonien pour la température + friction au sol

contains

  subroutine read_logic

    namelist /logic_nml/ purmats, fxyhypb, ysinus, ok_guide, iflag_phys

    !------------------------------------

    print *, "Enter namelist 'logic_nml'."
    read(unit=*, nml=logic_nml)
    write(unit=*, nml=logic_nml)

  end subroutine read_logic

end module logic
