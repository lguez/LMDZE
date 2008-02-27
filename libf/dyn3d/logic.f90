module logic

  implicit none

  LOGICAL:: purmats= .FALSE.
  ! Help = Choix du schema d'integration temporel.
  ! y = pure Matsuno sinon c'est du Matsuno-leapfrog

  LOGICAL, save:: forward, leapf, apphys, statcl, conser
  logical, save:: apdiss, apdelq, saison

  logical, save:: fxyhypb
  ! fonction f(y) hyperbolique, sinon sinuso�dale

  logical, save:: ysinus

  logical:: ok_guide= .FALSE.
  ! Help = Guidage

  INTEGER:: iflag_phys = 1
  ! contr�le l'appel � la physique :
  ! 0 : pas de physique
  ! 1 : physique normale (appel � phylmd, phymars...) (default)
  ! 2 : rappel Newtonien pour la temp�rature + friction au sol

end module logic
