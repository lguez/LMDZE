module comdissnew

  ! The variables declared here are read from the file "run.def"

  implicit none

  LOGICAL:: lstardis= .TRUE.
  ! Help = choix de l'operateur de dissipation
  ! 'y' si on veut star et 'n' si on veut non-star
  ! Moi y en a pas comprendre ! 

  INTEGER:: nitergdiv= 1
  ! Help = nombre d'iterations de l'operateur de dissipation 
  ! gradiv

  INTEGER:: nitergrot= 2
  ! Help = nombre d'iterations de l'operateur de dissipation  
  ! nxgradrot

  INTEGER:: niterh= 2
  ! Help = nombre d'iterations de l'operateur de dissipation
  ! divgrad

  REAL::     tetagdiv= 7200.
  ! Help = temps de dissipation des plus petites longeur 
  ! d'ondes pour u,v (gradiv)

  REAL:: tetagrot= 7200.
  ! Help = temps de dissipation des plus petites longeur 
  ! d'ondes pour u,v (nxgradrot)

  REAL::  tetatemp= 7200.
  ! Help =  temps de dissipation des plus petites longeur 
  ! d'ondes pour h (divgrad)   

  REAL:: coefdis= 0.
  ! Help = coefficient pour gamdissip  

contains

  subroutine read_comdissnew

    namelist /comdissnew_nml/lstardis, nitergdiv, nitergrot, niterh, &
         tetagdiv, tetagrot, tetatemp, coefdis

    !-------------------------------------------------

    print *, "Call sequence information: read_comdissnew"
    print *, "Enter namelist 'comdissnew_nml'."
    read(unit=*, nml=comdissnew_nml)
    write(unit=*, nml=comdissnew_nml)

  end subroutine read_comdissnew

end module comdissnew
