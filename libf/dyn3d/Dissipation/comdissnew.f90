module comdissnew

  implicit none

  LOGICAL:: lstardis = .TRUE. ! choix de l'opérateur de dissipation star

  INTEGER:: nitergdiv = 1
  ! nombre d'itérations de l'opérateur de dissipation gradiv

  INTEGER:: nitergrot = 2
  ! nombre d'itérations de l'opérateur de dissipation nxgradrot

  INTEGER:: niterh = 2
  ! nombre d'itérations de l'opérateur de dissipation divgrad

  REAL:: tetagdiv = 7200. ! in s
  ! temps de dissipation des plus petites longueurs d'ondes pour u, v (grad div)

  REAL:: tetagrot = 7200. ! in s
  ! temps de dissipation des plus petites longueurs d'ondes pour u, v
  ! (nxgradrot)

  REAL:: tetatemp = 7200. ! in s
  ! temps de dissipation des plus petites longueurs d'ondes pour h (divgrad) 

  REAL:: coefdis = 0.

contains

  subroutine read_comdissnew

    use unit_nml_m, only: unit_nml
    
    namelist /comdissnew_nml/lstardis, nitergdiv, nitergrot, niterh, &
         tetagdiv, tetagrot, tetatemp, coefdis

    !-------------------------------------------------

    print *, "Enter namelist 'comdissnew_nml'."
    read(unit=*, nml=comdissnew_nml)
    write(unit_nml, nml=comdissnew_nml)

  end subroutine read_comdissnew

end module comdissnew
