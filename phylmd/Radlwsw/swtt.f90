SUBROUTINE swtt(knu, ka, pu, ptr)
  USE dimensions
  USE dimphy
  use conf_phys_m, only: kdlon
  IMPLICIT NONE

  ! -----------------------------------------------------------------------
  ! PURPOSE.
  ! --------
  ! THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
  ! ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN THE TWO SPECTRAL
  ! INTERVALS.

  ! METHOD.
  ! -------

  ! TRANSMISSION FUNCTION ARE COMPUTED USING PADE APPROXIMANTS
  ! AND HORNER'S ALGORITHM.

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
  ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 88-12-15
  ! -----------------------------------------------------------------------

  ! * ARGUMENTS

  INTEGER knu ! INDEX OF THE SPECTRAL INTERVAL
  INTEGER ka ! INDEX OF THE ABSORBER
  DOUBLE PRECISION pu(kdlon) ! ABSORBER AMOUNT

  DOUBLE PRECISION ptr(kdlon) ! TRANSMISSION FUNCTION

  ! * LOCAL VARIABLES:

  DOUBLE PRECISION zr1(kdlon), zr2(kdlon)
  INTEGER jl, i, j

  ! * Prescribed Data:

  DOUBLE PRECISION apad(2, 3, 7), bpad(2, 3, 7), d(2, 3)
  SAVE apad, bpad, d
  DATA ((apad(1,i,j),i=1,3), j=1, 7)/0.912418292D+05, 0.000000000D-00, &
    0.925887084D-04, 0.723613782D+05, 0.000000000D-00, 0.129353723D-01, &
    0.596037057D+04, 0.000000000D-00, 0.800821928D+00, 0.000000000D-00, &
    0.000000000D-00, 0.242715973D+02, 0.000000000D-00, 0.000000000D-00, &
    0.878331486D+02, 0.000000000D-00, 0.000000000D-00, 0.191559725D+02, &
    0.000000000D-00, 0.000000000D-00, 0.000000000D+00/
  DATA ((apad(2,i,j),i=1,3), j=1, 7)/0.376655383D-08, 0.739646016D-08, &
    0.410177786D+03, 0.978576773D-04, 0.131849595D-03, 0.672595424D+02, &
    0.387714006D+00, 0.437772681D+00, 0.000000000D-00, 0.118461660D+03, &
    0.151345118D+03, 0.000000000D-00, 0.119079797D+04, 0.233628890D+04, &
    0.000000000D-00, 0.293353397D+03, 0.797219934D+03, 0.000000000D-00, &
    0.000000000D+00, 0.000000000D+00, 0.000000000D+00/

  DATA ((bpad(1,i,j),i=1,3), j=1, 7)/0.912418292D+05, 0.000000000D-00, &
    0.925887084D-04, 0.724555318D+05, 0.000000000D-00, 0.131812683D-01, &
    0.602593328D+04, 0.000000000D-00, 0.812706117D+00, 0.100000000D+01, &
    0.000000000D-00, 0.249863591D+02, 0.000000000D-00, 0.000000000D-00, &
    0.931071925D+02, 0.000000000D-00, 0.000000000D-00, 0.252233437D+02, &
    0.000000000D-00, 0.000000000D-00, 0.100000000D+01/
  DATA ((bpad(2,i,j),i=1,3), j=1, 7)/0.376655383D-08, 0.739646016D-08, &
    0.410177786D+03, 0.979023421D-04, 0.131861712D-03, 0.731185438D+02, &
    0.388611139D+00, 0.437949001D+00, 0.100000000D+01, 0.120291383D+03, &
    0.151692730D+03, 0.000000000D+00, 0.130531005D+04, 0.237071130D+04, &
    0.000000000D+00, 0.415049409D+03, 0.867914360D+03, 0.000000000D+00, &
    0.100000000D+01, 0.100000000D+01, 0.000000000D+00/

  DATA (d(1,i), i=1, 3)/0.00d0, 0.00d0, 0.00d0/
  DATA (d(2,i), i=1, 3)/0.000000000d0, 0.000000000d0, 0.800000000d0/

  ! -----------------------------------------------------------------------

  ! *         1.      HORNER'S ALGORITHM TO COMPUTE TRANSMISSION FUNCTION


  DO jl = 1, kdlon
    zr1(jl) = apad(knu, ka, 1) + pu(jl)*(apad(knu,ka,2)+pu(jl)*(apad(knu,ka, &
      3)+pu(jl)*(apad(knu,ka,4)+pu(jl)*(apad(knu,ka,5)+pu(jl)*(apad(knu,ka,6) &
      +pu(jl)*(apad(knu,ka,7)))))))

    zr2(jl) = bpad(knu, ka, 1) + pu(jl)*(bpad(knu,ka,2)+pu(jl)*(bpad(knu,ka, &
      3)+pu(jl)*(bpad(knu,ka,4)+pu(jl)*(bpad(knu,ka,5)+pu(jl)*(bpad(knu,ka,6) &
      +pu(jl)*(bpad(knu,ka,7)))))))


    ! *         2.      ADD THE BACKGROUND TRANSMISSION



    ptr(jl) = (zr1(jl)/zr2(jl))*(1.-d(knu,ka)) + d(knu, ka)
  END DO

  RETURN
END SUBROUTINE swtt
