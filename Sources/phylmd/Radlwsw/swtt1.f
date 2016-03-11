SUBROUTINE swtt1(knu, kabs, kind, pu, ptr)
  USE dimens_m
  USE dimphy
  USE raddim
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
  ! ORIGINAL : 95-01-20
  ! -----------------------------------------------------------------------
  ! * ARGUMENTS:

  INTEGER knu ! INDEX OF THE SPECTRAL INTERVAL
  INTEGER kabs ! NUMBER OF ABSORBERS
  INTEGER kind(kabs) ! INDICES OF THE ABSORBERS
  DOUBLE PRECISION pu(kdlon, kabs) ! ABSORBER AMOUNT

  DOUBLE PRECISION ptr(kdlon, kabs) ! TRANSMISSION FUNCTION

  ! * LOCAL VARIABLES:

  DOUBLE PRECISION zr1(kdlon)
  DOUBLE PRECISION zr2(kdlon)
  DOUBLE PRECISION zu(kdlon)
  INTEGER jl, ja, i, j, ia

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

  DATA (d(1,i), i=1, 3)/0.00d0, 0.00d0, 0d0/
  DATA (d(2,i), i=1, 3)/0.000000000d0, 0.000000000d0, 0.800000000d0/
  ! -----------------------------------------------------------------------

  ! *         1.      HORNER'S ALGORITHM TO COMPUTE TRANSMISSION FUNCTION


  DO ja = 1, kabs
    ia = kind(ja)
    DO jl = 1, kdlon
      zu(jl) = pu(jl, ja)
      zr1(jl) = apad(knu, ia, 1) + zu(jl)*(apad(knu,ia,2)+zu(jl)*(apad(knu, &
        ia,3)+zu(jl)*(apad(knu,ia,4)+zu(jl)*(apad(knu,ia,5)+zu(jl)*(apad(knu, &
        ia,6)+zu(jl)*(apad(knu,ia,7)))))))

      zr2(jl) = bpad(knu, ia, 1) + zu(jl)*(bpad(knu,ia,2)+zu(jl)*(bpad(knu, &
        ia,3)+zu(jl)*(bpad(knu,ia,4)+zu(jl)*(bpad(knu,ia,5)+zu(jl)*(bpad(knu, &
        ia,6)+zu(jl)*(bpad(knu,ia,7)))))))


      ! *         2.      ADD THE BACKGROUND TRANSMISSION


      ptr(jl, ja) = (zr1(jl)/zr2(jl))*(1.-d(knu,ia)) + d(knu, ia)
    END DO
  END DO

  RETURN
END SUBROUTINE swtt1
