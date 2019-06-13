module yoegwd

  ! From phylmd/YOEGWD.h, version 1.1.1.1 2004/05/19 12:53:08
  ! Parameters for gravity wave drag calculations

  implicit none

  integer NKTOPG, NSTRA
  real, parameter:: GFRCRIT = 1., GKWAKE = 0.5, GRCRIT = 0.01, GVCRIT = 0.
  real, parameter:: GKDRAG = 0.2, GKLIFT = 0.5, GRAHILO = 1., GSIGCR = 0.8

  ! SECURITY PARAMETERS:
  real, parameter:: GVSEC = 0.1, GSSEC = 1E-12, GTSEC = 1E-7

contains

  SUBROUTINE sugwd(paprs, pplay)

    ! Initialize yoegwd, the common that controls the gravity wave
    ! drag parametrization.

    ! REFERENCE: ECMWF Research Department documentation of the IFS
    ! AUTHOR: MARTIN MILLER *ECMWF*
    ! ORIGINAL : 90-01-01

    use nr_util, only: assert_eq, ifirstloc

    REAL, INTENT(IN):: paprs(:, :) ! (klon, llm + 1)
    REAL, INTENT(IN):: pplay(:, :) ! (klon, llm)

    ! Local:
    INTEGER klon, llm
    real zpm1r(size(pplay, 2)) ! (llm)

    !------------------------------------------------------------

    print *, "Call sequence information: sugwd"
    klon = assert_eq(size(paprs, 1), size(pplay, 1), "sugwd klon")
    llm = assert_eq(size(paprs, 2) - 1, size(pplay, 2), "sugwd llm")

    ! 1. SET THE VALUES OF THE PARAMETERS

    zpm1r = pplay(klon / 2, llm:1:- 1) / paprs(klon / 2, 1)
    ! inversion car dans orodrag on compte les niveaux \`a l'envers

    nktopg = ifirstloc(zpm1r >= 0.94)
    nstra = ifirstloc(zpm1r >= 0.1) - 1
    PRINT *, 'nktopg = ', nktopg
    PRINT *, 'nstra = ', nstra

  END SUBROUTINE sugwd

end module yoegwd
