module sugwd_m

  IMPLICIT NONE

contains

  SUBROUTINE sugwd(paprs, pplay)

    ! Initialize yoegwd, the common that controls the gravity wave
    ! drag parametrization.

    ! REFERENCE.
    ! ECMWF Research Department documentation of the IFS

    ! AUTHOR.
    ! MARTIN MILLER *ECMWF*

    ! ORIGINAL : 90-01-01

    USE yoegwd, ONLY : gfrcrit, ghmax, gkdrag, gklift, gkwake, grahilo, &
         grcrit, gsigcr, gssec, gtsec, gvcrit, gvsec, nktopg, nstra
    use nr_util, only: assert_eq

    REAL, INTENT(IN):: paprs(:, :) ! (nlon, nlev+1)
    REAL, INTENT(IN):: pplay(:, :) ! (nlon, nlev)

    ! Local:
    INTEGER nlon, nlev
    integer jk
    REAL zpr, zstra, zsigt, zpm1r

    !------------------------------------------------------------

    print *, "Call sequence information: sugwd"
    nlon = assert_eq(size(paprs, 1), size(pplay, 1), "sugwd nlon")
    nlev = assert_eq(size(paprs, 2) - 1, size(pplay, 2), "sugwd nlon")

    ! 1. SET THE VALUES OF THE PARAMETERS

    ghmax = 10000.

    zpr = 100000.
    zstra = 0.1
    zsigt = 0.94

    DO jk = 1, nlev
       zpm1r = pplay(nlon / 2, jk) / paprs(nlon / 2, 1)
       IF (zpm1r >= zsigt) nktopg = jk
       IF (zpm1r >= zstra) nstra = jk
    end DO

    ! inversion car dans orodrag on compte les niveaux a l'envers
    nktopg = nlev - nktopg + 1
    nstra = nlev - nstra
    PRINT *, 'nktopg=', nktopg
    PRINT *, 'nstra=', nstra

    gsigcr = 0.8

    gkdrag = 0.2
    grahilo = 1.
    grcrit = 0.01
    gfrcrit = 1.
    gkwake = 0.5

    gklift = 0.5
    gvcrit = 0.

    ! 2. SET VALUES OF SECURITY PARAMETERS
    gvsec = 0.1
    gssec = 1E-12
    gtsec = 1E-7

  END SUBROUTINE sugwd

end module sugwd_m
