    SUBROUTINE sugwd(nlon,nlev,paprs,pplay)

!**** *SUGWD* INITIALIZE COMMON YOEGWD CONTROLLING GRAVITY WAVE DRAG

!     PURPOSE.
!     --------
!           INITIALIZE YOEGWD, THE COMMON THAT CONTROLS THE
!           GRAVITY WAVE DRAG PARAMETRIZATION.

!**   INTERFACE.
!     ----------
!        CALL *SUGWD* FROM *SUPHEC*
!              -----        ------

!        EXPLICIT ARGUMENTS :
!        --------------------
!        PSIG        : VERTICAL COORDINATE TABLE
!        NLEV        : NUMBER OF MODEL LEVELS

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOEGWD

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!        NONE

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        MARTIN MILLER             *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 90-01-01
!     ------------------------------------------------------------------
      USE yoegwd
      IMPLICIT NONE

!     -----------------------------------------------------------------
!      ----------------------------------------------------------------

      INTEGER nlon, nlev, jk
      REAL, INTENT (IN) :: paprs(nlon,nlev+1)
      REAL, INTENT (IN) :: pplay(nlon,nlev)
      REAL zpr, zstra, zsigt, zpm1r

!*       1.    SET THE VALUES OF THE PARAMETERS
!              --------------------------------

100   CONTINUE

      PRINT *, ' DANS SUGWD NLEV=', nlev
      ghmax = 10000.

      zpr = 100000.
      zstra = 0.1
      zsigt = 0.94
!old  ZPR=80000.
!old  ZSIGT=0.85

      DO 110 jk = 1, nlev
        zpm1r = pplay(nlon/2,jk)/paprs(nlon/2,1)
        IF (zpm1r>=zsigt) THEN
          nktopg = jk
        END IF
        zpm1r = pplay(nlon/2,jk)/paprs(nlon/2,1)
        IF (zpm1r>=zstra) THEN
          nstra = jk
        END IF
110   CONTINUE

!  inversion car dans orodrag on compte les niveaux a l'envers
      nktopg = nlev - nktopg + 1
      nstra = nlev - nstra
      PRINT *, ' DANS SUGWD nktopg=', nktopg
      PRINT *, ' DANS SUGWD nstra=', nstra

      gsigcr = 0.80

      gkdrag = 0.2
      grahilo = 1.
      grcrit = 0.01
      gfrcrit = 1.0
      gkwake = 0.50

      gklift = 0.50
      gvcrit = 0.0


!      ----------------------------------------------------------------

!*       2.    SET VALUES OF SECURITY PARAMETERS
!              ---------------------------------

200   CONTINUE

      gvsec = 0.10
      gssec = 1.E-12

      gtsec = 1.E-07

!      ----------------------------------------------------------------

      RETURN
    END
