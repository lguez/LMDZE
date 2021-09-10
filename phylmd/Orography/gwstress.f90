module gwstress_m

  IMPLICIT NONE

contains

  SUBROUTINE gwstress(nlon,nlev,ktest,kkenvh,prho,pstab,pvph, &
       pstd,psig,pmea,ppic,ptau,pgeom1,pdmod)

    !**** *gwstress*

    !     purpose.
    !     --------

    !**   interface.
    !     ----------
    !     call *gwstress*  from *gwdrag*

    !        explicit arguments :
    !        --------------------
    !     ==== inputs ===
    !     ==== outputs ===

    !        implicit arguments :   none
    !        --------------------

    !     method.
    !     -------


    !     externals.
    !     ----------


    !     reference.
    !     ----------

    !        see ecmwf research department documentation of the "i.f.s."

    !     author.
    !     -------

    !     modifications.
    !     --------------
    !     f. lott put the new gwd on ifs      22/11/93

    !-----------------------------------------------------------------------
    USE dimensions
    USE dimphy
    USE suphec_m
    USE yoegwd

    !-----------------------------------------------------------------------

    !*       0.1   arguments
    !              ---------

    INTEGER nlon, nlev
    integer, intent(in):: ktest(nlon)
    integer, intent(in):: kkenvh(nlon)

    REAL prho(nlon,nlev+1), pstab(nlon,nlev+1), ptau(nlon,nlev+1), &
         pvph(nlon,nlev+1), pgeom1(nlon,nlev)
    REAL, INTENT(IN):: pstd(nlon)

    REAL, INTENT(IN):: psig(nlon)
    REAL, INTENT(IN):: pmea(nlon), ppic(nlon)
    REAL pdmod(nlon)

    !-----------------------------------------------------------------------

    !*       0.2   local arrays
    !              ------------
    INTEGER jl
    REAL zblock, zvar, zeff

    !-----------------------------------------------------------------------

    !*       0.3   functions
    !              ---------
    !     ------------------------------------------------------------------

    !*         1.    initialization
    !                --------------

    !*         3.1     gravity wave stress.

    DO jl = 1, klon
       IF (ktest(jl)==1) THEN

          !  effective mountain height above the blocked flow

          IF (kkenvh(jl)==klev) THEN
             zblock = 0.0
          ELSE
             zblock = (pgeom1(jl,kkenvh(jl))+pgeom1(jl,kkenvh(jl)+1))/2./rg
          END IF

          zvar = ppic(jl) - pmea(jl)
          zeff = amax1(0.,zvar-zblock)

          ptau(jl,klev+1) = prho(jl,klev+1)*gkdrag*psig(jl)*zeff**2/4./ &
               pstd(jl)*pvph(jl,klev+1)*pdmod(jl)*sqrt(pstab(jl,klev+1))

          !  too small value of stress or  low level flow include critical level
          !  or low level flow:  gravity wave stress nul.
       ELSE

          ptau(jl,klev+1) = 0.0

       END IF

    end DO

  END SUBROUTINE gwstress

end module gwstress_m
