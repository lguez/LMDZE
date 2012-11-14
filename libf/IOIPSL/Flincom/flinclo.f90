module flinclo_m

  IMPLICIT NONE

contains

  SUBROUTINE flinclo (fid_in)

    use flininfo_m, only: NCFILEOPEN, NCIDS
    USE netcdf, ONLY: nf90_close

    INTEGER:: fid_in

    INTEGER:: iret

    !---------------------------------------------------------------------

    iret = NF90_CLOSE (ncids(fid_in))
    ncfileopen(fid_in) = .FALSE.

  END SUBROUTINE flinclo

end module flinclo_m
