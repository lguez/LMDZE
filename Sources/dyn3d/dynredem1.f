module dynredem1_m

  IMPLICIT NONE

contains

  SUBROUTINE dynredem1(vcov, ucov, teta, q, masse, ps, itau)

    ! From dyn3d/dynredem.F, version 1.2, 2004/06/22 11:45:30
    ! Ecriture du fichier de red\'emarrage au format NetCDF

    USE dimens_m, ONLY: iim, jjm, llm, nqmx
    use dynredem0_m, only: ncid
    USE iniadvtrac_m, ONLY: tname
    use netcdf, only: nf90_write
    use netcdf95, only: nf95_close, nf95_inq_varid, nf95_open, nf95_put_var
    use nr_util, only: assert

    REAL, INTENT(IN):: vcov(:, :, :) ! (iim + 1, jjm, llm)
    REAL, INTENT(IN):: ucov(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, INTENT(IN):: teta(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, INTENT(IN):: q(:, :, :, :) ! (iim + 1, jjm + 1, llm, nqmx)
    REAL, INTENT(IN):: masse(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, INTENT(IN):: ps(:, :) ! (iim + 1, jjm + 1)
    INTEGER, INTENT(IN):: itau

    ! Local:
    INTEGER varid, iq

    !---------------------------------------------------------

    PRINT *, 'Call sequence information: dynredem1'

    call assert((/size(vcov, 1), size(ucov, 1), size(teta, 1), size(q, 1), &
         size(masse, 1), size(ps, 1)/) == iim + 1, "dynredem1 iim")
    call assert((/size(vcov, 2) + 1, size(ucov, 2), size(teta, 2), size(q, 2), &
         size(masse, 2), size(ps, 2)/) == jjm + 1, "dynredem1 jjm")
    call assert((/size(vcov, 3), size(ucov, 3), size(teta, 3), size(q, 3), &
         size(masse, 3)/) == llm, "dynredem1 llm")
    call assert(size(q, 4) == nqmx, "dynredem1 nqmx")

    ! \'Ecriture/extension de la coordonn\'ee temps
    call nf95_inq_varid(ncid, 'temps', varid)
    call nf95_put_var(ncid, varid, values = 0.)

    ! R\'ecriture du tableau de contr\^ole, "itaufin" n'est pas d\'efini quand
    ! on passe dans "dynredem0"
    call nf95_inq_varid(ncid, 'controle', varid)
    call nf95_put_var(ncid, varid, real(itau), start=(/31/))

    ! \'Ecriture des champs

    call nf95_inq_varid(ncid, 'ucov', varid)
    call nf95_put_var(ncid, varid, ucov)

    call nf95_inq_varid(ncid, 'vcov', varid)
    call nf95_put_var(ncid, varid, vcov)

    call nf95_inq_varid(ncid, 'teta', varid)
    call nf95_put_var(ncid, varid, teta)

    DO iq = 1, nqmx
       call nf95_inq_varid(ncid, tname(iq), varid)
       call nf95_put_var(ncid, varid, q(:, :, :, iq))
    END DO

    call nf95_inq_varid(ncid, 'masse', varid)
    call nf95_put_var(ncid, varid, masse)

    call nf95_inq_varid(ncid, 'ps', varid)
    call nf95_put_var(ncid, varid, ps)

    call nf95_close(ncid)

  END SUBROUTINE dynredem1

end module dynredem1_m
