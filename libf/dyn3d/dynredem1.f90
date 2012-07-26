module dynredem1_m

  IMPLICIT NONE

contains

  SUBROUTINE dynredem1(fichnom, vcov, ucov, teta, q, masse, ps, itau)

    ! From dyn3d/dynredem.F, v 1.2 2004/06/22 11:45:30
    ! Ecriture du fichier de red�marrage au format NetCDF

    USE dimens_m, ONLY : iim, jjm, llm, nqmx
    USE iniadvtrac_m, ONLY : tname
    use netcdf, only: nf90_write
    use netcdf95, only: nf95_close, nf95_inq_varid, nf95_open, nf95_put_var

    CHARACTER(len=*), INTENT (IN) :: fichnom
    REAL, INTENT (IN) :: vcov(iim + 1, jjm, llm), ucov(iim+1, jjm+1, llm)
    REAL, INTENT (IN) :: teta(iim+1, jjm+1, llm)
    REAL, INTENT (IN) :: q(iim+1, jjm+1, llm, nqmx)
    REAL, INTENT (IN) :: ps(iim+1, jjm+1), masse(iim+1, jjm+1, llm)
    INTEGER, INTENT (IN) :: itau

    ! Variables local to the procedure:
    INTEGER ncid, varid
    INTEGER iq
    INTEGER:: nb = 0

    !---------------------------------------------------------

    PRINT *, 'Call sequence information: dynredem1'

    call nf95_open(fichnom, nf90_write, ncid)

    ! �criture/extension de la coordonn�e temps
    nb = nb + 1
    call nf95_inq_varid(ncid, 'temps', varid)
    call nf95_put_var(ncid, varid, values=0., start=(/nb/))
    PRINT *, "Enregistrement pour nb = ", nb

    ! R�criture du tableau de contr�le, "itaufin" n'est pas d�fini quand
    ! on passe dans "dynredem0"
    call nf95_inq_varid(ncid, 'controle', varid)
    call nf95_put_var(ncid, varid, real(itau), start=(/31/))

    ! �criture des champs

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
