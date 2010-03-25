module dynredem1_m

  IMPLICIT NONE

contains

  SUBROUTINE dynredem1(fichnom, vcov, ucov, teta, q, masse, ps, itau)

    ! From dyn3d/dynredem.F, v 1.2 2004/06/22 11:45:30
    ! Ecriture du fichier de redémarrage au format NetCDF

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
    INTEGER nid, nvarid
    INTEGER iq
    INTEGER:: nb = 0

    !---------------------------------------------------------

    PRINT *, 'Call sequence information: dynredem1'

    call nf95_open(fichnom, nf90_write, nid)

    ! Écriture/extension de la coordonnée temps
    nb = nb + 1
    call nf95_inq_varid(nid, 'temps', nvarid)
    call nf95_put_var(nid, nvarid, values=0., start=(/nb/))
    PRINT *, "Enregistrement pour nb = ", nb

    ! Récriture du tableau de contrôle, "itaufin" n'est pas défini quand
    ! on passe dans "dynredem0"
    call nf95_inq_varid(nid, 'controle', nvarid)
    call nf95_put_var(nid, nvarid, real(itau), start=(/31/))

    ! Écriture des champs

    call nf95_inq_varid(nid, 'ucov', nvarid)
    call nf95_put_var(nid, nvarid, ucov)

    call nf95_inq_varid(nid, 'vcov', nvarid)
    call nf95_put_var(nid, nvarid, vcov)

    call nf95_inq_varid(nid, 'teta', nvarid)
    call nf95_put_var(nid, nvarid, teta)

    DO iq = 1, nqmx
       call nf95_inq_varid(nid, tname(iq), nvarid)
       call nf95_put_var(nid, nvarid, q(:, :, :, iq))
    END DO

    call nf95_inq_varid(nid, 'masse', nvarid)
    call nf95_put_var(nid, nvarid, masse)

    call nf95_inq_varid(nid, 'ps', nvarid)
    call nf95_put_var(nid, nvarid, ps)

    call nf95_close(nid)

  END SUBROUTINE dynredem1

end module dynredem1_m
