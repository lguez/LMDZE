module interfoce_lim_m

  implicit none

contains

  SUBROUTINE interfoce_lim(jour, pctsrf_new_oce, pctsrf_new_sic)

    ! lecture conditions limites
    ! Cette routine sert d'interface entre le modèle atmosphérique et
    ! un fichier de conditions aux limites.

    ! Laurent FAIRHEAD, February 2000

    USE netcdf, ONLY: nf90_nowrite
    use netcdf95, only: NF95_CLOSE, nf95_get_var, NF95_INQ_VARID, nf95_open

    integer, intent(IN):: jour ! jour \`a lire dans l'ann\'ee

    real, intent(out):: pctsrf_new_oce(:), pctsrf_new_sic(:) ! (klon)
    ! sous-maille fractionnelle

    ! Local:
    integer ncid, varid ! pour NetCDF

    ! --------------------------------------------------

    call NF95_OPEN ('limit.nc', NF90_NOWRITE, ncid)

    ! Fraction "ocean"
    call NF95_INQ_VARID(ncid, 'FOCE', varid)
    call NF95_GET_VAR(ncid, varid, pctsrf_new_oce, start = (/1, jour/))

    ! Fraction "glace de mer"
    call NF95_INQ_VARID(ncid, 'FSIC', varid)
    call NF95_GET_VAR(ncid, varid, pctsrf_new_sic, start = (/1, jour/))

    call NF95_CLOSE(ncid)

  END SUBROUTINE interfoce_lim

end module interfoce_lim_m
