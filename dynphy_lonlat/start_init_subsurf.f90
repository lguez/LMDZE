module start_init_subsurf_m

  implicit none

contains

  subroutine start_init_subsurf(pctsrf)
    
    ! From "etat0_netcdf.F", version 1.3, 2005/05/25 13:10:09

    ! On initialise les sous-surfaces. Lecture du fichier glace de
    ! terre pour fixer la fraction de terre et de glace de terre.

    use jumble, only: pi, deg_to_rad
    use netcdf, only: nf90_nowrite
    use netcdf95, only: nf95_close, nf95_gw_var, nf95_inq_varid, nf95_open, &
         nf95_find_coord

    use dimensions, only: iim, jjm, llm
    use dynetat0_chosen_m, only: day_ref
    use dynetat0_m, only: rlatu, rlonv
    use grid_change, only: dyn_phy
    use grille_m_m, only: grille_m
    use indicesol, only: is_oce, is_ter, is_lic, epsfra
    use phyetat0_m, only: masque

    REAL, intent(out):: pctsrf(:, :) ! (klon, nbsrf)
    ! "pctsrf(i, :)" is the composition of the surface at horizontal
    ! position "i".

    ! Local:
    INTEGER ncid, varid
    REAL, ALLOCATABLE:: dlon_lic(:), dlat_lic(:)
    REAL, ALLOCATABLE:: landice(:, :) ! fraction land ice
    REAL flic_tmp(iim + 1, jjm + 1) ! fraction land ice temporary
    INTEGER ji

    !---------------------------------

    print *, "Call sequence information: start_init_subsurf"
    call nf95_open("landiceref.nc", nf90_nowrite, ncid)
    call nf95_find_coord(ncid, std_name = "longitude", varid = varid)
    call nf95_gw_var(ncid, varid, dlon_lic)

    call nf95_find_coord(ncid, std_name = "latitude", varid = varid)
    call nf95_gw_var(ncid, varid, dlat_lic)

    call nf95_inq_varid(ncid, 'landice', varid)
    call nf95_gw_var(ncid, varid, landice)

    call nf95_close(ncid)

    ! Interpolation sur la grille T du mod\`ele :

    ! Si les coordonn\'ees sont en degr\'es, on les transforme :
    IF (MAXVAL(dlon_lic) > pi) THEN
       dlon_lic = dlon_lic * deg_to_rad
    ENDIF
    IF (maxval(dlat_lic) > pi) THEN 
       dlat_lic = dlat_lic * deg_to_rad
    ENDIF

    flic_tmp(:iim, :) = grille_m(dlon_lic, dlat_lic, landice, rlonv(:iim), &
         rlatu)
    flic_tmp(iim + 1, :) = flic_tmp(1, :)

    ! Passage sur la grille physique :
    pctsrf = 0.
    pctsrf(:, is_lic) = pack(flic_tmp, dyn_phy)

    ! Ad\'equation avec le maque terre-mer :
    WHERE (pctsrf(:, is_lic) < EPSFRA) pctsrf(:, is_lic) = 0.
    WHERE (masque < EPSFRA) pctsrf(:, is_lic) = 0.
    where (masque <= EPSFRA) pctsrf(:, is_ter) = masque
    where (masque > EPSFRA)
       where (pctsrf(:, is_lic) >= masque)
          pctsrf(:, is_lic) = masque
          pctsrf(:, is_ter) = 0.
       elsewhere
          pctsrf(:, is_ter) = masque - pctsrf(:, is_lic)
          where (pctsrf(:, is_ter) < EPSFRA)
             pctsrf(:, is_ter) = 0.
             pctsrf(:, is_lic) = masque
          end where
       end where
    end where

    ! Sous-surface oc\'ean et glace de mer (pour d\'emarrer on met glace
    ! de mer \`a 0) :
    pctsrf(:, is_oce) = 1. - masque
    WHERE (pctsrf(:, is_oce) < EPSFRA) pctsrf(:, is_oce) = 0.

    ! V\'erification que la somme des sous-surfaces vaut 1 :
    ji = count(abs(sum(pctsrf, dim = 2) - 1.) > EPSFRA)
    IF (ji /= 0) then
       PRINT *, 'Bad surface percentages for ', ji, 'points'
    end IF

  end subroutine start_init_subsurf

end module start_init_subsurf_m
