module start_inter_3d_m

  IMPLICIT NONE

contains

  subroutine start_inter_3d(varname, lon_in2, lat_in2, pls, var3d)

    ! This procedure gets a 3D variable from a file and interpolates it.

    use netcdf, only: nf90_nowrite
    use netcdf95, only: nf95_open, nf95_close, nf95_get_var, nf95_inq_varid
    use jumble, only: assert_eq
    use numer_rec_95, only: spline, splint

    use inter_barxy_m, only: inter_barxy
    use gr_int_dyn_m, only: gr_int_dyn
    use conf_dat3d_m, only: conf_dat3d
    use start_init_dyn_m, only: IML_DYN, JML_DYN, LLM_DYN, LAT_INI, &
         LEVDYN_INI, LON_INI

    CHARACTER(len=*), intent(in):: varname
    REAL, intent(in):: lon_in2(:) ! (iml) longitude, in rad
    REAL, intent(in):: lat_in2(:) ! latitude, in rad
    REAL, intent(in):: pls(:, :, :) ! (iml, jml, lml)
    REAL, intent(out):: var3d(:, :, :) ! (iml, jml, lml)

    ! LOCAL:
    INTEGER iml, jml, lml, ncid, varid
    INTEGER ii, ij, il
    REAL lon_rad(iml_dyn), lat_rad(jml_dyn)
    REAL lev_dyn(llm_dyn)
    REAL var_tmp2d(size(lon_in2)-1, size(pls, 2))
    real var_tmp3d(size(lon_in2), size(pls, 2), llm_dyn)
    REAL ax(llm_dyn), ay(llm_dyn), yder(llm_dyn)
    real var_ana3d(iml_dyn, jml_dyn, llm_dyn)

    !--------------------------------

    print *, "Call sequence information: start_inter_3d"

    iml = assert_eq(size(pls, 1), size(lon_in2), size(var3d, 1), &
         "start_inter_3d iml")
    jml = assert_eq(size(pls, 2), size(var3d, 2), "start_inter_3d jml")
    lml = assert_eq(size(pls, 3), size(var3d, 3), "start_inter_3d lml")

    print *, "iml = ", iml, ", jml = ", jml
    print *, "varname = ", varname
    call nf95_open('ECDYN.nc', nf90_nowrite, ncid)
    call nf95_inq_varid(ncid, varname, varid)
    call nf95_get_var(ncid, varid, var_ana3d)
    call nf95_close(ncid)
    CALL conf_dat3d(lon_ini, lat_ini, levdyn_ini, lon_rad, lat_rad, lev_dyn, &
         var_ana3d)

    DO il = 1, llm_dyn
       CALL inter_barxy(lon_rad, lat_rad(:jml_dyn-1), var_ana3d(:, :, il), &
            lon_in2(:iml-1), lat_in2, var_tmp2d) 
       var_tmp3d(:, :, il) = gr_int_dyn(var_tmp2d)
    ENDDO

    ! Pour l'interpolation verticale, on interpole du haut de l'atmosphère
    ! vers le sol :
    ax = lev_dyn(llm_dyn:1:-1)
    DO ij=1, jml
       DO ii=1, iml-1
          ay = var_tmp3d(ii, ij, llm_dyn:1:-1)
          yder = SPLINE(ax, ay)
          do il=1, lml
             var3d(ii, ij, il) = SPLINT(ax, ay, yder, pls(ii, ij, il))
          END do
       ENDDO
    ENDDO
    var3d(iml, :, :) = var3d(1, :, :) 

  END subroutine start_inter_3d

end module start_inter_3d_m
