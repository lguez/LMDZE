module trans_buff_m

  implicit none

contains

  SUBROUTINE trans_buff (ox, sx, oy, sy, oz, sz, xsz, ysz, zsz, v3d, sl, v1d)
    !- This subroutine extracts from the full 3D variable the slab of
    !- data that will be used later. Perhaps there are hardware routines
    !- for this task on some computers. This routine will be obsolete in
    !- a F90 environnement

    !- INPUT
    !- ox  : Origin of slab of data in X
    !- sx  : Size of slab in X
    !- oy  : Origin of slab of data in Y
    !- sy  : Size of slab in Y
    !- oz  : Origin of slab of data in Z
    !- sz  : Size of slab in Z
    !- xsz, ysz, zsz : 3 sizes of full variable v3d
    !- v3d : The full 3D variable
    !- sl  : size of variable v1d
    !- v1d : The 1D variable containing the slab

    INTEGER :: ox, sx, oy, sy, oz, sz
    INTEGER :: xsz, ysz, zsz
    INTEGER :: sl
    REAL :: v3d(xsz, ysz, zsz)
    REAL :: v1d(sl)

    !---------------------------------------------------------------------

    V1d(:sx*sy*sz) = reshape(V3d(ox:ox-1+sx, oy:oy-1+sy, oz:oz-1+sz), &
         SHAPE=(/sx*sy*sz/))

  END SUBROUTINE trans_buff

end module trans_buff_m
