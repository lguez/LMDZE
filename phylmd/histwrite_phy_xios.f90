module histwrite_phy_xios_m

  ! Libraries:
  use xios, only: xios_send_field

  use gr_phy_write_m, only: gr_phy_write

  implicit none

  INTERFACE histwrite_phy_xios
     MODULE PROCEDURE histwrite2d_phy, histwrite3d_phy
  end INTERFACE histwrite_phy_xios

  private
  public histwrite_phy_xios

contains

  subroutine histwrite2d_phy(var, field)

    character(len=*), intent(in):: var
    real, intent(in):: field(:)

    !-----------------------------------------------------------

    CALL xios_send_field(var, gr_phy_write(field))

  end subroutine histwrite2d_phy

  !*************************************************************

  subroutine histwrite3d_phy(var, field)

    character(len=*), intent(in):: var
    real, intent(in):: field(:, :)

    !-----------------------------------------------------------

    CALL xios_send_field(var, gr_phy_write(field))

  end subroutine histwrite3d_phy

end module histwrite_phy_xios_m
