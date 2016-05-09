module histwrite_phy_m

  use clesphys, only: ok_instan
  use gr_phy_write_m, only: gr_phy_write
  use histwrite_m, only: histwrite
  use ini_histins_m, only: nid_ins
  use time_phylmdz, only: itau_w

  implicit none

  INTEGER:: itap = 0 ! number of calls to "physiq"

  INTERFACE histwrite_phy
     MODULE PROCEDURE histwrite2d_phy, histwrite3d_phy
  end INTERFACE histwrite_phy

  private
  public itap, histwrite_phy

contains

  subroutine histwrite2d_phy(var, field)

    character(len=*), intent(in):: var
    real, intent(in):: field(:)

    !-----------------------------------------------------------

    IF (ok_instan) CALL histwrite(nid_ins, var, itau_w, gr_phy_write(field))

  end subroutine histwrite2d_phy

  !*************************************************************

  subroutine histwrite3d_phy(var, field)

    character(len=*), intent(in):: var
    real, intent(in):: field(:, :)

    !-----------------------------------------------------------

    IF (ok_instan) CALL histwrite(nid_ins, var, itau_w, gr_phy_write(field))

  end subroutine histwrite3d_phy

end module histwrite_phy_m
