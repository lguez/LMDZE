module ctherm

  implicit none

  integer, save:: iflag_thermals,nsplit_thermals

  !  A inclure eventuellement dans les fichiers de configuration

  real:: r_aspect_thermals = 4.
  real:: l_mix_thermals = 10.
  real:: tho_thermals = 0.
  integer:: w2di_thermals = 0

end module ctherm
