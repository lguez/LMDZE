module cvparam

  ! From LMDZ4/libf/phylmd/cvparam.h,v 1.1.1.1 2004/05/19 12:53:08

  ! Parameters for convectL:
  ! (includes - microphysical parameters, 
  !			- parameters that control the rate of approach 
  !               to quasi-equilibrium)
  !			- noff & minorig (previously in input of convect1)

  implicit none

  integer noff, minorig, nl, nlp, nlm
  real elcrit, tlcrit
  real entp
  real sigs, sigd
  real omtrain, omtsnow, coeffr, coeffs
  real dtmax
  real cu
  real betad
  real alpha, damp
  real delta

end module cvparam
