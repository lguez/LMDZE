module cvparam3

  ! From LMDZ4/libf/phylmd/cvparam3.h,v 1.1.1.1 2004/05/19 12:53:09
  ! Parameters for convectL, iflag_con=3:
  ! (includes - microphysical parameters, 
  !			- parameters that control the rate of approach 
  !               to quasi-equilibrium)
  !			- noff & minorig (previously in input of convect1)

  implicit none

  integer noff, minorig, nl, nlp, nlm
  real sigd, spfac
  real pbcrit, ptcrit
  real omtrain
  real dtovsh, dpbase, dttrig
  real dtcrit, tau, beta, alpha
  real delta
  real betad

end module cvparam3
