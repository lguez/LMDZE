!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/cvparam.h,v 1.1.1.1 2004/05/19 12:53:08 lmdzadmin Exp $
!
c------------------------------------------------------------
c Parameters for convectL:
c (includes - microphysical parameters, 
c			- parameters that control the rate of approach 
c               to quasi-equilibrium)
c			- noff & minorig (previously in input of convect1)
c------------------------------------------------------------

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

      COMMON /cvparam/ noff, minorig, nl, nlp, nlm
     :                ,elcrit, tlcrit
     :                ,entp, sigs, sigd
     :                ,omtrain, omtsnow, coeffr, coeffs
     :                ,dtmax, cu, betad, alpha, damp, delta

