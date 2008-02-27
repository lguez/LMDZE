!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/cvparam3.h,v 1.1.1.1 2004/05/19 12:53:09 lmdzadmin Exp $
!
c------------------------------------------------------------
c Parameters for convectL, iflag_con=3:
c (includes - microphysical parameters, 
c			- parameters that control the rate of approach 
c               to quasi-equilibrium)
c			- noff & minorig (previously in input of convect1)
c------------------------------------------------------------

      integer noff, minorig, nl, nlp, nlm
      real sigd, spfac
cIM cf. FH : pour compatibilite avec conema3 TEMPORAIRE   real pbcrit, ptcrit, epmax
      real pbcrit, ptcrit
      real omtrain
      real dtovsh, dpbase, dttrig
      real dtcrit, tau, beta, alpha
      real delta
      real betad

      COMMON /cvparam3/  noff, minorig, nl, nlp, nlm
     :                ,  sigd, spfac
cIM cf. FH : pour compatibilite avec conema3 TEMPORAIRE  :                ,pbcrit, ptcrit, epmax
     :                ,pbcrit, ptcrit
     :                ,omtrain
     :                ,dtovsh, dpbase, dttrig
     :                ,dtcrit, tau, beta, alpha, delta, betad

