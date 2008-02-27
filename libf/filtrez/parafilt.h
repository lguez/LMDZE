!
! $Header: /home/cvsroot/LMDZ4/libf/filtrez/parafilt.h,v 1.1.1.1 2004/05/19 12:53:09 lmdzadmin Exp $
!
        INTEGER nfilun, nfilus, nfilvn, nfilvs
      PARAMETER (nfilun=6, nfilus=5, nfilvn=5, nfilvs=5)

c      La premiere fois que  le Gcm  rentrera  dans le Filtre ,
c      il indiquera  les bonnes valeurs  de  nfilun , nflius, nfilvn  et 
c      nfilvs  a  mettre .  Il suffira alors de changer ces valeurs dans
c      Parameter  ci-dessus  et de relancer  le  run.
