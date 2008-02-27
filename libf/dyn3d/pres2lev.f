!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/pres2lev.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
c******************************************************
      SUBROUTINE   pres2lev(varo,varn,lmo,lmn,po,pn,
     &                      ni,nj)
c
c interpolation lineaire pour passer
c a une nouvelle discretisation verticale pour
c les variables de GCM
c Francois Forget (01/1995)
c 
c MOdif remy roca 12/97 pour passer de pres2sig
c**********************************************************

      IMPLICIT NONE

c   Declarations:
c ==============
c
c  ARGUMENTS
c  """""""""

       INTEGER lmo ! dimensions ancienne couches (input)
       INTEGER lmn ! dimensions nouvelle couches (input)
       INTEGER lmomx ! dimensions ancienne couches (input)
       INTEGER lmnmx ! dimensions nouvelle couches (input)

       parameter(lmomx=10000,lmnmx=10000)

        real po(lmo)! niveau de pression en millibars
        integer ni,nj
        real pn(ni,nj,lmn) ! niveau de pression en pascals

       INTEGER i,j,Nhoriz ! nombre de point horizontale (input)

       REAL varo(ni,nj,lmo) ! var dans l'ancienne grille (input)
       REAL varn(ni,nj,lmn) ! var dans la nouvelle grille (output)

       real zvaro(lmomx),zpo(lmomx)

c Autres variables
c """"""""""""""""
       INTEGER n, ln ,lo 
       REAL coef

c run
c ====
        do i=1,ni
        do j=1,nj
c a chaque point de grille correspond un nouveau sigma old
c qui vaut pres(l)/ps(i,j)
           do lo=1,lmo
              zpo(lo)=po(lmo+1-lo)
              zvaro(lo)=varo(i,j,lmo+1-lo)
           enddo
        
           do ln=1,lmn
              if (pn(i,j,ln).ge.zpo(1))then
                 varn(i,j,ln) =  zvaro(1)
              else if (pn(i,j,ln).le.zpo(lmo)) then
                 varn(i,j,ln) =  zvaro(lmo)
              else
                 do lo=1,lmo-1 
                    if ( (pn(i,j,ln).le.zpo(lo)).and.
     &                 (pn(i,j,ln).gt.zpo(lo+1)) )then
                       coef=(pn(i,j,ln)-zpo(lo))
     &                 /(zpo(lo+1)-zpo(lo))
                       varn(i,j,ln)=zvaro(lo)
     &                 +coef*(zvaro(lo+1)-zvaro(lo))
c       print*,'pn(',ln,')=',pn(i,j,ln),varn(i,j,ln)
                    end if
                 enddo           
              endif
           enddo

        enddo
        enddo
      return
      end    
