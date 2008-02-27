!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/interpre.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
       subroutine interpre(q,qppm,w,fluxwppm,masse,
     s            apppm,bpppm,massebx,masseby,pbaru,pbarv,
     s            unatppm,vnatppm,psppm)

       use dimens_m
      use paramet_m
      use comconst
      use comvert
      use conf_gcm_m
      use logic
      use comgeom
      use temps
      use ener
      use comdissip
       implicit none

c---------------------------------------------------
c Arguments     
      real   apppm(llm+1),bpppm(llm+1)
      real   q(iip1,jjp1,llm),qppm(iim,jjp1,llm)
c---------------------------------------------------
      real   masse(iip1,jjp1,llm) 
      real   massebx(iip1,jjp1,llm),masseby(iip1,jjm,llm)      
      real   w(iip1,jjp1,llm+1)
      real   fluxwppm(iim,jjp1,llm)
      real   pbaru(iip1,jjp1,llm )
      real   pbarv(iip1,jjm,llm)
      real   unatppm(iim,jjp1,llm)
      real   vnatppm(iim,jjp1,llm)
      real   psppm(iim,jjp1)
c---------------------------------------------------
c Local
      real   vnat(iip1,jjp1,llm)
      real   unat(iip1,jjp1,llm)
      real   fluxw(iip1,jjp1,llm)
      real   smass(iip1,jjp1)
c----------------------------------------------------
      integer l,ij,i,j

c       CALCUL DE LA PRESSION DE SURFACE
c       Les coefficients ap et bp sont passés en common 
c       Calcul de la pression au sol en mb optimisée pour 
c       la vectorialisation
                   
         do j=1,jjp1
             do i=1,iip1
                smass(i,j)=0.
             enddo
         enddo

         do l=1,llm
             do j=1,jjp1
                 do i=1,iip1
                    smass(i,j)=smass(i,j)+masse(i,j,l)
                 enddo
             enddo
         enddo
      
         do j=1,jjp1
             do i=1,iim
                 psppm(i,j)=smass(i,j)/aire_2d(i,j)*g*0.01
             end do
         end do                        
       
c RECONSTRUCTION DES CHAMPS CONTRAVARIANTS
c Le programme ppm3d travaille avec les composantes
c de vitesse et pas les flux, on doit donc passer de l'un à l'autre
c Dans le même temps, on fait le changement d'orientation du vent en v
      do l=1,llm
          do j=1,jjm
              do i=1,iip1
                  vnat(i,j,l)=-pbarv(i,j,l)/masseby(i,j,l)*cv_2d(i,j)             
              enddo
          enddo
          do  i=1,iim
          vnat(i,jjp1,l)=0.
          enddo
          do j=1,jjp1
              do i=1,iip1
                  unat(i,j,l)=pbaru(i,j,l)/massebx(i,j,l)*cu_2d(i,j)
              enddo
          enddo
      enddo
              
c CALCUL DU FLUX MASSIQUE VERTICAL
c Flux en l=1 (sol) nul
      fluxw=0.        
      do l=1,llm
           do j=1,jjp1
              do i=1,iip1              
               fluxw(i,j,l)=w(i,j,l)*g*0.01/aire_2d(i,j)
C               print*,i,j,l,'fluxw(i,j,l)=',fluxw(i,j,l),
C     c                      'w(i,j,l)=',w(i,j,l)
              enddo
           enddo
      enddo
      
c INVERSION DES NIVEAUX
c le programme ppm3d travaille avec une 3ème coordonnée inversée par rapport
c de celle du LMDZ: z=1<=>niveau max, z=llm+1<=>surface
c On passe donc des niveaux du LMDZ à ceux de Lin
     
      do l=1,llm+1
          apppm(l)=ap(llm+2-l)
          bpppm(l)=bp(llm+2-l)         
      enddo 
     
      do l=1,llm
          do j=1,jjp1
             do i=1,iim     
                 unatppm(i,j,l)=unat(i,j,llm-l+1)
                 vnatppm(i,j,l)=vnat(i,j,llm-l+1)
                 fluxwppm(i,j,l)=fluxw(i,j,llm-l+1)
                 qppm(i,j,l)=q(i,j,llm-l+1)                              
             enddo
          enddo                                
      enddo
   
      return
      end






