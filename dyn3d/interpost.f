!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/interpost.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
        subroutine interpost(q,qppm)

       use dimens_m
      use paramet_m
      use comconst
      use disvert_m
      use comgeom
       implicit none



c Arguments   
      real   q(iip1,jjp1,llm)
      real   qppm(iim,jjp1,llm)
c Local
      integer l,i,j
  
c RE-INVERSION DES NIVEAUX
c le programme ppm3d travaille avec une 3ème coordonnée inversée par rapport
c de celle du LMDZ: z=1<=>niveau max, z=llm+1<=>surface
c On passe donc des niveaux de Lin à ceux du LMDZ
           
        do l=1,llm
          do j=1,jjp1
             do i=1,iim
                 q(i,j,l)=qppm(i,j,llm-l+1)
             enddo
          enddo
         enddo
            
c BOUCLAGE EN LONGITUDE PAS EFFECTUE DANS PPM3D

         do l=1,llm
           do j=1,jjp1
            q(iip1,j,l)=q(1,j,l)
           enddo
         enddo
  
      
       return

       end
