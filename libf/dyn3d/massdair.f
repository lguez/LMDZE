!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/massdair.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
      SUBROUTINE massdair( p, masse )
c
c *********************************************************************
c       ....  Calcule la masse d'air  dans chaque maille   ....
c *********************************************************************
c
c    Auteurs : P. Le Van , Fr. Hourdin  .
c   ..........
c
c  ..    p                      est  un argum. d'entree pour le s-pg ...
c  ..  masse                    est un  argum.de sortie pour le s-pg ...
c     
c  ....  p est defini aux interfaces des llm couches   .....
c
      use dimens_m
      use paramet_m
      use comconst
      use comgeom, only: airesurg
      IMPLICIT NONE
c
c
c  .....   arguments  ....
c
      REAL,intent(in):: p(ip1jmp1,llmp1)
      real masse(ip1jmp1,llm)

c   ....  Variables locales  .....

      INTEGER l,ij
      REAL massemoyn, massemoys

      REAL SSUM
c
c
c   Methode pour calculer massebx et masseby .
c   ----------------------------------------
c
c    A chaque point scalaire P (i,j) est affecte 4 coefficients d'aires
c       alpha1(i,j)  calcule  au point ( i+1/4,j-1/4 )
c       alpha2(i,j)  calcule  au point ( i+1/4,j+1/4 )
c       alpha3(i,j)  calcule  au point ( i-1/4,j+1/4 )
c       alpha4(i,j)  calcule  au point ( i-1/4,j-1/4 )
c
c    Avec  alpha1(i,j) = aire(i+1/4,j-1/4)/ aire(i,j)        
c
c    N.B .  Pour plus de details, voir s-pg  ...  iniconst ...
c
c
c
c   alpha4 .         . alpha1    . alpha4
c    (i,j)             (i,j)       (i+1,j)
c
c             P .        U .          . P
c           (i,j)       (i,j)         (i+1,j)
c
c   alpha3 .         . alpha2    .alpha3 
c    (i,j)              (i,j)     (i+1,j)
c
c             V .        Z .          . V
c           (i,j)
c
c   alpha4 .         . alpha1    .alpha4
c   (i,j+1)            (i,j+1)   (i+1,j+1) 
c
c             P .        U .          . P
c          (i,j+1)                    (i+1,j+1)
c
c
c
c                       On  a :
c
c    massebx(i,j) = masse(i  ,j) * ( alpha1(i  ,j) + alpha2(i,j))   +
c                   masse(i+1,j) * ( alpha3(i+1,j) + alpha4(i+1,j) )
c     localise  au point  ... U (i,j) ...
c
c    masseby(i,j) = masse(i,j  ) * ( alpha2(i,j  ) + alpha3(i,j  )  +
c                   masse(i,j+1) * ( alpha1(i,j+1) + alpha4(i,j+1)  
c     localise  au point  ... V (i,j) ...
c
c
c=======================================================================

CC      print *, "Call sequence information: massdair"
      DO   100    l = 1 , llm
c
        DO    ij     = 1, ip1jmp1
         masse(ij,l) = airesurg(ij) * ( p(ij,l) - p(ij,l+1) )
        ENDDO
c
        DO   ij = 1, ip1jmp1,iip1
         masse(ij+ iim,l) = masse(ij,l)
        ENDDO
100   CONTINUE
c
      RETURN
      END
