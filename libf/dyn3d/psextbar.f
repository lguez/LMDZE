!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/psextbar.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE psextbar ( ps, psexbarxy )
      use dimens_m
      use paramet_m
      use comgeom
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c
c   Objet:
c   ------
c
c **********************************************************************
c calcul des moyennes en x et en y de (pression au sol*aire variable) ..
c **********************************************************************
c
c         ps          est un  argum. d'entree  pour le s-pg ..
c         psexbarxy   est un  argum. de sortie pour le s-pg ..
c
c   Methode:
c   --------
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
c
c                       On  a :
c
c    pbarx(i,j) = Pext(i  ,j) * ( alpha1(i  ,j) + alpha2(i,j))      +
c                 Pext(i+1,j) * ( alpha3(i+1,j) + alpha4(i+1,j) )
c     localise  au point  ... U (i,j) ...
c
c    pbary(i,j) = Pext(i,j  ) * ( alpha2(i,j  ) + alpha3(i,j  )     +
c                 Pext(i,j+1) * ( alpha1(i,j+1) + alpha4(i,j+1)  
c     localise  au point  ... V (i,j) ...
c
c  pbarxy(i,j)= Pext(i,j) *alpha2(i,j) + Pext(i+1,j) *alpha3(i+1,j) +
c               Pext(i,j+1)*alpha1(i,j+1)+ Pext(i+1,j+1)*alpha4(i+1,j+1)
c     localise  au point  ... Z (i,j) ...
c
c
c
c=======================================================================



      REAL, intent(in):: ps( ip1jmp1 )
      real psexbarxy ( ip1jm ), pext( ip1jmp1 )

      INTEGER  l, ij
c

      DO ij = 1, ip1jmp1
       pext(ij) = ps(ij) * aire(ij)
      ENDDO


      DO     5     ij = 1, ip1jm - 1
      psexbarxy( ij ) = pext(ij)*alpha2(ij) + pext(ij+1)*alpha3(ij+1) +
     *   pext(ij+iip1)*alpha1(ij+iip1) + pext(ij+iip2)*alpha4(ij+iip2)
   5  CONTINUE


c    ....  correction pour     psexbarxy( iip1,j )  ........

CDIR$ IVDEP

      DO 7 ij = iip1, ip1jm, iip1
      psexbarxy( ij ) = psexbarxy( ij - iim )
   7  CONTINUE


      RETURN
      END
