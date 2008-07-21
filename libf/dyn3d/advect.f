!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/advect.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE advect(ucov,vcov,teta,w,massebx,masseby,du,dv,dteta,
     $     conser)

      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use comgeom
      use ener
      IMPLICIT NONE
c=======================================================================
c
c   Auteurs:  P. Le Van , Fr. Hourdin  .
c   -------
c
c   Objet:
c   ------
c
c   *************************************************************
c   .... calcul des termes d'advection vertic.pour u,v,teta,q ...
c   *************************************************************
c        ces termes sont ajoutes a du,dv,dteta et dq .
c  Modif F.Forget 03/94 : on retire q de advect
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------


c   Arguments:
c   ----------

      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm),teta(ip1jmp1,llm)
      REAL massebx(ip1jmp1,llm),masseby(ip1jm,llm),w(ip1jmp1,llm)
      REAL dv(ip1jm,llm),du(ip1jmp1,llm),dteta(ip1jmp1,llm)
      logical, intent(in):: conser

c   Local:
c   ------

      REAL uav(ip1jmp1,llm),vav(ip1jm,llm),wsur2(ip1jmp1)
      REAL unsaire2(ip1jmp1), ge(ip1jmp1)
      REAL deuxjour, ww, gt, uu, vv

      INTEGER  ij,l

      REAL      SSUM

c-----------------------------------------------------------------------
c   2. Calculs preliminaires:
c   -------------------------

      IF (conser)  THEN
         deuxjour = 2. * daysec

         DO   1  ij   = 1, ip1jmp1
         unsaire2(ij) = unsaire(ij) * unsaire(ij)
   1     CONTINUE
      END IF


c------------------  -yy ----------------------------------------------
c   .  Calcul de     u

      DO  l=1,llm
         DO    ij     = iip2, ip1jmp1
            uav(ij,l) = 0.25 * ( ucov(ij,l) + ucov(ij-iip1,l) )
         ENDDO
         DO    ij     = iip2, ip1jm
            uav(ij,l) = uav(ij,l) + uav(ij+iip1,l)
         ENDDO
         DO      ij         = 1, iip1
            uav(ij      ,l) = 0.
            uav(ip1jm+ij,l) = 0.
         ENDDO
      ENDDO

c------------------  -xx ----------------------------------------------
c   .  Calcul de     v

      DO  l=1,llm
         DO    ij   = 2, ip1jm
          vav(ij,l) = 0.25 * ( vcov(ij,l) + vcov(ij-1,l) )
         ENDDO
         DO    ij   = 1,ip1jm,iip1
          vav(ij,l) = vav(ij+iim,l)
         ENDDO
         DO    ij   = 1, ip1jm-1
          vav(ij,l) = vav(ij,l) + vav(ij+1,l)
         ENDDO
         DO    ij       = 1, ip1jm, iip1
          vav(ij+iim,l) = vav(ij,l)
         ENDDO
      ENDDO

c-----------------------------------------------------------------------

c
      DO 20 l = 1, llmm1


c       ......   calcul de  - w/2.    au niveau  l+1   .......

      DO 5   ij   = 1, ip1jmp1
      wsur2( ij ) = - 0.5 * w( ij,l+1 )
   5  CONTINUE


c     .....................     calcul pour  du     ..................

      DO 6 ij = iip2 ,ip1jm-1
      ww        = wsur2 (  ij  )     + wsur2( ij+1 ) 
      uu        = 0.5 * ( ucov(ij,l) + ucov(ij,l+1) )
      du(ij,l)  = du(ij,l)   - ww * ( uu - uav(ij, l ) )/massebx(ij, l )
      du(ij,l+1)= du(ij,l+1) + ww * ( uu - uav(ij,l+1) )/massebx(ij,l+1)
   6  CONTINUE

c     .....  correction pour  du(iip1,j,l)  ........
c     .....     du(iip1,j,l)= du(1,j,l)   .....

CDIR$ IVDEP
      DO   7  ij   = iip1 +iip1, ip1jm, iip1
      du( ij, l  ) = du( ij -iim, l  )
      du( ij,l+1 ) = du( ij -iim,l+1 )
   7  CONTINUE

c     .................    calcul pour   dv      .....................

      DO 8 ij = 1, ip1jm
      ww        = wsur2( ij+iip1 )   + wsur2( ij )
      vv        = 0.5 * ( vcov(ij,l) + vcov(ij,l+1) )
      dv(ij,l)  = dv(ij, l ) - ww * (vv - vav(ij, l ) )/masseby(ij, l )
      dv(ij,l+1)= dv(ij,l+1) + ww * (vv - vav(ij,l+1) )/masseby(ij,l+1)
   8  CONTINUE

c

c     ............................................................
c     ...............    calcul pour   dh      ...................
c     ............................................................

c                       ---z
c       calcul de  - d( teta  * w )      qu'on ajoute a   dh
c                   ...............

        DO 15 ij = 1, ip1jmp1
         ww            = wsur2(ij) * (teta(ij,l) + teta(ij,l+1) )
         dteta(ij, l ) = dteta(ij, l )  -  ww
         dteta(ij,l+1) = dteta(ij,l+1)  +  ww
  15    CONTINUE

      IF( conser)  THEN
        DO 17 ij = 1,ip1jmp1
        ge(ij)   = wsur2(ij) * wsur2(ij) * unsaire2(ij)
  17    CONTINUE
        gt       = SSUM( ip1jmp1,ge,1 )
        gtot(l)  = deuxjour * SQRT( gt/ip1jmp1 )
      END IF

  20  CONTINUE
 
      RETURN
      END
