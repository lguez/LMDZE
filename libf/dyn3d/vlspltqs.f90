!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/vlspltqs.F,v 1.2 2005/02/24 12:16:57 fairhead Exp $
!
       SUBROUTINE vlspltqs ( q,pente_max,masse,w,pbaru,pbarv,pdt, &
                                        p,pk,teta                 )
!
!     Auteurs:   P.Le Van, F.Hourdin, F.Forget, F.Codron
!
!    ********************************************************************
!          Shema  d'advection " pseudo amont " .
!      + test sur humidite specifique: Q advecte< Qsat aval
!                   (F. Codron, 10/99)
!    ********************************************************************
!     q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....
!
!     pente_max facteur de limitation des pentes: 2 en general
!                                                0 pour un schema amont
!     pbaru,pbarv,w flux de masse en u ,v ,w
!     pdt pas de temps
!
!     teta temperature potentielle, p pression aux interfaces,
!     pk exner au milieu des couches necessaire pour calculer Qsat
!   --------------------------------------------------------------------
      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use logic
      IMPLICIT NONE
!

!
!   Arguments:
!   ----------
      REAL masse(ip1jmp1,llm),pente_max
      REAL, intent(in):: pbaru( ip1jmp1,llm ),pbarv( ip1jm,llm)
      REAL q(ip1jmp1,llm)
      REAL w(ip1jmp1,llm)
      real, intent(in):: pdt
      REAL, intent(in):: p(ip1jmp1,llmp1)
      real teta(ip1jmp1,llm),pk(ip1jmp1,llm)
!
!      Local
!   ---------
!
      INTEGER i,ij,l,j,ii
!
      REAL qsat(ip1jmp1,llm)
      REAL zm(ip1jmp1,llm)
      REAL mu(ip1jmp1,llm)
      REAL mv(ip1jm,llm)
      REAL mw(ip1jmp1,llm+1)
      REAL zq(ip1jmp1,llm)
      REAL temps1,temps2,temps3
      REAL zzpbar, zzw
      LOGICAL testcpu
      SAVE testcpu
      SAVE temps1,temps2,temps3

      REAL qmin,qmax
      DATA qmin,qmax/0.,1.e33/
      DATA testcpu/.false./
      DATA temps1,temps2,temps3/0.,0.,0./

!--pour rapport de melange saturant--

      REAL rtt,retv,r2es,r3les,r3ies,r4les,r4ies,play
      REAL ptarg,pdelarg,foeew,zdelta
      REAL tempe(ip1jmp1)

!    fonction psat(T)

       FOEEW ( PTARG,PDELARG ) = EXP ( &
                (R3LES*(1.-PDELARG)+R3IES*PDELARG) * (PTARG-RTT) &
       / (PTARG-(R4LES*(1.-PDELARG)+R4IES*PDELARG)) )

        r2es  = 380.11733
        r3les = 17.269
        r3ies = 21.875
        r4les = 35.86
        r4ies = 7.66
        retv = 0.6077667
        rtt  = 273.16

!-- Calcul de Qsat en chaque point
!-- approximation: au milieu des couches play(l)=(p(l)+p(l+1))/2
!   pour eviter une exponentielle.
        DO l = 1, llm
         DO ij = 1, ip1jmp1
          tempe(ij) = teta(ij,l) * pk(ij,l) /cpp
         ENDDO
         DO ij = 1, ip1jmp1
          zdelta = MAX( 0., SIGN(1., rtt - tempe(ij)) )
          play   = 0.5*(p(ij,l)+p(ij,l+1))
          qsat(ij,l) = MIN(0.5, r2es* FOEEW(tempe(ij),zdelta) / play )
          qsat(ij,l) = qsat(ij,l) / ( 1. - retv * qsat(ij,l) )
         ENDDO
        ENDDO

        zzpbar = 0.5 * pdt
        zzw    = pdt
      DO l=1,llm
        DO ij = iip2,ip1jm
            mu(ij,l)=pbaru(ij,l) * zzpbar
         ENDDO
         DO ij=1,ip1jm
            mv(ij,l)=pbarv(ij,l) * zzpbar
         ENDDO
         DO ij=1,ip1jmp1
            mw(ij,l)=w(ij,l) * zzw
         ENDDO
      ENDDO

      DO ij=1,ip1jmp1
         mw(ij,llm+1)=0.
      ENDDO

      CALL SCOPY(ijp1llm,q,1,zq,1)
      CALL SCOPY(ijp1llm,masse,1,zm,1)

!      call minmaxq(zq,qmin,qmax,'avant vlxqs     ')
      call vlxqs(zq,pente_max,zm,mu,qsat)


!     call minmaxq(zq,qmin,qmax,'avant vlyqs     ')

      call vlyqs(zq,pente_max,zm,mv,qsat)


!      call minmaxq(zq,qmin,qmax,'avant vlz     ')

      call vlz(zq,pente_max,zm,mw)


!     call minmaxq(zq,qmin,qmax,'avant vlyqs     ')
!     call minmaxq(zm,qmin,qmax,'M avant vlyqs     ')

      call vlyqs(zq,pente_max,zm,mv,qsat)


!     call minmaxq(zq,qmin,qmax,'avant vlxqs     ')
!     call minmaxq(zm,qmin,qmax,'M avant vlxqs     ')

      call vlxqs(zq,pente_max,zm,mu,qsat)

!     call minmaxq(zq,qmin,qmax,'apres vlxqs     ')
!     call minmaxq(zm,qmin,qmax,'M apres vlxqs     ')


      DO l=1,llm
         DO ij=1,ip1jmp1
           q(ij,l)=zq(ij,l)
         ENDDO
         DO ij=1,ip1jm+1,iip1
            q(ij+iim,l)=q(ij,l)
         ENDDO
      ENDDO

      RETURN
      END
