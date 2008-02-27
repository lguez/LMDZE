!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/addfi.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE addfi(nq, pdt, 
     S          pucov, pvcov, pteta, pq   , pps ,
     S          pdufi, pdvfi, pdhfi,pdqfi, pdpfi  )
      use dimens_m
      use paramet_m
      use comconst
      use comgeom
      use serre
      IMPLICIT NONE
c
c=======================================================================
c
c    Addition of the physical tendencies
c
c    Interface :
c    -----------
c
c      Input :
c      -------
c      pdt                    time step of integration
c      pucov(ip1jmp1,llm)     first component of the covariant velocity
c      pvcov(ip1ip1jm,llm)    second component of the covariant velocity
c      pteta(ip1jmp1,llm)     potential temperature
c      pts(ip1jmp1,llm)       surface temperature
c      pdufi(ip1jmp1,llm)     |
c      pdvfi(ip1jm,llm)       |   respective
c      pdhfi(ip1jmp1)         |      tendencies
c      pdtsfi(ip1jmp1)        |
c
c      Output :
c      --------
c      pucov
c      pvcov
c      ph
c      pts
c
c
c=======================================================================
c
c-----------------------------------------------------------------------
c
c    0.  Declarations :
c    ------------------
c
c
c    Arguments :
c    -----------
c
      INTEGER nq

      REAL pdt
c
      REAL pvcov(ip1jm,llm),pucov(ip1jmp1,llm)
      REAL pteta(ip1jmp1,llm),pq(ip1jmp1,llm,nq),pps(ip1jmp1)
c
      REAL pdvfi(ip1jm,llm),pdufi(ip1jmp1,llm)
      REAL pdqfi(ip1jmp1,llm,nq),pdhfi(ip1jmp1,llm),pdpfi(ip1jmp1)
c
c    Local variables :
c    -----------------
c
      REAL xpn(iim),xps(iim),tpn,tps
      INTEGER j,k,iq,ij
      REAL qtestw, qtestt
      PARAMETER ( qtestw = 1.0e-15 )
      PARAMETER ( qtestt = 1.0e-40 )

      REAL SSUM
c
c-----------------------------------------------------------------------

      !!print *, "Call sequence information: addfi"

      DO k = 1,llm
         DO j = 1,ip1jmp1
            pteta(j,k)= pteta(j,k) + pdhfi(j,k) * pdt
         ENDDO
      ENDDO

      DO  k    = 1, llm
       DO  ij   = 1, iim
         xpn(ij) = aire(   ij   ) * pteta(  ij    ,k)
         xps(ij) = aire(ij+ip1jm) * pteta(ij+ip1jm,k)
       ENDDO
       tpn      = SSUM(iim,xpn,1)/ apoln
       tps      = SSUM(iim,xps,1)/ apols

       DO ij   = 1, iip1
         pteta(   ij   ,k)  = tpn
         pteta(ij+ip1jm,k)  = tps
       ENDDO
      ENDDO
c

      DO k = 1,llm
         DO j = iip2,ip1jm
            pucov(j,k)= pucov(j,k) + pdufi(j,k) * pdt
         ENDDO
      ENDDO

      DO k = 1,llm
         DO j = 1,ip1jm
            pvcov(j,k)= pvcov(j,k) + pdvfi(j,k) * pdt
         ENDDO
      ENDDO

c
      DO j = 1,ip1jmp1
         pps(j) = pps(j) + pdpfi(j) * pdt
      ENDDO
 
      DO iq = 1, 2
         DO k = 1,llm
            DO j = 1,ip1jmp1
               pq(j,k,iq)= pq(j,k,iq) + pdqfi(j,k,iq) * pdt
               pq(j,k,iq)= AMAX1( pq(j,k,iq), qtestw )
            ENDDO
         ENDDO
      ENDDO

      DO iq = 3, nq
         DO k = 1,llm
            DO j = 1,ip1jmp1
               pq(j,k,iq)= pq(j,k,iq) + pdqfi(j,k,iq) * pdt
               pq(j,k,iq)= AMAX1( pq(j,k,iq), qtestt )
            ENDDO
         ENDDO
      ENDDO


      DO  ij   = 1, iim
        xpn(ij) = aire(   ij   ) * pps(  ij     )
        xps(ij) = aire(ij+ip1jm) * pps(ij+ip1jm )
      ENDDO
      tpn      = SSUM(iim,xpn,1)/apoln
      tps      = SSUM(iim,xps,1)/apols

      DO ij   = 1, iip1
        pps (   ij     )  = tpn
        pps ( ij+ip1jm )  = tps
      ENDDO


      DO iq = 1, nq
        DO  k    = 1, llm
          DO  ij   = 1, iim
            xpn(ij) = aire(   ij   ) * pq(  ij    ,k,iq)
            xps(ij) = aire(ij+ip1jm) * pq(ij+ip1jm,k,iq)
          ENDDO
          tpn      = SSUM(iim,xpn,1)/apoln
          tps      = SSUM(iim,xps,1)/apols

          DO ij   = 1, iip1
            pq (   ij   ,k,iq)  = tpn
            pq (ij+ip1jm,k,iq)  = tps
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END
