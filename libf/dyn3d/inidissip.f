!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/inidissip.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE inidissip ( lstardis,nitergdiv,nitergrot,niterh  ,
     *                       tetagdiv,tetagrot,tetatemp             )
c=======================================================================
c   initialisation de la dissipation horizontale
c=======================================================================
c-----------------------------------------------------------------------
c   declarations:
c   -------------

      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use conf_gcm_m
            use comdissipn
      IMPLICIT NONE

      LOGICAL lstardis
      INTEGER nitergdiv,nitergrot,niterh
      REAL    tetagdiv,tetagrot,tetatemp
      REAL fact,zvert(llm),zz
      REAL zh(ip1jmp1),zu(ip1jmp1),zv(ip1jm),deltap(ip1jmp1,llm)
      REAL ullm,vllm,umin,vmin,zhmin,zhmax
      REAL zllm,z1llm

      INTEGER l,ij,idum,ii
      REAL tetamin

      REAL ran1


c-----------------------------------------------------------------------

      print *, "Call sequence information: inidissip"
c
c   calcul des valeurs propres des operateurs par methode iterrative:
c   -----------------------------------------------------------------

      crot     = -1.
      cdivu    = -1.
      cdivh    = -1.

c   calcul de la valeur propre de divgrad:
c   --------------------------------------
      idum = 0
      DO l = 1, llm
       DO ij = 1, ip1jmp1
        deltap(ij,l) = 1.
       ENDDO
      ENDDO

      idum  = -1
      zh(1) = RAN1(idum)-.5
      idum  = 0
      DO ij = 2, ip1jmp1
        zh(ij) = RAN1(idum) -.5
      ENDDO

      CALL filtreg (zh,jjp1,1,2,1,.TRUE.,1)

      CALL minmax(iip1*jjp1,zh,zhmin,zhmax )

      IF ( zhmin .GE. zhmax  )     THEN
         PRINT*,'  Inidissip  zh min max  ',zhmin,zhmax
         STOP'probleme generateur alleatoire dans inidissip'
      ENDIF

      zllm = ABS( zhmax )
      DO l = 1,50
         IF(lstardis) THEN
            CALL divgrad2(1,zh,deltap,niterh,zh)
         ELSE
            CALL divgrad (1,zh,niterh,zh)
         ENDIF

        CALL minmax(iip1*jjp1,zh,zhmin,zhmax )

         zllm  = ABS( zhmax )
         z1llm = 1./zllm
         DO ij = 1,ip1jmp1
            zh(ij) = zh(ij)* z1llm
         ENDDO
      ENDDO

      IF(lstardis) THEN
         cdivh = 1./ zllm
      ELSE
         cdivh = zllm ** ( -1./niterh )
      ENDIF

c   calcul des valeurs propres de gradiv (ii =1) et  nxgrarot(ii=2)
c   -----------------------------------------------------------------
      print*,'calcul des valeurs propres'

      DO  20  ii = 1, 2
c
         DO ij = 1, ip1jmp1
           zu(ij)  = RAN1(idum) -.5
         ENDDO
         CALL filtreg (zu,jjp1,1,2,1,.TRUE.,1)
         DO ij = 1, ip1jm
            zv(ij) = RAN1(idum) -.5
         ENDDO
         CALL filtreg (zv,jjm,1,2,1,.FALSE.,1)

         CALL minmax(iip1*jjp1,zu,umin,ullm )
         CALL minmax(iip1*jjm, zv,vmin,vllm )

         ullm = ABS ( ullm )
         vllm = ABS ( vllm )

         DO  5  l = 1, 50
            IF(ii.EQ.1) THEN
ccccc             CALL covcont( 1,zu,zv,zu,zv )
               IF(lstardis) THEN
                  CALL gradiv2( 1,zu,zv,nitergdiv,zu,zv )
               ELSE
                  CALL gradiv ( 1,zu,zv,nitergdiv,zu,zv )
               ENDIF
            ELSE
               IF(lstardis) THEN
                  CALL nxgraro2( 1,zu,zv,nitergrot,zu,zv )
               ELSE
                  CALL nxgrarot( 1,zu,zv,nitergrot,zu,zv )
               ENDIF
            ENDIF

            CALL minmax(iip1*jjp1,zu,umin,ullm )
            CALL minmax(iip1*jjm, zv,vmin,vllm )

            ullm = ABS  ( ullm )
            vllm = ABS  ( vllm )

            zllm  = MAX( ullm,vllm )
            z1llm = 1./ zllm
            DO ij = 1, ip1jmp1
              zu(ij) = zu(ij)* z1llm
            ENDDO
            DO ij = 1, ip1jm
               zv(ij) = zv(ij)* z1llm
            ENDDO
 5       CONTINUE

         IF ( ii.EQ.1 ) THEN
            IF(lstardis) THEN
               cdivu  = 1./zllm
            ELSE
               cdivu  = zllm **( -1./nitergdiv )
            ENDIF
         ELSE
            IF(lstardis) THEN
               crot   = 1./ zllm
            ELSE
               crot   = zllm **( -1./nitergrot )
            ENDIF
         ENDIF

 20   CONTINUE

c   petit test pour les operateurs non star:
c   ----------------------------------------

c     IF(.NOT.lstardis) THEN
         fact    = rad*24./float(jjm)
         fact    = fact*fact
         PRINT*,'coef u ', fact/cdivu, 1./cdivu
         PRINT*,'coef r ', fact/crot , 1./crot
         PRINT*,'coef h ', fact/cdivh, 1./cdivh
c     ENDIF

c-----------------------------------------------------------------------
c   variation verticale du coefficient de dissipation:
c   --------------------------------------------------

      DO l=1,llm
         zvert(l)=1.
      ENDDO

      fact=2.
c
      DO l = 1, llm
         zz      = 1. - preff/presnivs(l)
         zvert(l)= fact -( fact-1.)/( 1.+zz*zz )
      ENDDO


      PRINT*,'Constantes de temps de la diffusion horizontale'

      tetamin =  1.e+6

      DO l=1,llm
        tetaudiv(l)   = zvert(l)/tetagdiv
        tetaurot(l)   = zvert(l)/tetagrot
        tetah(l)      = zvert(l)/tetatemp

        IF( tetamin.GT. (1./tetaudiv(l)) ) tetamin = 1./ tetaudiv(l)
        IF( tetamin.GT. (1./tetaurot(l)) ) tetamin = 1./ tetaurot(l)
        IF( tetamin.GT. (1./   tetah(l)) ) tetamin = 1./    tetah(l)
      ENDDO

      PRINT *,' INIDI tetamin dtvr ',tetamin,dtvr,iperiod
      idissip = INT( tetamin/( 2.*dtvr*iperiod) ) * iperiod
      PRINT *,' tetamin = ',tetamin
      idissip = MAX(iperiod,idissip)
      dtdiss  = idissip * dtvr
      PRINT *,' INIDI idissip dtdiss ',idissip,dtdiss

      DO l = 1,llm
         PRINT*,zvert(l),dtdiss*tetaudiv(l),dtdiss*tetaurot(l),
     *                   dtdiss*tetah(l)
      ENDDO

c
      RETURN
      END
