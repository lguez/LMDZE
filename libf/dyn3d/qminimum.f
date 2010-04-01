!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/qminimum.F,v 1.1.1.1 2004/05/19 12:53:05 lmdzadmin Exp $
!
      SUBROUTINE qminimum( q,nq,deltap )

      use dimens_m
      use paramet_m
      use comvert
      IMPLICIT none
c
c  -- Objet : Traiter les valeurs trop petites (meme negatives)
c             pour l'eau vapeur et l'eau liquide
c
c
      INTEGER nq
      REAL q(ip1jmp1,llm,nq), deltap(ip1jmp1,llm)
c
      INTEGER iq_vap, iq_liq
      PARAMETER ( iq_vap = 1 ) ! indice pour l'eau vapeur
      PARAMETER ( iq_liq = 2 ) ! indice pour l'eau liquide
      REAL seuil_vap, seuil_liq
      PARAMETER ( seuil_vap = 1.0e-10 ) ! seuil pour l'eau vapeur
      PARAMETER ( seuil_liq = 1.0e-11 ) ! seuil pour l'eau liquide
c
c  NB. ....( Il est souhaitable mais non obligatoire que les valeurs des
c            parametres seuil_vap, seuil_liq soient pareilles a celles 
c            qui  sont utilisees dans la routine    ADDFI       )
c     .................................................................
c
      INTEGER i, k, iq
      REAL zx_defau, zx_abc, zx_pump(ip1jmp1), pompe
c
      REAL SSUM
c
      INTEGER imprim
      SAVE imprim
      DATA imprim /0/
c
c Quand l'eau liquide est trop petite (ou negative), on prend
c l'eau vapeur de la meme couche et la convertit en eau liquide
c (sans changer la temperature !)
c
      DO 1000 k = 1, llm
      DO 1040 i = 1, ip1jmp1
            zx_defau      = AMAX1( seuil_liq - q(i,k,iq_liq), 0.0 )
            q(i,k,iq_vap) = q(i,k,iq_vap) - zx_defau
            q(i,k,iq_liq) = q(i,k,iq_liq) + zx_defau
 1040 CONTINUE
 1000 CONTINUE
c
c Quand l'eau vapeur est trop faible (ou negative), on complete
c le defaut en prennant de l'eau vapeur de la couche au-dessous.
c
      iq = iq_vap
c
      DO k = llm, 2, -1
ccc      zx_abc = dpres(k) / dpres(k-1)
      DO i = 1, ip1jmp1
         zx_abc = deltap(i,k)/deltap(i,k-1)
         zx_defau    = AMAX1( seuil_vap - q(i,k,iq), 0.0 )
         q(i,k-1,iq) =  q(i,k-1,iq) - zx_defau * zx_abc
         q(i,k,iq)   =  q(i,k,iq)   + zx_defau  
      ENDDO
      ENDDO
c
c Quand il s'agit de la premiere couche au-dessus du sol, on
c doit imprimer un message d'avertissement (saturation possible).
c
      DO i = 1, ip1jmp1
         zx_pump(i) = AMAX1( 0.0, seuil_vap - q(i,1,iq) )
         q(i,1,iq)  = AMAX1( q(i,1,iq), seuil_vap )
      ENDDO
      pompe = SSUM(ip1jmp1,zx_pump,1)
      IF (imprim.LE.500 .AND. pompe.GT.0.0) THEN
         WRITE(6,'(1x,"ATT!:on pompe de l eau au sol",e15.7)') pompe
         DO i = 1, ip1jmp1
            IF (zx_pump(i).GT.0.0) THEN
               imprim = imprim + 1
            ENDIF
         ENDDO
      ENDIF
c
      RETURN
      END
