!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/bernoui.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE bernoui (ngrid,nlay,pphi,pecin,pbern)
      use dimens_m
      use paramet_m
      use conf_gcm_m
      use filtreg_m, only: filtreg
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:   P. Le Van
c   -------
c
c   Objet:
c   ------
c     calcul de la fonction de Bernouilli aux niveaux s  .....
c     phi  et  ecin  sont des arguments d'entree pour le s-pg .......
c          bern       est un  argument de sortie pour le s-pg  ......
c
c    fonction de Bernouilli = bern = filtre de( geopotentiel + 
c                              energ.cinet.)
c
c=======================================================================
c
c-----------------------------------------------------------------------
c   Decalrations:
c   -------------
c
c
c   Arguments:
c   ----------
c
      INTEGER nlay,ngrid
      REAL pphi(ngrid*nlay),pecin(ngrid*nlay),pbern(ngrid*nlay)
c
c   Local:
c   ------
c
      INTEGER   ijl
c
c-----------------------------------------------------------------------
c   calcul de Bernouilli:
c   ---------------------
c
      DO 4 ijl = 1,ngrid*nlay
         pbern( ijl ) =  pphi( ijl ) + pecin( ijl )
   4  CONTINUE
c
c-----------------------------------------------------------------------
c   filtre:
c   -------
c
      CALL filtreg( pbern, jjp1, llm, 2,1, .true.)
c
c-----------------------------------------------------------------------
      RETURN
      END
