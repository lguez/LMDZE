!
! $Header: /home/cvsroot/LMDZ4/libf/bibio/writedynav.F,v 1.1.1.1 2004/05/19 12:53:05 lmdzadmin Exp $
!
      subroutine writedynav( histid, nq, time, vcov, 
     ,                          ucov,teta,ppk,phi,q,masse,ps,phis)

C
C   Ecriture du fichier histoire au format IOIPSL
C
C   Appels succesifs des routines: histwrite
C
C   Entree:
C      histid: ID du fichier histoire
C      time: temps de l'ecriture
C      vcov: vents v covariants
C      ucov: vents u covariants
C      teta: temperature potentielle
C      phi : geopotentiel instantane
C      q   : traceurs
C      masse: masse
C      ps   :pression au sol
C      phis : geopotentiel au sol
C      
C
C   Sortie:
C      fileid: ID du fichier netcdf cree
C
C   L. Fairhead, LMD, 03/99
C
C =====================================================================
C
C   Declarations
      USE histwrite_m
      use histcom
      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use logic
      use comgeom
      use serre
      use temps
      use ener
      use iniadvtrac_m
      implicit none


C
C   Arguments
C

      INTEGER histid, nq
      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm) 
      REAL teta(ip1jmp1*llm),phi(ip1jmp1,llm),ppk(ip1jmp1*llm)                  
      REAL ps(ip1jmp1),masse(ip1jmp1,llm)                   
      REAL phis(ip1jmp1)                  
      REAL q(ip1jmp1,llm,nq)
      integer, intent(in):: time


C   Variables locales
C
      integer ndex2d(iip1*jjp1),ndex3d(iip1*jjp1*llm),iq, ii, ll
      real us(ip1jmp1*llm), vs(ip1jmp1*llm)
      real tm(ip1jmp1*llm)
      REAL vnat(ip1jm,llm),unat(ip1jmp1,llm) 
      logical ok_sync
      integer itau_w
C
C  Initialisations
C
      ndex3d = 0
      ndex2d = 0
      ok_sync = .TRUE.
      us = 999.999
      vs = 999.999
      tm = 999.999
      vnat = 999.999
      unat = 999.999
      itau_w = itau_dyn + time

C Passage aux composantes naturelles du vent
      call covnat(llm, ucov, vcov, unat, vnat)

C
C  Appels a histwrite pour l'ecriture des variables a sauvegarder
C
C  Vents U scalaire
C
      call gr_u_scal(llm, unat, us)
      call histwrite(histid, 'u', itau_w, us)
C
C  Vents V scalaire
C
      call gr_v_scal(llm, vnat, vs)
      call histwrite(histid, 'v', itau_w, vs)
C
C  Temperature potentielle moyennee
C
      call histwrite(histid, 'theta', itau_w, teta)
C
C  Temperature moyennee
C
      do ii = 1, ijp1llm
        tm(ii) = teta(ii) * ppk(ii)/cpp
      enddo
      call histwrite(histid, 'temp', itau_w, tm)
C
C  Geopotentiel
C
      call histwrite(histid, 'phi', itau_w, phi)
C
C  Traceurs
C
        DO iq=1,nq
          call histwrite(histid, ttext(iq), itau_w, q(:,:,iq))
        enddo
C
C  Masse
C
       call histwrite(histid, 'masse', itau_w, masse)
C
C  Pression au sol
C
       call histwrite(histid, 'ps', itau_w, ps)
C
C  Geopotentiel au sol
C
       call histwrite(histid, 'phis', itau_w, phis)
C
C  Fin
C
      if (ok_sync) call histsync(histid)
      return
      end
