!
! $Id: writehist.F 1403 2010-07-01 09:02:53Z fairhead $
!
      subroutine writehist(time,vcov,ucov,teta,phi,q,masse,ps,phis)

      use dimens_m, only: nqmx, llm, jjm
      USE iniadvtrac_m, ONLY: ttext
      use com_io_dyn, only: histid,histvid,histuid
      use paramet_m, only: ip1jm, ip1jmp1, iip1, jjp1
      use temps, only: itau_dyn
      use histwrite_m, only: histwrite
      use histcom, only: histsync

      implicit none

C
C   Ecriture du fichier histoire au format IOIPSL
C
C   Appels succesifs des routines: histwrite
C
C   Entree:
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
C   L. Fairhead, LMD, 03/99
C
C =====================================================================
C
C   Declarations

C
C   Arguments
C

      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm) 
      REAL teta(ip1jmp1,llm),phi(ip1jmp1,llm)                   
      REAL ps(ip1jmp1),masse(ip1jmp1,llm)                   
      REAL phis(ip1jmp1)                  
      REAL q(ip1jmp1,llm,nqmx)
      integer time


! This routine needs IOIPSL to work
C   Variables locales
C
      integer iq, ii, ll
      integer ndexu(ip1jmp1*llm),ndexv(ip1jm*llm),ndex2d(ip1jmp1)
      logical ok_sync
      integer itau_w
      REAL vnat(ip1jm,llm),unat(ip1jmp1,llm)

C
C  Initialisations
C
      ndexu = 0
      ndexv = 0
      ndex2d = 0
      ok_sync =.TRUE.
      itau_w = itau_dyn + time
!  Passage aux composantes naturelles du vent
      call covnat(llm, ucov, vcov, unat, vnat)
C
C  Appels a histwrite pour l'ecriture des variables a sauvegarder
C
C  Vents U
C
      call histwrite(histuid, 'u', itau_w, unat)
C
C  Vents V
C
      call histwrite(histvid, 'v', itau_w, vnat)

C
C  Temperature potentielle
C
      call histwrite(histid, 'teta', itau_w, teta)
C
C  Geopotentiel
C
      call histwrite(histid, 'phi', itau_w, phi)
C
C  Traceurs
C
!        DO iq=1,nqmx
!          call histwrite(histid, ttext(iq), itau_w, q(:,:,iq), 
!     .                   iip1*jjp1*llm, ndexu)
!        enddo
!C
C  Masse
C
      call histwrite(histid,'masse',itau_w, masse)
C
C  Pression au sol
C
      call histwrite(histid, 'ps', itau_w, ps)
C
C  Geopotentiel au sol
C
!      call histwrite(histid, 'phis', itau_w, phis, iip1*jjp1, ndex2d)
C
C  Fin
C
      if (ok_sync) then
        call histsync(histid)
        call histsync(histvid)
        call histsync(histuid)
      endif
      return
      end
