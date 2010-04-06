      SUBROUTINE vly(q,pente_max,masse,masse_adv_v)
!
!     Auteurs:   P.Le Van, F.Hourdin, F.Forget
!
!    ***********************************************************
!     Shema  d'advection " pseudo amont " .
!    *************************************************************
!     q,masse_adv_v,w sont des arguments d'entree  pour le s-pg ....
!     dq              sont des arguments de sortie pour le s-pg ....
!
!
!   ----------------------------------------------------------------
      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use logic
      use comgeom
      IMPLICIT NONE
!
!
!
!   Arguments:
!   ----------
      REAL masse(ip1jmp1,llm),pente_max
      REAL masse_adv_v( ip1jm,llm)
      REAL q(ip1jmp1,llm), dq( ip1jmp1,llm)
!
!      Local
!   ---------
!
      INTEGER i,ij,l
!
      REAL airej2,airejjm,airescb(iim),airesch(iim)
      REAL dyq(ip1jmp1,llm),dyqv(ip1jm),zdvm(ip1jmp1,llm)
      REAL adyqv(ip1jm),dyqmax(ip1jmp1)
      REAL qbyv(ip1jm,llm)

      REAL qpns,qpsn,apn,aps,dyn1,dys1,dyn2,dys2,newmasse,fn,fs
!     REAL newq,oldmasse
      Logical extremum,first,testcpu
      REAL temps0,temps1,temps2,temps3,temps4,temps5,second
      SAVE temps0,temps1,temps2,temps3,temps4,temps5
      SAVE first,testcpu

      REAL convpn,convps,convmpn,convmps
      real massepn,masseps,qpn,qps
      REAL sinlon(iip1),sinlondlon(iip1)
      REAL coslon(iip1),coslondlon(iip1)
      SAVE sinlon,coslon,sinlondlon,coslondlon
      SAVE airej2,airejjm
!
!
      REAL      SSUM

      DATA first,testcpu/.true.,.false./
      DATA temps0,temps1,temps2,temps3,temps4,temps5/0.,0.,0.,0.,0.,0./

      IF(first) THEN
         PRINT*,'Shema  Amont nouveau  appele dans  Vanleer   '
         first=.false.
         do i=2,iip1
            coslon(i)=cos(rlonv(i))
            sinlon(i)=sin(rlonv(i))
            coslondlon(i)=coslon(i)*(rlonu(i)-rlonu(i-1))/pi
            sinlondlon(i)=sinlon(i)*(rlonu(i)-rlonu(i-1))/pi
         ENDDO
         coslon(1)=coslon(iip1)
         coslondlon(1)=coslondlon(iip1)
         sinlon(1)=sinlon(iip1)
         sinlondlon(1)=sinlondlon(iip1)
         airej2 = SSUM( iim, aire(iip2), 1 )
         airejjm= SSUM( iim, aire(ip1jm -iim), 1 )
      ENDIF

!
      DO l = 1, llm
!
!   --------------------------------
!      CALCUL EN LATITUDE
!   --------------------------------

!   On commence par calculer la valeur du traceur moyenne sur le
!     premier cercle
!   de latitude autour du pole (qpns pour le pole nord et qpsn pour
!    le pole nord) qui sera utilisee pour evaluer les pentes au pole.

      DO i = 1, iim
      airescb(i) = aire(i+ iip1) * q(i+ iip1,l)
      airesch(i) = aire(i+ ip1jm- iip1) * q(i+ ip1jm- iip1,l)
      ENDDO
      qpns   = SSUM( iim,  airescb ,1 ) / airej2
      qpsn   = SSUM( iim,  airesch ,1 ) / airejjm

!   calcul des pentes aux points v

      DO ij=1,ip1jm
         dyqv(ij)=q(ij,l)-q(ij+iip1,l)
         adyqv(ij)=abs(dyqv(ij))
      ENDDO

!   calcul des pentes aux points scalaires

      DO ij=iip2,ip1jm
         dyq(ij,l)=.5*(dyqv(ij-iip1)+dyqv(ij))
         dyqmax(ij)=min(adyqv(ij-iip1),adyqv(ij))
         dyqmax(ij)=pente_max*dyqmax(ij)
      ENDDO

!   calcul des pentes aux poles

      DO ij=1,iip1
         dyq(ij,l)=qpns-q(ij+iip1,l)
         dyq(ip1jm+ij,l)=q(ip1jm+ij-iip1,l)-qpsn
      ENDDO

!   filtrage de la derivee
      dyn1=0.
      dys1=0.
      dyn2=0.
      dys2=0.
      DO ij=1,iim
         dyn1=dyn1+sinlondlon(ij)*dyq(ij,l)
         dys1=dys1+sinlondlon(ij)*dyq(ip1jm+ij,l)
         dyn2=dyn2+coslondlon(ij)*dyq(ij,l)
         dys2=dys2+coslondlon(ij)*dyq(ip1jm+ij,l)
      ENDDO
      DO ij=1,iip1
         dyq(ij,l)=dyn1*sinlon(ij)+dyn2*coslon(ij)
         dyq(ip1jm+ij,l)=dys1*sinlon(ij)+dys2*coslon(ij)
      ENDDO

!   calcul des pentes limites aux poles

      goto 8888
      fn=1.
      fs=1.
      DO ij=1,iim
         IF(pente_max*adyqv(ij).lt.abs(dyq(ij,l))) THEN
            fn=min(pente_max*adyqv(ij)/abs(dyq(ij,l)),fn)
         ENDIF
      IF(pente_max*adyqv(ij+ip1jm-iip1).lt.abs(dyq(ij+ip1jm,l))) THEN
         fs=min(pente_max*adyqv(ij+ip1jm-iip1)/abs(dyq(ij+ip1jm,l)),fs)
         ENDIF
      ENDDO
      DO ij=1,iip1
         dyq(ij,l)=fn*dyq(ij,l)
         dyq(ip1jm+ij,l)=fs*dyq(ip1jm+ij,l)
      ENDDO
 8888   continue
      DO ij=1,iip1
         dyq(ij,l)=0.
         dyq(ip1jm+ij,l)=0.
      ENDDO

!   calcul des pentes limitees

      DO ij=iip2,ip1jm
         IF(dyqv(ij)*dyqv(ij-iip1).gt.0.) THEN
            dyq(ij,l)=sign(min(abs(dyq(ij,l)),dyqmax(ij)),dyq(ij,l))
         ELSE
            dyq(ij,l)=0.
         ENDIF
      ENDDO

      ENDDO

      DO l=1,llm
       DO ij=1,ip1jm
          IF(masse_adv_v(ij,l).gt.0) THEN
              qbyv(ij,l)=q(ij+iip1,l)+dyq(ij+iip1,l)*                   &
     &                   0.5*(1.-masse_adv_v(ij,l)/masse(ij+iip1,l))
          ELSE
              qbyv(ij,l)=q(ij,l)-dyq(ij,l)*                             &
     &                   0.5*(1.+masse_adv_v(ij,l)/masse(ij,l))
          ENDIF
          qbyv(ij,l)=masse_adv_v(ij,l)*qbyv(ij,l)
       ENDDO
      ENDDO


      DO l=1,llm
         DO ij=iip2,ip1jm
            newmasse=masse(ij,l)                                        &
     &      +masse_adv_v(ij,l)-masse_adv_v(ij-iip1,l)
            q(ij,l)=(q(ij,l)*masse(ij,l)+qbyv(ij,l)-qbyv(ij-iip1,l))    &
     &         /newmasse
            masse(ij,l)=newmasse
         ENDDO
!.-. ancienne version
!        convpn=SSUM(iim,qbyv(1,l),1)/apoln
!        convmpn=ssum(iim,masse_adv_v(1,l),1)/apoln

         convpn=SSUM(iim,qbyv(1,l),1)
         convmpn=ssum(iim,masse_adv_v(1,l),1)
         massepn=ssum(iim,masse(1,l),1)
         qpn=0.
         do ij=1,iim
            qpn=qpn+masse(ij,l)*q(ij,l)
         enddo
         qpn=(qpn+convpn)/(massepn+convmpn)
         do ij=1,iip1
            q(ij,l)=qpn
         enddo

!        convps=-SSUM(iim,qbyv(ip1jm-iim,l),1)/apols
!        convmps=-ssum(iim,masse_adv_v(ip1jm-iim,l),1)/apols

         convps=-SSUM(iim,qbyv(ip1jm-iim,l),1)
         convmps=-ssum(iim,masse_adv_v(ip1jm-iim,l),1)
         masseps=ssum(iim, masse(ip1jm+1,l),1)
         qps=0.
         do ij = ip1jm+1,ip1jmp1-1
            qps=qps+masse(ij,l)*q(ij,l)
         enddo
         qps=(qps+convps)/(masseps+convmps)
         do ij=ip1jm+1,ip1jmp1
            q(ij,l)=qps
         enddo

!.-. fin ancienne version

!._. nouvelle version
!        convpn=SSUM(iim,qbyv(1,l),1)
!        convmpn=ssum(iim,masse_adv_v(1,l),1)
!        oldmasse=ssum(iim,masse(1,l),1)
!        newmasse=oldmasse+convmpn
!        newq=(q(1,l)*oldmasse+convpn)/newmasse
!        newmasse=newmasse/apoln
!        DO ij = 1,iip1
!           q(ij,l)=newq
!           masse(ij,l)=newmasse*aire(ij)
!        ENDDO
!        convps=-SSUM(iim,qbyv(ip1jm-iim,l),1)
!        convmps=-ssum(iim,masse_adv_v(ip1jm-iim,l),1)
!        oldmasse=ssum(iim,masse(ip1jm-iim,l),1)
!        newmasse=oldmasse+convmps
!        newq=(q(ip1jmp1,l)*oldmasse+convps)/newmasse
!        newmasse=newmasse/apols
!        DO ij = ip1jm+1,ip1jmp1
!           q(ij,l)=newq
!           masse(ij,l)=newmasse*aire(ij)
!        ENDDO
!._. fin nouvelle version
      ENDDO

      RETURN
      END
