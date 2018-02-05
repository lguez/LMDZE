      SUBROUTINE vlxqs(q,pente_max,masse,u_m,qsat)
!
!     Auteurs:   P.Le Van, F.Hourdin, F.Forget
!
!    ********************************************************************
!     Shema  d'advection " pseudo amont " .
!    ********************************************************************
!
!   --------------------------------------------------------------------
      use dimens_m
      use paramet_m
      use comconst
      use disvert_m
      use conf_gcm_m
      IMPLICIT NONE
!
!
!
!   Arguments:
!   ----------
      REAL masse(ip1jmp1,llm),pente_max
      REAL u_m( ip1jmp1,llm )
      REAL q(ip1jmp1,llm)
      REAL qsat(ip1jmp1,llm)
!
!      Local
!   ---------
!
      INTEGER ij,l,j,i,iju,ijq,indu(ip1jmp1),niju
      INTEGER n0,iadvplus(ip1jmp1,llm),nl(llm)
!
      REAL new_m,zu_m,zdum(ip1jmp1,llm)
      REAL dxq(ip1jmp1,llm),dxqu(ip1jmp1)
      REAL zz(ip1jmp1)
      REAL adxqu(ip1jmp1),dxqmax(ip1jmp1,llm)
      REAL u_mq(ip1jmp1,llm)

!   calcul de la pente a droite et a gauche de la maille


      IF (pente_max.gt.-1.e-5) THEN
!     IF (pente_max.gt.10) THEN

!   calcul des pentes avec limitation, Van Leer scheme I:
!   -----------------------------------------------------

!   calcul de la pente aux points u
         DO l = 1, llm
            DO ij=iip2,ip1jm-1
               dxqu(ij)=q(ij+1,l)-q(ij,l)
!              IF(u_m(ij,l).lt.0.) stop'limx n admet pas les U<0'
!              sigu(ij)=u_m(ij,l)/masse(ij,l)
            ENDDO
            DO ij=iip1+iip1,ip1jm,iip1
               dxqu(ij)=dxqu(ij-iim)
!              sigu(ij)=sigu(ij-iim)
            ENDDO

            DO ij=iip2,ip1jm
               adxqu(ij)=abs(dxqu(ij))
            ENDDO

!   calcul de la pente maximum dans la maille en valeur absolue

            DO ij=iip2+1,ip1jm
               dxqmax(ij,l)=pente_max* &
            min(adxqu(ij-1),adxqu(ij))
! limitation subtile
!    ,      min(adxqu(ij-1)/sigu(ij-1),adxqu(ij)/(1.-sigu(ij)))


            ENDDO

            DO ij=iip1+iip1,ip1jm,iip1
               dxqmax(ij-iim,l)=dxqmax(ij,l)
            ENDDO

            DO ij=iip2+1,ip1jm
               IF(dxqu(ij-1)*dxqu(ij).gt.0) THEN
                  dxq(ij,l)=dxqu(ij-1)+dxqu(ij)
               ELSE
!   extremum local
                  dxq(ij,l)=0.
               ENDIF
               dxq(ij,l)=0.5*dxq(ij,l)
               dxq(ij,l)= &
               sign(min(abs(dxq(ij,l)),dxqmax(ij,l)),dxq(ij,l))
            ENDDO

         ENDDO ! l=1,llm

      ELSE ! (pente_max.lt.-1.e-5)

!   Pentes produits:
!   ----------------

         DO l = 1, llm
            DO ij=iip2,ip1jm-1
               dxqu(ij)=q(ij+1,l)-q(ij,l)
            ENDDO
            DO ij=iip1+iip1,ip1jm,iip1
               dxqu(ij)=dxqu(ij-iim)
            ENDDO

            DO ij=iip2+1,ip1jm
               zz(ij)=dxqu(ij-1)*dxqu(ij)
               zz(ij)=zz(ij)+zz(ij)
               IF(zz(ij).gt.0) THEN
                  dxq(ij,l)=zz(ij)/(dxqu(ij-1)+dxqu(ij))
               ELSE
!   extremum local
                  dxq(ij,l)=0.
               ENDIF
            ENDDO

         ENDDO

      ENDIF ! (pente_max.lt.-1.e-5)

!   bouclage de la pente en iip1:
!   -----------------------------

      DO l=1,llm
         DO ij=iip1+iip1,ip1jm,iip1
            dxq(ij-iim,l)=dxq(ij,l)
         ENDDO

         DO ij=1,ip1jmp1
            iadvplus(ij,l)=0
         ENDDO

      ENDDO


!   calcul des flux a gauche et a droite

!   on cumule le flux correspondant a toutes les mailles dont la masse
!   au travers de la paroi pENDant le pas de temps.
!   le rapport de melange de l'air advecte est min(q_vanleer, Qsat_downwind)
      DO l=1,llm
       DO ij=iip2,ip1jm-1
          IF (u_m(ij,l).gt.0.) THEN
             zdum(ij,l)=1.-u_m(ij,l)/masse(ij,l)
             u_mq(ij,l)=u_m(ij,l)* &
               min(q(ij,l)+0.5*zdum(ij,l)*dxq(ij,l),qsat(ij+1,l))
          ELSE
             zdum(ij,l)=1.+u_m(ij,l)/masse(ij+1,l)
             u_mq(ij,l)=u_m(ij,l)* &
               min(q(ij+1,l)-0.5*zdum(ij,l)*dxq(ij+1,l),qsat(ij,l))
          ENDIF
       ENDDO
      ENDDO


!   detection des points ou on advecte plus que la masse de la
!   maille
      DO l=1,llm
         DO ij=iip2,ip1jm-1
            IF(zdum(ij,l).lt.0) THEN
               iadvplus(ij,l)=1
               u_mq(ij,l)=0.
            ENDIF
         ENDDO
      ENDDO
      DO l=1,llm
       DO ij=iip1+iip1,ip1jm,iip1
          iadvplus(ij,l)=iadvplus(ij-iim,l)
       ENDDO
      ENDDO



!   traitement special pour le cas ou on advecte en longitude plus que le
!   contenu de la maille.
!   cette partie est mal vectorisee.

!   pas d'influence de la pression saturante (pour l'instant)

!  calcul du nombre de maille sur lequel on advecte plus que la maille.

      n0=0
      DO l=1,llm
         nl(l)=0
         DO ij=iip2,ip1jm
            nl(l)=nl(l)+iadvplus(ij,l)
         ENDDO
         n0=n0+nl(l)
      ENDDO

      IF(n0.gt.0) THEN
         DO l=1,llm
            IF(nl(l).gt.0) THEN
               iju=0
!   indicage des mailles concernees par le traitement special
               DO ij=iip2,ip1jm
                  IF(iadvplus(ij,l).eq.1.and.mod(ij,iip1).ne.0) THEN
                     iju=iju+1
                     indu(iju)=ij
                  ENDIF
               ENDDO
               niju=iju

!  traitement des mailles
               DO iju=1,niju
                  ij=indu(iju)
                  j=(ij-1)/iip1+1
                  zu_m=u_m(ij,l)
                  u_mq(ij,l)=0.
                  IF(zu_m.gt.0.) THEN
                     ijq=ij
                     i=ijq-(j-1)*iip1
!   accumulation pour les mailles completements advectees
                     do while(zu_m.gt.masse(ijq,l))
                        u_mq(ij,l)=u_mq(ij,l)+q(ijq,l)*masse(ijq,l)
                        zu_m=zu_m-masse(ijq,l)
                        i=mod(i-2+iim,iim)+1
                        ijq=(j-1)*iip1+i
                     ENDDO
!   ajout de la maille non completement advectee
                     u_mq(ij,l)=u_mq(ij,l)+zu_m* &
                     (q(ijq,l)+0.5*(1.-zu_m/masse(ijq,l))*dxq(ijq,l))
                  ELSE
                     ijq=ij+1
                     i=ijq-(j-1)*iip1
!   accumulation pour les mailles completements advectees
                     do while(-zu_m.gt.masse(ijq,l))
                        u_mq(ij,l)=u_mq(ij,l)-q(ijq,l)*masse(ijq,l)
                        zu_m=zu_m+masse(ijq,l)
                        i=mod(i,iim)+1
                        ijq=(j-1)*iip1+i
                     ENDDO
!   ajout de la maille non completement advectee
                     u_mq(ij,l)=u_mq(ij,l)+zu_m*(q(ijq,l)- &
                     0.5*(1.+zu_m/masse(ijq,l))*dxq(ijq,l))
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF  ! n0.gt.0



!   bouclage en latitude

      DO l=1,llm
        DO ij=iip1+iip1,ip1jm,iip1
           u_mq(ij,l)=u_mq(ij-iim,l)
        ENDDO
      ENDDO


!   calcul des tendances

      DO l=1,llm
         DO ij=iip2+1,ip1jm
            new_m=masse(ij,l)+u_m(ij-1,l)-u_m(ij,l)
            q(ij,l)=(q(ij,l)*masse(ij,l)+ &
            u_mq(ij-1,l)-u_mq(ij,l)) &
            /new_m
            masse(ij,l)=new_m
         ENDDO
!   Modif Fred 22 03 96 correction d'un bug (les scopy ci-dessous)
         DO ij=iip1+iip1,ip1jm,iip1
            q(ij-iim,l)=q(ij,l)
            masse(ij-iim,l)=masse(ij,l)
         ENDDO
      ENDDO

!     CALL SCOPY((jjm-1)*llm,q(iip1+iip1,1),iip1,q(iip2,1),iip1)
!     CALL SCOPY((jjm-1)*llm,masse(iip1+iip1,1),iip1,masse(iip2,1),iip1)


      RETURN
      END
