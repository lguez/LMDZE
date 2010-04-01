      SUBROUTINE vlx(q,pente_max,masse,u_m)

!     Auteurs:   P.Le Van, F.Hourdin, F.Forget
!
!   *******************************************************
!     Shema  d'advection " pseudo amont " .
!    ************************************************************
!     nq,iq,q,pbaru,pbarv,w sont des arguments d'entree  pour le
!     s-pg ....
!
!
!   --------------------------------------------------------------
      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use logic
      IMPLICIT NONE
!
!
!
!   Arguments:
!   ----------
      REAL masse(ip1jmp1,llm),pente_max
      REAL u_m( ip1jmp1,llm ),pbarv( iip1,jjm,llm)
      REAL q(ip1jmp1,llm)
      REAL w(ip1jmp1,llm)
!
!      Local
!   ---------
!
      INTEGER ij,l,j,i,iju,ijq,indu(ip1jmp1),niju
      INTEGER n0,iadvplus(ip1jmp1,llm),nl(llm)
!
      REAL new_m,zu_m,zdum(ip1jmp1,llm)
      REAL sigu(ip1jmp1),dxq(ip1jmp1,llm),dxqu(ip1jmp1)
      REAL zz(ip1jmp1)
      REAL adxqu(ip1jmp1),dxqmax(ip1jmp1,llm)
      REAL u_mq(ip1jmp1,llm)

      Logical extremum,first,testcpu
      SAVE first,testcpu

      REAL      SSUM
      REAL temps0,temps1,temps2,temps3,temps4,temps5,second
      SAVE temps0,temps1,temps2,temps3,temps4,temps5

      REAL z1,z2,z3

      DATA first,testcpu/.true.,.false./

      IF(first) THEN
         temps1=0.
         temps2=0.
         temps3=0.
         temps4=0.
         temps5=0.
         first=.false.
      ENDIF

!   calcul de la pente a droite et a gauche de la maille


      IF (pente_max.gt.-1.e-5) THEN
!       IF (pente_max.gt.10) THEN

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
               dxqmax(ij,l)=pente_max*                                  &
     &      min(adxqu(ij-1),adxqu(ij))
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
               dxq(ij,l)=                                               &
     &         sign(min(abs(dxq(ij,l)),dxqmax(ij,l)),dxq(ij,l))
            ENDDO

               ! l=1,llm
         ENDDO
!print*,'Ok calcul des pentes'

           ! (pente_max.lt.-1.e-5)
      ELSE

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

            ! (pente_max.lt.-1.e-5)
      ENDIF

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

! print*,'Bouclage en iip1'

!   calcul des flux a gauche et a droite

!   on cumule le flux correspondant a toutes les mailles dont la masse
!   au travers de la paroi pENDant le pas de temps.
!print*,'Cumule ....'

      DO l=1,llm
       DO ij=iip2,ip1jm-1
!      print*,'masse(',ij,')=',masse(ij,l)
          IF (u_m(ij,l).gt.0.) THEN
             zdum(ij,l)=1.-u_m(ij,l)/masse(ij,l)
             u_mq(ij,l)=u_m(ij,l)*(q(ij,l)+0.5*zdum(ij,l)*dxq(ij,l))
          ELSE
             zdum(ij,l)=1.+u_m(ij,l)/masse(ij+1,l)
             u_mq(ij,l)=u_m(ij,l)                                       &
     &            *(q(ij+1,l)-0.5*zdum(ij,l)*dxq(ij+1,l))
          ENDIF
       ENDDO
      ENDDO
!      stop

!      go to 9999
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
!print*,'Ok test 1'
      DO l=1,llm
       DO ij=iip1+iip1,ip1jm,iip1
          iadvplus(ij,l)=iadvplus(ij-iim,l)
       ENDDO
      ENDDO
! print*,'Ok test 2'


!   traitement special pour le cas ou on advecte en longitude plus
!     que le
!   contenu de la maille.
!   cette partie est mal vectorisee.

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
!C      PRINT*,'Nombre de points pour lesquels on advect plus que le'
!C     &       ,'contenu de la maille : ',n0

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
!              PRINT*,'niju,nl',niju,nl(l)

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
                     u_mq(ij,l)=u_mq(ij,l)+zu_m*                        &
     &               (q(ijq,l)+0.5*(1.-zu_m/masse(ijq,l))*dxq(ijq,l))
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
                     u_mq(ij,l)=u_mq(ij,l)+zu_m*(q(ijq,l)-              &
     &               0.5*(1.+zu_m/masse(ijq,l))*dxq(ijq,l))
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
             ! n0.gt.0
      ENDIF
 9999   continue


!   bouclage en latitude
!print*,'cvant bouclage en latitude'
      DO l=1,llm
        DO ij=iip1+iip1,ip1jm,iip1
           u_mq(ij,l)=u_mq(ij-iim,l)
        ENDDO
      ENDDO


!   calcul des tENDances

      DO l=1,llm
         DO ij=iip2+1,ip1jm
            new_m=masse(ij,l)+u_m(ij-1,l)-u_m(ij,l)
            q(ij,l)=(q(ij,l)*masse(ij,l)+                               &
     &      u_mq(ij-1,l)-u_mq(ij,l))                                    &
     &      /new_m
            masse(ij,l)=new_m
         ENDDO
!   ModIF Fred 22 03 96 correction d'un bug (les scopy ci-dessous)
         DO ij=iip1+iip1,ip1jm,iip1
            q(ij-iim,l)=q(ij,l)
            masse(ij-iim,l)=masse(ij,l)
         ENDDO
      ENDDO


      RETURN
      END
