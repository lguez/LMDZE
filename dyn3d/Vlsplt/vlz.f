      SUBROUTINE vlz(q,pente_max,masse,w)
!
!     Auteurs:   P.Le Van, F.Hourdin, F.Forget
!
!    *************************************************************
!     Shema  d'advection " pseudo amont " .
!    ****************************************************************
!    q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....
!     dq              sont des arguments de sortie pour le s-pg ....
!
!
!   ----------------------------------------------------------------
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
      REAL q(ip1jmp1,llm)
      REAL w(ip1jmp1,llm+1)
!
!      Local
!   ---------
!
      INTEGER i,ij,l,j,ii
!
      REAL wq(ip1jmp1,llm+1),newmasse

      REAL dzq(ip1jmp1,llm),dzqw(ip1jmp1,llm),adzqw(ip1jmp1,llm),dzqmax
      REAL sigw

      LOGICAL testcpu
      SAVE testcpu

      REAL temps0,temps1,temps2,temps3,temps4,temps5,second
      SAVE temps0,temps1,temps2,temps3,temps4,temps5
      REAL      SSUM

      DATA testcpu/.false./
      DATA temps0,temps1,temps2,temps3,temps4,temps5/0.,0.,0.,0.,0.,0./

!    On oriente tout dans le sens de la pression c'est a dire dans le
!    sens de W

      DO l=2,llm
         DO ij=1,ip1jmp1
            dzqw(ij,l)=q(ij,l-1)-q(ij,l)
            adzqw(ij,l)=abs(dzqw(ij,l))
         ENDDO
      ENDDO

      DO l=2,llm-1
         DO ij=1,ip1jmp1
            IF(dzqw(ij,l)*dzqw(ij,l+1).gt.0.) THEN
                dzq(ij,l)=0.5*(dzqw(ij,l)+dzqw(ij,l+1))
            ELSE
                dzq(ij,l)=0.
            ENDIF
            dzqmax=pente_max*min(adzqw(ij,l),adzqw(ij,l+1))
            dzq(ij,l)=sign(min(abs(dzq(ij,l)),dzqmax),dzq(ij,l))
         ENDDO
      ENDDO

      DO ij=1,ip1jmp1
         dzq(ij,1)=0.
         dzq(ij,llm)=0.
      ENDDO

! ---------------------------------------------------------------
!   .... calcul des termes d'advection verticale  .......
! ---------------------------------------------------------------

! calcul de  - d( q   * w )/ d(sigma)    qu'on ajoute a  dq pour
!     calculer dq

       DO l = 1,llm-1
         do  ij = 1,ip1jmp1
          IF(w(ij,l+1).gt.0.) THEN
             sigw=w(ij,l+1)/masse(ij,l+1)
             wq(ij,l+1)=w(ij,l+1)*(q(ij,l+1)+0.5*(1.-sigw)*dzq(ij,l+1))
          ELSE
             sigw=w(ij,l+1)/masse(ij,l)
             wq(ij,l+1)=w(ij,l+1)*(q(ij,l)-0.5*(1.+sigw)*dzq(ij,l))
          ENDIF
         ENDDO
       ENDDO

       DO ij=1,ip1jmp1
          wq(ij,llm+1)=0.
          wq(ij,1)=0.
       ENDDO

      DO l=1,llm
         DO ij=1,ip1jmp1
            newmasse=masse(ij,l)+w(ij,l+1)-w(ij,l)
            q(ij,l)=(q(ij,l)*masse(ij,l)+wq(ij,l+1)-wq(ij,l))           &
     &         /newmasse
            masse(ij,l)=newmasse
         ENDDO
      ENDDO

      END
