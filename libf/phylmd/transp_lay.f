      SUBROUTINE transp_lay (paprs,tsol,
     e                   t, q, u, v, geom,
     s                   vtran_e, vtran_q, utran_e, utran_q)
c
      use dimens_m
      use dimphy
      use SUPHEC_M
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X.Li (LMD/CNRS)
c Date: le 25 avril 1994
c Objet: Calculer le transport de l'energie et de la vapeur d'eau
c======================================================================
c
c
      REAL, intent(in):: paprs(klon,klev+1)
      real tsol(klon)
      REAL t(klon,klev), q(klon,klev), u(klon,klev), v(klon,klev)
      REAL utran_e(klon,klev), utran_q(klon,klev)
      REAL vtran_e(klon,klev), vtran_q(klon,klev)
c
      INTEGER i, l
c     ------------------------------------------------------------------
      REAL geom(klon,klev), esh
c     ------------------------------------------------------------------
      DO l = 1, klev
      DO i = 1, klon
         utran_e(i,l) = 0.0
         utran_q(i,l) = 0.0
         vtran_e(i,l) = 0.0
         vtran_q(i,l) = 0.0
      ENDDO
      ENDDO
c
      DO l = 1, klev
      DO i = 1, klon
         esh = RCPD*t(i,l) + RLVTT*q(i,l) + geom(i,l)
         utran_e(i,l)=utran_e(i,l)+ u(i,l)*esh*
     .                (paprs(i,l)-paprs(i,l+1))/RG
         utran_q(i,l)=utran_q(i,l)+ u(i,l)*q(i,l)
     .                *(paprs(i,l)-paprs(i,l+1))/RG
         vtran_e(i,l)=vtran_e(i,l)+ v(i,l)*esh*
     .                (paprs(i,l)-paprs(i,l+1))/RG
         vtran_q(i,l)=vtran_q(i,l)+ v(i,l)*q(i,l)
     .                *(paprs(i,l)-paprs(i,l+1))/RG
      ENDDO
      ENDDO
c
      RETURN
      END
