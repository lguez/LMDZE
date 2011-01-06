!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/transp.F,v 1.1.1.1 2004/05/19 12:53:09 lmdzadmin Exp $
!
      SUBROUTINE transp (paprs,tsol,
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
      REAL utran_e(klon), utran_q(klon), vtran_e(klon), vtran_q(klon)
c
      INTEGER i, l
c     ------------------------------------------------------------------
      REAL geom(klon,klev), e
c     ------------------------------------------------------------------
      DO i = 1, klon
         utran_e(i) = 0.0
         utran_q(i) = 0.0
         vtran_e(i) = 0.0
         vtran_q(i) = 0.0
      ENDDO
c
      DO l = 1, klev
      DO i = 1, klon
         e = RCPD*t(i,l) + RLVTT*q(i,l) + geom(i,l)
         utran_e(i)=utran_e(i)+ u(i,l)*e*(paprs(i,l)-paprs(i,l+1))/RG
         utran_q(i)=utran_q(i)+ u(i,l)*q(i,l)
     .                         *(paprs(i,l)-paprs(i,l+1))/RG
         vtran_e(i)=vtran_e(i)+ v(i,l)*e*(paprs(i,l)-paprs(i,l+1))/RG
         vtran_q(i)=vtran_q(i)+ v(i,l)*q(i,l)
     .                         *(paprs(i,l)-paprs(i,l+1))/RG
      ENDDO
      ENDDO
c
      RETURN
      END
