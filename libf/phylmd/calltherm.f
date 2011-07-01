!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/calltherm.F,v 1.2 2004/12/10 11:27:46 lmdzadmin Exp $
!
      subroutine calltherm(dtime
     s      ,pplay,paprs,pphi
     s      ,u_seri,v_seri,t_seri,q_seri
     s      ,d_u_ajs,d_v_ajs,d_t_ajs,d_q_ajs
     s      ,fm_therm,entr_therm)

      use dimens_m
      use dimphy
      use ctherm
      implicit none

      REAL, intent(in):: dtime

      REAL u_seri(klon,klev),v_seri(klon,klev)
      REAL t_seri(klon,klev),q_seri(klon,klev)
      REAL, intent(in):: paprs(klon,klev+1)
      REAL, intent(in):: pplay(klon,klev)
      REAL, intent(in):: pphi(klon,klev)

CFH Update Thermiques
      REAL d_t_ajs(klon,klev), d_q_ajs(klon,klev)
      REAL d_u_ajs(klon,klev),d_v_ajs(klon,klev)
      real fm_therm(klon,klev+1),entr_therm(klon,klev)


c variables locales
      REAL d_t_the(klon,klev), d_q_the(klon,klev)
      REAL d_u_the(klon,klev),d_v_the(klon,klev)
c
      real zfm_therm(klon,klev+1),zentr_therm(klon,klev),zdt
      save zentr_therm,zfm_therm

      integer i,k, isplit

*********************************************************

c  Modele du thermique
c  ===================
c         print*,'thermiques: WARNING on passe t au lieu de t_seri'
       print*,'avant isplit ',nsplit_thermals


         fm_therm(:,:)=0.
         entr_therm(:,:)=0.

c   tests sur les valeurs negatives de l'eau
         do k=1,klev
            do i=1,klon
               if (.not.q_seri(i,k).ge.0.) then
                   print*,'WARN eau<0 avant therm i=',i,'  k=',k
     s         ,' dq,q',d_q_the(i,k),q_seri(i,k)
                  q_seri(i,k)=1.e-15
               endif
            enddo
         enddo


         zdt=dtime/float(nsplit_thermals)
         do isplit=1,nsplit_thermals

            CALL thermcell(klon,klev,zdt
     s      ,pplay,paprs,pphi
     s      ,u_seri,v_seri,t_seri,q_seri
     s      ,d_u_the,d_v_the,d_t_the,d_q_the
     s      ,zfm_therm,zentr_therm
     s      ,r_aspect_thermals,l_mix_thermals,w2di_thermals
     s      ,tho_thermals)

c  transformation de la derivee en tendance
            d_t_the(:,:)=d_t_the(:,:)*dtime/float(nsplit_thermals)
            d_u_the(:,:)=d_u_the(:,:)*dtime/float(nsplit_thermals)
            d_v_the(:,:)=d_v_the(:,:)*dtime/float(nsplit_thermals)
            d_q_the(:,:)=d_q_the(:,:)*dtime/float(nsplit_thermals)
            fm_therm(:,:)=fm_therm(:,:)
     s      +zfm_therm(:,:)/float(nsplit_thermals)
            entr_therm(:,:)=entr_therm(:,:)
     s       +zentr_therm(:,:)/float(nsplit_thermals)
            fm_therm(:,klev+1)=0.



c  accumulation de la tendance
            d_t_ajs(:,:)=d_t_ajs(:,:)+d_t_the(:,:)
            d_u_ajs(:,:)=d_u_ajs(:,:)+d_u_the(:,:)
            d_v_ajs(:,:)=d_v_ajs(:,:)+d_v_the(:,:)
            d_q_ajs(:,:)=d_q_ajs(:,:)+d_q_the(:,:)

c  incrementation des variables meteo
            t_seri(:,:) = t_seri(:,:) + d_t_the(:,:)
            u_seri(:,:) = u_seri(:,:) + d_u_the(:,:)
            v_seri(:,:) = v_seri(:,:) + d_v_the(:,:)
            q_seri(:,:) = q_seri(:,:) + d_q_the(:,:)

c   tests sur les valeurs negatives de l'eau
            DO k = 1, klev
            DO i = 1, klon
               if (.not.q_seri(i,k).ge.0.) then
                   print*,'WARN eau<0 apres therm i=',i,'  k=',k
     s         ,' dq,q',d_q_the(i,k),q_seri(i,k)
                  q_seri(i,k)=1.e-15
               endif
            ENDDO
            ENDDO

         enddo ! isplit

      return

      end
