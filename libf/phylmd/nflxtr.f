!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/nflxtr.F,v 1.1.1.1 2004/05/19 12:53:08 lmdzadmin Exp $
!
      SUBROUTINE nflxtr(pdtime,pmfu,pmfd,pen_u,pde_u,pen_d,pde_d,
     .                 paprs,x,dx) 
      use dimens_m
      use dimphy
      use SUPHEC_M
            use yoecumf 
      IMPLICIT NONE 
c=====================================================================
c Objet : Melange convectif de traceurs a partir des flux de masse 
c Date : 13/12/1996 -- 13/01/97
c Auteur: O. Boucher (LOA) sur inspiration de Z. X. Li (LMD),
c         Brinkop et Sausen (1996) et Boucher et al. (1996).
c ATTENTION : meme si cette routine se veut la plus generale possible, 
c             elle a herite de certaines notations et conventions du 
c             schema de Tiedtke (1993). 
c --En particulier, les couches sont numerotees de haut en bas !!!
c   Ceci est valable pour les flux
c   mais pas pour les entrees x, paprs !!!!
c --pmfu est positif, pmfd est negatif 
c --Tous les flux d'entrainements et de detrainements sont positifs 
c   contrairement au schema de Tiedtke d'ou les changements de signe!!!! 
c=====================================================================
c
c
      REAL, intent(in):: pdtime
c--les flux sont definis au 1/2 niveaux
c--pmfu(klev+1) et pmfd(klev+1) sont implicitement nuls
      REAL pmfu(klon,klev)  ! flux de masse dans le panache montant 
      REAL pmfd(klon,klev)  ! flux de masse dans le panache descendant
      REAL pen_u(klon,klev) ! flux entraine dans le panache montant
      REAL pde_u(klon,klev) ! flux detraine dans le panache montant
      REAL pen_d(klon,klev) ! flux entraine dans le panache descendant
      REAL pde_d(klon,klev) ! flux detraine dans le panache descendant

      REAL, intent(in):: paprs(klon,klev+1) ! pression aux 1/2 couches (bas en haut)
      REAL, intent(in):: x(klon,klev)        ! q de traceur (bas en haut) 
      REAL dx(klon,klev)     ! tendance de traceur  (bas en haut)
c
c--flux convectifs mais en variables locales
      REAL zmfu(klon,klev+1) 
      REAL zmfd(klon,klev+1) 
      REAL zen_u(klon,klev) 
      REAL zde_u(klon,klev)
      REAL zen_d(klon,klev) 
      REAL zde_d(klon,klev)
      real zmfe
c
c--variables locales      
c--les flux de x sont definis aux 1/2 niveaux 
c--xu et xd sont definis aux niveaux complets
      REAL xu(klon,klev)        ! q de traceurs dans le panache montant
      REAL xd(klon,klev)        ! q de traceurs dans le panache descendant
      REAL zmfux(klon,klev+1)   ! flux de x dans le panache montant
      REAL zmfdx(klon,klev+1)   ! flux de x dans le panache descendant
      REAL zmfex(klon,klev+1)   ! flux de x dans l'environnement 
      INTEGER i, k 
      REAL zmfmin
      PARAMETER (zmfmin=1.E-10)

c =========================================
c
c
c   Extension des flux UP et DN sur klev+1 niveaux
c =========================================
      do k=1,klev
         do i=1,klon
            zmfu(i,k)=pmfu(i,k)
            zmfd(i,k)=pmfd(i,k)
         enddo
      enddo
      do i=1,klon
         zmfu(i,klev+1)=0.
         zmfd(i,klev+1)=0.
      enddo

c--modif pour diagnostiquer les detrainements
c =========================================
c   on privilegie l'ajustement de l'entrainement dans l'ascendance.

      do k=1, klev
         do i=1, klon
            zen_d(i,k)=pen_d(i,k)
            zde_u(i,k)=pde_u(i,k)
            zde_d(i,k) =-zmfd(i,k+1)+zmfd(i,k)+zen_d(i,k)
            zen_u(i,k) = zmfu(i,k+1)-zmfu(i,k)+zde_u(i,k)
         enddo 
      enddo 
c
c--calcul des flux dans le panache montant
c =========================================
c
c Dans la premiere couche, on prend q comme valeur de qu
c
      do i=1, klon
         zmfux(i,1)=0.0 
      enddo
c
c Autres couches
      do k=1,klev
         do i=1, klon
            if ((zmfu(i,k+1)+zde_u(i,k)).lt.zmfmin) THEN
               xu(i,k)=x(i,k)
            else
               xu(i,k)=(zmfux(i,k)+zen_u(i,k)*x(i,k))
     s               /(zmfu(i,k+1)+zde_u(i,k))
            endif
            zmfux(i,k+1)=zmfu(i,k+1)*xu(i,k)
         enddo
      enddo
c
c--calcul des flux dans le panache descendant
c =========================================
c   
      do i=1, klon
         zmfdx(i,klev+1)=0.0 
      enddo
c
      do k=klev,1,-1
         do i=1, klon
            if ((zde_d(i,k)-zmfd(i,k)).lt.zmfmin) THEN
               xd(i,k)=x(i,k)
            else
               xd(i,k)=(zmfdx(i,k+1)-zen_d(i,k)*x(i,k)) /
     .               (zmfd(i,k)-zde_d(i,k))
            endif
            zmfdx(i,k)=zmfd(i,k)*xd(i,k)
         enddo
      enddo
c
c--introduction du flux de retour dans l'environnement
c =========================================
c
      do k=2, klev
         do i=1, klon
            zmfe=-zmfu(i,k)-zmfd(i,k)
            if (zmfe.le.0.) then
               zmfex(i,k)= zmfe*x(i,k)
            else
               zmfex(i,k)= zmfe*x(i,k-1)
            endif
         enddo
      enddo

      do i=1, klon
         zmfex(i,1)=0.
         zmfex(i,klev+1)=0.
      enddo
c
c--calcul final des tendances
c
      do k=1, klev
         do i=1, klon
            dx(i,k)=RG/(paprs(i,k)-paprs(i,k+1))*pdtime*
     .                      ( zmfux(i,k) - zmfux(i,k+1) +
     .                        zmfdx(i,k) - zmfdx(i,k+1) +
     .                        zmfex(i,k) - zmfex(i,k+1) )
         enddo
      enddo
c
      return 
      end
