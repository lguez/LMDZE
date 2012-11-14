module aeropt_m

  IMPLICIT none

contains

  SUBROUTINE aeropt(pplay, paprs, t_seri, msulfate, RHcl, tau_ae, piz_ae, &
       cg_ae, ai)

    ! From LMDZ4/libf/phylmd/aeropt.F, v 1.1.1.1 2004/05/19 12:53:09

    ! Author: Olivier Boucher
    ! Calculate aerosol optical properties.

    USE dimphy, ONLY: klev, klon
    USE suphec_m, ONLY: rd, rg

    REAL, intent(in):: pplay(klon, klev), paprs(klon, klev + 1)
    REAL, intent(in):: t_seri(klon, klev)

    REAL, intent(in):: msulfate(klon, klev)
    ! masse sulfate ug SO4 / m3 (ug / m^3)

    REAL, intent(in):: RHcl(klon, klev) ! humidité relative ciel clair

    REAL, intent(out):: tau_ae(klon, klev, 2) ! épaisseur optique aérosols
    REAL, intent(out):: piz_ae(klon, klev, 2) ! single scattering albedo aerosol
    REAL, intent(out):: cg_ae(klon, klev, 2) ! asymmetry parameter aerosol
    REAL, intent(out):: ai(klon) ! POLDER aerosol index

    ! Local:

    INTEGER i, k, inu
    INTEGER RH_num
    INTEGER, PARAMETER:: nbre_RH = 12

    REAL:: RH_tab(nbre_RH) = (/0., 10., 20., 30., 40., 50., 60., 70., 80., &
         85., 90., 95./)

    REAL DELTA, rh
    REAL, PARAMETER:: RH_MAX = 95.
    REAL zrho, zdz
    REAL taue670(klon) ! épaisseur optique aerosol absorption 550 nm
    REAL taue865(klon) ! épaisseur optique aerosol extinction 865 nm
    REAL alpha_aer_sulfate(nbre_RH, 5) ! unit m2 / g SO4
    REAL alphasulfate

    ! Propriétés optiques

    REAL alpha_aer(nbre_RH, 2) ! unit m2 / g SO4
    REAL cg_aer(nbre_RH, 2)
    DATA alpha_aer/.500130E+01, .500130E+01, .500130E+01, &
         .500130E+01, .500130E+01, .616710E+01, &
         .826850E+01, .107687E+02, .136976E+02, &
         .162972E+02, .211690E+02, .354833E+02, &
         .139460E+01, .139460E+01, .139460E+01, &
         .139460E+01, .139460E+01, .173910E+01, &
         .244380E+01, .332320E+01, .440120E+01, &
         .539570E+01, .734580E+01, .136038E+02 /
    DATA cg_aer/.619800E+00, .619800E+00, .619800E+00, &
         .619800E+00, .619800E+00, .662700E+00, &
         .682100E+00, .698500E+00, .712500E+00, &
         .721800E+00, .734600E+00, .755800E+00, &
         .545600E+00, .545600E+00, .545600E+00, &
         .545600E+00, .545600E+00, .583700E+00, &
         .607100E+00, .627700E+00, .645800E+00, &
         .658400E+00, .676500E+00, .708500E+00 /
    DATA alpha_aer_sulfate/ &
         4.910, 4.910, 4.910, 4.910, 6.547, 7.373, &
         8.373, 9.788, 12.167, 14.256, 17.924, 28.433, &
         1.453, 1.453, 1.453, 1.453, 2.003, 2.321, &
         2.711, 3.282, 4.287, 5.210, 6.914, 12.305, &
         4.308, 4.308, 4.308, 4.308, 5.753, 6.521, &
         7.449, 8.772, 11.014, 12.999, 16.518, 26.772, &
         3.265, 3.265, 3.265, 3.265, 4.388, 5.016, &
         5.775, 6.868, 8.745, 10.429, 13.457, 22.538, &
         2.116, 2.116, 2.116, 2.116, 2.882, 3.330, &
         3.876, 4.670, 6.059, 7.327, 9.650, 16.883/

    !----------------------------------------------------------------------

    taue670 = 0.
    taue865 = 0.

    DO k = 1, klev
       DO i = 1, klon
          if (t_seri(i, k).eq.0) write (*, *) 'aeropt T ', i, k, t_seri(i, k)
          if (pplay(i, k).eq.0) write (*, *) 'aeropt p ', i, k, pplay(i, k)
          zrho = pplay(i, k) / t_seri(i, k) / RD ! kg / m3
          zdz = (paprs(i, k) - paprs(i, k + 1)) / zrho / RG ! m
          rh = MIN(RHcl(i, k) * 100., RH_MAX)
          RH_num = INT(rh / 10. + 1.)
          IF (rh < 0.) then
             print *, 'aeropt: RH < 0 not possible'
             STOP 1
          end IF
          IF (rh > 85.) RH_num = 10
          IF (rh > 90.) RH_num = 11
          DELTA = (rh - RH_tab(RH_num)) / (RH_tab(RH_num + 1) - RH_tab(RH_num))

          do inu = 1, 2
             tau_ae(i, k, inu) = alpha_aer(RH_num, inu) + DELTA &
                  * (alpha_aer(RH_num + 1, inu) - alpha_aer(RH_num, inu))
             tau_ae(i, k, inu) = tau_ae(i, k, inu) * msulfate(i, k) * zdz * 1e-6
             piz_ae(i, k, inu) = 1.
             cg_ae(i, k, inu) = cg_aer(RH_num, inu) &
                  + DELTA * (cg_aer(RH_num + 1, inu) - cg_aer(RH_num, inu))
          end do

          alphasulfate = alpha_aer_sulfate(RH_num, 4) + DELTA &
               * (alpha_aer_sulfate(RH_num + 1, 4) &
               - alpha_aer_sulfate(RH_num, 4)) ! m2 / g

          taue670(i) = taue670(i) + alphasulfate * msulfate(i, k) * zdz * 1e-6

          alphasulfate = alpha_aer_sulfate(RH_num, 5) + DELTA &
               * (alpha_aer_sulfate(RH_num + 1, 5) &
               - alpha_aer_sulfate(RH_num, 5)) ! m2 / g

          taue865(i) = taue865(i) + alphasulfate * msulfate(i, k) * zdz * 1e-6
       ENDDO
    ENDDO

    DO i = 1, klon
       ai(i) = (- log(MAX(taue670(i), 0.0001) / MAX(taue865(i), 0.0001)) &
            / log(670. / 865.)) * taue865(i)
    ENDDO

  END SUBROUTINE aeropt

end module aeropt_m
