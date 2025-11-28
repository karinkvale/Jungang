! source file: /Users/dkeller/Desktop/UVic_ESCM/2.9/source/mom/timeavgs.h
!====================== include file "timeavgs.h" ======================

!     imtav =  # of longitudes used for the time averages grid
!     jmtav =  # of latitudes used for the time averages grid
!     kmav  =  # of levels used for the time averages grid

      integer imtav, jmtav, kmav
      parameter (imtav=imt, jmtav=jmt-2, kmav=km)
      real ta_vetiso, ta_vntiso, ta_vbtiso
      common /ta_gm_r/ ta_vetiso(imt,km,jmt), ta_vntiso(imt,km,jmt)
      common /ta_gm_r/ ta_vbtiso(imt,km,jmt)
      real ta_dc14
      common /ta_dc14/ ta_dc14(imt,km,jmt)
      real ta_rnpp, ta_rgraz, ta_rmorp, ta_rmorpt, ta_rmorz, ta_rexcr
      real ta_rremi, ta_rexpo, ta_rnpp_D, ta_rgraz_D, ta_rmorp_D
      real ta_rnfix, ta_rdeni, ta_rexpocal, ta_rprocal, ta_rgraz_Z
      real ta_rgraz_Det, ta_ravej, ta_ravej_D, ta_rgmax, ta_rno3P
      real ta_rpo4P, ta_rpo4_D
#if defined O_zoop_det
      real ta_rremiz, ta_rexpoz
#endif
#if defined O_npzd_iron
      real ta_rfeorgads, ta_rfecol
      real ta_rchl, ta_rchl_D
      real ta_rdeffe, ta_rremife, ta_rexpofe, ta_rfeprime
      real ta_rfesed, ta_rbfe
# if defined O_kk_caco3tr
      real ta_rfeorgads_ca
# endif
#endif
#if defined O_plastic
      real ta_rexpomp_a, ta_rexpomp_p, ta_ragg_upt, ta_rzoop_upt
      real ta_ragg_rel, ta_rzoop_rel
#endif
#if defined O_kk_ballast
      real ta_rremi_B, ta_rexpo_B, ta_rgraz_Det_B
#endif
#if defined O_kk_caco3tr
      real ta_romca, ta_rco3, ta_rco3_sat, ta_rdel_sat
      real ta_rexpocaco3, ta_rdissl, ta_rimpocaco3
      real ta_rold_diss, ta_rdissl_new
      real ta_rcalatt
#endif
#if defined O_cal_kk
      real ta_rmorpt_C
      real ta_rnpp_C, ta_rgraz_C, ta_rmorp_C
#endif
#if defined O_kk_diat
      real ta_rmorpt_Diat
      real ta_rnpp_Diat, ta_rgraz_Diat, ta_rmorp_Diat
#endif
#if defined O_kk_si
      real ta_rproopl, ta_rrivsil, ta_rexpoopl, ta_rdisopl
      real ta_rbur_sil
      common /ta_npzd_r/ ta_rproopl(imt,jmt)
      common /ta_npzd_r/ ta_rrivsil(imt,jmt)
      common /ta_npzd_r/ ta_rexpoopl(imt,km,jmt)
      common /ta_npzd_r/ ta_rdisopl(imt,km,jmt)
      common /ta_npzd_r/ ta_rbur_sil(imt,km,jmt)
#endif
#if defined O_kk_caco3tr
      common /ta_npzd_r/ ta_romca(imt,km,jmt)
      common /ta_npzd_r/ ta_rco3(imt,km,jmt)
      common /ta_npzd_r/ ta_rco3_sat(imt,km,jmt)
      common /ta_npzd_r/ ta_rdel_sat(imt,km,jmt)
#endif
#if defined O_cal_kk
      common /ta_npzd_r/ ta_rnpp_C(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rgraz_C(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorp_C(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorpt_C(imt,kpzd,jmt)
#endif
#if defined O_kk_diat
      common /ta_npzd_r/ ta_rnpp_Diat(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rgraz_Diat(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorp_Diat(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorpt_Diat(imt,kpzd,jmt)
#endif
#if defined O_kk_caco3tr
      common /ta_npzd_r/ ta_rexpocaco3(imt,km,jmt)
      common /ta_npzd_r/ ta_rdissl(imt,km,jmt)
      common /ta_npzd_r/ ta_rdissl_new(imt,km,jmt)
      common /ta_npzd_r/ ta_rold_diss(imt,km,jmt)
      common /ta_npzd_r/ ta_rcalatt(imt,kpzd,jmt)
#endif
#if defined O_kk_ballast
      common /ta_npzd_r/ ta_rremi_B(imt,km,jmt)
      common /ta_npzd_r/ ta_rexpo_B(imt,km,jmt)
      common /ta_npzd_r/ ta_rgraz_Det_B(imt,kpzd,jmt)
#endif
      common /ta_npzd_r/ ta_rnpp(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rgraz(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorp(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorpt(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorz(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rexcr(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rremi(imt,km,jmt)
      common /ta_npzd_r/ ta_rexpo(imt,km,jmt)
      common /ta_npzd_r/ ta_rexpocal(imt,km,jmt)
#if defined O_kk_caco3tr
      common /ta_npzd_r/ ta_rprocal(imt,km,jmt)
#else
      common /ta_npzd_r/ ta_rprocal(imt,jmt)
#endif
#if defined O_zoop_det
      common /ta_npzd_r/ ta_rremiz(imt,km,jmt)
      common /ta_npzd_r/ ta_rexpoz(imt,km,jmt)
#endif
      common /ta_npzd_r/ ta_rnpp_D(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rgraz_D(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorp_D(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rnfix(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rdeni(imt,km,jmt)
      common /ta_npzd_r/ ta_rgraz_Det(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rgraz_Z(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_ravej(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_ravej_D(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rgmax(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rno3P(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rpo4P(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rpo4_D(imt,kpzd,jmt)
#if defined O_npzd_iron
      common /ta_npzd_r/ ta_rfeorgads(imt,kpzd,jmt)
# if defined O_kk_caco3tr
      common /ta_npzd_r/ ta_rfeorgads_ca(imt,kpzd,jmt)
# endif
      common /ta_npzd_r/ ta_rfecol(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rdeffe(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rremife(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rexpofe(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rfeprime(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rfesed(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rbfe(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rchl(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rchl_D(imt,kpzd,jmt)
#endif
#if defined O_plastic
      common /ta_npzd_r/ ta_rexpomp_a(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rexpomp_p(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_ragg_upt(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rzoop_upt(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_ragg_rel(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rzoop_rel(imt,kpzd,jmt)
#endif
      integer nta_conv
      common /ta_conv_i/ nta_conv

      real ta_totalk, ta_vdepth, ta_pe
      common /ta_conv_r/ ta_totalk(imt,jmt), ta_vdepth(imt,jmt)
      common /ta_conv_r/ ta_pe(imt,jmt)
      integer nta_sscar
      common /ta_car_i/ nta_sscar

      real ta_sspH, ta_ssCO3, ta_ssOc, ta_ssOa, ta_sspCO2
      common /ta_car_r/ ta_sspH(imt,jmt), ta_ssCO3(imt,jmt)
      common /ta_car_r/ ta_ssOc(imt,jmt), ta_ssOa(imt,jmt)
      common /ta_car_r/ ta_sspCO2(imt,jmt)
