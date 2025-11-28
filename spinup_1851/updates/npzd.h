!====================== include file "npzd.h" =========================

!   variables for npzd model

!   ntnpzd   = number of npzd tracers (moved to size.h)
!   nbio     = number of npzd timesteps per ocean timestep
!   trcmin   = minimum tracer for flux calculations
!   alpha    = initial slope P-I curve [(W/m^2)^(-1)/day]
!   kw       = light attenuation due to water [1/m]
!   kc       = light attenuation by phytoplankton [1/(m*mmol m-3)]
!   ki       = light attenuation through sea ice & snow
!   abio     = maximum growth rate parameter [1/day]
!   bbio     = b
!   cbio     = [1/deg_C]
!   k1n      = half saturation constant for N uptake [mmol m-3]
!   nup      = specific mortality rate (Phytoplankton) [day-1]
!   gamma1   = assimilation efficiency (zpk)
!   gbio0     = maximum grazing rate at 0 deg C [day-1]
!   nuz      = quadratic mortality (zpk)
!   nud0     = remineralization rate [day-1]
!   LFe      = Iron limitation
!   wd       = sinking speed of detritus [m day-1]
!   ztt      = depth to top of grid cell [cm]
!   rkwz     = reciprical of light attenuation times grid depth
!   par      = fraction of photosythetically active radiation
!   dtnpzd   = time step of biology
!   capr     = carbonate to carbon production ratio
!   dcaco3   = remineralisation depth of calcite [cm]
!   rcak     = array used in calculating calcite remineralization
!   rcab     = array used in calculating bottom calcite remineralization
!   nupt0    = specific mortality rate (Phytoplankton) [1/day]
!   wd0      = sinking speed of detritus at surface [m/day]
!   wdd      = inc. sinking speed of detritus at surface [m/day]
!   k1p      = half saturation constant for P uptake
!   jdiar    = factor reducing the growth rate of diazotrophs
!   redctn   = C/N Redfield ratio (includes mol to mmol conversion)
!   redctp   = C/P Redfield ratio (includes mol to mmol conversion)
!   redptn   = P/N Redfield ratio
!   redntp   = N/P Redfield ratio
!   redotn   = O/N Redfield ratio (includes mol to mmol conversion)
!   redotp   = O/P Redfield ratio (includes mol to mmol conversion)
!   rnbio    = reciprical of nbio
!   rdtts    = reciprical of dtts [s-1]
!   dtbio    = npzd time step [s]
!   rnpp     = rate of net primary production [nmol cm-3 s-1]
!   rgraz    = rate of grazing [nmol cm-3 s-1]
!   rmorp    = rate of mortality of phytoplankton [nmol cm-3 s-1]
!   rmorz    = rate of mortality of zooplankton [nmol cm-3 s-1]
!   rremi    = rate of remineralization [nmol cm-3 s-1]
!   rexcr    = rate of excretion [nmol cm-3 s-1]
!   rexpo    = rate of export through the bottom [nmol cm-3 s-1]
!   rnpp_D   = npp for diazotraphs [nmol cm-3 s-1]
!   rgraz_D  = rgraz for diazotraphs [nmol cm-3 s-1]
!   rmorp_D  = rmorp for diazotraphs [nmol cm-3 s-1]
!   rnfix    = rate of nitrogen fixation [nmol cm-3 s-1]
!   rdeni    = rate of denitrification [nmol cm-3 s-1]
!   kzoo     = half saturation constant for Z grazing
!   zprefP   = Z preference for grazing on P
!   zprefMP   = Z preference for grazing on MP
!   zprefDiaz   = Z preference for grazing on Diaz
!   zprefZ   = Z preference for grazing on other Z
!   zprefDet = Z preference for grazing on Detritus
!   rgraz_Det = rate of grazing on Detritus [nmol cm-3 s-1]
!   rgraz_Z   = rate of grazing on other Zooplankton [nmol cm-3 s-1]
!   geZ      = growth efficiency of zooplankton
!   silwflx  = silica weathering flux
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! New diagnostic output
!   ravej      = light-dependant growth rate of P
!   ravej_D    = light-dependant growth rate of Diaz
!   rgmax      = temp-dependant growth rate of zoo
!   rno3P      = nitrate-dependant growth rate of P
!   rpo4P       = phosphate-dependant growth rate of P
!   rpo4_D     = phosphate-dependant growth rate of D
!
!   fe_dissolved = dissolved iron concentration
!   kfe = Fe limitation half saturation parameter
!   kfe_D = Fe limitation half sat. param. for diaz.
!   kmfe = number of depth levels for the iron mask
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Dynamic iron cycle
!   kfeleq = Fe-ligand stability constant [m^3/(mmol ligand)]
!   lig = Ligand concentration  [mmol/m^3]
!   thetamaxhi = Maximum Chl:C ratio, abundant iron [gChl/(gC)]
!   thetamaxlo = Maximum Chl:C ratio, extreme iron limitation [gChl/(gC)]
!   alphamax = Maximum initial slope in PI-curve [day^-1 (W/m^2)^-1 (gChl/(gC))^-1]
!   alphamin = Minimum intital slope in PI-curve [day^-1 (W/m^2)^-1 (gChl/(gC))^-1]
!   mc = Molar mass of carbon [g/mol]
!   fetopsed = Fe:P ratio for sedimentary iron source [molFe/molP]
!   o2min = Minimum O2 concentration for aerobic respiration [mmolO_2/m^3]
!   kfeorg = Organic-matter dependent scavenging rate [(m^3/(gC s))^0.58]
!   kfecol = Colloidal Fe production rate (inorganic scavenging) [s^-1]


      integer nbio, kmfe
      common /npzd_i/ nbio(km)

      real trcmin
      parameter (trcmin=5e-12)

      real alpha, kw, kc, ki, abio, bbio, cbio, k1n, nup, gamma1, gbio0
      real epsbio, nuz, nud0, LFe, wd, ztt, rkwz, par, dtnpzd
      real capr, dcaco3, rcak, rcab, nupt0, wd0, k1p, jdiar, redctn
      real redctp, redptn, redntp, redotn, redotp, rnbio, rdtts, dtbio
      real rnpp, rgraz, rmorp, rmorpt, rmorz, rremi, rexcr, rexpo
      real rnpp_D, rgraz_D, rmorp_D, rnfix, rdeni, kzoo, zprefP
      real zprefDiaz, zprefZ, zprefDet, rgraz_Det, rgraz_Z, geZ, kfe
      real ravej, ravej_D, rgmax, rno3P, rpo4P, rpo4_D, kfe_D
      real kfemax, kfemin, pmax, wdd, alpha_D, zprefC, zprefDiat
      real zprefMP, ro2n, o2t
#if defined O_phyt_albedo
      real phin, caco3in
      common /npzd_r/ phin(imt,jmt), caco3in(imt,jmt)
#endif
# if defined O_kk_ballast
      real rgraz_Det_B, rexpo_B,rremi_B,bapr
# endif
# if defined O_cal_kk
      real abioc , k1n_C, k1p_C, rnpp_C, nuc, nuct0, alpha_C
      real rgraz_C, rmorp_C, rmorpt_C
      real kfe_C, wdc
# endif
# if defined O_kk_caco3tr
      real rcalpro, wc, wc0, dissk0, rdissl, kcal
      real rexpocaco3, rcaldiss, rcalatt, rimpocaco3
      real kc_c, romca, rco3, rco3_sat, rdel_sat, rdissl_new
# endif
#if defined O_npzd_iron
      real kfeleq, lig, thetamaxhi, alphamax, alphamin
      real thetamaxlo, mc, fetopsed, o2min, kfeorg, rfeton
      real kfecol
      real rfeorgads, rdeffe, rremife, rexpofe, rfeprime
      real rfesed, rbfe, rfecol
      real fe_hydr
# if defined O_kk_caco3tr
      real kfeorg_ca, rfeorgads_ca
# endif
# if defined O_npzd_chl
      real rchl, rchl_D 
# endif
#endif
#if defined O_kk_diat
      real abiodiat, k1n_Diat, k1p_Diat, rnpp_Diat, nu_diat, nudt0
      real alpha_Diat, k1si, si_msk
      real rgraz_Diat, rmorp_Diat, rmorpt_Diat
      real kfe_Diat
# if defined O_kk_si
      real sipr, sipr0, ropk, dopal, si_sol, si_h_sol, si_dis
      real silwflx, globalsilwflx, ws0
      real si_hydr
# endif
#endif
# if defined O_plastic
      real wmp0, wmp, pinp, pbc, pagg, prise, kp, rexpomp_a
      real rexpomp_p, ragg_upt, rzoop_upt, ping
      real ragg_rel, rzoop_rel
# endif
# if defined O_zoop_det
      real rremiz, rexpoz, wdz, wdz0, rpoo
# endif
#if defined O_npzd_fe_limitation
      real fe_dissolved
      common /fe_dissolved_r/ fe_dissolved(imt,jmt,km,12)
      common /fe_dissolved_i/ kmfe
#endif
      common /npzd_r/ alpha, kw, kc, ki, abio, bbio, cbio, k1n, nup
      common /npzd_r/ gamma1, gbio0, epsbio, nuz, nud0, LFe
      common /npzd_r/ wd(km), ztt(km), rkwz(km), par, dtnpzd, capr
      common /npzd_r/ dcaco3, rcak(km), rcab(km), nupt0, wd0, k1p
      common /npzd_r/ jdiar, redctn, redctp, redptn, redntp, redotn
      common /npzd_r/ redotp, rnbio(km), rdtts(km), dtbio(km), geZ
      common /npzd_r/ kzoo, zprefP, zprefDiaz, zprefZ, zprefDet
      common /npzd_r/ kfe, kfe_D, wdd, alpha_D, zprefMP, zprefDiat
      common /npzd_r/ zprefC, ro2n, o2t
# if defined O_kk_si
      common /npzd_r/ sipr, sipr0, ropk(km), dopal, ws0
      common /npzd_r/ silwflx, globalsilwflx, si_sol, si_h_sol, si_dis
# endif
# if defined O_plastic
      common /npzd_r/ wmp0, wmp(km), pinp, pbc, pagg, prise, kp, ping
      common /npzd_r/ rexpomp_a(imt,km,jsmw:jemw)
      common /npzd_r/ rexpomp_p(imt,km,jsmw:jemw)
      common /npzd_r/ ragg_upt(imt,km,jsmw:jemw)
      common /npzd_r/ rzoop_upt(imt,km,jsmw:jemw)
      common /npzd_r/ ragg_rel(imt,km,jsmw:jemw)
      common /npzd_r/ rzoop_rel(imt,km,jsmw:jemw)
# endif
# if defined O_zoop_det
      common /npzd_r/ wdz0, wdz(km), rpoo
      common /npzd_r/ rremiz(imt,km,jsmw:jemw)
      common /npzd_r/ rexpoz(imt,km,jsmw:jemw)
# endif
# if defined O_cal_kk
      common /npzd_r/ abioc, k1n_C, k1p_C, nuc, nuct0, alpha_C
      common /npzd_r/ kfe_C
# endif
# if defined O_kk_ballast
      common /npzd_r/ bapr
# endif
# if defined O_kk_caco3tr
      common /npzd_r/ wc0, dissk0, wc(km)
      common /npzd_r/ kc_c, kcal, wdc
# endif
#if defined O_npzd_iron
      common /npzd_r/ kfeleq, alphamax, alphamin
      common /npzd_r/ thetamaxhi, thetamaxlo, lig, fetopsed, o2min
      common /npzd_r/ mc, kfeorg, rfeton, kfecol
      common /npzd_r/ kfemax, kfemin, pmax, fe_hydr(imt,jmt,km)
#endif
#if defined O_kk_diat
      common /npzd_r/ abiodiat, k1n_Diat, k1p_Diat, nu_diat, nudt0
      common /npzd_r/ alpha_Diat, k1si
      common /npzd_r/ kfe_Diat, si_msk(imt,jmt,km)
#endif
#if defined O_kk_si
      common /npzd_r/ si_hydr(imt,jmt,km)
#endif
#if defined O_save_npzd
      common /npzd_r/ rnpp(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rgraz(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rmorp(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rmorpt(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rmorz(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rexcr(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rremi(imt,km,jsmw:jemw)
      common /npzd_r/ rexpo(imt,km,jsmw:jemw)
      common /npzd_r/ rgraz_Det(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rgraz_Z(imt,kpzd,jsmw:jemw)
# if defined O_kk_ballast
      common /npzd_r/ rgraz_Det_B(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rremi_B(imt,km,jsmw:jemw)
      common /npzd_r/ rexpo_B(imt,km,jsmw:jemw)
# endif
#if defined O_kk_caco3tr
      common /npzd_r/ rcalpro(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rcaldiss(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rcalatt(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rdissl(imt,km,jsmw:jemw)
      common /npzd_r/ rdissl_new(imt,km,jsmw:jemw)
      common /npzd_r/ rexpocaco3(imt,km,jsmw:jemw)
      common /npzd_r/ rimpocaco3(imt,km,jsmw:jemw)
      common /npzd_r/ romca(imt,km,jsmw:jemw)
      common /npzd_r/ rco3(imt,km,jsmw:jemw)
      common /npzd_r/ rco3_sat(imt,km,jsmw:jemw)
      common /npzd_r/ rdel_sat(imt,km,jsmw:jemw)
#endif
#if defined O_cal_kk
      common /npzd_r/ rnpp_C(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rgraz_C(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rmorp_C(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rmorpt_C(imt,kpzd,jsmw:jemw)
#endif
# if defined O_kk_diat
      common /npzd_r/ rnpp_Diat(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rgraz_Diat(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rmorp_Diat(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rmorpt_Diat(imt,kpzd,jsmw:jemw)
# endif
# if defined O_npzd_nitrogen
      common /npzd_r/ rnpp_D(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rgraz_D(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rmorp_D(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rnfix(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rdeni(imt,km,jsmw:jemw)
# endif
# if defined O_npzd_iron
      common /npzd_r/ rremife(imt,km,jsmw:jemw)
      common /npzd_r/ rexpofe(imt,km,jsmw:jemw)
#  if defined O_npzd_iron_diagnostics
      common /npzd_r/ rfeorgads(imt,km,jsmw:jemw)
      common /npzd_r/ rdeffe(imt,km,jsmw:jemw)
      common /npzd_r/ rfeprime(imt,km,jsmw:jemw)
      common /npzd_r/ rfesed(imt,km,jsmw:jemw)
      common /npzd_r/ rbfe(imt,km,jsmw:jemw)
      common /npzd_r/ rfecol(imt,km,jsmw:jemw)
#  endif
# if defined O_kk_caco3tr
      common /npzd_r/ kfeorg_ca, rfeorgads_ca(imt,km,jsmw:jemw)
# endif
#  if defined O_npzd_chl
      common /npzd_r/ rchl(imt,km,jsmw:jemw)
#   if defined O_npzd_nitrogen
      common /npzd_r/ rchl_D(imt,km,jsmw:jemw)
#   endif
#  endif
# endif
# if defined O_npzd_extra_diagnostics
      common /npzd_r/ ravej(imt,kpzd,jsmw:jemw)
      common /npzd_r/ ravej_D(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rgmax(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rno3P(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rpo4P(imt,kpzd,jsmw:jemw)
      common /npzd_r/ rpo4_D(imt,kpzd,jsmw:jemw)
# endif
#endif


