!======================= include file "mw.h" ===========================

!                     M E M O R Y    W I N D O W

!                         Data layout on disk

!     On disk, a latitude consists of a row of prognostic quantities.
!     Latitudes are ordered from south to north on two logical disk
!     units: one for time level "tau" and one for time level "tau-1".
!     Newly computed "tau+1" quantities are written to the "tau-1" disk
!     which serves double duty as a "tau+1" disk.

!                         Data layout in memory

!     A memory window "MW" is constructed to hold a group of "jmw"
!     adjacent latitude rows of prognostic quantities from disks "tau-1"
!     and "tau".  Parameter "jmw" controls the size of the MW and can be
!     set anywhere from a minimum of three latitude rows to all "jmt"
!     latitudes.  In addition, the MW holds diagnostic quantities
!     (density, hydrostatic pressure gradients, and advective
!     velocities) along with work arrays for constructing intermediate
!     fluxes used in solving the tracer and momentum equations. A
!     latitude row in the MW is referred to by index "j" and corresponds
!     to the physical latitude row "jrow" on disk.

!                   Data flow between disk and memory

!     The MW is loaded with prognostic data from the first "jww"
!     latitude rows on disks "tau" and "tau-1". As the tracer and
!     momentum equations are solved for rows j=2 through jmw-1 in the
!     MW, the solutions are written to disk "tau+1". When finished, the
!     MW is moved northward by moving quantities from the northernmost
!     two rows into the southernmost two rows within the MW. The
!     remaining MW rows are then loaded with more latitude rows from
!     disk. The process continues until all latitude rows 2 through
!     jmt-1 on disk "tau+1" have been updated.

!=======================================================================

!     taum1 = tau-1 time level for variables in MW
!     tau   = tau   time level for variables in MW
!     taup1 = tau+1 time level for variables in MW

      integer taum1, tau, taup1
      common /mw_i/ taum1, tau, taup1

!-----------------------------------------------------------------------
!     MW arrays for prognostic equations:
!-----------------------------------------------------------------------

!     u(i,k,j,n,tau) = total velocity where:
!      i   = index for longitude
!      k   = index for depth
!      j   = index for latitude row within MW
!      n   = component (1 = zonal, 2 = meridional)
!      tau = time level (tau-1, tau, tau+1)

!      (only internal modes are on disk and at tau+1 in the MW)

!     t(i,k,j,n,tau) = tracer where:
!      i   = index for longitude
!      k   = index for depth
!      j   = index for latitude row within MW
!      n   = component (1 = temperature, 2 = salinity)
!            if nt > 2 then other tracers are allowed.
!      tau = time level (tau-1, tau, tau+1)

!     note: temperature is potential temperature in degrees Celsius and
!           salinity is in "model units", the deviation from 0.035 grams
!           of salt/cm**3 of water, or, assuming a water density of
!           1 gram/cm**3, the deviation from 0.035 g of salt/g of water.
!           one can convert model units to the more common units of parts
!           per thousand (ppt) by adding 0.035 grams/cm**3 to the model
!           units and then multiplying by 1000.

      real u,t
      common /mw_r/ u(imt,km,jmw,2,-1:1), t(imt,km,jmw,nt,-1:1)

!     indicies for ocean tracer array
!       iidealage = index for ideal age tracer
!       itemp  = index for temperature
!       isalt  = index for salinity
!       idic   = index for dissolved inorganic carbon
!       imp    = index for microplastic
!       impa   = index for microplastic in phytoplankton aggregates
!       impp   = index for microplastic in zooplankton poo
!       ic14   = index for carbon 14
!       icfc11 = index for cfc11
!       icfc12 = index for cfc12
!       ialk   = index for alkalinity
!       ipo4   = index for phosphate
!       iphyt  = index for phytoplankton
!       izoop  = index for zooplankton
!       idetr  = index for detritus
#if defined O_npzd_iron
!       idfe   = index for iron
!       idetrfe= index for particulate iron
#endif
!       io2    = index for oxygen
!       ino3   = index for nitrate
!       idiaz  = index for diazotrophs
!       mapt   = map for ocean tracer names

      character(10) :: mapt
      common /mw_c/ mapt(nt)

      integer itemp, isalt, idic, ic14, icfc11, icfc12, ialk, ipo4
      integer iphyt, izoop, idetr, io2, ino3, iidealage
#if !defined O_npzd_simple
      integer idiaz
#endif
#if defined O_npzd_iron
      integer idfe, idetrfe
#endif
#if defined O_plastic
      integer imp, impa, impp
#endif
#if defined O_zoop_det
      integer idetrz
#endif
#if defined O_kk_ballast
      integer idetr_B
#endif
#if defined O_cal_kk
      integer icocc
#endif
#if defined O_kk_caco3tr
      integer icaco3
#endif
#if defined O_kk_diat
      integer idiat
#endif
#if defined O_kk_si
      integer isil
#endif
      common /mw_i/ itemp, isalt, idic, ic14, icfc11, icfc12, ialk
      common /mw_i/ ipo4, iphyt, izoop, idetr, io2, ino3
#if !defined O_npzd_simple
      common /mw_i/ idiaz
#endif
      common /mw_i/ iidealage
#if defined O_npzd_iron
      common /mw_i/ idfe, idetrfe
#endif
#if defined O_plastic
      common /mw_i/ imp, impa, impp
#endif
#if defined O_zoop_det
      common /mw_i/ idetrz
#endif
#if defined O_kk_ballast
      common /mw_i/ idetr_B
#endif
#if defined O_cal_kk
      common /mw_i/ icocc
#endif
#if defined O_kk_caco3tr
      common /mw_i/ icaco3
#endif
#if defined O_kk_diat
      common /mw_i/ idiat
#endif
#if defined O_kk_si
      common /mw_i/ isil
#endif
#if defined O_isotopes
      integer idic13
      integer iphytc13,izoopc13,idetrc13
      integer icoccc13,icoccn15,icaco3c13,idiazc13
      integer idiatc13

      common /mw_i/ idic13
      common /mw_i/ iphytc13,izoopc13,idetrc13
      common /mw_i/ icoccc13,icoccn15,icaco3c13,idiazc13
      common /mw_i/ idiatc13

      integer idin15
      integer iphytn15,izoopn15,idetrn15
      integer idiazn15
      integer idiatn15

      common /mw_i/ idin15
      common /mw_i/ iphytn15,izoopn15,idetrn15
      common /mw_i/ idiazn15
      common /mw_i/ idiatn15

      integer isil30, idiats30
      common /mw_i/ isil30, idiats30
#endif

!     indicies for ocean tracer source array
!       isidealage = index for ideal age source
!       istemp  = index for temperature
!       issalt  = index for salinity
!       isdic   = index for carbon
!       isc14   = index for carbon 14
!       isalk   = index for alkalinity
!       ispo4   = index for phosphate
!       isphyt  = index for phytoplankton
!       iszoop  = index for zooplankton
!       isdetr  = index for detritus
#if defined O_npzd_iron
!       isdfe   = index for iron
!       isdetrfe= index for particulate iron
#endif
!       itrc    = index of tracer sources for all tracers (0 = no source)
!       mapst   = map for ocean tracer source names

      character(10) :: mapst
      common /mw_c/ mapst(nt)

      integer itrc
      common /mw_i/ itrc(nt)

      integer istemp, issalt, isdic, isc14, isalk, ispo4, isphyt
      integer iszoop, isdetr, iso2, isno3
#if !defined O_npzd_simple
      integer isdiaz
#endif
#if defined O_idealage
      integer isidealage
#endif
#if defined O_kk_ballast
      integer isdetr_B
#endif
#if defined O_cal_kk
      integer iscocc
#endif
#if defined O_kk_caco3tr
      integer iscaco3
#endif
#if defined O_kk_diat
      integer isdiat
#endif
#if defined O_kk_si
      integer issil
#endif
#if defined O_plastic
      integer ismp, ismpa, ismpp
#endif
#if defined O_zoop_det
      integer isdetrz
#endif
#if defined O_npzd_iron
      integer isdfe, isdetrfe
#endif

      common /mw_i/ istemp, issalt, isdic, isc14,isalk, ispo4, isphyt
      common /mw_i/ iszoop, isdetr, iso2, isno3
#if !defined O_npzd_simple
      common /mw_i/ isdiaz
#endif
#if defined O_plastic
      common /mw_i/ ismp, ismpa, ismpp
#endif
#if defined O_zoop_det
      common /mw_i/ isdetrz
#endif
#if defined O_npzd_iron
      common /mw_i/ isdfe, isdetrfe
#endif
#if defined O_idealage
      common /mw_i/ isidealage
#endif
#if defined O_kk_ballast
      common /mw_i/ isdetr_B
#endif
#if defined O_cal_kk
      common /mw_i/ iscocc
#endif
#if defined O_kk_caco3tr
      common /mw_i/ iscaco3
#endif
#if defined O_kk_diat
      common /mw_i/ isdiat
#endif
#if defined O_kk_si
      common /mw_i/ issil
#endif
#if defined O_isotopes
      integer isdic13
      integer isphytc13,iszoopc13,isdetrc13
      integer iscoccc13,iscaco3c13,isdiazc13
      integer isdiatc13

      common /mw_i/ isdic13
      common /mw_i/ isphytc13,iszoopc13,isdetrc13
      common /mw_i/ iscoccc13,iscaco3c13,isdiazc13
      common /mw_i/ isdiatc13

      integer isdin15, issil30
      integer isphytn15,iszoopn15,isdetrn15
      integer iscoccn15,isdiazn15
      integer isdiatn15, isdiats30

      common /mw_i/ isdin15, issil30
      common /mw_i/ isphytn15,iszoopn15,isdetrn15
      common /mw_i/ iscoccn15,isdiazn15
      common /mw_i/ isdiatn15, isdiats30
#endif
!-----------------------------------------------------------------------
!     MW arrays for diagnostic equations and workspace:
!-----------------------------------------------------------------------

!     diagnostic advective velocities are in units of cm/sec

!     adv_vet  = advective velocity on the eastern face of a "T" cell
!     adv_vnt  = advective velocity on the northern face of a "T" cell
!     adv_veu  = advective velocity on the eastern face of a "u" cell
!     adv_vnu  = advective velocity on the northern face of a "u" cell
!     adv_vbt  = advective velocity on the bottom face of a "T" cell
!     adv_vbu  = advective velocity on the bottom face of a "u" cell

!     rho      = density at centre of a "T" cell in units of gm/cm**3
!                note: there is an arbitrary constant which is only a
!                function of depth in "rho". It is related to
!                subtracting a reference level density for purposes of
!                accuracy.

!     grad_p   = hydrostatic pressure gradient for "u" cell. There are
!                two components: (1,2) is for (dp/dx, dp/dy)

      real adv_vet, adv_vnt, adv_veu, adv_vnu, adv_vbt, adv_vbu, rho
      real rhotaum1, rhotaup1, rhotilde, grad_p
      common /mw_r/ adv_vet(imt,km,jsmw:jmw), adv_vnt(imt,km,1:jmw)
      common /mw_r/ adv_veu(imt,km,jsmw:jemw)
      common /mw_r/ adv_vnu(imt,km,1:jemw)
      common /mw_r/ adv_vbt(imt,0:km,jsmw:jmw)
      common /mw_r/ adv_vbu(imt,0:km,jsmw:jemw)
      common /mw_r/ rho(imt,km,jsmw:jmw)
#if defined O_pressure_gradient_average
      common /mw_r/ rhotaum1(imt,km,jsmw:jmw)
      common /mw_r/ rhotaup1(imt,km,jsmw:jmw)
      common /mw_r/ rhotilde(imt,km,jsmw:jmw)
#endif
      common /mw_r/ grad_p(imt,km,jsmw:jemw,2)

!     tmask = tracer cell land/sea mask   = (0.0, 1.0) on (land, sea)
!     umask = velocity cell land/sea mask = (0.0, 1.0) on (land, sea)

      real tmask, umask
      common /mw_r/ tmask(imt,km,1:jmw), umask(imt,km,1:jmw)

#if defined O_fourth_order_tracer_advection || defined O_quicker

!     adv_f4n = 4th order advective flux

      real adv_f4n
      common /mw_r/ adv_f4n(imt,km,1:jemw,nt)
#endif

#if defined O_biharmonic

!     delsq = del**2 of prognostic variables

      real del2
      common /mw_r/ del2(imt,km,1:jmw,nvarbh)
#endif

!     these workspace arrays are recalculated for each component of the
!     equations so do not have to be moved as the MW moves northward.

!     adv_fe   = advective flux across the eastern face of a cell
!     adv_fn   = advective flux across the northern face of a cell
!                (removed in most cases and put directly into the
!                 statement functions for speed optimization.)
!     adv_fb   = advective flux across the bottom face of a cell

!     diff_fe  = diffusive flux across the eastern face of a cell
!     diff_fn  = diffusive flux across the northern face of a cell
!     diff_fb  = diffusive flux across the bottom face of a cell
!     source   = source term

      real adv_fe
      common /mw_r/ adv_fe(imt,km,jsmw:jemw)
#if defined O_fct
      real adv_fn
      common /mw_r/ adv_fn(imt,km,1:jemw)
#endif
      real adv_fb, diff_fe
      common /mw_r/ adv_fb(imt,0:km,jsmw:jemw)
      common /mw_r/ diff_fe(imt,km,jsmw:jemw)
#if !defined O_consthmix || defined O_biharmonic || defined O_isopycmix
      real diff_fn
      common /mw_r/ diff_fn(imt,km,1:jemw)
#endif
      real diff_fb
      common /mw_r/ diff_fb(imt,0:km,jsmw:jemw)
#if defined O_isopycmix
      real diff_fbiso
      common /mw_r/ diff_fbiso(imt,0:km,jsmw:jemw)
#endif
#if defined O_source_term || defined O_npzd || defined O_carbon_14
      real source
      common /mw_r/ source(imt,km,jsmw:jemw)
#endif

#if defined O_implicitvmix || defined O_isopycmix
      real zzi
      common /mw_r/ zzi(imt,km,jsmw:jemw)
#endif

!     these grid factors are for optimizations and are recalculated as
!     the MW moves northward so they do not have to be moved.

      real cstdxtr, cstdxur, cstdxt2r, ah_cstdxur, csudxur, csudxtr
      real csudxu2r, am_csudxtr

      common /mw_r/ cstdxtr(imt,jsmw:jmw), cstdxur(imt,jsmw:jmw)
      common /mw_r/ cstdxt2r(imt,jsmw:jmw), ah_cstdxur(imt,jsmw:jmw)
      common /mw_r/ csudxur(imt,jsmw:jmw)
      common /mw_r/ csudxtr(imt,jsmw:jmw)
      common /mw_r/ csudxu2r(imt,jsmw:jmw)
#if defined O_anisotropic_viscosity
      common /mw_r/ am_csudxtr(imt,km,jsmw:jmw)
#else
      common /mw_r/ am_csudxtr(imt,jsmw:jmw)
#endif

!     these variables are either constant or globally dimensioned by
!     "jmt", so they do not need to be moved as the MW moves northward

!     advmet = coeff for metric advection.
!     cori   = Coriolis parameter for velocity component "n"

      real advmet, cori, cf
      common /advec_r/ advmet(jmt,2)
      common /coriol_r/ cori(imt,jmt,2)
      common /extwrk_r/ cf(imt,jmt,-1:1,-1:1)

#if defined O_linearized_advection

!     initial temperature profile independent of latitude and longitude

      real tbarz
      common /profile_r/ tbarz(km,nt)
#endif

!     smf = surface momentum flux
!          1 => zonal wind stress (dynes/cm**2)
!          2 => meridional wind stress (dynes/cm**2)
!     bmf = bottom momentum flux
!          1 => zonal bottom drag (dynes/cm**2)
!          2 => meridional bottom drag (dynes/cm**2)
!     stf = surface tracer flux
!          1 => surface heat flux (cal/cm**2/sec = cm*degC/sec = ly/sec)
!          2 => surface salt flux (grams of salt/cm**2/sec)
!     btf = bottom tracer flux (for consistency but normally zero!)
!          1 => bottom heat flux (cal/cm**2/sec = cm*degC/sec = ly/sec)
!          2 => bottom salt flux (grams of salt/cm**2/sec)

      real smf,  bmf, stf, btf
      common /mw_r/ smf(imt,1:jmw,2),  bmf(imt,1:jmw,2)
      common /mw_r/ stf(imt,1:jmw,nt), btf(imt,1:jmw,nt)

#if defined O_fct
!               "antidiffusive" flux as in Zalesak, 1989( see FCTstm.h)
!               same for R+, R-
!     anti_fe  = antidiffusive flux across the eastern face of a T cell
!     anti_fn  = antidiffusive flux across the northern face of a T cell
!     anti_fb  = antidiffusive flux across the bottom face of a T cell
!     R_plusY  = ratio of maximal feasible to maximal possible change
!                of tracer T in subroutine tracer.F, N-S dimension delimiter
!     R_minusY = ratio of minimal feasible to minimal possible change
!               of tracer T in subroutine tracer.F, N-S dimension delimiter

      real anti_fe, anti_fn, anti_fb, R_plusY, R_minusY
      common /mw_r/ anti_fe(imt,km,jsmw:jmw,nt)
      common /mw_r/ anti_fn(imt,km,1:jmw-1+jmw/jmt,nt)
      common /mw_r/ anti_fb(imt,0:km,jsmw:jmw,nt)
      common /mw_r/ R_plusY(imt,km,1:jmw-1+jmw/jmt,nt)
      common /mw_r/ R_minusY(imt,km,1:jmw-1+jmw/jmt,nt)
# if defined O_fct_dlm2 && !defined O_fct_dlm1

!     t_lo     = low order solution for tracer

      real t_lo
      common /mw_r/ t_lo(imt,km,jsmw:jmw,nt)
# endif
# if defined O_fct_3d

!     R_plus3 = ratio of maximal feasible to maximal possible change
!               of tracer T in subroutine tracer.F, 3-d delimiter
!     R_minus3= ratio of minimal feasible to minimal possible change
!               of tracer T in subroutine tracer.F, 3-d delimiter

      real R_plus3, R_minus3
      common /mw_r/ R_plus3(imt,km,1:jmw-1+jmw/jmt,nt)
      common /mw_r/ R_minus3(imt,km,1:jmw-1+jmw/jmt,nt)
# endif
#endif
