#if defined O_mom

      subroutine tracer (joff, js, je, is, ie)

!#if defined O_mom
!=======================================================================
!     compute tracers at "tau+1" for rows js through je in the MW.

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!=======================================================================

      implicit none

      character(120) :: fname, new_file_name

      integer istrt, iend, i, k, j, ip, kr, jq, n, jp, jrow, iou, js
      integer je, limit, joff, is, ie, kmx, m, kb, idiag, index
      integer it(10), iu(10), ib(10), ic(10), nfnpzd, mfnpzd, mxfnpzd
      integer id_time, id_xt, id_yt, id_zt, fe_jlo, fe_m, fe_k, fe_n

      parameter (fe_n = 14)

      logical inqvardef, exists

      real rctheta, declin, gl, impo, expo, npp, time
      real remi, excr, graz, morp, morpt, morz, temp, swr, dayfrac
      real graz_Det, graz_Z, avej, avej_D, gmax, no3P, po4P, po4_D
      real phin_loc, dz, prca, dprca, nud, bct, tap, fo2, so2, ai
      real hi, hs
      real npp_D, graz_D, morp_D, no3flag, deni, nfix, felimit
      real t_i, t_j, dz_t2r, dz_tr, dz_wtr, dx_t2r, dx_tr, dy_t2r
      real dy_tr, adv_tx, adv_ty, adv_tz, adv_txiso, adv_tyiso
      real adv_tziso, diff_tx, diff_ty, diff_tz, zmax, cont, drho
      real drhom1, wt, ahbi_cstr, ahbi_csu_dyur, gamma, rrstd, fy, fyz
      real fe_dy, fe_conc, fe_x(fe_n), fe_y(fe_n), bctz, felimit_D
      real Paulmier_a, Paulmier_z, Paulmier_R0

      real expofe, impofe, feorgads, remife, thetamax, deffe, fecol
      real thetachl, chl, chl_D, feprime, fesed, bfe, sgb
# if defined O_kk_caco3tr || defined O_kk_si || defined O_kk_diat
      real sil_in, feorgads_ca
# endif
# if defined O_kk_ballast
      real remi_B,expo_B, impo_B, graz_Det_B
# endif
# if defined O_kk_diat
      real npp_Diat, graz_Diat, morp_Diat, morpt_Diat,felimit_Diat
#  if defined O_npzd_iron
      real deffe_Diat, thetamax_Diat
#  endif
# endif
# if defined O_kk_si
      real dprop, prop, expo_op, op_diss, bsi
# endif
# if defined O_cal_kk
      real npp_C, graz_C, morp_C, morpt_C,felimit_C
#  if defined O_npzd_iron
      real deffe_C, thetamax_C
#  endif
# endif
# if defined O_benthic_denitrification
      real bdeni, din15flag, lno3, sg_bdeni
# endif
# if defined O_kk_caco3tr
      real temp_in, d_in,alk_in, sal_in, expocaco3,old_diss
      real dissl,calpro,impocaco3,calatt
      real expo_out, impo_out, caco3in_loc, dissk1
      real atmpres1, pHlo1, pHhi1, pH1, p_in
      real co2star1, dco2star1, pCO21, dpco21, CO31
      real omegaca, omegaar, c_in,depth1
#include "cembm.h"
# endif
#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "accel.h"
#include "coord.h"
#include "cregin.h"
#include "csbc.h"
#include "emode.h"
#include "grdvar.h"
#include "hmixc.h"
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "timeavgs.h"
#include "tmngr.h"
#include "vmixc.h"
# if defined O_save_convection || defined O_carbon_14
#include "diaga.h"
# endif
# if defined O_ice
#  if defined O_ice_cpts
#include "cpts.h"
#  endif
#include "ice.h"
# endif
# if defined O_npzd
#  if defined O_embm
#include "atm.h"
#  endif
#include "npzd.h"
# endif
# if defined O_carbon_fnpzd
#include "calendar.h"
# endif
# if defined O_matrix
#include "matrix.h"
      integer ntile
      logical isTile
# endif

# ifndef O_TMM
      parameter (istrt=2, iend=imt-1)
# else
      parameter (istrt=1, iend=imt)
# endif
      real twodt(km)
# if defined O_plume
      real  work(imt,km,jsmw:jemw)
# endif
# if defined O_npzd
      real snpzd(ntnpzd), tnpzd(ntnpzd), gbio
# endif
# if defined O_kk_si
      real settle, constant, temp_2(km)
      real bctd(imt,km,jsmw:jemw)
# endif
# if defined O_isotopes
      real c13impo, c13expo, ac13_DIC_aq, fcnpp, bc13npp
      real ac13_aq_POC, ac13b, prca13, rdic13, rtdic13(km)
      real caco3c13impo, caco3c13expo, rtcaco3c13(km)
      real n15impo, n15expo, rno3, eps_bdeni, eps_bdeni0
      real bbdeni, uno3, bwcdeni, eps_wcdeni
      real eps_assim, eps_excr, eps_nfix, eps_recy
      real rtsil30(km), sil30_in, rtdiats30(km)
      real eps_hydr, eps_dust, eps_diss, eps_riv
      real b2si, fcsi
# endif
#ifndef O_TMM
# if defined O_npzd || defined O_carbon_14
      real src(imt,km,jsmw:jemw,nsrc)
# endif
#else
#include "uvok_comm.h"
#endif

# if defined O_carbon_fnpzd
      real tmpijk(imtm2,jmtm2,km), tmpi(imtm2), tmpj(jmtm2)
      real, allocatable :: fnpzd(:,:,:,:)
      save fnpzd, nfnpzd, mfnpzd, mxfnpzd
# endif
# if defined O_isopycmix
#include "isopyc.h"
# endif
#include "fdift.h"
# if defined O_plastic
      real zoop1_r, phyt_s, zoop1_s, agg_upt, agg_rel, zoop_rel
      real phyt_r, zoop_r, agg_con, expomp, impomp, zoop_upt
      real poo_con, expomp_a, impomp_a, expomp_p, impomp_p
# endif
# if defined O_zoop_det
      real remiz, expoz, impoz
# endif
# if defined O_carbon_fnpzd
!-----------------------------------------------------------------------
!     create file for writing npzd 3d fluxes of carbon and alkalinity
!-----------------------------------------------------------------------
      if (.not. allocated (fnpzd)) then
        allocate ( fnpzd(imt,jmt,km,2) )
        fnpzd(:,:,:,:) = 0.
!       initialize time record counter
        mfnpzd = 0
!       initialize records written counter
        nfnpzd = 0
!       set the repeating time interval to one year
        mxfnpzd = nint(daylen*yrlen/dtts)
        fname = new_file_name ("O_fnpzd.nc")
        inquire (file=trim(fname), exist=exists)
        if (exists) then
!         limit time interval to what exists
          call openfile (fname, iou)
          call getdimlen ('time', iou, nfnpzd)
          mxfnpzd = min(mxfnpzd, nfnpzd)
#  if defined O_carbon
          if (inqvardef('O_fnpzd1', iou)) nfnpzd = mxfnpzd
#  endif
#  if defined O_npzd_alk
          if (inqvardef('O_fnpzd2', iou)) nfnpzd = mxfnpzd
#  endif
        else
!         define file if not there
          call openfile (fname, iou)
          call redef (iou)
          call defdim ('time', iou, 0, id_time)
          call defdim ('longitude', iou, imtm2, id_xt)
          call defdim ('latitude', iou, jmtm2, id_yt)
          call defdim ('depth', iou, km, id_zt)
          it(1) = id_time
          call defvar ('time', iou, 1, it, c0, c0, 'T', 'D'
     &,     'time', 'time', 'years since 0-1-1')
          call putatttext (iou, 'time', 'calendar', calendar)
          it(1) = id_xt
          call defvar ('longitude', iou, 1, it, c0, c0, 'X', 'D'
     &,     'longitude', 'longitude', 'degrees')
          it(1) = id_yt
          call defvar ('latitude', iou, 1, it, c0, c0, 'Y', 'D'
     &,     'latitude', 'latitude', 'degrees')
          it(1) = id_zt
          call defvar ('depth', iou, 1, it, c0, c0, 'Z', 'D'
     &,     'depth', 'depth', 'm')
          it(1) = id_xt
          it(2) = id_yt
          it(3) = id_zt
          it(4) = id_time
#  if defined O_carbon
          call defvar ('O_fnpzd1', iou, 4, it, c0, c0, ' ', 'D'
     &,       ' ', ' ', ' ')
#  endif
#  if defined O_npzd_alk
          call defvar ('O_fnpzd2', iou, 4, it, c0, c0, ' ', 'D'
     &,       ' ', ' ', ' ')
#  endif
          call enddef (iou)
          ib(:) = 1
          ic(:) = imtm2
          tmpi(1:imtm2) = xt(2:imtm1)
          call putvara ('longitude', iou, imtm2, ib, ic, xt, c1, c0)
          ib(:) = 1
          ic(:) = jmtm2
          tmpj(1:jmtm2) = yt(2:jmtm1)
          call putvara ('latitude', iou, jmtm2, ib, ic, yt, c1, c0)
          ib(:) = 1
          ic(:) = km
          call putvara ('depth', iou, km, ib, ic, zt, c1, c0)
        endif
      endif
# endif

# ifndef O_TMM
!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------

      if (js .gt. je) return

!-----------------------------------------------------------------------
!     limit the longitude indices based on those from the argument list
!     Note: this is currently bypassed. istrt and iend are set as
!           parameters to optimize performance
!-----------------------------------------------------------------------

!      istrt = max(2,is)
!      iend  = min(imt-1,ie)

!-----------------------------------------------------------------------
!     build coefficients to minimize advection and diffusion computation
!-----------------------------------------------------------------------

#  if defined O_fourth_order_tracer_advection || defined O_fct || defined O_quicker || defined O_pressure_gradient_average || defined O_biharmonic || defined O_isopycmix
      limit = min(je+1+joff,jmt) - joff
      do j=js,limit
#  else
      do j=js,je
#  endif
        jrow = j + joff
        do i=istrt-1,iend
          cstdxtr(i,j)    = cstr(jrow)*dxtr(i)
          cstdxt2r(i,j)   = cstr(jrow)*dxtr(i)*p5
          cstdxur(i,j)    = cstr(jrow)*dxur(i)
#  if defined O_consthmix && !defined O_bryan_lewis_horizontal && !defined O_biharmonic
          ah_cstdxur(i,j) = diff_cet*cstr(jrow)*dxur(i)
#  endif
        enddo
      enddo
#  if defined O_plume

!-----------------------------------------------------------------------
!     find depth of penetration for the llnl simple plume model
!-----------------------------------------------------------------------

!     set parameters
      zmax = 160.0e2 ! use constant depth penetration (if cont=0)
!      cont = 0.0    ! use rho contrast for penetration (if cont>0)
      cont = 0.4e-3  ! 0.4 kg/m3 = 0.4e-3g/cm3

      do j=js,je
        jrow = j + joff
        do i=is,ie
          subz(i,jrow) = zmax
        enddo
        if (cont .gt. c0) then
          call state_ref (t(1,1,1,1,tau), t(1,1,1,2,tau)
     &,                   work(1,1,jsmw), j, j, 1, imt, 1)
          do i=is,ie
            kmx = kmt(i,jrow)
            subz(i,jrow) = zw(1)
            do k=2,kmx
              drho = work(i,k,j) - work(i,1,j)
              drhom1 = work(i,k-1,j) - work(i,1,j)
              if (drho .ge. cont .and. drhom1 .lt. cont)
     &          subz(i,jrow) = zt(k-1) + (zt(k) - zt(k-1))
     &                      *(cont - drhom1)/(drho - drhom1)
            enddo
            if (drho .le. cont) subz(i,jrow) = zw(kmx)
          enddo
        endif

!-----------------------------------------------------------------------
!     calculate "3d stf" multiplier for the llnl simple plume model
!-----------------------------------------------------------------------

        do i=is,ie
          if (kmt(i,jrow) .gt. 0) then
            subz(i,jrow) = max(min(subz(i,jrow),zw(kmt(i,jrow))),zw(1))
            work(i,1,j) = c1/subz(i,jrow)
            do k=2,km
              wt =  min(max(c0,(subz(i,jrow) - zw(k-1))),dzt(k))*dztr(k)
              work(i,k,j) = wt*work(i,1,j)
            enddo
          else
            do k=1,km
              work(i,k,j) = c0
            enddo
          endif
        enddo
      enddo
#  endif

# endif ! not O_TMM

# if defined O_npzd
!-----------------------------------------------------------------------
!     calculation of biological interactions
!-----------------------------------------------------------------------

      declin = sin((mod(relyr,1.) - 0.22)*2.*pi)*0.4   ! declination

      do k=1,km
        twodt(k) = c2dtts*dtxcel(k)
        nbio(k) = twodt(k)/dtnpzd
        dtbio(k) = twodt(k)/nbio(k)
        rdtts(k) = 1./twodt(k)
        rnbio(k) = 1./nbio(k)
      enddo
 

      tap = 2.*par   ! alpha is added in npzd_src

#  if defined O_kk_si
      rivsil = 0.0        !river Si deposition mol si/s??
      silwflx = 0.0       !Opal flux through bottom umol si/s
      bctd(:,:,:) = 0.0     !temperature scaling of microbial dissolution
#  endif  
      do j=js,je
        jrow = j + joff

        do i=is,ie

          if (kmt(i,jrow) .gt. 0) then

#  if defined O_ice
#   if defined O_ice_cpts
            ai = 0.
            hi = 0.
            hs = 0.
            do n=1,ncat
              ai =  ai + A(i,jrow,2,n)
              hi =  hi + heff(i,jrow,2,n)
              hs =  hs + hseff(i,jrow,2,n)
            enddo
#   else
            ai = aice(i,jrow,2)
            hi = hice(i,jrow,2)
            hs = hsno(i,jrow,2)
#   endif
#  else
            ai = 0.
            hi = 0.
            hs = 0.
#  endif
!           calculate day fraction and incoming solar
!           angle of incidence = lat - declin, refraction index = 1.33
            rctheta = max(-1.5, min(1.5, tlat(i,jrow)/radian - declin))
#  if defined O_npzd_cdom_attenuation
            rctheta = (1.2*kw)/sqrt(1. - (1. - cos(rctheta)**2.)/1.33
     &                **2.)
#  else
            rctheta = kw/sqrt(1. - (1. - cos(rctheta)**2.)/1.33**2.)
#  endif
            dayfrac = min( 1., -tan(tlat(i,jrow)/radian)*tan(declin))
            dayfrac = max(1e-12, acos(max(-1., dayfrac))/pi)
#  if defined O_embm
            swr = dnswr(i,jrow)*1e-3*(1. + ai*(exp(-ki*(hi + hs)) - 1.))
#  else
            swr = 200.
#  endif
            expo = 0.0
            impo = 0.0
# if defined O_kk_ballast
            expo_B = 0.0 ! detritus from coccs and zoop
            impo_B = 0.0 ! detritus from coccs and zoop
            remi_B = 0.0 ! detritus shifting ballasted and free pools
# endif
#  if defined O_kk_si
            expo_op = 0.0 ! export of opal umol Si/cm^2/s
            dprop = 0.0   ! umol Si/cm^3/s
            prop = 0.0          ! integrated production of opal mol Si/m^2/s
            op_diss = 0.0 ! opal dissolution rate umol/cm^2/s ??       
#  endif
# if defined O_kk_caco3tr
            caco3in_loc = 0.0
            old_diss = 0.0
            calpro = 0.0
            calatt = 0.0
            dissl = 0.0
            impocaco3 = 0.0
            expocaco3 = 0.0
            expo_out = 0.0
            impo_out = 0.0
            dissk1 = 0.0
# endif
#  if defined O_npzd_iron
            expofe = 0.0
            impofe = 0.0
#  endif
            phin_loc = 0.0 ! integrated phytoplankton
            prca = 0.0 ! integrated production of calcite
# if defined O_phyt_albedo
            phin(:,:) = 0.
            caco3in(:,:) = 0.
# endif
# if defined O_isotopes
            c13impo = 0.0
            c13expo = 0.0
            n15impo = 0.0
            n15expo = 0.0
            prca13 = 0.0
            caco3c13impo = 0.0
            caco3c13expo = 0.0

            eps_assim = 6.
            eps_excr = 4.
            eps_nfix = 1.
            eps_wcdeni = 20.
            eps_bdeni0 = 6.
            eps_recy = 1.3

            eps_hydr = -0.3 !Frings et al 2016
            eps_dust = -0.67 !Frings et al 2016
            eps_diss = -0.55 !Demarest et al 2018
            eps_riv = 1.28 !Frings et al 2016
# endif
# if defined O_zoop_det
            impoz = 0.0
            expoz = 0.0
# endif
#  if defined O_benthic_denitrification
            sg_bdeni = 0.
            bdeni = 0.
#  endif
# if defined O_plastic
            impomp_a = 0.0
            expomp_a = 0.0
            impomp_p = 0.0
            expomp_p = 0.0
            phyt_r = 0.0
            zoop_r = 0.0
            zoop1_r = 0.0
            phyt_s = 0.0
            zoop1_s = 0.0
            agg_upt = 0.0
            zoop_upt = 0.0
            agg_rel = 0.0
            zoop_rel = 0.0
!convert a portion of detritus to marine snow mmolN=>aggs
! 12gC/molC x 6.6molC/mmolN x 1e6 ugC/gC x 1 agg/8.8ugC
            agg_con = 12.011*redctn*1e6/8.8*pagg 
! particles/cm^3 <=> mmol N/m^3 (equivalent to nmol N/cm^3)
!236000*1e6/51.2e12/588g food *16 mol*1e3=> mmol N
            poo_con = 73.75/588.0*ping!523416.0 !convert plastic to poo
! this value is produced from 355.4e8 g MP=485e10 particles from Eriksen et al 2014. or Van Sebille 236000 tons MP = 51.2e12 particles
! It is assumed that 1 g MP is going to replace 1 G of FOOD in the zooplankton's diet, and so
! MP particles/cm^3 are converted to an equivalent mmol N/m3 (redfield ratio)

#   if defined O_plastic_conc
      call getplastic !call the plastic into the model
#   endif
# endif
            kmx = min(kmt(i,jrow), kpzd)
            do k=1,kmx

!-----------------------------------------------------------------------
!     Initialize tracers.       
!-----------------------------------------------------------------------
              src(i,k,j,:)=0.0
              tnpzd(ibion) = t(i,k,j,ipo4,taum1)
              tnpzd(ibiop) = t(i,k,j,iphyt,taum1)
              tnpzd(ibioz) = t(i,k,j,izoop,taum1)
              tnpzd(ibiod) = t(i,k,j,idetr,taum1)
#  if defined O_plastic
              tnpzd(ibiomp) = t(i,k,j,imp,taum1)
              tnpzd(ibiompa) = t(i,k,j,impa,taum1)
              tnpzd(ibiompp) = t(i,k,j,impp,taum1)
#  endif
#  if defined O_zoop_det
              tnpzd(ibiodz) = t(i,k,j,idetrz,taum1)
#  endif
#  if defined O_npzd_nitrogen
              tnpzd(ibiono3) = t(i,k,j,ino3,taum1)
#   if !defined O_npzd_simple
              tnpzd(ibiodiaz) = t(i,k,j,idiaz,taum1)
#   endif
#  endif
#  if defined O_npzd_iron
              tnpzd(ibiodfe) = t(i,k,j,idfe,taum1)
              tnpzd(ibiodetrfe) = t(i,k,j,idetrfe,taum1)
#  endif
#  if defined O_cal_kk
              tnpzd(ibioc) = t(i,k,j,icocc,taum1)
#  endif
#  if defined O_kk_diat
              tnpzd(ibiodiat) = t(i,k,j,idiat,taum1)
#  endif
#  if defined O_kk_caco3tr
              tnpzd(ibiocaco3) = t(i,k,j,icaco3,taum1)
! mmol/m3
#  endif
#  if defined O_kk_ballast
              tnpzd(ibiod_B) = t(i,k,j,idetr_B,taum1)
#  endif
#  if defined O_isotopes
              tnpzd(ibiocaco3c13) = max(trcmin,t(i,k,j,icaco3c13,taum1))
              tnpzd(ibiodic13) = max(trcmin,t(i,k,j,idic13,taum1))
              tnpzd(ibiophytc13) = max(trcmin,t(i,k,j,iphytc13,taum1))
              tnpzd(ibiodiazc13) = max(trcmin,t(i,k,j,idiazc13,taum1))
              tnpzd(ibiozoopc13) = max(trcmin,t(i,k,j,izoopc13,taum1))
              tnpzd(ibiococcc13) = max(trcmin,t(i,k,j,icoccc13,taum1)) 
              tnpzd(ibiodiatc13) = max(trcmin,t(i,k,j,idiatc13,taum1))
              tnpzd(ibiodetrc13) = max(trcmin,t(i,k,j,idetrc13,taum1))

              tnpzd(ibiodin15) = max(trcmin,t(i,k,j,idin15,taum1))
              tnpzd(ibiophytn15) = max(trcmin,t(i,k,j,iphytn15,taum1))
              tnpzd(ibiodiazn15) = max(trcmin,t(i,k,j,idiazn15,taum1))
              tnpzd(ibiozoopn15) = max(trcmin,t(i,k,j,izoopn15,taum1))
              tnpzd(ibiococcn15) = max(trcmin,t(i,k,j,icoccn15,taum1))
              tnpzd(ibiodiatn15) = max(trcmin,t(i,k,j,idiatn15,taum1))
              tnpzd(ibiodetrn15) = max(trcmin,t(i,k,j,idetrn15,taum1))

!              tnpzd(ibiosil30) = t(i,k,j,isil30,taum1)
              tnpzd(ibiodiats30) = max(trcmin,t(i,k,j,idiats30,taum1))
#  endif
#  if defined O_kk_caco3tr
              swr = swr*exp(-(kc*phin_loc+kc_c*caco3in_loc))
#  else
              swr = swr*exp(-kc*phin_loc)
#  endif
! phin_loc has been removed from the following code to fix a bug
              phin_loc = (max(tnpzd(ibiop),trcmin)
#  if defined O_npzd_nitrogen && !defined O_npzd_simple
     &             + max(tnpzd(ibiodiaz),trcmin) 
#  endif
#  if defined O_cal_kk
     &             + max(tnpzd(ibioc),trcmin) 
#  endif
#  if defined O_kk_diat
     &             + max(tnpzd(ibiodiat),trcmin) 
#  endif
     &              )*dzt(k)
#  if defined O_kk_caco3tr
              caco3in_loc = max(tnpzd(ibiocaco3),trcmin)*dzt(k)
#  endif
#  if defined O_phyt_albedo
              phin(i,j) = phin(i,j) + phin_loc 
              caco3in(i,j) = caco3in(i,j) + caco3in_loc
#  endif
# if defined O_kk_caco3tr
              impocaco3 = expocaco3*dztr(k)
# endif
# if defined O_kk_ballast
              impo_B = expo_B*dztr(k) !ballasted detritus
# endif
              gl = tap*swr*exp(ztt(k)*rctheta)
              impo = expo*dztr(k)
#  if defined O_plastic
              impomp_a = expomp_a*dztr(k)
              impomp_p = expomp_p*dztr(k)
#  endif
#  if defined O_zoop_det
              impoz = expoz*dztr(k)
#  endif
#  if defined O_npzd_iron
              impofe = expofe*dztr(k)
#  endif
              bct = bbio**(cbio*t(i,k,j,itemp,taum1))
#  if defined O_kk_si
              bctd(i,k,j) = bct !store it for using below
#  endif
#  if defined O_zoop_graz_upper_temp_limit
              if (t(i,k,j,itemp,taum1).gt.20) then
                 bctz = (0.5*(tanh(t(i,k,j,io2,taum1)*1000. - o2t)+1))
     &                  *bbio**(cbio*20)
              else
                 bctz = (0.5*(tanh(t(i,k,j,io2,taum1)*1000. - o2t)+1))
     &                  *bct
              endif
#  else
                 bctz = (0.5*(tanh(t(i,k,j,io2,taum1)*1000. - o2t)+1))
     &                  *bct
#  endif
#  if defined O_npzd_fe_limitation
              if (k.le.kmfe) then
#   ifndef O_TMM
! create x (time) and y (data) arrays for fe interpolation
! first zero them out just in case
                 fe_y(:) = 0
                 fe_x(:) = 0
! write in values for days 0 and 365 for interp. bounds
                 fe_y(1) = fe_dissolved(i,j,k,1)
                 fe_x(1) = 0
                 fe_y(14) = fe_dissolved(i,j,k,12)
                 fe_x(14) = 365
! now create the rest of the array based on the BLING
! data which is from day 16 of each month
                 do m=2,13
                    fe_y(m) = fe_dissolved(i,j,k,m-1)
                    if (m.eq.2) then
                       fe_x(2) = 16
                    else if (m.eq.3) then
                       fe_x(3) = 47
                    else if (m.eq.4) then
                       fe_x(4) = 75
                    else if (m.eq.5) then
                       fe_x(5) = 106
                    else if (m.eq.6) then
                       fe_x(6) = 136
                    else if (m.eq.7)  then
                       fe_x(7) = 167
                    else if (m.eq.8) then
                       fe_x(8) = 197
                    else if (m.eq.9) then
                       fe_x(9) = 228
                    else if (m.eq.10) then
                       fe_x(10) = 259
                    else if (m.eq.11) then
                       fe_x(11) = 289
                    else if (m.eq.12) then
                       fe_x(12) = 320
                    else if (m.eq.13) then
                       fe_x(13) = 350
                    endif
                 enddo
                    fe_jlo = 2
                    fe_m = 4
! find fe data day index
                    call hunt (fe_x,fe_n,dayoyr,fe_jlo)
! initialize the fe data array at the right day 
                    fe_k = min(max(fe_jlo-(fe_m-1)/2,1),fe_n+1-fe_m)
! interpolate the fe data, note this does not use the whole array
                    call polint (fe_x(fe_k),fe_y(fe_k)
     &,                          fe_m,dayoyr,fe_conc,fe_dy)
! calculate the fe limitation term                
#   else
               fe_conc = fe_dissolved(i,j,k,1) !this is time-interpolated value
#   endif !not O_TMM
                    felimit = fe_conc/(kfe + fe_conc)
                    felimit_D = fe_conc/(kfe_D + fe_conc)
#   if defined O_cal_kk
                    felimit_C = fe_conc/(kfe_C + fe_conc)
#   endif
#   if defined O_kk_diat
                    felimit_Diat = fe_conc/(kfe_Diat + fe_conc)
#   endif
                 else
                    felimit = 1
                    felimit_D = 1
#   if defined O_cal_kk
                    felimit_C = 1
#   endif
#   if defined O_kk_diat
                    felimit_Diat = 1
#   endif
                 endif
#  endif
#  if defined O_npzd_o2
! Decrease remineralisation rate in oxygen minimum zone and when there is no
! nitrate.
              if (t(i,k,j,ino3,taum1).le.0.) then
                 nud = 0
              else
                 nud = nud0*(ro2n+(1.-ro2n)*tanh(t(i,k,j,io2,taum1)
     &                 *1000.-o2t))
              endif
#  else
              nud = nud0
#  endif

#  if defined O_kk_caco3tr
#     if defined O_co2ccn_data || defined O_co2ccn_data_transient || defined O_co2emit_track_co2
!-----------------------------------------------------------------------
!     set co2 concentration or emissions by tracking average co2
!-----------------------------------------------------------------------
       call co2ccndata
#     endif
!-----------------------------------------------------------------------
!             calculate delCO3
!-----------------------------------------------------------------------
             sal_in =  1000.0*t(i,k,j,isalt,taum1) + 35.0
             sal_in = min(45.,max(sal_in,0.))
             temp_in = t(i,k,j,itemp,taum1)
             temp_in = min(35.,max(temp_in,-2.))
             d_in = max(t(i,k,j,idic,taum1),trcmin)
             alk_in = max(t(i,k,j,ialk,taum1),trcmin)
             pHlo1 = 6.
             pHhi1 = 10.
#    if defined O_kk_si
             sil_in = max(t(i,k,j,isil,taum1),trcmin) !silicate mol/m^3
#    else
             sil_in = 7.6875e-03 !mol/m^3
#    endif
             p_in = 0.5125e-3   !phos mol/m^3
             atmpres1 = 1.0       !atm
             depth1 = dzt(k)/100    !depth in meters
#     if defined O_carbon_co2_2d
            c_in = at(imt,jmw,2,ico2)
#     else
            c_in = co2ccn
#     endif
            call co2calc_SWS (temp_in, sal_in,d_in,alk_in,c_in, p_in 
     &,                       sil_in, atmpres1, depth1, pHlo1, pHhi1,pH1
     &,                       co2star1, dco2star1, pCO21, dpco21, CO31
     &,                       omegaca, omegaar)

            romca(i,k,j) = omegaca
! assume KK that Ca ion is constant (S&G 2006 pg 366) 
            rco3(i,k,j) = CO31 !mol/m3
            rco3_sat(i,k,j) = CO31/omegaca
            rdel_sat(i,k,j) = rco3(i,k,j)- rco3_sat(i,k,j) !mol/m3
#  if defined O_isotopes
!        c13 biological fractionation 
         ac13_DIC_aq = -1.0512994e-4*temp_in +1.011765
!        Popp et al. (1989) Am. J. Sci.
         ac13_aq_POC = -0.017*log10(min(max(co2star1*1000.,2.),74.))
     &                 +1.0034
         ac13b = ac13_aq_POC/ac13_DIC_aq
         rdic13 = max(t(i,k,j,idic13,taum1),trcmin*rc13std/(1+rc13std))/
     &          max(max(t(i,k,j,idic,taum1),trcmin)
     &          - max(t(i,k,j,idic13,taum1),trcmin*rc13std/(1+rc13std))
     &          ,trcmin*rc13std/(1+rc13std))
         rdic13 = min(rdic13, 2.*rc13std)
         rdic13 = max(rdic13, 0.5*rc13std)
         bc13npp = ac13b*rdic13
         fcnpp = bc13npp/(1+bc13npp)

         c13impo = c13expo*dztr(k) !N units
         n15impo = n15expo*dztr(k) !N units
         caco3c13impo = caco3c13expo*dztr(k)
#  endif
            dissk0 = (1-(CO31-CO31/omegaca))/(kcal+(CO31-CO31/omegaca))
! dissk0 and dissk1 dimensionless
          if(dissk0.ge.0) then
            dissk1 =  min(c1,dissk0)
          else
            dissk1 = 0.0
          endif
!           if(omegaca-1.lt.0) then
!             dissk1 = dissk0*(1-omegaca) !Andreas S. 11 d^-1 from Gehlen 2006
!           else
!             dissk1 = 0.
!           endif
#  endif
#  if defined O_kk_si
              sil_in =  max(t(i,k,j,isil,taum1)*1e3,trcmin) !convert to mmol m-3 for npzd_src!
#  endif
#  if defined O_kk_variable_sipr
              k1si = 0.8 + 7.2 * (sil_in/(30.+sil_in)) !Aumont et al (2003) in mmol m-3
#  endif
#  if defined O_kk_diat && !defined O_kk_si
              sil_in = si_msk(i,j,k)*1e3  !convert to mmol m-3 for npzd_src!
#  endif
#  if defined O_isotopes
! KK different placement from MOBI
         rtdic13(k) = max(t(i,k,j,idic13,taum1)
     &             ,trcmin*rc13std/(1+rc13std))
     &             / max(t(i,k,j,idic,taum1),trcmin)
         rtdic13(k) = min(rtdic13(k), 2.*rc13std/(1+rc13std))
         rtdic13(k) = max(rtdic13(k), 0.5*rc13std/(1+rc13std))
         rtcaco3c13(k) = max(t(i,k,j,icaco3c13,taum1)
     &                     ,trcmin*rc13std/(1+rc13std))
     &             / max(t(i,k,j,icaco3,taum1),trcmin)
         rtcaco3c13(k) = min(rtcaco3c13(k), 2.*rc13std/(1+rc13std))
         rtcaco3c13(k) = max(rtcaco3c13(k), 0.5*rc13std/(1+rc13std))

         rtsil30(k) = max(t(i,k,j,isil30,taum1)
     &              ,trcmin*rs30std/(1+rs30std))
     &              /max(t(i,k,j,isil,taum1),trcmin)
         rtsil30(k) = min(rtsil30(k), 2.*rs30std/(1+rs30std))
         rtsil30(k) = max(rtsil30(k), 0.5*rs30std/(1+rs30std))
         sil30_in = max(t(i,k,j,isil30,taum1)*1e3,trcmin)

#  endif
#  if defined O_kk_si
! repeated down below unfortunately, needs tidying up
              sipr = 0.
#   if defined O_kk_variable_sipr && (defined O_npzd_iron || defined O_fe_limitation)
!     From HAMOCC5 (Aumont et al, GBC, 2003) values 0-0.5, mostly ~0.1
              sipr = sipr0 * min(1.,sil_in/k1si)
#    if defined O_fe_limitation || defined O_npzd_iron
     &             * (4.-3.*min(1.,
#     if !defined O_npzd_iron
     & fe_conc
#     else
     & max(t(i,k,j,idfe,taum1),trcmin)
#     endif
#    if defined O_kk_diat
     &/kfe_Diat))
#    elif defined O_npzd_simple
     &/kfe))
#    endif
#    endif
#   else
              sipr = sipr0
#   endif
#  endif
#  if defined O_plastic

      if(k.eq.1) then
! calculate a plastic inventory from atmospheric CO2 concentration
! plastic emissions forecast to increase linearly with increasing CO2 emission
! Borelle et al 2017 !*(co2ccn/338.3)
! plastic concentrations in world ocean forecast to triple by 2100 Koelmans et al 2015
! huge uncertainties here- Isobe et al 2019 forecast 2x by 2030 for NPac
! tied to CO2 concentrations may be conservative (< 3x in 2100  wrt 2010)
#  if defined O_plastic_conc
! formulation reflects plastic reported in particles/m^2 can really be considered an average in m^3 in the upper ~50 or so meters. Kooi et al. (2018)
        t(i,1,j,imp,taum1) = max(trcmin,plastic(i,j,2)*1.e-2) !set surface mask
#  endif

! now calculate the number of organic aggregates that fill up with plastic
! check what proportion the number of affected particles is in relation to the total detritus
! convert plastic affected aggregates  to detritus units
! conversion of 30e4 aggs per 0.285 g POC taken from Shanks and Trent 1980
! is probably conservative. 
! find the ratio of total detritus holding plastic in that grid cell
        phyt_r = min(1.
     &         ,max(trcmin,t(i,k,j,imp,taum1)*agg_con)
     &         /max(trcmin,t(i,k,j,idetr,taum1))) !proportion of full detritus

! now do the same for poo, which are bigger
! let 1 mm^3 fecal pellet = 3.4 ug N as unpublished estimate from Hofmann et al., 1981 
! assume an average fecal pellet volume of 1.13 um^3 from Cole et al., 2016
! this converts to 1 fecal pellet containing 2.74e-10 mol N
! No averaged data available- 3 photos of poo in Cole et al 2016 show at least 20 plastic
! particles per poo pellet
        zoop1_r = min(1.,max(trcmin,t(i,k,j,imp,taum1)*poo_con)
     &          /max(trcmin,t(i,k,j,idetrz,taum1)))
         zoop_r = max(trcmin,t(i,k,j,imp,taum1)*11.58*0.094/14)
     &          /max(trcmin,t(i,k,j,izoop,taum1))

! now scale the effect on detrial sinking rates
        phyt_s = wd0 + wd0*0.72*phyt_r !Moehlenkamp et al 2018
!        phyt_s = wd0 - wd0*0.65*phyt_r !Long et al 2015
!        phyt_s = wd0 + wd0*0.28*phyt_r !Porter et al 2018
        zoop1_s = wdz0 - wdz0*0.56*zoop1_r !Cole et al 2016

      endif !end k condition

!scale parameters based on plastic concentrations
!           gbio = gbio0 - gbio0*0.4*zoop_r
!           if(gbio.lt.0) gbio = 0.
!       linear increase wd0-200m with depth
!           wd(k) = (phyt_s+4.0e-2*zt(k))/86400.0/dzt(k)  ! [s-1]
!           wdz(k) = (zoop1_s+4.0e-2*zt(k))/86400.0/dzt(k)  ! [s-1]
!KK causing mass loss in microplastic!!!VVVVV
!           wd(k) = (wd0+4.0e-2*zt(k))/86400.0/dzt(k)  ! [s-1]
!           wdz(k) = (wdz0+4.0e-2*zt(k))/86400.0/dzt(k)  ! [s-1]
           gbio = gbio0
#  endif
#  if !defined O_plastic
           gbio = gbio0
#  endif
!-----------------------------------------------------------------------
!             call the npzd model
!-----------------------------------------------------------------------

              call npzd_src (tnpzd, nbio(k), dtbio(k), gl, bct, impo
     &,                      dzt(k), dayfrac, wd(k), gbio
#  if defined O_plastic
     &,                      agg_con, poo_con, expomp_a, impomp_a
     &,                      expomp_p, impomp_p,agg_upt,zoop_upt
     &,                      agg_rel, zoop_rel
#  endif
#  if defined O_zoop_det
     &,                      wdz(k), remiz, expoz, impoz
#  endif
     &,                      rkwz(k), nud
     &,                      snpzd, expo, graz, morp, morz, graz_Det
     &,                      graz_Z
#  if defined O_save_npzd
     &,                      npp, morpt, remi, excr
#   if defined O_kk_ballast
     &,                      impo_B, expo_B, remi_B, graz_Det_B
#   endif
#   if defined O_npzd_nitrogen
#    if !defined O_npzd_simple
     &,                      npp_D, graz_D, morp_D
#    endif
     &,                      nfix
#   endif
#   if defined O_cal_kk 
     &,                      npp_C, morpt_C, graz_C, morp_C
#    if defined O_npzd_fe_limitation
     &,                      felimit_C
#    endif
#    if defined O_npzd_iron
     &,                      thetamax_C, deffe_C
#    endif
#   endif
#   if defined O_kk_diat
     &,                      npp_Diat, morpt_Diat, graz_Diat, morp_Diat
#    if defined O_npzd_fe_limitation
     &,                      felimit_Diat
#    endif
#    if defined O_npzd_iron
     &,                      thetamax_Diat, deffe_Diat
#    endif
     &,                      sil_in
#  endif
#   if defined O_kk_caco3tr
     &,                      impocaco3,wc(k),expocaco3,dissl
     &,                      calpro
     &,                      calatt,dissk1
#   endif
#   if defined O_npzd_extra_diagnostics
     &,                      avej, avej_D, gmax, no3P, po4P, po4_D
#   endif
#  endif
#  if defined O_npzd_fe_limitation
     &,                      felimit, felimit_D
#  endif
     &,                      bctz
#  if defined O_npzd_iron
     &,                      expofe, impofe, remife
#   if defined O_npzd_iron_diagnostics
     &,                      feorgads, thetamax, deffe, feprime
     &,                      fecol
#    if defined O_kk_caco3tr
     &,                      feorgads_ca
#    endif
#   endif
     &,                      t(i,k,j,io2,taum1)
#   if defined O_npzd_chl
     &,                      chl, thetachl
#    if defined O_npzd_nitrogen
     &,                      chl_D
#    endif
#   endif
#  endif
#  if defined O_isotopes
     &,                      rtdic13(k), rtcaco3c13(k), caco3c13impo
     &,                      caco3c13expo, c13expo, c13impo, fcnpp
     &,                      n15expo, n15impo, sil30_in
     &,                      eps_assim, eps_excr, eps_nfix
     &,                      rtdiats30(k)
#  endif
     &                       )
! These are source/sink terms 
              snpzd(ibion) = snpzd(ibion)*rdtts(k)
              snpzd(ibiop) = snpzd(ibiop)*rdtts(k)
              snpzd(ibioz) = snpzd(ibioz)*rdtts(k)
              snpzd(ibiod) = snpzd(ibiod)*rdtts(k)
#  if defined O_plastic
              snpzd(ibiomp) = snpzd(ibiomp)*rdtts(k)
              snpzd(ibiompa) = snpzd(ibiompa)*rdtts(k)
              snpzd(ibiompp) = snpzd(ibiompp)*rdtts(k)
#  endif
#  if defined O_zoop_det
              snpzd(ibiodz) = snpzd(ibiodz)*rdtts(k)
#  endif
#  if defined O_npzd_nitrogen
              snpzd(ibiono3) = snpzd(ibiono3)*rdtts(k)
#   if !defined O_npzd_simple
              snpzd(ibiodiaz) = snpzd(ibiodiaz)*rdtts(k)
#   endif
#  endif
#  if defined O_cal_kk
              snpzd(ibioc) = snpzd(ibioc)*rdtts(k)
#  endif
#  if defined O_kk_diat
              snpzd(ibiodiat) = snpzd(ibiodiat)*rdtts(k)
#  endif
#  if defined O_kk_caco3tr
              expocaco3 = expocaco3*rnbio(k)
              snpzd(ibiocaco3) = snpzd(ibiocaco3)*rdtts(k)
#  endif
#  if defined O_kk_ballast
              snpzd(ibiod_B) = snpzd(ibiod_B)*rdtts(k)
              expo_B = expo_B*rnbio(k) !ballast detritus
#  endif
#  if defined O_npzd_iron
              snpzd(ibiodfe) = snpzd(ibiodfe)*rdtts(k)
              snpzd(ibiodetrfe) = snpzd(ibiodetrfe)*rdtts(k)
              expofe=expofe*rnbio(k)
#  endif
#  if defined O_isotopes
             c13expo = c13expo*rnbio(k) !N units
             n15expo = n15expo*rnbio(k) !N units
             caco3c13expo = caco3c13expo*rnbio(k)
             snpzd(ibiocaco3c13) = snpzd(ibiocaco3c13)*rdtts(k)
             snpzd(ibiodic13) = snpzd(ibiodic13)*rdtts(k)
             snpzd(ibiophytc13) = snpzd(ibiophytc13)*rdtts(k)
             snpzd(ibiodiazc13) = snpzd(ibiodiazc13)*rdtts(k)
             snpzd(ibiozoopc13) = snpzd(ibiozoopc13)*rdtts(k)
             snpzd(ibiococcc13) = snpzd(ibiococcc13)*rdtts(k)
             snpzd(ibiodiatc13) = snpzd(ibiodiatc13)*rdtts(k)
             snpzd(ibiodetrc13) =  snpzd(ibiodetrc13)*rdtts(k)

             snpzd(ibiodin15) = snpzd(ibiodin15)*rdtts(k)
             snpzd(ibiophytn15) = snpzd(ibiophytn15)*rdtts(k)
             snpzd(ibiodiazn15) = snpzd(ibiodiazn15)*rdtts(k)
             snpzd(ibiozoopn15) = snpzd(ibiozoopn15)*rdtts(k)
             snpzd(ibiococcn15) = snpzd(ibiococcn15)*rdtts(k)
             snpzd(ibiodiatn15) = snpzd(ibiodiatn15)*rdtts(k)
             snpzd(ibiodetrn15) =  snpzd(ibiodetrn15)*rdtts(k)

!             snpzd(ibiosil30) = snpzd(ibiosil30)*rdtts(k)
             snpzd(ibiodiats30) = snpzd(ibiodiats30)*rdtts(k)
#  endif
              expo = expo*rnbio(k)
#  if defined O_plastic
              expomp_a = expomp_a*rnbio(k)
              expomp_p = expomp_p*rnbio(k)
              rexpomp_a(i,k,j) = expomp_a
              rexpomp_p(i,k,j) = expomp_p
              ragg_upt(i,k,j) = agg_upt*rnbio(k)
              rzoop_upt(i,k,j) = zoop_upt*rnbio(k)
              ragg_rel(i,k,j) = agg_rel*rnbio(k)
              rzoop_rel(i,k,j) = zoop_rel*rnbio(k)
#  endif
#  if defined O_zoop_det
              expoz = expoz*rnbio(k)
#  endif

#  if defined O_save_npzd
              rgraz(i,k,j) = graz*rnbio(k) 
              rgraz_Det(i,k,j) = graz_Det*rnbio(k)
              rgraz_Z(i,k,j) = graz_Z*rnbio(k)
              rmorp(i,k,j) = morp*rnbio(k)
              rmorz(i,k,j) = morz*rnbio(k)
              rnpp(i,k,j) = npp*rnbio(k)
              rmorpt(i,k,j) = morpt*rnbio(k)
              rremi(i,k,j) = remi*rnbio(k)
              rexcr(i,k,j) = excr*rnbio(k)
#   if defined O_kk_ballast
              rgraz_Det_B(i,k,j) = graz_Det_B*rnbio(k)
              rremi_B(i,k,j) = remi_B*rdtts(k)
#   endif
#   if defined O_npzd_nitrogen
#    if !defined O_npzd_simple
              rnpp_D(i,k,j) = npp_D*rnbio(k)
              rgraz_D(i,k,j) = graz_D*rnbio(k)
              rmorp_D(i,k,j) = morp_D*rnbio(k)
#    endif
              rnfix(i,k,j) = nfix*rnbio(k)
#   endif
#   if defined O_kk_caco3tr
              rcalpro(i,k,j) = calpro*rnbio(k)
              rcalatt(i,k,j) = calatt*rnbio(k)
              rdissl(i,k,j) = dissl*rdtts(k)
#   endif
#   if defined O_kk_diat
              rnpp_Diat(i,k,j) = npp_Diat*rnbio(k)
              rgraz_Diat(i,k,j) = graz_Diat*rnbio(k)
              rmorp_Diat(i,k,j) = morp_Diat*rnbio(k)
              rmorpt_Diat(i,k,j) = morpt_Diat*rnbio(k)
#   endif
#   if defined O_cal_kk
              rnpp_C(i,k,j) = npp_C*rnbio(k)
              rgraz_C(i,k,j) = graz_C*rnbio(k)
              rmorp_C(i,k,j) = morp_C*rnbio(k)
              rmorpt_C(i,k,j) = morpt_C*rnbio(k)
#   endif
#   if defined O_npzd_iron
              rremife(i,k,j) = remife*rnbio(k)
#    if defined O_npzd_iron_diagnostics
              rfeorgads(i,k,j) = feorgads*rnbio(k)
#     if defined O_kk_caco3tr
              rfeorgads_ca(i,k,j) = feorgads_ca*rnbio(k)
#     endif
              rdeffe(i,k,j) = deffe*rnbio(k)
              rfeprime(i,k,j) = feprime*rnbio(k)
              rfecol(i,k,j) = fecol*rnbio(k)
#    endif
#    if defined O_npzd_chl
              rchl(i,k,j) = chl*rnbio(k)
#     if defined O_npzd_nitrogen
              rchl_D(i,k,j) = chl_D*rnbio(k)
#     endif
#    endif
#   endif
#   if defined O_zoop_det
              rremiz(i,k,j) = remiz*rnbio(k)
#   endif
#   if defined O_npzd_extra_diagnostics
              ravej(i,k,j) = avej*rnbio(k)
              ravej_D(i,k,j) = avej_D*rnbio(k)
              rgmax(i,k,j) = gmax*rnbio(k)
              rno3P(i,k,j) = no3P*rnbio(k)
              rpo4P(i,k,j) = po4P*rnbio(k)
              rpo4_D(i,k,j) = po4_D*rnbio(k)
#   endif
#  endif

!-----------------------------------------------------------------------
!             calculate detritus at the bottom and remineralize
!-----------------------------------------------------------------------

#  if defined O_npzd_subgridbathy
              sgb = sg_bathy(i,j,k)
              if (sgb .gt. 0) then
#   if defined O_sed
#    if defined O_kk_ballast
                if (addflxo .and. eots)
     &            sbc(i,jrow,irorg) = sbc(i,jrow,irorg)
     &               + (expo+expo_B
#    if defined O_zoop_det
     &                         + expoz
#    endif
     &               )*sgb*redctn*twodt(1)*dzt(k)
#    else 
                if (addflxo .and. eots)
     &            sbc(i,jrow,irorg) = sbc(i,jrow,irorg)
     &                              + (expo
#    if defined O_zoop_det
     &                              + expoz
#    endif
     &                              )*sgb*redctn*twodt(1)*dzt(k)
#    endif
#   endif
#   if defined O_save_npzd
#    if defined O_kk_ballast
                rremi(i,k,j) = rremi(i,k,j) + (expo+expo_B)*sgb 
#    else
                rremi(i,k,j) = rremi(i,k,j) + expo*sgb !free N
#    endif
#     if defined O_zoop_det
                rremiz(i,k,j) = rremiz(i,k,j) + expoz*sgb
#     endif
#   endif
#   if defined O_plastic
!                ragg_rel(i,k,j) = ragg_rel(i,k,j) + expomp_a*pbc*sgb
!                rzoop_rel(i,k,j) = rzoop_rel(i,k,j) + expomp_p*pbc*sgb
#   endif
#   if defined O_kk_ballast
                snpzd(ibion) = snpzd(ibion) + redptn*(expo+expo_B
#    if defined O_zoop_det
     &                         + expoz
#    endif
     &                         )*sgb 
#   else
                snpzd(ibion) = snpzd(ibion) + redptn*(expo
#    if defined O_zoop_det
     &                         + expoz
#    endif
     &                        )*sgb 
#   endif
#   if defined O_npzd_nitrogen
#    if defined O_kk_ballast
                snpzd(ibiono3) = snpzd(ibiono3) + (expo+expo_B
#     if defined O_zoop_det
     &                         + expoz
#     endif
     &                         )*sgb 
#    else
                snpzd(ibiono3) = snpzd(ibiono3) + (expo
#     if defined O_zoop_det
     &                         + expoz
#     endif
     &                         )*sgb 
#    endif
#    if defined O_benthic_denitrification
!-----------------------------------------------------------------------
!             benthic denitrification model of Bohlen et al., 2012, GBC 
!             NO3 is removed out of bottom water nitrate.
!             See Somes et al., 2012, BGS for additional details/results
!-----------------------------------------------------------------------
!     limit denitrification as nitrate approaches 0 uM
                no3flag = 0.5+sign(0.5,t(i,k,j,ino3,taum1)-trcmin)
                din15flag = 1.0
                lno3 = 0.5*tanh(t(i,k,j,ino3,taum1)*10 - 5.0)

                sg_bdeni = (0.06 + 0.19*0.99
     &                        **(max(t(i,k,j,io2,taum1)*1000,trcmin)
     &                        - max(t(i,k,j,ino3,taum1),trcmin)))
     &                        *max((expo
#     if defined O_zoop_det
     &                         + expoz
#     endif
#     if defined O_kk_ballast
     &                         + expo_B
#     endif
     &                        )*sgb,trcmin)*6.625

                sg_bdeni = min(sg_bdeni, sgb*(
#     if defined O_zoop_det
     &                     expoz +
#     endif
#     if defined O_kk_ballast
     &                     expo_B +
#     endif
     &                     expo))
                sg_bdeni = max(sg_bdeni, 0.)
                sg_bdeni = sg_bdeni*(0.5 + lno3)*no3flag*din15flag
                snpzd(ibiono3) = snpzd(ibiono3) - sg_bdeni
#    endif
#   endif
#   if defined O_kk_caco3tr
         src(i,k,j,isdic) = src(i,k,j,isdic)+expocaco3*1.e-3*sgb
         src(i,k,j,isalk) = src(i,k,j,isalk)
     &                         + 2.*expocaco3*1.e-3*sgb
#   endif
!----------------------------
! sediment transfer function for iron
!----------------------------
#  if defined O_npzd_iron
#   if defined O_kk_ballast
                fesed=fetopsed*bct*redptn*(expo+expo_B
#    if defined O_zoop_det
     &                         + expoz
#    endif
     &               )*sgb
#   else
                fesed=fetopsed*bct*redptn*(expo
#    if defined O_zoop_det
     &                         + expoz
#    endif
     &                )*sgb
#   endif   
                bfe=fesed
                snpzd(ibiodfe) = snpzd(ibiodfe) + fesed
                if (t(i,k,j,io2,taum1)*1e3 .lt. o2t) then
                  snpzd(ibiodfe) = snpzd(ibiodfe) + expofe*sgb
                  bfe=bfe+expofe*sgb
                  expofe = expofe - sgb * expofe
                endif
                rremife(i,k,j) = rremife(i,k,j) + bfe
#   if defined O_npzd_iron_diagnostics
                rfesed=fesed
                rbfe=bfe
#   endif
#  endif
                expo = expo - sgb * expo
#  if defined O_zoop_det
                expoz = expoz - sgb * expoz
#  endif
#  if defined O_kk_ballast
                expo_B = expo_B - sgb * expo_B
#  endif
#  if defined O_isotopes
         snpzd(ibiodic13) = snpzd(ibiodic13) 
     &                     + c13expo*sgb*redctn
     &                     + caco3c13expo*1.e-3*sgb !added KK
! KK missing 'expo' or not? MOBI not consistent between C13 and N15!
         snpzd(ibiodin15) = snpzd(ibiodin15) + n15expo*sgb
         c13expo = c13expo - sgb*c13expo
         n15expo = n15expo - sgb*n15expo
         caco3c13expo = caco3c13expo - sgb*caco3c13expo
#  endif
#  if defined O_kk_caco3tr
!----------------------------------------------------------------------
!             calculate calcite at the bottom and dissolve
!---------------------------------------------------------------------
#      if defined O_sed
                if (addflxo .and. eots)
     &           sbc(i,jrow,ircal) = sbc(i,jrow,ircal) 
     &           + expocaco3*sgb*twodt(1)*dzt(k)*1.e-3 !mol C
#      endif
               expocaco3 = expocaco3 - sgb*expocaco3
#  endif
              endif !end of bathy condition
# endif
# if defined O_save_npzd
              rexpo(i,k,j) = expo
# if defined O_zoop_det
              rexpoz(i,k,j) = expoz
# endif
# if defined O_kk_caco3tr
#  if defined O_kk_ballast
              rexpo_B(i,k,j) = expo_B
#  endif
              rexpocaco3(i,k,j) = expocaco3
# endif
# if defined O_npzd_iron
              rexpofe(i,k,j) = expofe
# endif
# endif
              if (k .eq. kmt(i,jrow)) then
#   if defined O_sed
#    if defined O_kk_ballast
                if (addflxo .and. eots)
     &            sbc(i,jrow,irorg) = sbc(i,jrow,irorg)
     &               + (expo+expo_B
#    if defined O_zoop_det
     &                         + expoz
#    endif
     &               )*redctn*twodt(1)*dzt(k)
#    else 
                if (addflxo .and. eots)
     &            sbc(i,jrow,irorg) = sbc(i,jrow,irorg)
     &                              + (expo
#    if defined O_zoop_det
     &                         + expoz
#    endif
     &                              )*redctn*twodt(1)*dzt(k)
#    endif
#   endif
#   if defined O_save_npzd
#    if defined O_kk_ballast
                rremi(i,k,j) = rremi(i,k,j) + expo + expo_B 
#    else
                rremi(i,k,j) = rremi(i,k,j) + expo !free nitrogen
#    endif
#    if defined O_zoop_det
                rremiz(i,k,j) = rremiz(i,k,j) + expoz
#    endif
#   endif
#   if defined O_plastic
                ragg_rel(i,k,j) = ragg_rel(i,k,j) + expomp_a*pbc
                rzoop_rel(i,k,j) = rzoop_rel(i,k,j) + expomp_p*pbc
#   endif
#   if defined O_kk_ballast
                snpzd(ibion) = snpzd(ibion) + redptn*(expo+expo_B) !nitrogen
#   else
                snpzd(ibion) = snpzd(ibion) + redptn*expo !nitrogen
#   endif
#   if defined O_zoop_det
                snpzd(ibion) = snpzd(ibion) + redptn*expoz
#   endif
#   if defined O_npzd_nitrogen
#    if defined O_kk_ballast
                snpzd(ibiono3) = snpzd(ibiono3) + (expo+expo_B)
#    else
                snpzd(ibiono3) = snpzd(ibiono3) + expo
#    endif
#    if defined O_zoop_det
                snpzd(ibiono3) = snpzd(ibiono3) + expoz
#    endif
#    if defined O_benthic_denitrification
!-----------------------------------------------------------------------
!             benthic denitrification model of Bohlen et al., 2012, GBC 
!             NO3 is removed out of bottom water nitrate.
!             See Somes et al., 2012, BGS for additional details/results
!-----------------------------------------------------------------------
!     limit denitrification as nitrate approaches 0 uM
                no3flag = 0.5+sign(0.5,t(i,k,j,ino3,taum1)-trcmin)
                din15flag = 1.0
                lno3 = 0.5*tanh(t(i,k,j,ino3,taum1)*10 - 5.0)

                bdeni = (0.06 + 0.19*0.99
     &                        **(max(t(i,k,j,io2,taum1)*1000,trcmin)
     &                        - max(t(i,k,j,ino3,taum1),trcmin)))
     &                        *max((expo
#     if defined O_zoop_det
     &                         + expoz
#     endif
#     if defined O_kk_ballast
     &                         + expo_B
#     endif
     &                        ),trcmin)*6.625

                bdeni = min(bdeni,(
#     if defined O_zoop_det
     &                     expoz +
#     endif
#     if defined O_kk_ballast
     &                     expo_B +
#     endif
     &                     expo))
                bdeni = max(bdeni, 0.)
                bdeni = bdeni*(0.5 + lno3)*no3flag*din15flag
                snpzd(ibiono3) = snpzd(ibiono3) - bdeni
#    endif
#   endif
#   if defined O_plastic
                snpzd(ibiomp) = snpzd(ibiomp) 
     &                        + (expomp_a
     &                                   +expomp_p
     &                           )*pbc
#   endif
#   if defined O_npzd_iron
#    if defined O_kk_ballast
                fesed=fetopsed*bct*redptn*(expo+expo_B
#    if defined O_zoop_det
     &               + expoz
#    endif
     &                )
#    else
                fesed=fetopsed*bct*redptn*(expo
#    if defined O_zoop_det
     &               + expoz
#    endif
     &               )
#    endif   
                bfe=fesed
                snpzd(ibiodfe) = snpzd(ibiodfe) + fesed
                if (t(i,k,j,io2,taum1)*1e3 .lt. o2t) then
                  snpzd(ibiodfe) = snpzd(ibiodfe) + expofe
                  bfe=bfe+expofe
                endif
                rremife(i,k,j) = rremife(i,k,j) + bfe
#    if defined O_npzd_iron_diagnostics
                rfesed=fesed
                rbfe=bfe
#    endif
#   endif
#   if defined O_kk_caco3tr
!----------------------------------------------------------------------
!             calculate calcite at the bottom and dissolve
!---------------------------------------------------------------------
#      if defined O_sed
                if (addflxo .and. eots)
     &           sbc(i,jrow,ircal) = sbc(i,jrow,ircal) 
     &           + expocaco3*twodt(1)*dzt(k)*1.e-3 !mol C
#      endif
#   endif
              endif !end of bottom box loop

!-----------------------------------------------------------------------
!             set source/sink terms
!-----------------------------------------------------------------------
              src(i,k,j,ispo4) = snpzd(ibion)
              src(i,k,j,isphyt) = snpzd(ibiop)
              src(i,k,j,iszoop) = snpzd(ibioz)
              src(i,k,j,isdetr) = snpzd(ibiod)
#  if defined O_plastic
              src(i,k,j,ismp) = snpzd(ibiomp)
              src(i,k,j,ismpa) = snpzd(ibiompa)
              src(i,k,j,ismpp) = snpzd(ibiompp)
#  endif
#  if defined O_zoop_det
              src(i,k,j,isdetrz) = snpzd(ibiodz)
#  endif
#  if defined O_npzd_nitrogen
              src(i,k,j,isno3) = snpzd(ibiono3)
#   if !defined O_npzd_simple
              src(i,k,j,isdiaz) = snpzd(ibiodiaz)
#   endif
#  endif
#  if defined O_cal_kk
              src(i,k,j,iscocc) = snpzd(ibioc)
#  endif
#  if defined O_kk_diat
              src(i,k,j,isdiat) = snpzd(ibiodiat)
#  endif
#  if defined O_kk_caco3tr
              src(i,k,j,iscaco3) = snpzd(ibiocaco3)
#  endif
#  if defined O_kk_ballast
              src(i,k,j,isdetr_B) = snpzd(ibiod_B)
#  endif
#  if defined O_npzd_iron             
              src(i,k,j,isdfe) = snpzd(ibiodfe) + fe_hydr(i,j,k)
              src(i,k,j,isdetrfe) = snpzd(ibiodetrfe)
#  endif
#  if defined O_isotopes
              src(i,k,j,iscaco3c13) = snpzd(ibiocaco3c13)
              src(i,k,j,isdic13) = snpzd(ibiodic13)
              src(i,k,j,isphytc13) = snpzd(ibiophytc13)
              src(i,k,j,isdiazc13) = snpzd(ibiodiazc13)
              src(i,k,j,iszoopc13) = snpzd(ibiozoopc13)
              src(i,k,j,iscoccc13) = snpzd(ibiococcc13)
              src(i,k,j,isdiatc13) = snpzd(ibiodiatc13)
              src(i,k,j,isdetrc13) = snpzd(ibiodetrc13)

              src(i,k,j,isdin15) = snpzd(ibiodin15)
              src(i,k,j,isphytn15) = snpzd(ibiophytn15)
              src(i,k,j,isdiazn15) = snpzd(ibiodiazn15)
              src(i,k,j,iszoopn15) = snpzd(ibiozoopn15)
              src(i,k,j,iscoccn15) = snpzd(ibiococcn15)
              src(i,k,j,isdiatn15) = snpzd(ibiodiatn15)
              src(i,k,j,isdetrn15) = snpzd(ibiodetrn15)

              src(i,k,j,isdiats30) = snpzd(ibiodiats30)
#  endif 
!-----------------------------------------------------------------------
!             production, export, and dissolution of opal
!-----------------------------------------------------------------------
#  if defined O_kk_si

!              sipr = 0.
!#   if defined O_kk_variable_sipr && (defined O_npzd_iron || defined O_fe_limitation)
!     From HAMOCC5 (Aumont et al, GBC, 2003) values 0-0.5, mostly ~0.1
!              sipr = sipr0 * min(1.,sil_in/k1si)
!#    if defined O_fe_limitation || defined O_npzd_iron
!     &             * (4.-3.*min(1.,
!#     if !defined O_npzd_iron
!     & fe_conc
!#     else
!     & max(t(i,k,j,idfe,taum1),trcmin)
!#     endif
!     &/kfe_Diat))
!#    endif
!#   else
!              sipr = sipr0
!#   endif
!             production of opal.
#  if defined O_kk_diat 
              dprop = (morp_Diat + graz_Diat*(1.-gamma1))
     &              *sipr*redctn*rnbio(k) !umol Si cm-3 s-1
#  elif defined O_npzd_simple
              dprop = (morp + graz*(1.-gamma1))
     &              *sipr*redctn*rnbio(k) !umol Si cm-3 s-1
#  endif
              prop = prop + dprop*dzt(k) !umol Si cm-2 s-1
              src(i,k,j,issil) = -dprop + si_hydr(i,j,k)*si_h_sol !umol Si cm-3 s-1
#   if defined O_isotopes
! KK check ratio for inputs 
! Kleine et al 2018 says hydrothermal enrichment higly variable
              b2si = rs30std - eps_hydr*rs30std/1000.
              fcsi = b2si/(1+b2si)
              src(i,k,j,issil30) = -rtsil30(k)*dprop 
     &                         + fcsi*si_hydr(i,j,k)*si_h_sol !umol Si cm-3 s-1
#   endif
              silwflx = silwflx - si_hydr(i,j,k)*si_h_sol
     &        *dzt(k)*dxt(i)*dyt(j)*cst(j) !umol Si/s
#  endif
!-----------------------------------------------------------------------
!             production of calcite
!-----------------------------------------------------------------------
!     Still calculated but now only used if kk_caco3tr is not implemented!
#  if defined O_cal_kk
          dprca = (capr*morp_C+capr*morz
     &          +(capr*graz_C+capr*graz_Z)*(1.-gamma1))
     &              *redctn*rnbio(k)
#  else
          dprca = (morp+morz+(graz+graz_Z)*(1.-gamma1))  !molC/m3/s
     &            *capr*redctn*rnbio(k)
#  endif
              prca = prca + dprca*dzt(k) !mol C/m3/s
#  if defined O_isotopes
         prca13 = prca13 + dprca*dzt(k)*rtdic13(k)       
! sources and sinks calculated below
#  endif
!             These are sources and sinks of DIC (i.e. remin - pp)
!             all are based on po4 uptake and remineralization
!             dprca is a correction term
#  if defined O_carbon
#   if defined O_kk_caco3tr
!              calculated below
#   else
              src(i,k,j,isdic) = snpzd(ibion)*redctp - dprca
#   endif
#  endif
#  if defined O_npzd_alk
#   if defined O_kk_caco3tr
!              calculated below
#   else
              src(i,k,j,isalk) = -snpzd(ibion)*redntp*1.e-3 - 2.*dprca
#   endif
#  endif
!       print*,'KK3 tr', src(i,k,j,:)

#  if (defined O_time_averages || defined O_TMM) && defined O_save_npzd

!-----------------------------------------------------------------------
!             accumulate time averages
!-----------------------------------------------------------------------

              if (timavgperts .and. .not. euler2) then
                ta_rnpp(i,k,jrow) = ta_rnpp(i,k,jrow) + rnpp(i,k,j)
                ta_rgraz(i,k,jrow) = ta_rgraz(i,k,jrow) + rgraz(i,k,j)
                ta_rgraz_Z(i,k,jrow) = ta_rgraz_Z(i,k,jrow) 
     &                               + rgraz_Z(i,k,j)
                ta_rgraz_Det(i,k,jrow) = ta_rgraz_Det(i,k,jrow)
     &                                 + rgraz_Det(i,k,j)
                ta_rmorp(i,k,jrow) = ta_rmorp(i,k,jrow) + rmorp(i,k,j)
                ta_rmorpt(i,k,jrow)= ta_rmorpt(i,k,jrow) + rmorpt(i,k,j)
                ta_rmorz(i,k,jrow) = ta_rmorz(i,k,jrow) + rmorz(i,k,j)
                ta_rexcr(i,k,jrow) = ta_rexcr(i,k,jrow) + rexcr(i,k,j)
#   if defined O_zoop_det
                ta_rremiz(i,k,jrow) = ta_rremiz(i,k,jrow) +rremiz(i,k,j)
               ta_rexpoz(i,k,jrow) = ta_rexpoz(i,k,jrow) + rexpoz(i,k,j) !detrz
#   endif
#   if defined O_cal_kk
                ta_rnpp_C(i,k,jrow) = ta_rnpp_C(i,k,jrow)
     &                              + rnpp_C(i,k,j)
                ta_rgraz_C(i,k,jrow) = ta_rgraz_C(i,k,jrow) 
     &                               + rgraz_C(i,k,j)
                ta_rmorp_C(i,k,jrow) = ta_rmorp_C(i,k,jrow) 
     &                               + rmorp_C(i,k,j)
                ta_rmorpt_C(i,k,jrow)= ta_rmorpt_C(i,k,jrow) 
     &                               + rmorpt_C(i,k,j)
#   endif
#   if defined O_kk_diat
                ta_rnpp_Diat(i,k,jrow) = ta_rnpp_Diat(i,k,jrow)
     &                              + rnpp_Diat(i,k,j)
                ta_rgraz_Diat(i,k,jrow) = ta_rgraz_Diat(i,k,jrow)
     &                               + rgraz_Diat(i,k,j)
                ta_rmorp_Diat(i,k,jrow) = ta_rmorp_Diat(i,k,jrow)
     &                               + rmorp_Diat(i,k,j)
                ta_rmorpt_Diat(i,k,jrow)= ta_rmorpt_Diat(i,k,jrow)
     &                               + rmorpt_Diat(i,k,j)
#   endif
#  if defined O_kk_caco3tr
                ta_romca(i,k,jrow) =  ta_romca(i,k,jrow) + romca(i,k,j) 
                ta_rco3(i,k,jrow) = ta_rco3(i,k,jrow) + rco3(i,k,j)
                ta_rco3_sat(i,k,jrow) = ta_rco3_sat(i,k,jrow) 
     &                                + rco3_sat(i,k,j)
                ta_rdel_sat(i,k,jrow) = ta_rdel_sat(i,k,jrow) 
     &                                + rdel_sat(i,k,j)
                ta_rcalatt(i,k,jrow) = ta_rcalatt(i,k,jrow)
     &                               + rcalatt(i,k,j)
                ta_rprocal(i,k,jrow) = ta_rprocal(i,k,jrow)
     &                               + rcalpro(i,k,j)
#  endif
#  if defined O_kk_ballast
                ta_rgraz_Det_B(i,k,jrow) = ta_rgraz_Det_B(i,k,jrow)
     &                                 + rgraz_Det_B(i,k,j)
#  endif
#   if defined O_npzd_nitrogen
#    if !defined O_npzd_simple
                ta_rnpp_D(i,k,jrow) = ta_rnpp_D(i,k,jrow)
     &                              + rnpp_D(i,k,j)
                ta_rgraz_D(i,k,jrow) = ta_rgraz_D(i,k,jrow)
     &                               + rgraz_D(i,k,j)
                ta_rmorp_D(i,k,jrow) = ta_rmorp_D(i,k,jrow)
     &                               + rmorp_D(i,k,j)
#    endif
                ta_rnfix(i,k,jrow) = ta_rnfix(i,k,jrow) + rnfix(i,k,j)
#   endif
#   if defined O_npzd_iron
#    if defined O_npzd_iron_diagnostics
                ta_rfeorgads(i,k,jrow) = ta_rfeorgads(i,k,jrow)
     &                              + rfeorgads(i,k,j)
#    if defined O_kk_caco3tr
                ta_rfeorgads_ca(i,k,jrow) = ta_rfeorgads_ca(i,k,jrow)
     &                              + rfeorgads_ca(i,k,j)
#    endif
                ta_rdeffe(i,k,jrow) = ta_rdeffe(i,k,jrow)
     &                              + rdeffe(i,k,j)
                ta_rfeprime(i,k,jrow) = ta_rfeprime(i,k,jrow)
     &                              + rfeprime(i,k,j)
                ta_rfesed(i,k,jrow) = ta_rfesed(i,k,jrow)
     &                              + rfesed(i,k,j)
                ta_rbfe(i,k,jrow) = ta_rbfe(i,k,jrow)
     &                              + rbfe(i,k,j)
                ta_rfecol(i,k,jrow) = ta_rfecol(i,k,jrow)
     &                              + rfecol(i,k,j)
#    endif
#    if defined O_npzd_chl
                ta_rchl(i,k,jrow) = ta_rchl(i,k,jrow)
     &                              + rchl(i,k,j)
#     if defined O_npzd_nitrogen
                ta_rchl_D(i,k,jrow) = ta_rchl_D(i,k,jrow)
     &                              + rchl_D(i,k,j)
#     endif
#    endif
#   endif
#   if defined O_npzd_extra_diagnostics
                ta_ravej(i,k,jrow) = ta_ravej(i,k,jrow)
     &                              + ravej(i,k,j)
                ta_ravej_D(i,k,jrow) = ta_ravej_D(i,k,jrow)
     &                              + ravej_D(i,k,j)
                ta_rgmax(i,k,jrow) = ta_rgmax(i,k,jrow)
     &                              + rgmax(i,k,j)
                ta_rno3P(i,k,jrow) = ta_rno3P(i,k,jrow)
     &                              + rno3P(i,k,j)
                ta_rpo4P(i,k,jrow) = ta_rpo4P(i,k,jrow)
     &                              + rpo4P(i,k,j)
                ta_rpo4_D(i,k,jrow) = ta_rpo4_D(i,k,jrow)
     &                              + rpo4_D(i,k,j)
#   endif
              endif
#  endif
!             calculate total export to get total import for next layer
              expo = expo*dzt(k)
#  if defined O_plastic
              expomp_a = expomp_a*dzt(k)
              expomp_p = expomp_p*dzt(k)
              ta_rexpomp_a(i,k,jrow) = ta_rexpomp_a(i,k,jrow) 
     &                               + rexpomp_a(i,k,j)
              ta_rexpomp_p(i,k,jrow) = ta_rexpomp_p(i,k,jrow) 
     &                               + rexpomp_p(i,k,j)
              ta_ragg_upt(i,k,jrow) = ta_ragg_upt(i,k,jrow)
     &                              + ragg_upt(i,k,j)
              ta_rzoop_upt(i,k,jrow) = ta_rzoop_upt(i,k,jrow)
     &                              + rzoop_upt(i,k,j)
              ta_ragg_rel(i,k,jrow) = ta_ragg_rel(i,k,jrow)
     &                              + ragg_rel(i,k,j)
              ta_rzoop_rel(i,k,jrow) = ta_rzoop_rel(i,k,jrow)
     &                              + rzoop_rel(i,k,j)
#  endif
#  if defined O_zoop_det
              expoz = expoz*dzt(k)
#  endif
#  if defined O_kk_ballast
              expo_B = expo_B*dzt(k)  !ballast detritus
#  endif
#  if defined O_kk_caco3tr
              expocaco3 = expocaco3*dzt(k)
#  endif
#  if defined O_npzd_iron
              expofe = expofe*dzt(k)
#  endif
#  if defined O_isotopes
              c13expo = c13expo*dzt(k)
              n15expo = n15expo*dzt(k)
              caco3c13expo = caco3c13expo*dzt(k)
#  endif  
            enddo
#  if defined O_phyt_albedo
! sum up all column biomass and convert to concentration for surface layer
! these numbers go to fluxes.f
! this creates an assumption of positive buoyancy and should be explored
! a better approach may be to use surface biomass values only
            caco3in(i,j) = caco3in(i,j)*dztr(1) !mmolC/m3
            phin(i,j) = phin(i,j)*dztr(1)*redctn*1.e3 !mmolC/m3            
#  endif
#  if defined O_npzd_o2
            kmx = kmt(i,jrow)
            do k=1,kmx
!             limit oxygen consumption below concentrations of
!             5umol/kg as recommended in OCMIP
! KK try optimising this
              fo2 = 0.5*tanh(t(i,k,j,io2,taum1)*1000. - o2t)
!             sink of oxygen
! O2 is needed to generate the equivalent of NO3 from N2 during N2 fixation
! 0.5 H2O + 0.5 N2+1.25O2 -> HNO3
! note that so2 is -dO2/dt
              so2 = src(i,k,j,ispo4)*redotp
#   if defined O_npzd_nitrogen
     &            + rnfix(i,k,j)*1.25e-3
#   endif
              src(i,k,j,iso2) = -so2*(0.5 + fo2)
#   if defined O_npzd_nitrogen
!             add denitrification as source term for NO3
              no3flag = 0.5+sign(0.5,t(i,k,j,ino3,taum1)-trcmin)
!             800 = 0.8*1000 = (elec/mol O2)/(elec/mol NO3)*(mmol/mol)
!KK switch to namelist parameter for sensitivity testing
              deni = max(0., 800.*no3flag*so2*(0.5 - fo2))
              src(i,k,j,isno3) = src(i,k,j,isno3) - deni
! Correct the ALK stoichiometry to account for N2 fixation
#     if defined O_npzd_nitrogen_paulmier_stoich
! New stoichiometric model parameters formulated as in Paulmier et al. 2009 BG
              Paulmier_a = 1.e3*redctn
              Paulmier_z = 4*1.e3*redotp - 4.*Paulmier_a - 8.*redntp
              Paulmier_R0 = Paulmier_a + 0.25*Paulmier_z
              src(i,k,j,isalk) = src(i,k,j,isalk)
     &        + no3flag*src(i,k,j,ispo4)*(0.5-fo2)
     &        *(4./5.*Paulmier_R0 + (3./5. +1.)*redntp) * 1.e-3
#     else
! Now account for N2 fixation (ALK production is tied to PO4 change and
! thus in the case of N2 fixation was not correct without this fix).
              src(i,k,j,isalk) = src(i,k,j,isalk)
     &         + deni*1.e-3
#     endif
     &         - rnfix(i,k,j)*1.e-3 
#    if defined O_benthic_denitrification
     &                        + (sg_bdeni + bdeni)*1.e-3
#    endif
#    if defined O_save_npzd
              rdeni(i,k,jrow) = deni
#     if defined O_benthic_denitrification
     &                        + sg_bdeni + bdeni
#     endif
#    endif
#    if defined O_isotopes
!        calculate isotope effect of water column denitrification
         uno3 = deni*twodt(k)/max(trcmin,t(i,k,j,ino3,taum1))
         uno3 = min(uno3, 0.999)
         uno3 = max(uno3, trcmin)
         rno3 = max(t(i,k,j,idin15,taum1),trcmin*rn15std/(1+rn15std))
     &             / max(t(i,k,j,ino3,taum1)-t(i,k,j,idin15,taum1),
     &                   trcmin*rn15std/(1+rn15std))
         rno3 = min(rno3, 2.*rn15std)
         rno3 = max(rno3, rn15std/2.)
! benthic denitrification !skip for now KK
!         eps_bdeni = eps_bdeni
!         bbdeni = rno3 - eps_bdeni*rno3/1000.
! KK check this
         bwcdeni = rno3 + eps_wcdeni*(1-uno3)/uno3*log(1-uno3)*rno3/1000.
         src(i,k,j,isdin15) = src(i,k,j,isdin15) 
     &                      - (bwcdeni/(1+bwcdeni))*deni
#    endif
#   endif
            enddo
#  endif 

!-----------------------------------------------------------------------
!           remineralize implicit calcite and opal
!-----------------------------------------------------------------------
            kmx = kmt(i,jrow)
            do k=1,kmx-1
#  if defined O_carbon
#   if defined O_kk_caco3tr 
              src(i,k,j,isdic) = src(i,k,j,isdic)
     &                         + src(i,k,j,ispo4)*redctp !mmol/m3/s
     &                         + rdissl(i,k,j)*1.e-3
     &                         - rcalpro(i,k,j)*1.e-3 
     &                         - rcalatt(i,k,j)*1.e-3 
              old_diss = prca*rcak(k)
#   else
              src(i,k,j,isdic) = src(i,k,j,isdic) + prca*rcak(k)
#   endif
#  endif
#  if defined O_isotopes
         src(i,k,j,isdic13) = src(i,k,j,isdic13)
     &        + rdissl(i,k,j)*1.e-3*rtcaco3c13(k)
     &        - rcalpro(i,k,j)*1.e-3*rtdic13(k)
     &        - rcalatt(i,k,j)*1.e-3*rtdic13(k)
#  endif
#  if defined O_npzd_alk
#   if defined O_kk_caco3tr
              src(i,k,j,isalk) = src(i,k,j,isalk) 
     &                         - src(i,k,j,ispo4)*redntp*1.e-3 
     &                         + 2.*rdissl(i,k,j)*1.e-3
     &                         - 2.*rcalpro(i,k,j)*1.e-3
     &                         - 2.*rcalatt(i,k,j)*1.e-3 
#   else
              src(i,k,j,isalk) = src(i,k,j,isalk)+2.*prca*rcak(k)
#   endif
#  endif
            enddo
#  if defined O_sed
#   if defined O_kk_caco3tr
!      calculated earlier
#   else
            if (addflxo .and. eots)
     &        sbc(i,jrow,ircal) = sbc(i,jrow,ircal) + (prca*rcab(kmx)
     &                          - prca*rcak(kmx))*twodt(1)*dzt(kmx)
#   endif
#  endif
#  if defined O_carbon
#   if defined O_kk_caco3tr
              src(i,kmx,j,isdic) = src(i,kmx,j,isdic)
     &                         + src(i,kmx,j,ispo4)*redctp 
     &                         + rdissl(i,kmx,j)*1.e-3
     &                         - rcalpro(i,kmx,j)*1.e-3
     &                         - rcalatt(i,kmx,j)*1.e-3
     &                         + rexpocaco3(i,kmx,j)*1.e-3
            old_diss = prca*rcab(kmx)
#   else
            src(i,kmx,j,isdic) = src(i,kmx,j,isdic)+prca*rcab(kmx)
#   endif
#  endif
#  if defined O_isotopes
         src(i,kmx,j,isdic13) = src(i,kmx,j,isdic13)
     &        + rdissl(i,kmx,j)*1.e-3*rtcaco3c13(kmx)
     &        - rcalpro(i,kmx,j)*1.e-3*rtdic13(kmx)
     &        - rcalatt(i,kmx,j)*1.e-3*rtdic13(kmx)
     &        + rexpocaco3(i,kmx,j)*1.e-3*rtcaco3c13(kmx)
#  endif
#  if defined O_npzd_alk
#   if defined O_kk_caco3tr
              src(i,kmx,j,isalk) = src(i,kmx,j,isalk)
     &                         - src(i,kmx,j,ispo4)*redntp*1.e-3 
     &                         + 2.*rdissl(i,kmx,j)*1.e-3
     &                         - 2.*rcalpro(i,kmx,j)*1.e-3  
     &                         - 2.*rcalatt(i,kmx,j)*1.e-3
     &                         + 2.*rexpocaco3(i,kmx,j)*1.e-3
#   else
            src(i,kmx,j,isalk) = src(i,kmx,j,isalk)+2.*prca*rcab(kmx)
#   endif
#  endif
#  if defined O_kk_si
!  Test some dissolution parameters, a default exponential from Aumont et al. (2003)
!  can be found in setmom.F    
      expo_op = prop  !flux

      do k=1,kmx
       temp_2(k) = t(i,k,j,itemp,taum1)+273.15 !convert to deg K 
!    (1) From Gnanadesikan 1999 !lower rates than HAMOCC
!     might over-estimate dissl in deep ocean (discussed in Ridgwell et al 2002)       
!       ropk(k) = dzt(k)*si_dis/ws0*exp(-11481./temp_2(k)) !unitless theoretical dissl
!     (2) From PlankTOM10, includes oxygen. Similar to HAMOCC except it always takes the minimum=1 and you end up with high dissolution.Plus there is no evidence of an oxygen dependency! 
!       ropk(k) = min(1.,si_dis*exp(-11200./temp_2(k)))
 !    &           *avail_oxy(k)/ws0*dzt(k)      
!     (3) Let's include the microbial loop instead. I made this up.
       ropk(k) = dzt(k)*si_dis/ws0*bctd(i,k,j) !exp(-11481./temp_2(k))
!      if(rmorpt_Diat(i,k,j).gt.0)then
!       ropk(k) = ropk(k)+ropk(k)*bctd(i,k,j)!*0.5!0.1     
!       ropk(k)=ropk(k)*bctd(i,k,j)
!      endif
       if(ropk(k).gt.1) ropk(k) = 1.
       op_diss = expo_op*ropk(k) !umol Si cm-2 s-1
       expo_op = expo_op - op_diss
       src(i,k,j,issil) = src(i,k,j,issil) + op_diss*dztr(k) !umol Si cm-3 s-1
#   if defined O_isotopes
! KK change value to one from Dermarest et al 2009
              b2si = rs30std - eps_diss*rs30std/1000.
              fcsi = b2si/(1+b2si)
       src(i,k,j,issil30) = src(i,k,j,issil30)
     &                    +fcsi*op_diss*dztr(k)
#   endif
#   if defined O_npzd_subgridbathy
              sgb = sg_bathy(i,j,k)
              if (sgb .gt. 0) then
! Use a limit (2 mmol Si /m^2/day) and ratios from Sarmiento & Gruber (2006)
! Buried Si is lost from the system!!!
                if (expo_op .le. 0.2/24/60/60) then
                  bsi = 0.05*expo_op*sgb !buried sil umol/cm^2/s
                else
                  bsi = 0.3*expo_op*sgb
                endif
                expo_op = expo_op - bsi
                silwflx = silwflx + bsi*dxt(i)*dyt(j)*cst(j)
              endif
#   else
           if (k .eq. kmx) then
                if (expo_op .le. 0.2/24/60/60) then
                  bsi = 0.05*expo_op !buried sil
                else
                  bsi = 0.3*expo_op
                endif
                expo_op = expo_op - bsi
                silwflx = silwflx + bsi*dxt(i)*dyt(j)*cst(j)
           endif
#   endif
#   if (defined O_time_averages || defined O_TMM) && defined O_save_npzd
           if (timavgperts .and. .not. euler2) then
            ta_rbur_sil(i,k,jrow) = ta_rbur_sil(i,k,jrow) + bsi*dztr(k) !umol/cm3/s
           endif
#   endif 
      enddo ! K loop
      src(i,kmx,j,issil) = src(i,kmx,j,issil) + expo_op*dztr(kmx) !put the leftovers back into ocean
#   if defined O_isotopes
      src(i,kmx,j,issil30) = src(i,kmx,j,issil30) 
     &               + expo_op*dztr(kmx)*rtdiats30(k)
#   endif
#  endif
#  if defined O_plastic
! Reverse k loop for rising free plastic particles
       expomp = 0.
       impomp = 0.
       kmx = min(kmt(i,jrow), kpzd)
      do k = kmx,1,-1
        impomp = expomp*dztr(k) !conc
! assume 30% rises
        expomp = wmp(k)*max(0.,t(i,k,j,imp,taum1)*prise) !conc/s
       if (k.gt.1) then
        src(i,k,j,ismp) = src(i,k,j,ismp) + impomp - expomp
       else
         src(i,1,j,ismp) = src(i,1,j,ismp) + impomp
       endif
        expomp = expomp*dzt(k) !flux
      enddo          
#  endif
#  if defined O_idealage
!     Here I should add the age change of the ideal age tracer
            kmx = kmt(i,jrow)
            do k=1,kmx
              src(i,k,j,isidealage) = 1.0 / 3.6e3 /24.
            enddo
#  endif
#  if (defined O_time_averages || defined O_TMM) && defined O_save_npzd

!-----------------------------------------------------------------------
!           accumulate time averages for full depth variables
!-----------------------------------------------------------------------
            if (timavgperts .and. .not. euler2) then
              kmx = kmt(i,jrow)
              expo = prca
#   if defined O_kk_si
              ta_rproopl(i,jrow) =  ta_rproopl(i,jrow) + prop
#   endif
# if !defined O_kk_caco3tr
              ta_rprocal(i,jrow) =  ta_rprocal(i,jrow) + prca
# endif
#  if defined O_plastic_conc
!    not time dependent for now
              ta_rplastic(i,jrow) = plastic(i,j,2) !only meaningful for concentration-forced
              ta_wd(i,jrow) = phyt_r-wd0 !anomaly !wd(1)
              ta_wdz(i,jrow) = zoop_r-wdz0!wdz(1)
#  endif
              do k=1,kmx
                expo = expo*dztr(k)
#   if defined O_kk_caco3tr
                old_diss = old_diss*dztr(k)
                ta_rold_diss(i,k,jrow) = ta_rold_diss(i,k,jrow)
     &                                 + old_diss
                old_diss = old_diss*dzt(k)
#   endif
#   if defined O_kk_ballast
                ta_rremi_B(i,k,jrow) = ta_rremi_B(i,k,jrow) 
     &                                 + rremi_B(i,k,j)
                ta_rexpo_B(i,k,jrow) = ta_rexpo_B(i,k,jrow) 
     &                                 + rexpo_B(i,k,j) !ballast detritus
#   endif
                ta_rremi(i,k,jrow) = ta_rremi(i,k,jrow) + rremi(i,k,j)
                ta_rexpo(i,k,jrow) = ta_rexpo(i,k,jrow) + rexpo(i,k,j)
                expo = expo - prca*rcak(k)
                ta_rexpocal(i,k,jrow) = ta_rexpocal(i,k,jrow) + expo
                expo = expo*dzt(k)
#   if defined O_kk_caco3tr
              ta_rdissl(i,k,jrow) = ta_rdissl(i,k,jrow) 
     &                              + rdissl(i,k,j)
              ta_rexpocaco3(i,k,jrow) = 
     &          ta_rexpocaco3(i,k,jrow) + rexpocaco3(i,k,j)
#   endif
#   if defined O_kk_si
          ta_rexpoopl(i,k,jrow) = ta_rexpoopl(i,k,jrow) +expo_op*dztr(k) !umol/cm3/s
          ta_rdisopl(i,k,jrow) = ta_rdisopl(i,k,jrow) + op_diss*dztr(k)
!          ta_rbur_sil(i,k,jrow) = ta_rbur_sil(i,k,jrow) + bsi*dztr(k) !umol/cm3/s
#   endif
#   if defined O_npzd_nitrogen && defined O_npzd_o2
                ta_rdeni(i,k,jrow) = ta_rdeni(i,k,jrow) + rdeni(i,k,j)
#   endif
#   if defined O_npzd_iron
                ta_rremife(i,k,jrow) = ta_rremife(i,k,jrow) 
     &               + rremife(i,k,j)
                ta_rexpofe(i,k,jrow) = ta_rexpofe(i,k,jrow)
     &               + rexpofe(i,k,j)
#   endif
              enddo
            endif
#  endif

          endif

        enddo

      enddo
# if defined O_embm      
#  if defined O_kk_si 
#    if !defined O_TMM && defined O_kk_si_compensating_sources
      rivsil = (silwflx-globalsilwflx)/globaldisch   !scale Si over total discharge
#    else
      rivsil = globalsilwflx/globaldisch
#    endif
      do j=js,je
        jrow = j + joff
        do i=is,ie
          if (kmt(i,jrow) .gt. 0) then      
          src(i,1,j,issil) = src(i,1,j,issil)
     &                        + rivsil*disch(i,j)*dztr(1)
#    if defined O_isotopes
              b2si = rs30std - eps_riv*rs30std/1000.
              fcsi = b2si/(1+b2si)
          src(i,1,j,issil30) = src(i,1,j,issil30)
     &                        + rivsil*disch(i,j)*dztr(1)*fcsi
#    endif
           if (timavgperts .and. .not. euler2) then
              ta_rrivsil(i,jrow) = ta_rrivsil(i,jrow)
     &                           + rivsil*disch(i,j)
           endif
          endif
       enddo
      enddo 
#  endif
# endif ! O_embm
# endif  !NPZD
# if defined O_carbon_fnpzd

!-----------------------------------------------------------------------
!    option to write and read fixed biological fluxes for carbon and
!    alkalinity. code assumes that restarts line up with the first
!    record. the repeating time period and restarts should probably
!    line up with the start of the year.
!-----------------------------------------------------------------------

      do j=js,je
        jrow = j + joff

!       adjust counters and read fluxes, if first row of slab
        if (jrow .eq. 2) then
          mfnpzd = mfnpzd + 1
          if (mfnpzd .gt. mxfnpzd) mfnpzd = 1
          if (nfnpzd .le. mxfnpzd) nfnpzd = nfnpzd + 1
!         read fluxes, if enough records written
          if (nfnpzd .gt. mxfnpzd) then
            fname = new_file_name ("O_fnpzd.nc")
            inquire (file=trim(fname), exist=exists)
            if (exists) then
              ib(:) = 1
              ic(:) = 1
              ic(1) = imtm2
              ic(2) = jmtm2
              ic(3) = km
              ib(4) = mfnpzd
              print*, "Reading fnpzd from: ", trim(fname),
     &          " record:", mfnpzd
              call openfile (fname, iou)
#  if defined O_carbon
              call getvara ('O_fnpzd1', iou, imtm2*jmtm2*km, ib, ic
     &,         tmpijk, c1, c0)
              fnpzd(2:imtm1,2:jmtm1,1:km,1)=tmpijk(1:imtm2,1:jmtm2,1:km)
#  endif
#  if defined O_npzd_alk
              call getvara ('O_fnpzd2', iou, imtm2*jmtm2*km, ib, ic
     &,         tmpijk, c1, c0)
              fnpzd(2:imtm1,2:jmtm1,1:km,2)=tmpijk(1:imtm2,1:jmtm2,1:km)
#  endif
            endif
          endif
        endif

        if (nfnpzd .le. mxfnpzd) then

!         store fluxes, if enough records written
          do i=is,ie
            do k=1,km
#  if defined O_carbon
              fnpzd(i,jrow,k,1) = src(i,k,j,isdic)
#  endif
#  if defined O_npzd_alk
              fnpzd(i,jrow,k,2) = src(i,k,j,isalk)
#  endif
            enddo
          enddo
!         write fluxes, if last row of slab
          if (jrow .eq. jmt-1) then
            fname = new_file_name ("O_fnpzd.nc")
            call openfile (fname, iou)
            print*, "Writing fnpzd to: ", trim(fname),
     &          " record:", mfnpzd
            ib(:) = 1
            ic(:) = 1
            ic(1) = imtm2
            ic(2) = jmtm2
            ic(3) = km
            ib(4) = mfnpzd
            time = relyr
            call putvars ('time', iou, mfnpzd, time, c1, c0)
#  if defined O_carbon
            tmpijk(1:imtm2,1:jmtm2,1:km) = fnpzd(2:imtm1,2:jmtm1,1:km,1)
            call putvara ('O_fnpzd1', iou, imtm2*jmtm2*km, ib, ic
     &,       tmpijk, c1, c0)
#  endif
#  if defined O_npzd_alk
            tmpijk(1:imtm2,1:jmtm2,1:km) = fnpzd(2:imtm1,2:jmtm1,1:km,2)
            call putvara ('O_fnpzd2', iou, imtm2*jmtm2*km, ib, ic
     &,       tmpijk, c1, c0)
#  endif
          endif

        else

!         set fluxes to source terms
          do i=is,ie
            do k=1,km
#  if defined O_carbon
              src(i,k,j,isdic) = fnpzd(i,jrow,k,1)
#  endif
#  if defined O_npzd_alk
              src(i,k,j,isalk) = fnpzd(i,jrow,k,2)
#  endif
            enddo
          enddo

        endif

      enddo
# endif
# if defined O_carbon && defined O_carbon_14

!-----------------------------------------------------------------------
!     set source for c14
!-----------------------------------------------------------------------
      do j=js,je
        jrow = j + joff
        do i=is,ie
          if (kmt(i,jrow) .gt. 0) then
            do k=1,kmt(i,jrow)
#  if defined O_npzd
              src(i,k,j,isc14) = src(i,k,j,isdic)*rstd
     &                         - 3.836e-12*t(i,k,j,ic14,taum1)
#  else
              src(i,k,j,isc14) = - 3.836e-12*t(i,k,j,ic14,taum1)
#  endif
            enddo
          endif
        enddo
      enddo
# endif
#ifndef O_TMM
# if defined O_sed && !defined O_sed_uncoupled

!-----------------------------------------------------------------------
!     set source terms from sediment model
!-----------------------------------------------------------------------
      do j=js,je
        jrow = j + joff
        do i=is,ie
          if (kmt(i,jrow) .gt. 0) then
            k = kmt(i,jrow)
#  if defined O_carbon && defined O_npzd
            src(i,k,j,isdic) = src(i,k,j,isdic) + sbc(i,j,ibdicfx)
#   if defined O_global_sums
            if (addflxo .and. eots) dtoic = dtoic - sbc(i,j,ibdicfx)
     &        *c2dtts*dtxcel(k)*dxt(i)*dyt(jrow)*cst(jrow)*dzt(k)
#   endif
#  endif
#  if defined O_npzd_alk
            src(i,k,j,isalk) = src(i,k,j,isalk) + sbc(i,j,ibalkfx)
#  endif
          endif
        enddo
      enddo
# endif
#endif ! NOT O_TMM

# ifndef O_TMM
!-----------------------------------------------------------------------
!     solve for one tracer at a time
!     n = 1 => temperature
!     n = 2 => salinity
!     n > 2 => other tracers (if applicable)
!-----------------------------------------------------------------------

      do n=1,nt
#  if defined O_matrix
        isTile=.false.
        if ((n.ge.ntilestart).and.(n.le.ntileend)) then
          isTile=.true.
		  ntile=n-ntilestart+1  ! index for tiles
!         reset tracer to tiles/initial condition
		  t(:,:,:,n,tau)=tile(:,:,:,ntile)
		  t(:,:,:,n,taum1)=tile(:,:,:,ntile)         
        endif
#  endif

!-----------------------------------------------------------------------
!       calculate advective tracer flux
!-----------------------------------------------------------------------

        call adv_flux (joff, js, je, is, ie, n)

!-----------------------------------------------------------------------
!       calculate diffusive flux across eastern and northern faces
!       of "T" cells due to various parameterizations for diffusion.
!-----------------------------------------------------------------------

#  if defined O_consthmix
#   if !defined O_biharmonic || defined O_bryan_lewis_horizontal

!       diffusive flux on eastern face of "T" cells

        do j=js,je
          do k=1,km
            do i=istrt-1,iend
              diff_fe(i,k,j) =
#    if defined O_bryan_lewis_horizontal
     &                         diff_cet(k)*cstdxur(i,j)*
#    else
     &                         ah_cstdxur(i,j)*
#    endif
     &                         (t(i+1,k,j,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
#    if defined O_isopycmix

!       diffusive flux on northern face of "T" cells
!       (background for isopycnal mixing)

        do j=js-1,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              diff_fn(i,k,j) =
#     if defined O_bryan_lewis_horizontal
     &                         diff_cnt(k)*
#     else
     &                         diff_cnt*
#     endif
     &           csu_dyur(jrow)*(t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
#    else

!       diffusive flux on northern face of "T" cells is not calculated.
!       instead, it is incorporated into DIFF_Ty for performance issues.

#    endif
#   else

!-----------------------------------------------------------------------
!       diffusive flux across eastern face of "T" cell
!       diffusive flux across northern face of "T" cell
!-----------------------------------------------------------------------

        m = 2 + n
        do j=js,je
          jrow = j + joff
          ahbi_cstr = diff_cet*cstr(jrow)
          do k=1,km
            do i=istrt-1,iend
              diff_fe(i,k,j) = ahbi_cstr*dxur(i)*
     &                         (del2(i+1,k,j,m)-del2(i,k,j,m))
            enddo
          enddo
        enddo
        do j=js-1,je
          jrow = j + joff
          ahbi_csu_dyur = diff_cnt*csu(jrow)*dyur(jrow)
          do k=1,km
            do i=istrt,iend
              diff_fn(i,k,j) = ahbi_csu_dyur*
     &                        (del2(i,k,j+1,m) - del2(i,k,j,m))
            enddo
          enddo
        enddo
#   endif
#  else
#   if defined O_smagnlmix

!       diffusive flux on eastern and northern faces of "T" cells

        do j=js,je
          do k=1,km
            do i=istrt-1,iend
              diff_fe(i,k,j) = diff_cet(i,k,j)*cstdxur(i,j)*
     &                         (t(i+1,k,j,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
        do j=js-1,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              diff_fn(i,k,j) = diff_cnt(i,k,j)*csu_dyur(jrow)*
     &                     (t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
#   endif
#  endif
!-----------------------------------------------------------------------
!       calculate diffusive flux across bottom face of "T" cells
!-----------------------------------------------------------------------

        do j=js,je
          do k=1,km-1
            do i=istrt,iend
              diff_fb(i,k,j) = diff_cbt(i,k,j)*dzwr(k)*
     &                         (t(i,k,j,n,taum1) - t(i,k+1,j,n,taum1))
            enddo
          enddo
        enddo

#  if defined O_isopycmix

!-----------------------------------------------------------------------
!       compute isopycnal diffusive flux through east, north,
!       and bottom faces of T cells.
!-----------------------------------------------------------------------

        call isoflux (joff, js, je, is, ie, n)
#  endif

!-----------------------------------------------------------------------
!       set surface and bottom vert b.c. on "T" cells for diffusion
!       and advection. for isopycnal diffusion, set adiabatic boundary
!       conditions.
!       note: the b.c. at adv_fb(i,k=bottom,j) is set by the above code.
!             However, it is not set when k=km so it is set below.
!             adv_fb(i,km,j) is always zero (to within roundoff).
!-----------------------------------------------------------------------

        do j=js,je
          jrow   = j + joff
          do i=istrt,iend
            kb              = kmt(i,jrow)
#  if defined O_replacst
            diff_fb(i,0,j)  = c0
#  else
            diff_fb(i,0,j)  = stf(i,j,n)
#  endif
            diff_fb(i,kb,j) = btf(i,j,n)
            adv_fb(i,0,j)   = adv_vbt(i,0,j)*(t(i,1,j,n,tau) +
     &                                        t(i,1,j,n,tau))
            adv_fb(i,km,j)  = adv_vbt(i,km,j)*t(i,km,j,n,tau)
          enddo
        enddo

#  if defined O_source_term || defined O_npzd || defined O_carbon_14

!-----------------------------------------------------------------------
!       set source term for "T" cells
!-----------------------------------------------------------------------

        source(:,:,:) = c0

#   if defined O_npzd || defined O_carbon_14
        if (itrc(n) .ne. 0) then
          do j=js,je
            do k=1,km
              do i=istrt,iend
                source(i,k,j) = src(i,k,j,itrc(n))
              enddo
            enddo
          enddo
        endif
#   endif
#   if defined O_shortwave

!-----------------------------------------------------------------------
!       incorporate short wave penetration into source
!-----------------------------------------------------------------------

        if (n .eq. 1) then
          call swflux0 (joff, js, je, istrt, iend, source)
        endif
#   endif
#  endif

!-----------------------------------------------------------------------
!       solve for "tau+1" tracer using statement functions to represent
!       each component of the calculation
!-----------------------------------------------------------------------

!       1st: solve using all components which are treated explicitly

        do j=js,je
          jrow   = j + joff
          do k=1,km
            twodt(k) = c2dtts*dtxcel(k)
            do i=istrt,iend
              t(i,k,j,n,taup1) = t(i,k,j,n,taum1) + twodt(k)*(
     &          DIFF_Tx(i,k,j) + DIFF_Ty(i,k,j,jrow,n) + DIFF_Tz(i,k,j)
     &          - ADV_Tx(i,k,j) -  ADV_Ty(i,k,j,jrow,n) -  ADV_Tz(i,k,j)
#  if defined O_isopycmix && defined O_gent_mcwilliams && !defined O_fct && !defined O_quicker
     &          - ADV_Txiso(i,k,j,n) - ADV_Tyiso(i,k,j,jrow,n)
     &          - ADV_Tziso(i,k,j)
#  endif
#  if defined O_source_term || defined O_npzd || defined O_carbon_14
     &          + source(i,k,j)
#  endif
#  if defined O_plume
     &          + work(i,k,j)*subflux(i,jrow,n)
#  endif
     &          )*tmask(i,k,j)
            enddo
          enddo
        enddo

#  if defined O_matrix
!     accumulate explicit matrix and reset to tiles/initial condition
	  if (isTile) then
	    call MATRIX_STORE_EXPLICIT(ntile,t(1,1,1,n,taup1),twodt)
      endif      
#  endif        
#  if defined O_implicitvmix || defined O_isopycmix || defined O_redi_diffusion

!       2nd: add in portion of vertical diffusion handled implicitly
#   if defined O_matrix
		if (isTile) then
          call ivdift (joff, js, je, istrt, iend, n, dtmatrix)
        else
          call ivdift (joff, js, je, istrt, iend, n, twodt)        
		endif          
#   else
        call ivdift (joff, js, je, istrt, iend, n, twodt)
#   endif

#  endif

        do j=js,je
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo

!-----------------------------------------------------------------------
!       construct diagnostics associated with tracer "n"
!-----------------------------------------------------------------------

        call diagt1 (joff, js, je, istrt, iend, n, twodt)

!-----------------------------------------------------------------------
!       end of tracer component "n" loop
!-----------------------------------------------------------------------

      enddo

#  if defined O_replacst

      gamma = zw(1)/twodt(1)
      do j=js,je
        jrow = j + joff
        do i=istrt,iend
          stf(i,j,1) = gamma*(sbc(i,jrow,isst) - t(i,1,j,1,taup1))
          stf(i,j,2) = gamma*(sbc(i,jrow,isss) - t(i,1,j,2,taup1))
          t(i,1,j,1,taup1) = sbc(i,jrow,isst)
          t(i,1,j,2,taup1) = sbc(i,jrow,isss)
        enddo
      enddo
#  endif
#  if defined O_convect_brine

!-----------------------------------------------------------------------
!     explicit convection: adjust column if gravitationally unstable
!     convect brine rejected under all ice categories
!-----------------------------------------------------------------------

      call convect_brine (joff, js, je, is, ie)
#  else
#   if !defined O_implicitvmix || defined O_isopycmix

!-----------------------------------------------------------------------
!     explicit convection: adjust column if gravitationally unstable
!-----------------------------------------------------------------------

#    if defined O_fullconvect
      call convct2 (t(1,1,1,1,taup1), joff, js, je, is, ie, kmt)
      do j=js,je
        do n=1,nt
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo
      enddo
#    else
      call convct (t(1,1,1,1,taup1), ncon, joff, js, je, istrt, iend
     &,            kmt)
#    endif
#   endif
#  endif
#  if defined O_save_convection
      if (timavgperts .and. eots) then
        if (joff .eq. 0) nta_conv = nta_conv + 1
        do j=js,je
          jrow = j + joff
          do i=istrt,iend
            ta_totalk(i,jrow) = ta_totalk(i,jrow) + totalk(i,j)
            ta_vdepth(i,jrow) = ta_vdepth(i,jrow) + vdepth(i,j)
            ta_pe(i,jrow) = ta_pe(i,jrow) + pe(i,j)
          enddo
        enddo
      endif
#  endif

!-----------------------------------------------------------------------
!     construct diagnostics after convection
!-----------------------------------------------------------------------

      idiag = 10
      call diagt2 (joff, js, je, istrt, iend, idiag)

#  if defined O_matrix
!     accumulate implicit matrix
      do n=ntilestart,ntileend
        ntile=n-ntilestart+1
        call MATRIX_STORE_IMPLICIT(ntile,t(1,1,1,n,taup1))
      enddo
#  endif

#  if defined O_fourfil || defined O_firfil

!-----------------------------------------------------------------------
!     filter tracers at high latitudes
!-----------------------------------------------------------------------

      if (istrt .eq. 2 .and. iend .eq. imt-1) then
        call filt (joff, js, je)
      else
        write (stdout,'(a)')
     &  'Error: filtering requires is=2 and ie=imt-1 in tracer'
        stop '=>tracer'
      endif
#  endif
      do n=1,nt
        do j=js,je
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo
      enddo

!-----------------------------------------------------------------------
!     construct diagnostics after filtering (for total dT/dt)
!-----------------------------------------------------------------------

      idiag = 1
      call diagt2 (joff, js, je, istrt, iend, idiag)

#  if defined O_idealage
! ----------------------------------------------------------------------
! restore t(iidealage) to 0.0 in the surface layer
! ----------------------------------------------------------------------
       do j=js,je
          jrow = j + joff
          do i=istrt,iend
             t(i,1,j,iidealage,taup1)= 0.0
          enddo
       enddo
#  endif
#  if defined O_plastic_conc
! ----------------------------------------------------------------------
! restore t(imp) to plastic mask in the surface layer
! ----------------------------------------------------------------------
       do j=js,je
          jrow = j + joff
          do i=istrt,iend
             t(i,1,j,imp,taup1)= max(trcmin,plastic(i,j,2)*1.e-2)
          enddo
       enddo
#  endif
#  if !defined O_replacst

!-----------------------------------------------------------------------
!     if needed, construct the Atmos S.B.C.(surface boundary conditions)
!     averaged over this segment
!     eg: SST and possibly SSS
!-----------------------------------------------------------------------

      call asbct (joff, js, je, istrt, iend, isst, itemp)
      call asbct (joff, js, je, istrt, iend, isss, isalt)
#  endif
#  if defined O_carbon
      call asbct (joff, js, je, istrt, iend, issdic, idic)
#   if defined O_carbon_14
      call asbct (joff, js, je, istrt, iend, issc14, ic14)
#   endif
#  endif
#  if defined O_plastic
      call asbct (joff, js, je, istrt, iend, issmp, imp)
      call asbct (joff, js, je, istrt, iend, issmpa, impa)
      call asbct (joff, js, je, istrt, iend, issmpp, impp)
#  endif
#  if defined O_npzd_alk
      call asbct (joff, js, je, istrt, iend, issalk, ialk)
#  endif
#  if defined O_npzd_o2
      call asbct (joff, js, je, istrt, iend, isso2, io2)
#  endif
# if defined O_kk_si
      call asbct (joff, js, je, istrt, iend, isssil, isil)
# endif
#  if defined O_npzd
      call asbct (joff, js, je, istrt, iend, isspo4, ipo4)
#   if !defined O_npzd_no_vflux
      call asbct (joff, js, je, istrt, iend, issphyt, iphyt)
      call asbct (joff, js, je, istrt, iend, isszoop, izoop)
      call asbct (joff, js, je, istrt, iend, issdetr, idetr)
#   if defined O_zoop_det
      call asbct (joff, js, je, istrt, iend, issdetrz, idetrz)
#   endif
#   if defined O_kk_ballast
      call asbct (joff, js, je, istrt, iend, issdetr_B, idetr_B)
#   endif
#   if defined O_cal_kk
      call asbct (joff, js, je, istrt, iend, isscocc, icocc)
#   endif
#   if defined O_kk_diat
      call asbct (joff, js, je, istrt, iend, issdiat, idiat)
#   endif
#   if defined O_kk_caco3tr
      call asbct (joff, js, je, istrt, iend, isscaco3, icaco3)
#   endif
#   endif
#   if defined O_npzd_nitrogen
      call asbct (joff, js, je, istrt, iend, issno3, ino3)
#    if !defined O_npzd_no_vflux && !defined O_npzd_simple
      call asbct (joff, js, je, istrt, iend, issdiaz, idiaz)
#    endif
#   endif
#   if defined O_npzd_iron
      call asbct (joff, js, je, istrt, iend, issdfe, idfe)
#    if !defined O_npzd_no_vflux
      call asbct (joff, js, je, istrt, iend, issdetrfe, idetrfe)
#    endif
#   endif
#  endif
#  if defined O_isotopes
      call asbct (joff, js, je, istrt, iend, issdic13, idic13)
      call asbct (joff, js, je, istrt, iend, issdin15, idin15)
      call asbct (joff, js, je, istrt, iend, isssil30, isil30)
#    if !defined O_npzd_no_vflux
      call asbct (joff, js, je, istrt, iend, isscaco3c13, icaco3c13)
      call asbct (joff, js, je, istrt, iend, issphytc13, iphytc13)
      call asbct (joff, js, je, istrt, iend, issdiazc13, idiazc13)
      call asbct (joff, js, je, istrt, iend, isszoopc13, izoopc13)
      call asbct (joff, js, je, istrt, iend, isscoccc13, icoccc13)
      call asbct (joff, js, je, istrt, iend, issdiatc13, idiatc13)
      call asbct (joff, js, je, istrt, iend, issdetrc13, idetrc13)

      call asbct (joff, js, je, istrt, iend, issphytn15, iphytn15)
      call asbct (joff, js, je, istrt, iend, issdiazn15, idiazn15)
      call asbct (joff, js, je, istrt, iend, isszoopn15, izoopn15)
      call asbct (joff, js, je, istrt, iend, isscoccn15, icoccn15)
      call asbct (joff, js, je, istrt, iend, issdiatn15, idiatn15)
      call asbct (joff, js, je, istrt, iend, issdetrn15, idetrn15)
      call asbct (joff, js, je, istrt, iend, issdiats30, idiats30)
#    endif
#  endif
#  if defined O_cfcs_data || defined O_cfcs_data_transient
      call asbct (joff, js, je, istrt, iend, isscfc11, icfc11)
      call asbct (joff, js, je, istrt, iend, isscfc12, icfc12)
#  endif

#endif ! NOT O_TMM

#  if defined O_sed
!-----------------------------------------------------------------------
!     accumulate bottom tracer values for the sediment model
!-----------------------------------------------------------------------

      if (addflxo .and. eots) then
        if (joff .eq. 0) atsed = atsed + twodt(1)
        do j=js,je
          jrow = j + joff
          do i=is,ie
            if (kmt(i,jrow) .gt. 0) then
              k = kmt(i,jrow)
              sbc(i,jrow,ibtemp) = sbc(i,jrow,ibtemp)
     &                           + t(i,k,j,itemp,taup1)*twodt(1)
              sbc(i,jrow,ibsalt) = sbc(i,jrow,ibsalt)
     &                           + t(i,k,j,isalt,taup1)*twodt(1)
#   if defined O_carbon
              sbc(i,jrow,ibdic) = sbc(i,jrow,ibdic)
     &                          + t(i,k,j,idic,taup1)*twodt(1)
#   endif
#   if defined O_npzd_alk
              sbc(i,jrow,ibalk) = sbc(i,jrow,ibalk)
     &                          + t(i,k,j,ialk,taup1)*twodt(1)
#   endif
#   if defined O_npzd_o2
              sbc(i,jrow,ibo2) = sbc(i,jrow,ibo2)
     &                         + t(i,k,j,io2,taup1)*twodt(1)
#   endif
            endif
          enddo
        enddo
      endif
#  endif

#ifndef O_TMM

#  if defined O_carbon && defined O_carbon_14

!-----------------------------------------------------------------------
!     calculate diagnostic delta carbon 14
!-----------------------------------------------------------------------

      if (tsiperts .or. timavgperts) then
        rrstd = 1000./rstd
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              dc14(i,k,j) = (rrstd*t(i,k,j,ic14,taup1)
     &                      /(t(i,k,j,idic,taup1) + epsln) - 1000.)
     &                      *tmask(i,k,j)
            enddo
          enddo
        enddo
      endif
      if (tsiperts .and. eots) then
        if (js+joff .eq. 2) dc14bar = 0.
        do j=js,je
          jrow = j + joff
          fy = cst(jrow)*dyt(jrow)
          do k=1,km
            fyz = fy*dzt(k)
            do i=istrt,iend
              dc14bar = dc14bar + dc14(i,k,j)*dxt(i)*fyz*tmask(i,k,j)
            enddo
          enddo
        enddo
      endif
      if (timavgperts .and. .not. euler2) then
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              ta_dc14(i,k,jrow) = ta_dc14(i,k,jrow) + dc14(i,k,j)
            enddo
          enddo
        enddo
      endif
#  endif

# endif ! not O_TMM

      return
      end

# ifndef O_TMM

      subroutine diagt1 (joff, js, je, is, ie, n, twodt)

!-----------------------------------------------------------------------
!     construct diagnostics associated with tracer component "n"

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       n     = (1,2) = (u,v) velocity component
!       twodt = (2*dtts,dtts) on (leapfrog,mixing) time steps
!-----------------------------------------------------------------------

      implicit none

      integer i, k, j, ip, kr, jq, n, jp, jrow, js, je, joff, is, ie
      integer mask, m

      real t_i, t_j, dz_t2r, dz_tr, dz_wtr, dx_t2r, dx_tr, dy_t2r
      real dy_tr, adv_tx, adv_ty, adv_tz, adv_txiso, adv_tyiso
      real adv_tziso, diff_tx, diff_ty, diff_tz, dtdx, dtdy, dtdz
      real r2dt, cosdyt, fx, darea, boxar, rtwodt, sumdx, delx
      real sumdxr, dxdy, dxdydz

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "accel.h"
#include "coord.h"
#include "cregin.h"
#include "csbc.h"
#  if defined O_tracer_averages
#include "ctavg.h"
#  endif
# if defined O_kk_caco3tr
#include "npzd.h"
# endif
#include "diag.h"
#include "diaga.h"
#include "emode.h"
#include "grdvar.h"
#include "hmixc.h"
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "vmixc.h"

#  if defined O_meridional_tracer_budget
#include "ctmb.h"
      real stor(imt,km), div(imt,km), sorc(imt,km), dif(imt,km)
      real t1(imt), t2(imt), t3(imt), t4(imt)
#  endif
#  if defined O_time_step_monitor
      real temp1(imt,km), temp2(imt,km), temp3(imt,km)
#  endif
      real twodt(km)

#  if defined O_isopycmix
#include "isopyc.h"
#  endif
#include "fdift.h"

#  if defined O_save_mixing_coeff

!-----------------------------------------------------------------------
!     diagnostic: estimate mixing coefficients on east, north, and
!                 bottom face of T cells from the flux
!-----------------------------------------------------------------------

      if (cmixts .and. n .eq. 1 .and. eots) then
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=2,imt-1
              dtdx = (t(i+1,k,j,1,taum1)-t(i,k,j,1,taum1))
     &                *tmask(i+1,k,j)*tmask(i,k,j)
     &                *cstr(jrow)*dxur(i) + epsln
              ce(i,k,j,2) = diff_fe(i,k,j)/dtdx
#   if !defined O_consthmix || defined O_biharmonic || defined O_isopycmix
              dtdy = (t(i,k,j+1,1,taum1)-t(i,k,j,1,taum1))
     &                *tmask(i,k,j+1)*tmask(i,k,j)
     &                *dyur(jrow) + epsln
              cn(i,k,j,2) = csur(jrow)*diff_fn(i,k,j)/dtdy
#   else
              cn(i,k,j,2) = csur(jrow)*ah
#   endif
            enddo
          enddo
        enddo
        do j=js,je
          jrow = j + joff
          do k=1,km-1
            do i=2,imt-1
              dtdz = (t(i,k,j,1,taum1)-t(i,k+1,j,1,taum1))
     &                *tmask(i,k,j)*tmask(i,k+1,j)
     &                *dzwr(k) + epsln
              cb(i,k,j,2) = diff_fb(i,k,j)/dtdz
#   if defined O_isopycmix
     &                      + diff_fbiso(i,k,j)/dtdz
#   endif
            enddo
          enddo
          do i=2,imt-1
            cb(i,km,j,2) = 0.0
          enddo
        enddo

        do j=js,je
          call setbcx (ce(1,1,j,2), imt, km)
          call setbcx (cn(1,1,j,2), imt, km)
          call setbcx (cb(1,1,j,2), imt, km)
        enddo
      endif
#  endif

#  if defined O_save_convection_full

!-----------------------------------------------------------------------
!     diagnostic: initialize temperature before convection
!-----------------------------------------------------------------------

      if (exconvts .and. n .eq. 1 .and. eots) then
        do j=js,je
          do k=1,km
            do i=1,imt
              excnv0(i,k,j) = t(i,k,j,1,taup1)
            enddo
          enddo
        enddo
      endif

#  endif

#  if defined O_time_step_monitor

!-----------------------------------------------------------------------
!     diagnostic: integrate |d(tracer)/dt|  and tracer variance on "tau"
!                 globally
!-----------------------------------------------------------------------

      if (tsiperts .and. eots) then
        do j=js,je
          jrow = j + joff
          r2dt    = c1/c2dtts
          cosdyt  = cst(jrow)*dyt(jrow)
          do k=1,km
            fx = r2dt/dtxcel(k)
            do i=is,ie
              darea      = dzt(k)*dxt(i)*cosdyt*tmask(i,k,j)
              temp3(i,k) = t(i,k,j,n,tau)*darea
              temp1(i,k) = t(i,k,j,n,tau)**2*darea
              temp2(i,k) = abs(t(i,k,j,n,taup1)-t(i,k,j,n,taum1))*
     &                     darea*fx
            enddo
            do i=is,ie
              tbar(k,n,jrow)   = tbar(k,n,jrow) + temp3(i,k)
              travar(k,n,jrow) = travar(k,n,jrow) + temp1(i,k)
              dtabs(k,n,jrow)  = dtabs(k,n,jrow) + temp2(i,k)
            enddo
          enddo
        enddo
      endif
#  endif

#  if defined O_tracer_averages

!-----------------------------------------------------------------------
!     diagnostic: accumulate tracers for averages under horizontal
!                 regions (use units of meters, rather than cm)
!-----------------------------------------------------------------------

      if (tavgts .and. eots) then
        do j=js,je
          jrow = j + joff
          do i=is,ie
            mask = mskhr(i,jrow)
            if (mask .ne. 0) then
              boxar = cst(jrow)*dxt(i)*dyt(jrow)*tmask(i,1,j)*0.0001
              sumbf(mask,n) = sumbf(mask,n) + stf(i,j,n)*boxar
              do k=1,km
                sumbk(mask,k,n) = sumbk(mask,k,n) + t(i,k,j,n,tau)
     &                             *boxar*dzt(k)*tmask(i,k,j)*0.01
             enddo
            endif
          enddo
        enddo
      endif
#  endif

#  if defined O_tracer_yz

!----------------------------------------------------------------------
!     diagnostic: integrate tracer equations in longitude
!----------------------------------------------------------------------

      if (tyzts .and. eots) then
        do j=js,je
          jrow = j + joff
          do k=1,km
            rtwodt = c1/twodt(k)
            sumdx = c0
            do m=1,5
              tyz(jrow,k,n,m) = c0
            enddo
            do i=is,ie
              delx  = dxt(i)*tmask(i,k,j)
              sumdx = sumdx + delx
              tyz(jrow,k,n,1) = tyz(jrow,k,n,1)+ delx*t(i,k,j,n,tau)
              tyz(jrow,k,n,2) = tyz(jrow,k,n,2) + rtwodt*delx*
     &          (t(i,k,j,n,taup1) - t(i,k,j,n,taum1))
              tyz(jrow,k,n,3) =  tyz(jrow,k,n,3)
     &          - delx*(ADV_Ty(i,k,j,jrow,n) + ADV_Tz(i,k,j))
              tyz(jrow,k,n,4) =  tyz(jrow,k,n,4)
     &          + delx*(DIFF_Ty(i,k,j,jrow,n) + DIFF_Tz(i,k,j))
#   if defined O_source_term || defined O_npzd || defined O_carbon_14
              tyz(jrow,k,n,5) = tyz(jrow,k,n,5) + delx*source(i,k,j)
#   endif
            enddo
            sumdxr = c1/(sumdx + epsln)
            do m=1,5
              tyz(jrow,k,n,m) = tyz(jrow,k,n,m)*sumdxr
            enddo
          enddo
        enddo
        if (js+joff .eq. 2) then
          do m=1,5
            do k=1,km
              tyz(1,k,n,m)   = c0
              tyz(jmt,k,n,m) = c0
            enddo
          enddo
        endif
      endif
#  endif

#  if defined O_meridional_tracer_budget

!----------------------------------------------------------------------
!     diagnostic: integrate equations in depth, lon, and time for
!                 each basin and jrow. basin #  0 is used to catch
!                 values in land areas
!----------------------------------------------------------------------

      if (tmbperts .and. eots) then
        do j=js,je
          jrow = j + joff
          if (jrow .eq. 2 .and. n .eq. 1) numtmb = numtmb + 1
          cosdyt = cst(jrow)*dyt(jrow)
          do k=1,km
            do i=is,ie
              stor(i,k) = 0.0
              div(i,k)  = 0.0
              sorc(i,k) = 0.0
              dif(i,k)  = 0.0
            enddo
          enddo
          do k=1,km
            rtwodt = c1/twodt(k)
            do i=is,ie
              dxdy            = cosdyt*dxt(i)*tmask(i,k,j)
              dxdydz          = dxdy*dzt(k)
              stor(i,k) = stor(i,k) + rtwodt*dxdydz*
     &                            (t(i,k,j,n,taup1) - t(i,k,j,n,taum1))
              div(i,k)  = div(i,k)  - dxdydz*ADV_Ty(i,k,j,jrow,n)
#   if defined O_source_term || defined O_npzd || defined O_carbon_14
              sorc(i,k) = sorc(i,k) + dxdydz*source(i,k,j)
#   endif
              dif(i,k)  = dif(i,k)  + dxdydz*DIFF_Ty(i,k,j,jrow,n)
            enddo
          enddo

!         the accumulation is done this way for issues of speed

          do i=is,ie
            t1(i) = stor(i,1)
            t2(i) = div(i,1)
            t3(i) = sorc(i,1)
            t4(i) = dif(i,1)
          enddo

          do k=2,km
            do i=is,ie
              t1(i) = t1(i) + stor(i,k)
              t2(i) = t2(i)  + div(i,k)
              t3(i) = t3(i) + sorc(i,k)
              t4(i) = t4(i)  + dif(i,k)
            enddo
          enddo

          do i=is,ie
            m               = msktmb(i,jrow)
            tstor(jrow,n,m) = tstor(jrow,n,m) + t1(i)
            tdiv(jrow,n,m)  = tdiv(jrow,n,m)  + t2(i)
            tsorc(jrow,n,m) = tsorc(jrow,n,m) + t3(i)
            tdif(jrow,n,m)  = tdif(jrow,n,m)  + t4(i)
          enddo

          k = 1
          do i=is,ie
            m               = msktmb(i,jrow)
            dxdy            = cosdyt*dxt(i)*tmask(i,k,j)
            tflux(jrow,n,m) = tflux(jrow,n,m) + dxdy*stf(i,j,n)
          enddo
          if (n .eq. 1) then
            do k=1,km
              do i=is,ie
                m               = msktmb(i,jrow)
                dxdy            = cosdyt*dxt(i)*tmask(i,k,j)
                dxdydz          = dxdy*dzt(k)
                smdvol(jrow,m)  = smdvol(jrow,m)  + dxdydz
              enddo
            enddo
          endif
        enddo
      endif
#  endif

#  if defined O_gyre_components

!-----------------------------------------------------------------------
!     diagnostic: compute the northward transport components of
!                 each tracer
!-----------------------------------------------------------------------

      if (gyrets .and. eots)  call gyre (joff, js, je, is, ie, n)
#  endif

#  if defined O_term_balances

!-----------------------------------------------------------------------
!     diagnostic: integrate r.h.s. terms in the tracer equations
!                 over specified regional volumes.
!-----------------------------------------------------------------------

      if (trmbts .and. eots)  call ttb1 (joff, js, je, is, ie, n)
#  endif

#  if defined O_xbts

!-----------------------------------------------------------------------
!     diagnostic: accumulate r.h.s. terms in the tracer equations
!-----------------------------------------------------------------------

      if (xbtperts .and. eots) call txbt1 (joff, js, je, n)
#  endif
#  if defined O_mom_tbt

!-----------------------------------------------------------------------
!     diagnostic: accumulate r.h.s. terms in the tracer equations
!-----------------------------------------------------------------------

      if (tbtperts .and. eots) call tbt1 (joff, js, je, n)
#  endif

      return
      end

      subroutine diagt2 (joff, js, je, is, ie, idiag)

!-----------------------------------------------------------------------
!     construct d(tracer)/dt diagnostics

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       idiag = 1  => total tracer change
!       idiag = 10 => change of tracer due to filtering(also convection)
!-----------------------------------------------------------------------

      implicit none

      integer idiag, j, js, je, k, i, joff, iocv, jrow, is, ie

      real rdt, reltim, period

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "coord.h"
#include "diaga.h"
#include "iounit.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "tmngr.h"
#include "timeavgs.h"

#  if defined O_save_convection_full

!-----------------------------------------------------------------------
!     diagnostic: save tracer time change due to explicit convection
!     idiag = 10 signifies diagnostics immediately after convection
!-----------------------------------------------------------------------

      if (exconvts .and. idiag .eq. 10 .and. eots) then
        rdt = c1/c2dtts

!       excnv1 = epsln over land cells and d(convect)/dt over ocean

        do j=js,je
          do k=1,km
            do i=1,imt
              excnv1(i,k,j) = tmask(i,k,j)*
     &         (t(i,k,j,1,taup1)-excnv0(i,k,j))*rdt
     &                       + (c1-tmask(i,k,j))*epsln
            enddo
          enddo
        enddo

        reltim = relyr
        if (joff + js .eq. 2) then
          write (stdout,*) ' =>Writing explicit convection at ts=',itt
     &  , ' ',stamp
          call getunit (iocv, 'cvct.dta'
     &,                'unformatted sequential append ieee')

          period = 0.0
          iotext = 'read(iocv) reltim, period, imt, jmt, km, flag'
          write (iocv) stamp, iotext, expnam
          write (iocv) reltim, period, imt, jmt, km, epsln

          iotext = 'read(iocv) (xt(i),i=1,imt)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, xt, imt)

          iotext = 'read(iocv) (yt(j),j=1,jmt)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, yt, jmt)

          iotext = 'read(iocv) (zt(k),k=1,km)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, zt, km)

          call relunit (iocv)
        endif

        call getunit (iocv, 'cvct.dta'
     &,                'unformatted sequential append ieee')

        do j=js,je
          jrow = j+joff

          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocv) (convect(i,k),i=1,imt),k=1,km)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, excnv1(1,1,j), imt*km)

        enddo
        call relunit(iocv)
      endif
#  endif

#  if defined O_term_balances

!-----------------------------------------------------------------------
!     diagnostic: integrate d/dt(tracer) over specified regional volumes
!                  after convection and filtering
!-----------------------------------------------------------------------

      if (trmbts .and. eots) call ttb2 (joff, js, je, is, ie, idiag)
#  endif

#  if defined O_xbts

!-----------------------------------------------------------------------
!     diagnostic: accumulate d/dt(tracer) after convection
!                 and filtering
!-----------------------------------------------------------------------

      if (xbtperts .and. eots) call txbt2 (joff, js, je, idiag)
#  endif
#  if defined O_mom_tbt

!-----------------------------------------------------------------------
!     diagnostic: accumulate d/dt(tracer) after convection
!                 and filtering

!-----------------------------------------------------------------------

      if (tbtperts .and. eots) call tbt2 (joff, js, je, idiag)
#  endif

      return
      end

      subroutine asbct (joff, js, je, is, ie, isbc, itr)

!-----------------------------------------------------------------------
!     construct the Atmos S.B.C. (surface boundary conditions)

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       isbc  = index for sbc
!       itr   = index for tracer
!-----------------------------------------------------------------------

      implicit none

      integer isbc, itr, j, js, je, jrow, joff, i, is, ie

      real rts

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "csbc.h"
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"

!     initialize the Atmos S.B.C. at the start of each ocean segment
!     (do not alter values in land)

      if (isbc .le. 0 .or. itr .le. 0) return

      if (eots .and. osegs) then
        do j=js,je
          jrow   = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0) sbc(i,jrow,isbc) = c0
          enddo
        enddo
      endif

!     accumulate surface tracers for the Atmos S.B.C. every time step

      if (eots) then
        do j=js,je
          jrow = j + joff
          do i=is,ie
            sbc(i,jrow,isbc) = sbc(i,jrow,isbc)+t(i,1,j,itr,taup1)
          enddo
        enddo
      endif

!     average the surface tracers for the Atmos S.B.C. at the end of
!     each ocean segment. (do not alter values in land)

      if (eots .and. osege) then
        rts = c1/ntspos
        do j=js,je
          jrow   = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0)
     &        sbc(i,jrow,isbc) = rts*sbc(i,jrow,isbc)
          enddo
        enddo
      endif

      return
      end

      subroutine ivdift (joff, js, je, is, ie, n, twodt)

#  if defined O_implicitvmix || defined O_isopycmix || defined O_redi_diffusion
!-----------------------------------------------------------------------
!     solve vertical diffusion of tracers implicitly

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       n     = tracer component
!       twodt = (2*dtts, dtts) on (leapfrog, mixing) time steps
!-----------------------------------------------------------------------

      implicit none

      integer j, js, je, k, i, is, ie, n, joff

      real rc2dt

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "levind.h"
#include "mw.h"
#include "switch.h"
#include "vmixc.h"

      real twodt(km)

!     store terms to compute implicit vertical mixing on
!     diagnostic time steps

#   if defined O_xbts || defined O_mom_tbt
      if (eots) then
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = t(i,k,j,n,taup1)
            enddo
          enddo
        enddo
      endif
#   else
#    if defined O_term_balances
      if (trmbts .and. eots) then
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = t(i,k,j,n,taup1)
            enddo
          enddo
        enddo
      endif
#    endif
#   endif

      call invtri (t(1,1,1,n,taup1), stf(1,1,n), btf(1,1,n)
     &, diff_cbt(1,1,jsmw), twodt, kmt, tmask(1,1,1), is, ie
     &, joff, js, je)

!     compute residual implicit vertical mixing

#   if defined O_xbts || defined O_mom_tbt
      if (eots) then
        do j=js,je
          do k=1,km
            rc2dt = c1/twodt(k)
            do i=is,ie
              zzi(i,k,j) = rc2dt*(t(i,k,j,n,taup1) - zzi(i,k,j))
            enddo
          enddo
        enddo
      endif
#   else
#    if defined O_term_balances
      if (trmbts .and. eots) then
        do j=js,je
          do k=1,km
            rc2dt = c1/twodt(k)
            do i=is,ie
              zzi(i,k,j) = rc2dt*(t(i,k,j,n,taup1) - zzi(i,k,j))
            enddo
          enddo
        enddo
      endif
#    endif
#   endif
#  endif

      return
      end

      subroutine swflux0 (joff, js, je, is, ie, source)

#  if defined O_mom && defined O_shortwave
!=======================================================================
!       incorporate short wave penetration via the "source"
!       term. note that the divergence of shortwave for the first
!       level "divpen(1)" is compensating for the effect of having
!       the shortwave component already included in the total
!       surface tracer flux "stf(i,j,1)"

!       incorporating the penetrative shortwave radiative flux into
!       the vertical diffusive flux term directly is not correct when
!       vertical diffusion is solved implicitly.

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW

!     output:
!       source = source term penetrative short wave
!=======================================================================

      implicit none

      integer j, js, je, jrow, joff, k, i, is, ie

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "csbc.h"
#include "cshort.h"

      real source(imt,km,jsmw:jemw)

      if (iswr .ne. 0) then
        do j=js,je
          jrow   = j + joff
          do k=1,km
            do i=is,ie
              source(i,k,j) = source(i,k,j)
     &                       + sbc(i,jrow,ipsw)*divpen(k)
            enddo
          enddo
        enddo
      endif
#  endif
      return
      end
!
# endif !not O_TMM

# if defined O_npzd_fe_limitation
      subroutine polint (xa,ya,n,x,y,dy)
!
      implicit none
!
      integer n,i,m,ns, nmax
      parameter (nmax = 10) ! largest anticipated value of n
      real dy,x,y,xa(n),ya(n)
      real den,dif,dift,ho,hp,w,c(nmax),d(nmax)
!
! Given arrays xa and ya, each of length n, and a given value x, this routine
! returns a value y, and an error estimate dy. If P(x) is the polynomial of
! degree N-1 such that P(xai) = yai, i=1,....,n then the returned value y=P(x)
!
!
      ns=1
      dif=abs(x-xa(1))
!
      do i=1,n
         dift = abs(x-xa(i))
         IF (dift.lt.dif) THEN
            ns = i
            dif = dift
         END IF
!
       c(i) = ya(i)
       d(i) = ya(i)
!
      enddo
!
      y = ya (ns)        ! this is the initial approximation to y
      ns = ns-1
!
      do m=1,n-1
         do i=1,n-m
            ho = xa(i) - x
            hp = xa (i+m) -x
            w = c(i+1) - d(i)
            den = ho - hp
            IF (den.eq.0.) PAUSE 'failure in polint'
            den = w/den
            d(i) = hp*den
            c(i) = ho*den
         enddo
         IF (2*ns.lt.n-m) THEN
            dy = c(ns+1)
         ELSE
            dy = d(ns)
            ns = ns - 1
         END IF
         y = y + dy
      enddo
      return
      end
!
      subroutine hunt (xa,n,x,jlo)
!
      integer jlo,n,inc,jhi,jm
      real x, xa(n)
      logical ascnd
!
! Given an array xa(1:n), and given a value x, returns a value jlo
! such that x is between xa(jlo) and xa(jlo +1). xa(1:n) must be
! monotonic, either increasing or decreasing. jlo=0 or jlo=n is
! returned to indicate that x is out of range. jlo on input is taken
! as the initial guess for jlo on output
!
      ascnd = xa(n).ge.xa(1)    ! True if ascending order of table, false otherwise
      IF (jlo.le.0.or.jlo.gt.n) THEN ! Input guess not useful. Go to bisection
         jlo = 0
         jhi = n+1
         goto 3
      END IF
      inc = 1                   ! Set the hunting increment
      IF (x.ge.xa(jlo).eqv.ascnd) THEN ! Hunting up:
 1       jhi = jlo + inc
         IF (jhi.gt.n) THEN     ! Done hunting, since off end of table
            jhi = n + 1
         ELSE IF (x.ge.xa(jhi).eqv.ascnd) THEN ! Not done hunting
            jlo = jhi
            inc = inc + inc     ! so double the increment
      goto 1              ! and try again
      END IF                    ! Done hunting, value bracketed.
      ELSE                      ! Hunt down:
         jhi = jlo
 2       jlo = jhi-inc
         IF (jlo.lt.1) THEN     !Done hunting, since off end of table
            jlo = 0
         ELSE IF (x.lt.xa(jlo).eqv.ascnd) THEN ! Not done hunting
            jhi = jlo
            inc =  inc + inc    ! so double the increment
            goto 2              ! and try again
         END IF                 ! Done hunting, value bracketed
      END IF                    ! Hunt is done, so begin the final bisection phase:
!
 3    IF (jhi-jlo.eq.1) THEN
         IF (x.eq.xa(n)) jlo = n-1
         IF (x.eq.xa(1)) jlo = 1
         RETURN
      END IF
      jm = (jhi + jlo)/2
      IF (x.ge.xa(jm).eqv.ascnd) THEN
         jlo = jm
      ELSE
         jhi = jm
      END IF
!
      GOTO 3
!
      end
# endif
#endif !O_MOM
