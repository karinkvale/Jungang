!======================= include file "size.h" =========================

!-----------------------------------------------------------------------
!     USER INPUT:
!-----------------------------------------------------------------------
!   ntnpzd   = number of npzd tracers
!     imt    = number of grid points in the longitudinal direction
!              (calculated points are from 2 through imt-1. End points
!               are boundaries)
!     jmt    = number of grid points (latitude rows) in the latitudinal
!              direction (calculated points are from 2 through jmt-1.
!              End points are boundaries)
!     km     = number of grid points in the vertical direction
!              (calculated points are from 1 through km)
!     nt     = number of tracers (temperature, salinity, ...)
!     nsrc   = number of tracer with sources
!     kpzd   = depth for limited npzd model
!     jmz    = size for "unrotated" zonal averages
!     jmzm1  = jmz minus one
!     mnisle = maximum number of islands (unconnected land masses)
!     maxipp = maximum number of all island perimeter points
!-----------------------------------------------------------------------
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Index for npzd biology tracers to remove hard-wiring in tracer.F and npzd_src.F
! by K.Kvale
!   ibiomp = index for microplastic particles
!   ibiompa = index for microplastic in phyt aggregates
!   ibiompp = index for microplastic in poo
!   ibion = index for phosphate
!   ibiop = index for phyt
!   ibioz = index for zoop
!   ibiod = index for detritus
!   ibiono3 = index for no3
!   ibiodiaz = index for diaz
!   ibiodfe = index for dfe
!   ibiodetrfe = index for detrfe
!   ibiod_B = index for caco3-ballasted detritus
!   ibioc = index for calcifiers
!   ibiocaco3 = index for caco3
!   ibiocaco3c13 = index for caco3 c13
!   ibiodiat = index for diatoms
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      character(10) :: mapb

#if defined O_npzd
      integer ntnpzd
      parameter(ntnpzd=4
# if defined O_npzd_nitrogen
     $                +1
#  if !defined O_npzd_simple
     $                +1
#  endif
# endif
# if defined O_plastic
     $                +3
# endif
# if defined O_zoop_det
     $                +1
# endif
# if defined O_npzd_iron
     $                +2
# endif
#if defined O_cal_kk
     $                +1
#endif
#if defined O_kk_caco3tr
     $                +1
#endif
#if defined O_kk_ballast
     $                +1
#endif
#if defined O_kk_diat
     $                +1
#endif
#if defined O_isotopes
     $                +16
#endif
     $                  )
      common /npzd_c/ mapb(ntnpzd)
      integer ibion, ibiop, ibioz, ibiod
# if defined O_plastic
      integer ibiomp, ibiompa, ibiompp
# endif
# if defined O_zoop_det
      integer ibiodz
# endif
# if defined O_npzd_nitrogen
      integer ibiono3
#  if !defined O_npzd_simple
      integer ibiodiaz
#  endif
# endif
# if defined O_npzd_iron
      integer ibiodfe, ibiodetrfe
# endif
# if defined O_kk_diat
      integer ibiodiat
# endif
# if defined O_cal_kk
      integer ibioc
# endif
# if defined O_kk_caco3tr
      integer ibiocaco3
# endif
# if defined O_kk_ballast
      integer ibiod_B
# endif
# if defined O_isotopes
      integer ibiocaco3c13, ibiodic13, ibiophytc13, ibiodiazc13
      integer ibiozoopc13
      integer ibiococcc13, ibiodiatc13, ibiodetrc13
      common /npzd_i/ ibiocaco3c13, ibiodic13, ibiophytc13, ibiodiazc13
      common /npzd_i/ ibiozoopc13, ibiococcc13, ibiodiatc13, ibiodetrc13

      integer ibiodin15, ibiophytn15, ibiodiazn15
      integer ibiozoopn15
      integer ibiococcn15, ibiodiatn15, ibiodetrn15
      common /npzd_i/ ibiodin15, ibiophytn15, ibiodiazn15
      common /npzd_i/ ibiozoopn15, ibiococcn15, ibiodiatn15, ibiodetrn15

!      integer ibiosil30
      integer ibiodiats30
!      common /npzd_i/ ibiosil30
      common /npzd_i/ ibiodiats30
# endif
      common /npzd_i/ ibion, ibiop, ibioz, ibiod
# if defined O_plastic
      common /npzd_i/ ibiompa, ibiompp, ibiomp
# endif
# if defined O_zoop_det
      common /npzd_i/ ibiodz
# endif
# if defined O_npzd_nitrogen
      common /npzd_i/ ibiono3
#  if !defined O_npzd_simple
      common /npzd_i/ ibiodiaz
#  endif
# endif
# if defined O_npzd_iron
      common /npzd_i/ ibiodfe, ibiodetrfe
# endif
# if defined O_kk_diat
      common /npzd_i/ ibiodiat
# endif
# if defined O_cal_kk
      common /npzd_i/ ibioc
# endif
# if defined O_kk_caco3tr
      common /npzd_i/ ibiocaco3
# endif
# if defined O_kk_ballast
      common /npzd_i/ ibiod_B
# endif
#endif

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer imt, jmt, km, nt, nsrc, kpzd, nat, jmz, jmzm1, mnisle
      integer maxipp, jmw, jsmw, jemw

#ifndef O_TMM
      parameter (imt=  102, jmt=  102, km= 19)
#else
      parameter (imt= 1, jmt= 1, km= 19)
#endif
      parameter (nt=2
#if defined O_matrix
     $             +km
#endif            
#if defined O_idealage
     $             +1
#endif
#if defined O_plastic
     $             +3
#endif
#if defined O_zoop_det
     $             +1
#endif
#if defined O_carbon
     $             +1
# if defined O_carbon_14
     $             +1
# endif
#endif
#if defined O_cfcs_data || defined O_cfcs_data_transient
     $             +2
#endif
#if defined O_kk_si
     $             +1
#endif
#if defined O_npzd_alk
     $             +1
#endif
#if defined O_npzd_o2
     $             +1
#endif
#if defined O_npzd
     $             +4
# if defined O_kk_ballast
     $             +1
# endif
# if defined O_cal_kk
     $             +1
# endif
# if defined O_kk_caco3tr
     $             +1
# endif
# if defined O_npzd_nitrogen
     $             +1
#  if !defined O_npzd_simple
     $             +1
#  endif
# endif
# if defined O_npzd_iron
     $             +2
# endif
# if defined O_kk_diat
     $             +1
# endif
# if defined O_isotopes
     $             +17 
# endif
#endif
     $               )
      parameter (nsrc=0
#if defined O_idealage
     $               +1
#endif
#if defined O_plastic
     $               +3
#endif
#if defined O_zoop_det
     $               +1
#endif
#if defined O_carbon
     $               +1
# if defined O_carbon_14
     $               +1
# endif
#endif
#if defined O_kk_si
     $               +1
#endif
#if defined O_npzd_alk
     $               +1
#endif
#if defined O_npzd_o2
     $               +1
#endif
#if defined O_npzd
     $               +4
# if defined O_kk_ballast
     $               +1
# endif
# if defined O_cal_kk
     $               +1 
# endif
# if defined O_kk_caco3tr
     $               +1
# endif
# if defined O_npzd_nitrogen
     $             +1
#  if !defined O_npzd_simple
     $             +1
#  endif
# endif
# if defined O_npzd_iron
     $               +2 
# endif
# if defined O_kk_diat
     $               +1
# endif
#endif
#if defined O_isotopes
     $               +17
#endif
     $                 )
      parameter (kpzd=km)
      parameter (nat=2
#if defined O_carbon && defined O_carbon_co2_2d
     $              +1
#endif
     $, jmz=jmt, jmzm1=jmz-1)
      parameter (mnisle=50, maxipp=5000)

#if !defined O_min_window
      parameter (jmw=jmt)
#else

!     for UNI-TASKING: "jmw" is set to the minimum for each option class
!     "jmw" may be increased up to "jmt"

# if defined O_fourth_order_tracer_advection || defined O_fct || defined O_quicker || defined O_pressure_gradient_average || defined O_biharmonic
#  if defined O_pressure_gradient_average
#   if defined O_biharmonic  || defined O_fourth_order_tracer_advection || defined O_fct || defined O_quicker
      parameter (jmw=5)
#   else
      parameter (jmw=4)
#   endif
#  else
      parameter (jmw=4)
#  endif
# else
      parameter (jmw=3)
# endif
#endif

!-----------------------------------------------------------------------
!     set first and last calculated row within the MW. other rows
!     are used as buffers
!-----------------------------------------------------------------------

!     jsmw   = 1st calculated row within the MW
!     jemw   = last calculated row within the MW
#ifndef O_TMM
      parameter (jsmw=2, jemw=jmw-1)
#else
      parameter (jsmw=1, jemw=1)
#endif
! Moses-Triffid land model

! POINTS = Maximum number of points in grid.
! STEPSM = Maximum number of timesteps in a day.
! klmax = maximum ocean depth levels over which the land model can exist

      integer POINTS, STEPSM, klmax
      parameter (POINTS=14300, STEPSM=24, klmax=0)

! NNVG  = Number of non-vegetation surface types.
! NPFT  = Number of plant functional types.
! NTYPE = Number of surface types.
! SOIL  = Index of the surface type 'Soil'
! Land surface types :
!     1 - Broadleaf Tree
!     2 - Needleleaf Tree
!     3 - C3 Grass
!     4 - C4 Grass
!     5 - Shrub
!     6 - Soil

      integer NNVG, NPFT, NTYPE, SOIL
      parameter (NNVG=4, NPFT=5, NTYPE=6, SOIL=6)
