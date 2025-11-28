This directory contains updates to UVic ESCM version 2.9 (and subsequent modifications made 
at Kiel/GEOMAR) that allow for the extraction of transport matrices and coupling of the 
biogeochemical model to the Transport Matrix Method (TMM; https://github.com/samarkhatiwala/tmm/). 
For details see Kvale et al. (2017), Evaluation of the Transport Matrix Method for simulation 
of ocean biogeochemical tracers, GMD: 
https://thredds.geomar.de/thredds/catalog/open_access/kvale_et_al_2017_gmdd/catalog.html

This code must be used in conjunction with the full version 2.9 of UVic: 
UVic_ESCM.2.9.updated.tar.gz (2013-8-1) available from http://wikyonos.seos.uvic.ca/model/

To use this code you must set the paths in mk.in to the top level of the full UVic 2.9 
source tree (variable Version_Directory) and this directory (variable Source_Directory).
----------------------------------------------------------------------------------------
Feb 2019 update (K. Kvale):

Isotopes of 13C, 15N, and 28Si are added.
CaCO3 is allowed to scavenge Fe.
Mixed phytoplankton become "low latitude calcifiers" and zooplankton no longer calcify. This is done because the explicit calcifiers now only represent high latitude species, with no representation in the low latitudes, which is not the reality.

The Paulmier stoichiometric correction of isalk in tracer is removed because it caused alkalinity to be far too high.

----------------------------------------------------------------------------------------
Jan 2019 update (K. Kvale):

Bug fixes have been added to solve a PO4 leak as well as drifting alkalinity and ocean DIC concentrations.
An ideal age tracer (O_idealage) has also been included from Wolfgang Koeve.
And an option to include surface ocean albedo adjustment as a function of phytoplankton concentration (higher concentrations produce darker albedo) and calcite concentration (higher concentrations produce brighter albedo) has been added as O_phyt_albedo

----------------------------------------------------------------------------------------
October 2018 update (K.Kvale):

A new subgrid bathymetry file has been included. It was unclear how the old one was made, and the new one looks very different. But I have more confidence in it. Unforutnately the model is quite sensitive to this change and po4 and no3 must be re-tuned.

This version of the code works both online and offline. Models conserve carbon, po4, silica. Stable parameter sets are available. The code is also set up for offline parameter calibrations.

Offline sediments are still not working
----------------------------------------------------------------------------------------
October 2017 update (K. Kvale):

Hydrothermal iron and silicate inputs are added as masks.

Hardwired increase in detritus sinking speed away from the surface is now given a variable name wdd.
----------------------------------------------------------------------------------------
May 2017 update (K. Kvale):

This model is updated to include calcifiers (Kvale et al., Atmosphere-Ocean, 2015) and Diatoms
A index for biological tracers replaces hard-wired numbers in npzd_src.F and tracer.F:
ibion,ibiop,ibioz,ibiod,ibiono3,ibiodiaz,ibiodfe,ibiodetrfe

Silica weathering discharge is accounted for in a new variable globaldisch

A new subroutine to initialize the biology is added called npzd_ini

----------------------------------------------------------------------------------------
Feb 2017 update (K. Kvale):

This model has been modified to include L. Nickelsen et al. (2015) GMD iron model.

Capacity for the online model to read matrix output and idealized BGC input is also added. To use WOA or GLODAP BGC input, set initbgc=true in control.in and turn on O_idealized_bgc_ic in mk.in. To turn off, set initbgc.false and/or turn off option in mk.in. To use TMM input, set initbgc=true and turn on O_tmm_bgc_ic.

It is found O_correct_sea_ice_to_ocean_stress causes a slow salinity leak over 1000 year timescales.

Other additions to the code come from D. Keller (see below)

----------------------------------------------------------------------------------------
July, 2016 --- David Keller (dkeller@geomar)

The previous update merged the UVic code at Kiel with M. Eby's official
02 update for the UVic model.  Contained in this update is that code,
which includes bug fixes and improvements 1-12, plus new bug fixes and
now also options for a dynamic iron cycle.

1) A bug in the ocean pCO2 calculation was identified and patched by D.
Keller.  The offical 02 update incorporates this fix, but M. Eby has
re-written the patch and changed significant portions of the "co2calc.F"
file while doing so.  Therefore, in this merged version no mention or
option is available to turn on or off the pCO2 bug fix.

2) A bug in “setmom.F” that effected the amount of light in the upper
ocean has been corrected (see code line 967 for notes on the
correction). Since correcting this bug which affected phytoplankton
growth some NPZD parameters had to be changed (see “control.in”) to
achieve similar results to those reported in the Keller et al. 2012 GMD
article.


3) In “setmom.F” the equation that determines the sinking speed of
detritus has been changed (see code line 962).  We made these changes so
that the sinking speed of detritus better fits the observations made in:
Berelson, 2002. Particle settling rates increase with depth in the
ocean. Deep-Sea Research II. 49. pgs. 237-251.  The previous equation
had been based on the lowest estimates of sinking speed in Berelson
(2002) and the new one now reflects the mean values of this study.



4) The iron mask, O_fe_dissolved.nc, that is used to limit phytoplankton
growth has had a few very high outlying values removed.  These high
values were artifacts of the original regridding of the BLING iron data
and occurred mostly around Indonesia where there are no observations of
extremely high Fe. In addition, the dissolved Fe concentration at a few
grid points in the North Western Arabian sea has been increased to
better reproduce observed Fe concentrations.



5) Oxygen stoichiometry during N2 fixation has been corrected in
tracer.F (lines 572-583).  Previously, N2 fixation had not been properly
accounted for and 1.5 moles of O2 were consumed for each mole of N fixed
by diazotrophs.  With this correction now only 1.25 moles of O2 are
consumed during N2 fixation according to the following stoichiometry:



     On the O2 production side of N2 fixation within the nitrogenase
     complex:



N2 + (8 H+) + (8 e-) + 16ATP --> 2 NH3 + (H2) + 16ADP + 16PO4		(1)



assuming that the protons and electrons come form the splitting of water

(1) becomes the bulk formula:



N2 + 3H2O --> 2 NH3 + 1.5 O2								(2)



Then to oxidise NH3 to NO3 via nitirification:



NH3 + 3/2O2 --> NO2- + H2O +H+             	  				(3)

NO2- + 1/2O2 --> NO3-									(4)



Summing up (3) and (4) in the bulk formula :



NH3  +  2O2   --> NO3 + H2O 								(5)



Combining (2) and (5) for every mole of NO3 produced 1.25 of O2 are
needed.



6) The model now correctly accounts for sources and sinks of alkalinity
during N2 fixation in tracer.F (lines 584-595).  The modified equations
are based on stoichiometric equations and parameters from Paulmier et
al. 2009. Stoichiometries of remineralization and denitrification in
global biogeochemical ocean models.  Biogeosciences 6, pgs. 923-935.

7) The high Southern Ocean mixing scheme code that was described in
Keller et al., 2012 GMD is an option in our code, but not the official
02 update.

8) Following Getzlaff and Dietze, 2013, GRL Vol. 40. an option to
increase the tropical zonal isopycnal diffusivity has been added to
isopyc.F and isopyc.h.  The option is called O_increase_isopyc_diff and
can be applied in two ways.  The first option
O_increase_isopyc_diff_everywhere applies the scheme as in Julia's
paper.  The second option O_increase_isopyc_diff_smooth is a modified
version of this code that was added by C. Somes and it smooths the
transition between the area where the scheme is applied and other areas.
Note that the use of these options will improves ocean oxygen
distributions, but lowers N-fixation and denitrification to much lower
levels.  I have not had time to tune the parameters to get N-fixation
and denitrification back to the levels they were at in Keller et al.
2012 GMD.

9) In tracer.F on lines 340 and 342 the term "phin" has been removed
from the equations to correct another light bug.  This bug caused too
much light attenuation due to phytoplankton biomass to occur in the
third layer and below.
10) A max statement has been added to tracer.F, line 583, to prevent
"deni", the denitrification term, from going negative.

11) Note also that in setmom.F there is a bug that we fixed a long time
ago that has never been fixed by M. Eby.  This bug is on line 161 where
the incorrect code is "call getvar ('O_tidenrg', ...".  The correct code
is "call getvar ('tidenrg',...".  If you don't use the tidal mixing
scheme this bug is not a problem, which is why it has not been fixed at
UVic.

12) Tronje Kemena identified some code that affected momentum and
sea-ice when in a snowball Earth configuration.  Note this is not a bug as
perviously stated, the old formulation and this one are both correct. 
This code limitation was in the embm.F and evp.F files.  
The fix for this code limitation can be activated by using the "O_correct_ice_to_ocean_stress"
option. The details of the bug are:

Line 192 to 193 in the original  embm.F:

flux(i,j,nat+1) = flux(i,j,nat+1) +  dts*sbc(i,j,itaux) flux(i,j,nat+2)
= flux(i,j,nat+2) + dts*sbc(i,j,itauy)

defines the zonal and meridional stress exerted by the wind on the
ocean’s surface. Irrespective of sea ice-cover the winds drive the
ocean. This is unrealistic. A more realistic approach is to weight the
stress exerted by the winds on the ocean, sbc(i,j,itaux),  with the ice
cover, aice(i,j,2), such that a fully ice-covered water column is
shielded from the winds blowing over the ice:

flux(i,j,nat+1) = flux(i,j,nat+1) &                               +
dts*sbc(i,j,itaux)*(1-aice(i,j,2)) flux(i,j,nat+2) = flux(i,j,nat+2) &  
                            + dts*sbc(i,j,itauy)*(1-aice(i,j,2))

Line 194 to 195 in the original  embm.F:

flux(i,j,nat+1) = flux(i,j,nat+1) + dts*xint(i,j) flux(i,j,nat+2) =
flux(i,j,nat+2) + dts*yint(i,j)

adds the zonal and meridional stress (xint(i,j) and yint(i,j)) exerted
by drifting ice on the ocean’s surface to the total stress received by
the ocean. This stress should be scaled with the ice cover. Otherwise
the stress would be overestimated in situations with little ice cover.

Hence, a more realistic replacement for lines 194 to 195 is

flux(i,j,nat+1) = flux(i,j,nat+1) + dts*xint(i,j) &                     
         * aice(i,j,2) flux(i,j,nat+2) = flux(i,j,nat+2) + dts*yint(i,j)
&                               * aice(i,j,2).

Line 632 & 633, supposedly, assigns the stress exchanged between the
ocean’s surface and the bottom edge of sea ice:

xint(i,j) = s11 + s12 yint(i,j) = s21 + s22

However, this does not apply because (as suggested by a reconciliation
of the code with equations in Hunke and Dukowicz (1997))  s11, s12, s21
and s22 describe divergences of internal ice stresses rather than the
exchange of momentum between sea ice and the ocean’s surface.

The following formulation is more realistic in that it calculates the
respective momentum exchange as a function of the difference of ocean
surface currents and ice drift speeds:

xint(i,j) = – (vrel*waterx(i,j)-vrel*umsk(i,j) & *(uice(i,j)*costh –
vice(i,j)*sinth)) yint(i,j) = – (vrel*watery(i,j)-vrel*umsk(i,j) &
*(vice(i,j)*costh + uice(i,j)*sinth))

The terms waterx, watery and uice(i,j)*costh – vice(i,j)*sinth,
vice(i,j)*costh + uice(i,j)*sinth are the water and ocean surface
velocities acting at the boundary ocean-ice. This terms include a
horizontal rotation scheme to parametrize the unresolved Ekman boundary
layer.

In line 335, 336, 337, 338 of evp.f the water velocities at the boundary
are calculated and in line 617, 618 of evp.f the relative ocean- sea ice
velocities are calculated. In the original formulation the geostrophic
ocean velocities are used. More realistic is to use directly the ocean
velocities of the surface layer. The model code is modified as
following:

line 335, 336, 337, 338 waterx(i,j) = umsk(i,j)*(sbc(i,j,isu)*costh - &
sbc(i,j,igv)*sinth) watery(i,j) = umsk(i,j)*(sbc(i,j,isv)*costh + &
sbc(i,j,igu)*sinth)

line 617, 618 uorel = sbc(i,j,isu) – uice(i,j) vorel = sbc(i,j,isv) –
vice(i,j)

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

13) In setmom.F there was a problem with having a negative sign in the term that
determines rcab(1).  In the previous version calcite production was subtracted 
twice when the depth was only one grid box deep.  This neg. sign has been removed.

14) In tracer.F the tracers were limited to positive values before the npzd_src
routine is called.  This made the tracer limiting statements in the subroutine
useless and could have causes some minor numerical errors.  Therefore the tracers
are only limited to positive values in the npzd_src routine. To revert to this
method use the option O_limit_tracers_to_pos_in_tracers

15) The gas exchange parameters for CO2 and O2 have been updated in gasbc.F to 
match those in Wanninkhof 2014.  The old parameterizations are commented out, 
but still therefor those who are interested in the difference.

16) The code for a dynamic iron cycle, i.e., Nickelsen et al., 2014 GMD, has
been added as to the default version of UVic.  The old option to use a masking
approach is still available.

17)  An option has been added so that the attenuation of light due to CDOM can be
accounted for.  This new formulation, which increases the attenuation of pure water (kw)
by 20% uses values that are standard in the marine optics literature (see Tedetti et al.,
2007 GRL and Morel et al., 2007 LnO).



