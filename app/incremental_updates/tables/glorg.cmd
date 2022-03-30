* glorg command file for time series, short-term combinations and velocity solutions
* Set VEL  as a globk command-line option for velocity; since these commands follow
*          the commands for repeatabilities, they will take precedence when invoked.
* Add choice of plate for final reference frame (e.g. EURA14 for Eurasia defined by
* Altamimi et al.'s (2017) ITRF2014 plate motion model, NOAM08 for North America
* defined by Altamimi et al's (2012) ITRF2008 plate motion model, NA1214 for North
* America defined by Blewitt et al. (2013) using rotated ITRF2014 velocities, etc.

* Last edited by MAF 180519

* << column 1 must be blank if not comment >>

* Parameters to be estimated
 pos_org  xtran ytran ztran xrot yrot zrot
VEL rate_org  xtran ytran ztran xrot yrot zrot 
# Optionally, if estimated scale in GLOBK:
SCALE pos_org  xtran ytran ztran xrot yrot zrot scale
xSCALE rate_org  xtran ytran ztran xrot yrot zrot 

#   or if translation-only
x pos_org  xtran ytran ztran
xVEL rate_org  xtran ytran ztran

* Downweight of height relative to horizontal (default is 10 10)
#   Heavy downweight if reference frame robust and heights suspect
x cnd_hgt  1000 1000

* Controls for removing sites from the stabilization 
# Vary these to make the stabilization more robust or more precise
 stab_it  4 0.8 3.0
x stab_it  4 0.5 4.0

* A priori coordinates to define the analysis reference frame
 apr_file ~/gg/tables/igb14_comb.apr
ITRF08 apr_file ~/gg/tables/itrf08_comb.apr

# ITRF2014 plate motion model defined by Altamimi et al. (2017)
AMUR14 apr_file ~/gg/tables/igb14_comb_amur08.apr
ANTA14 apr_file ~/gg/tables/igb14_comb_anta.apr
ARAB14 apr_file ~/gg/tables/igb14_comb_arab.apr
AUST14 apr_file ~/gg/tables/igb14_comb_aust.apr
CARB14 apr_file ~/gg/tables/igb14_comb_carb08.apr
EURA14 apr_file ~/gg/tables/igb14_comb_eura.apr
INDI14 apr_file ~/gg/tables/igb14_comb_indi.apr
NAZC14 apr_file ~/gg/tables/igb14_comb_nazc.apr
NOAM14 apr_file ~/gg/tables/igb14_comb_noam.apr
NUBI14 apr_file ~/gg/tables/igb14_comb_nubi.apr
PCFC14 apr_file ~/gg/tables/igb14_comb_pcfc.apr
SOAM14 apr_file ~/gg/tables/igb14_comb_soam.apr
SOMA14 apr_file ~/gg/tables/igb14_comb_soma.apr
SUND14 apr_file ~/gg/tables/igb14_comb_sund08.apr
# North America defined by Blewitt et al. (2013)
NA1214 apr_file ~/gg/tables/igb14_comb_na12.apr
# North America defined by Kreemer et al. (2018)
NA1714 apr_file ~/gg/tables/igb14_comb_na17.apr
# Nubia and Somalia defined by Saria et al. (2013)
NU1314 apr_file ~/gg/tables/igb14_comb_nu13.apr
SO1314 apr_file ~/gg/tables/igb14_comb_so13.apr

# ITRF2008 plate motion model defined by Altamimi et al. (2012)
AMUR08 apr_file ~/gg/tables/itrf08_comb_amur.apr
ANTA08 apr_file ~/gg/tables/itrf08_comb_anta.apr
ARAB08 apr_file ~/gg/tables/itrf08_comb_arab.apr
AUST08 apr_file ~/gg/tables/itrf08_comb_aust.apr
CARB08 apr_file ~/gg/tables/itrf08_comb_carb.apr
EURA08 apr_file ~/gg/tables/itrf08_comb_eura.apr
INDI08 apr_file ~/gg/tables/itrf08_comb_indi.apr
NAZC08 apr_file ~/gg/tables/itrf08_comb_nazc.apr
NOAM08 apr_file ~/gg/tables/itrf08_comb_noam.apr
NUBI08 apr_file ~/gg/tables/itrf08_comb_nubi.apr
PCFC08 apr_file ~/gg/tables/itrf08_comb_pcfc.apr
SOAM08 apr_file ~/gg/tables/itrf08_comb_soam.apr
SOMA08 apr_file ~/gg/tables/itrf08_comb_soma.apr
SUND08 apr_file ~/gg/tables/itrf08_comb_sund.apr
# North America defined by Blewitt et al. (2013)
NA1208 apr_file ~/gg/tables/itrf08_comb_na12.apr
# North America defined by Kreemer et al. (2018)
NA1708 apr_file ~/gg/tables/itrf08_comb_na17.apr
# Nubia and Somalia defined by Saria et al. (2013)
NU1308 apr_file ~/gg/tables/itrf08_comb_nu13.apr
SO1308 apr_file ~/gg/tables/itrf08_comb_so13.apr
# Use a regional stablization if available from a prior solution (comment out the itrf08 file)
x apr_file ../tables/regional.apr

* List of stabilization sites
#   This should match the well-determined sites in the apr_file
 stab_site clear
x source ~/gg/tables/igb14_hierarchy.stab_site
 source ~/gg/tables/igb14_comb.stab_site
xITRF08 source ~/gg/tables/igb08_hierarchy.stab_site
ITRF08 source ~/gg/tables/itrf08.stab_site
# Use a regional stabililization if available from a prior solution
x source ../tables/regional_stab_site

* Estimate rotation (Euler) vectors to be used with sh_org2vel to
* to rotate the solution to a block- or region-specific reference frame
xVEL plate eurasia kosg_2ps onsa_2ps nyal_4ps graz_2ps tlse_2ps kit3_2ps
xVEL plate eurasia vill_3ps mars_3ps cbre
xVEL plate weura kosg_2ps tlse_2ps vill_3ps mars_3ps
xVEL plate aegean milo kyra xris dioa leon mkn2 bodr roml omal koun seva
# Constrain the center-of-mass to the apr-file in plate estimate (comment out for global solutions)
xVEL NOPLATETRAN

* Equate the velocities of co-located sites
VEL eq_dist 1000 ndot
VEL eq_dist 1000 edot
VEL eq_dist 1000 udot

* Equate a few horizontal velocities for sites farther apart
VEL equate trab_gps ndot akto_gps ndot
VEL equate trab_gps edot akto_gps edot

* Unequate velocities that are incompatible
VEL unequate mad2_gps ndot
VEL unequate mad2_gps edot
VEL unequate mad2_gps udot

