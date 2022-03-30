* globk command file for time series, short-term combinations and velocity solutions
* Set COMB as a globk command line option for short-term combinations; since these commands
*          follow the commands for repeatabilities, they will take precedence when invoked.
* Set VEL  as a globk command-line option for velocity; since these commands follow
*          the commands for repeatabilities, they will take precedence when invoked.
* Set KEEP to allow glorg to be run separately after run; useful for testing
*          as a globk command-line option.
* Any combination, e.g. COMB+KEEP or VEL+KEEP will use both (or more) options.

* Lasted edited 180705

* << column 1 must be blank if not comment >>

* Use wild-cards in file names so that parallel runs do not overwrite
* each other. The "@" symbol is replaced with the .gdl-file name root.

* This group of commands must appear before any others:

 eq_file ~/gg/tables/igb14_comb.eq
ITRF08 eq_file ~/gg/tables/itrf08_comb.eq
# Optionally add a second eq_file for analysis-specific renames 
x eq_file my_eq_file.eq

# For time-series glred runs, it's faster not to use the "com_file" command
# (com_file need not be created; only needed if glorg or glsave are to be
# run by themselves after globk completes).
* This group of commands must appear before any others:
KEEP com_file @.com    ! Allows parallel runs, provide h-file list name is different
 srt_file @.srt
 srt_dir +1

# Normally for combined files, make_svs is not needed. Only needed if
# orbits are to be estimated.
x make_svs @.svs A

* End of commands that must appear first

* Solution file pointed to by the com file, if used
VEL sol_file @.sol

* ITRF augmented by now-defunct sites and recent IGS solutions;
* matched to igb14_comb.eq (ITRF2014; default) or itrf08_comb.eq ("ITRF08" command-line option)
 apr_file ~/gg/tables/igb14_comb.apr
ITRF08 apr_file ~/gg/tables/itrf08_comb.apr
# Optionally add additional apr files for other sites

* Set maximum chi2, prefit coordinate difference (m), and rotation (mas) for an h-file to be used;
 max_chii 13 3 100
# increase tolerances to include all files for diagnostics
x max_chi  100  5.0 20000

# Do not used an a priori rotation file with multi-day H-files
x in_pmu ~/gg/tables/pmu.usno

* Correct the pole tide when not compatible with GAMIT.
* Use caution here is binary hfiles from SINEX are used since these
* file do not contain status of pole tide correction.
 app_ptid all

* Invoke glorg
 org_cmd glorg.cmd

* Print file options
 crt_opt NOPR
 prt_opt NOPR GDLF
 org_opt PSUM GDLF FIXA RNRP
COMB org_opt ERAS PSUM CMDS GDLF FIXA RNRP
VEL org_opt ERAS PSUM CMDS GDLF VSUM FIXA RNRP
# The output .org-file will be generated from the
# .prt-file name from the globk command line
# sh_glred will name the glorg print files 
# To set an explicit name for the glorg print file
# globk print file name from the globk command line
x org_out glred_------.org
xCOMB org_out globk_------.org
xVEL org_out globk_------.org

* Coordinate parameters to be estimated and a priori constraints
 apr_site  all  10 10 10  0 0 0
VEL apr_site  all  10 10 10  1 1 1

* Rotation parameters to be estimated and a priori constraints
 apr_wob  10 10 1 1
 apr_ut1  10 1
# EOP tight if translation-only stabilization in glorg (also comment out mar_wob, mar_ut1)
x apr_wob .25 .25 .1 .1
x apr ut1 .25 .1
* For multiday combinations allow EOP's to change between days.
* The rate terms here depend on nature of gamit solutions:
* For RELAX. solutions, EOP rates are estimated and when combined
* with daily MIT .GLX-files, the rates are allowed to change and
* will be jointly estimated between the MIT and local file. When
* baseline processing is used, the rates would only apply to MIT
* files.
COMB mar_wob  3650 3650  365 365 
COMB mar_ut1  365 365
VEL mar_wob  3650 3650  0 0 
VEL mar_ut1  365 0
* Allow rotation when using SINEX files that do not contain EOP parameters
x apr_rot  10 10 10  1 1 1
x apr_rot  .25 .25 .25  .1 .1 .1
xCOMB mar_rot  3650 3650 365  365 365 365  0 0 0
xVEL mar_rot  3650 3650 365  0 0 0  0 0 0

* Translation a priori constraints
 apr_tran 1 1 1 0 0 0
xVEL apr_tran 1 1 1 1 1 1
xVEL mar_tran   3.65  3.65  3.65     0  0  0

# Allow scale variations for pre-1995 data and for data analyses that feature a change
# in SV PCVs (e.g. operational, pre-repro1 MIT or SOPAC h-files before Week 1400, Nov 2006)
SCALE apr_scale 10 0
xVEL apr_scale 10 1.
xVEL mar_scale  365   0

* If orbits free in GAMIT (RELAX.) and you want them fixed, use:
x apr_svs all F F F F F F FR 
* but if you are combining with globk h-files, better to leave them
* on but, if the models are incompatible, turn off radiation-pressure parameters,
x apr_svs all 100 100 100 10 10 10 0R  

* When using MIT GLX files, which have satellite phase center positions 
* estimated, fix antenna offsets to IGS a priori values by using:
 apr_svan  all  F F F

* When combining multiple h-files from the same epoch,
* estimate atmospheric zenith delay at common sites
 apr_atm common 1
COMB mar_atm common 3.65
VEL mar_atm common 3.65

* Optionally put a uselist and/or sig_neu and mar_neu reweight in a source file
x source ../tables/uselist
x source ../tables/monthly_reweights

* Turn off quake log estimates if in the eq_file
 free_log -1
 
* Write out a combined H-file
# Can substitute your analysis name for 'COMB' in the file name below
COMB out_glb H------_COMB.GLX
VEL out_glb @.GLX

* Remove scratch files for repeatability and combination runs
 del_scra yes
KEEP del_scra no

