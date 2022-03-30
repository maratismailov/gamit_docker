#2345678901234567890123456789012345678901234567890123456789012345678901234567890
# This is an example and template of a generic track command file. It assumes
# that RINEX observation files for five sites, one fixed ("FFFF") and four
# kinematic ("AAAA", "BBBB", "CCCC" and "DDDD"), SP3 orbits and IONEX files
# exist in the current working directory for the day to be processed (2019-100
# in this example).
#
# Before running track using this template users must do the following, based
# on their own experiment needs:
# 1. Download or copy RINEX observation files for all sites, e.g.
#    "sh_get_rinex -archive cddis sopac unavco -yr 2019 -doy 100 -sites ffff aaaa bbbb cccc dddd"
# 2. Download precise orbits, e.g.
#    "sh_get_orbits -orbit igsf -yr 2019 -doy 100 -nofit"
# 3. Download IONEX files, e.g. "sh_get_ion -yr 2019 -doy 100"
# 4. Calculate the GPS week and day corresponding to the year and day to be
#    processed, e.g. "doy 2019 100" gives GPS week 2048, day 3.
# 5. Define (a) accurate a priori coordinates for each site under "site_pos",
#    (b) a priori constraints consistent with their accuracy under "site_stats"
#    and (c) antenna (and, optionally, receiver "PCN" code) information under
#    "ante_off".
# 6. Change the "interval" to match the data being processed.
# Users may also ultimately want or need to adjust solution parameters, e.g.
# "float_type", or output options, e.g. "out_type", "sum_file", etc.
# See ~/gg/help/track.hlp for more information.
#
# track may then be run with a command such as:
# track -f track.cmd -d 100 -week 20483 -s 2019 19 ffff aaaa bbbb cccc dddd >& ffff_2020-100.log


* INPUT FILES *

# RINEX observation files
@ OBS_FILE
@    Site   RX_file   Type
 obs_file
  <S03> <S03><day>0.<S02>o F
  <S04> <S04><day>0.<S02>o K
  <S05> <S05><day>0.<S02>o K
  <S06> <S06><day>0.<S02>o K
  <S07> <S07><day>0.<S02>o K

# Orbit file
@ NAV_FILE  <name>  <SP3/NAV>
 nav_file igs<week>.sp3 SP3

# Differential code biases file
@ DCB_FILE <file name>
 dcb_file ~/gg/tables/dcb.dat

# IONEX (ionosphere) file
@ IONEX_FILE <file name>
 ionex_file igsg<day>0.<S02>i


* SITE METADATA *

# A priori site coordinates
@ SITE_POS
@   Site   <X (m)>  <Y (m)>  <Z (m)>
 site_pos
  <S03>  1234567.89000  1234567.89000  1234567.89000
  <S04>  1234567.89000  1234567.89000  1234567.89000
  <S05>  1234567.89000  1234567.89000  1234567.89000
  <S06>  1234567.89000  1234567.89000  1234567.89000
  <S07>  1234567.89000  1234567.89000  1234567.89000

# A priori site constraints and process noise
@ SITE_STATS
@   Site  <Apriori Sigma in XYZ>  <RW noise in XYZ>
 site_stats
  all   0.05 0.05 0.05 1 1 1  ! Set values for all sites (except fixed site)
  <S04> 0.05 0.05 0.05 1 1 1
  <S05> 0.05 0.05 0.05 1 1 1
  <S06> 0.05 0.05 0.05 1 1 1
  <S07> 0.05 0.05 0.05 1 1 1

# Antenna (and receiver) metadata
# N.B. Antenna and radome name must match the entries in ~/gg/tables/antmod.dat.
# The antenna code is 20 characters long, the IGS standard, but note that the
# same field in station.info contains one extra space between the antenna and
# radome. The receiver tracking code may be found in ~/gg/tables/rcvant.dat for
# the receiver being used).
@ ANTE_OFF
@   Site   <ARP dN (m)> <ARP dE (m)> <ARP dU (m)> <Antenna Name> <rcv code>
 ante_off
  <S03>  0.0000  0.0000  0.0083 TRM29659.00     SCIT C
  <S04>  0.0000  0.0000  0.5003 LEIAR10         NONE C
  <S05>  0.0000  0.0000  1.4331 TRM57971.00     NONE C
  <S06>  0.0000  0.0000  1.0476 SEPALTUS_NR3    NONE N
  <S07>  0.0000  0.0000  1.1856 TRM41249.00     NONE C

# ANTEX (antenna calibrations) file
@ ANTMOD_FILE <file name>
 antmod_file ~/gg/tables/antmod.dat


* PROCESSING CONTROLS *

# Preset based on baseline length
@ MODE <Type>
 mode long

# This command is needed if interval is not specified, or is inconsistent, in
# the RINEX headers, or a longer interval than that in the RINEX file is wanted
# (e.g. use "interval 30" to test a 1 Hz data solution).
@ INTERVAL <seconds>
 interval 30

# Time unit so process noise values, above, will be in m^2/<time_unit>
@ TIME_UNIT <epoch/sec/min/hour/day>
 time_unit sec

# Atmospheric delay model
-use_gptgmf 0.5  ! Syntax prior to track version 1.40 (2020/02/27)
@ ATM_MODELC <MODEL> <Relative humidity (0-1)>
 atm_modelc GMF 0.5

# Ambiguity resolution controls
@ FLOAT_TYPE <Start> <Decimation> <Type> <Float sigma Limits(2)> <WL_Fact> <Ion_fact> <MAX_Fit> [RelRank]

# Kalman filter type
@ BACK_TYPE  <string>
 back_type smooth


* OUTPUT OPTIONS *

# Output position format(s)
@ OUT_TYPE <string>
 out_type XYZ+DHU

# Output position file prefix
@ POS_ROOT <string>
 pos_root <S03>_<S01>-<day>

# Output summary file
@ SUM_FILE <string>
 sum_file <S03>_<S01>-<day>.sum


* OTHER COMMANDS (see ~/gg/help/track.hlp for more information) *
@ TR_GNSS <RGEC>
@ REF_NEU <X (m)>  <Y (m)>  <Z (m)>
@ TIMEDEP_PROCNS
@  Site    Sig XYZ (m/sqrt(t))  Start YY MM DD MN Sec End YY MM DD MN Sec
@ ATM_STATS
@   Site  <Apriori Zenith delay sigma> <RW noise in Zenith delay> [<RW dH/dt noise> or SCALE]
@ ATM_BIAS
@    Site   <Atmospheric delay offset (m)>
@ ATM_FILE <File name>
@ RCV_TYPE
@   Site  <Receiver code N/P/C>
@ BF_SET  <Max gap>  <Min good>
@ DEBUG  <Start EP> <End EP>
@ DATA_NOISE <L1> <L2> <P1> <P2> <Elev Weight> [PRN]
@ DATA_TYPE <choice 1> <choice 2> ...
@ OUT_SIG_LIMIT <sigma (m)>
@ RMS_EDIT_TOL <n-sigma limit>
@ EDIT_SSV <site> <prn #> <start time> <stop time>
@ USR_ADDBF <site> <prn #> <time (ymdhms)>
@ USR_DELBF <site> <prn #> <time (ymdhms)>
@ AMBIN_FILE <file name>
@ RM_CSLIP <# MW WL> <# EX WL> <CSL_RR> <CSL_addChi> <CSL_minChi>
@ MIN_TOLS <min LC sig> <WL Tau> <Dynamic Tol>
@ AMB_CYCLE <Samples> <Relative Rank> <Max search>
@ ION_STATS <Jump> <ION PPM> <ION Weight> <ION height> <ION spatial>
@ RES_ROOT <string>
@ WLS_ROOT <string>
@ RWL_ROOT <string>
@ IONLOS_FILE <site name> <file name>
@ CUT_OFF <min elevation angle>
@ START_TIME  <Year Month day hour min sec>
@ NUM_EPOCHS <number>
@ EXCLUDE_SVS <list of PRN numbers to be excluded>
@ STOPGO_MODE <Variance reduction>
@ MWWL_JUMP <tol cycles>
@ USE_BLQ

