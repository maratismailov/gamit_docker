* FILE trackRT_pbo.cmd
* Example TrackRT command file to process some PBO sites. 
* Run BNC with rtgpsout.unavco.org as the TCP/IP port (login information
* can be obtained from rtgps@unavco.org.) and select sites P496, P497, 
* P498, P505 to be cast on specific port (port 3765 for example).
* Run trackRT as
* trackRT -p 3765 -r p497 -d p496 p498 p505 -f trackRT_pbo.cmd -n P497 &
* (if BNC is run on another machine used -m <machine name> as well.
* A relative directory ../sp3_files is assumed to contain up to date
* sp3 file.

 sp3_dir ../sp3_files
  
* These positions are post April 4, Baja earthquake (extracted from
* globk output with sh_exglk -f <globk org file> -apr <apr file> 
* entries from  apr_file
* Nore: Additional sites can be given here

 site_pos
 P498_GGU -2313650.15153 -4835550.43945  3444474.55150   -0.00658    0.02008    0.01478 2010.259  0.0033  0.0059  0.0042  -1.0000 -1.0000 -1.0000
 P497_GGU -2315937.74208 -4838712.09671  3438545.10828   -0.00945    0.01927    0.01884 2010.259  0.0035  0.0063  0.0044  -1.0000 -1.0000 -1.0000
 P496_GGU -2319722.10265 -4842511.39340  3430709.79068   -0.01697    0.02370    0.02587 2010.259  0.0050  0.0089  0.0062  -1.0000 -1.0000 -1.0000
 P503_GGU -2325056.31018 -4826758.67413  3449209.09091   -0.01167    0.01951    0.01660 2010.259  0.0033  0.0059  0.0042  -1.0000 -1.0000 -1.0000
 P505_AGU -2309739.55681 -4802072.50955  3493258.58085   -0.00377    0.00608    0.00160 2010.404  0.0033  0.0057  0.0040 

* PBO sites antenna and receiver information (all are NetRS and thus C type)
* Again more sites can given). (Note: If antenna information extracted from
* gamit station.info, one space must be removed before radome string).
* Antenna/Radom combination must be in antmod_file for model to be used).
   ante_off
     p475   0.00 0.00 0.0083   TRM29659.00     SCIT  C
     p066   0.00 0.00 0.0083   TRM29659.00     SCIT  C 
     p472   0.00 0.00 0.0083   TRM29659.00     SCIT  C 
     p478   0.00 0.00 0.0083   TRM29659.00     SCIT  C 
     p494   0.00 0.00 0.0083   TRM41249.00     SCIT  C 
     p496   0.00 0.00 0.0083   TRM29659.00     SCIT  C  
     p497   0.00 0.00 0.0083   TRM29659.00     SCIT  C  
     p498   0.00 0.00 0.0083   TRM41249.00     SCIT  C 
     p500   0.00 0.00 0.0083   TRM29659.00     SCIT  C 
     p503   0.00 0.00 0.0083   TRM29659.00     SCIT  C 
     p505   0.00 0.00 0.0083   TRM29659.00     SCIT  C 
     p510   0.00 0.00 0.0083   TRM29659.00     SCIT  C 

   antmod_file /home/tah/gg/tables/antmod.dat

   dcb_file    /home/tah/gg/tables/dcb.dat

   data_type LCPC

   site_stats
     all   0.1 0.1 0.1   0.025 0.025 0.025
     p497  0.0 0.0 0.0   0.0 0.0 0.0

   atm_stats
     all   0.20 0.00010 0.000   ! Unit m/sqrt(sec) -> 0.0001 = 0.03 m/sqrt(day)
     p497  0.00 0.00000 0.000

    pos_root ? 1h

# Used at MIT for web output. Results can be viewed at
# http://chandler.mit.edu/kmeduna/
#    csv_root /net/chandler/var/www/kmeduna/trackrt_view/data/P497
 
# Useful to have this option.  Note file should be removed before trackRT
# is run and created when needed.
    update_file upd_app.cmd

