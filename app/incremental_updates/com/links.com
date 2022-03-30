#!/bin/csh -f
if ($#argv == 0) then
   echo ' '
   echo '    LINKS.COM creates standard GAMIT links to gg/tables '
   echo "    where 'gg' is a link to your local gamit/globk path"
   echo ' '
   echo ' FORMAT' 
   echo ' '
   echo '    links.com frame year [day] [expt] [L-file] [UT1/pole]'
   echo ' '
   echo ' REQUIRED'   
   echo ' '
   echo '    frame    = B1950 or J2000    (for luntab, soltab if used ) '
   echo ' '
   echo '    year     = 1986-2000   (for nutabl, luntab, soltab if used) '
   echo ' '
   echo '    Default with 2 arguments:   ut1.bull_b   pole.bull_b'
   echo ' '
   echo ' OPTIONAL'
   echo ' '
   echo '    day, expt  = day-of-year, 4-char experiment name for L-file linkname'
   echo '                   set to l[expt][y].[day]'
   echo ' '                                       
   echo '    L-file       = L-file full path; default = ../tables/lfile.'
   echo ' '          
   echo '    ut1/pole = bull_b_      IERS Bulletin B      (1-day interval)   valid  1 Jan 1992 - current '
   echo '               bull_b_5day  IERS Bulletin B_5day (5-day interval)   valid  1 Sep 1987 - 25 Nov 1993 '
   echo '               usno         IERS Bulletin A      (rapid-service)    valid  1 Jan 1992 - current '  
   echo '               vlbi         IRIS                                    valid  5 Jan 1984 - 28 Jun 95'
   echo ' '
   echo ' EXAMPLES'          
   echo ' '                     
   echo '    links.com J2000 1992 '
   echo '    links.com J2000 1992 181 trex'
   echo '    links.com J2000 1992 181 trex lfile.new '
   echo '    links.com J2000 1992 181 trex lfile.new vlbi'   
   echo ' '                   
   echo " NOTE: If the link to 'nbody' is empty, indicating the absence of a PEP or JPL ephemeris, "
   echo "       then arc, model, and ngstot will revert to the old luntab., soltab., and nutabl. "
   echo "       Otherwise, the nbody ephemeris will be read and the nutation computed by calling "
   echo "       gamit/lib/MHB_2000 (copy of routine of same name in kf/gen_util/precess.f). "
   echo " " 
   echo ' NOTE: Now links to a read-only "core" file in gg/tables to'
   echo '        block HP machines from dumping core '
   echo '  '
   echo '  '
   echo ' CURRENT LINKS:'
   echo '  '

else      
# create standard links in a working directory
#

# Variable yr[1] = 4 char yr, yr[2] = 2 char yr, yr[3] = 1 char yr
   set yr = `sh_year -year $2`

   /bin/rm -f ut1.
   /bin/rm -f pole.
   /bin/rm -f nutabl.
   /bin/rm -f soltab.
   /bin/rm -f luntab.
# MOD TAH 200223: Don't remove local cooy or link                 
#  /bin/rm -f nbody
   /bin/rm -f gdetic.dat    
   /bin/rm -f guess_rcvant.dat    
   /bin/rm -f hi.dat
   /bin/rm -f leap.sec
# MOD TAH 190627: Only link svnav.dat if a local copy does
#  not exist.  Allows use of igs_metadata.snx file. 
#   /bin/rm -f svnav.dat
   /bin/rm -f svs_exclude.dat
   /bin/rm -f tform.dat  
   /bin/rm -f antmod.dat
   /bin/rm -f rcvant.dat 
   /bin/rm -f otl.grid 
   /bin/rm -f otl.list
   /bin/rm -f atl.grid 
   /bin/rm -f atl.list
   /bin/rm -f atml.grid 
   /bin/rm -f atml.list
   /bin/rm -f met.grid 
   /bin/rm -f met.list
   /bin/rm -f map.grid 
   /bin/rm -f map.list
# MOD TAH 190606: Fixed typo in line below:  gpt,grid
   /bin/rm -f gpt.grid 
   /bin/rm -f dcb.dat 
   /bin/rm -f otlcmc.dat 
   /bin/rm -f core

# MOD TAH 190606: Tested on nbody and only link old files
# if nbody does not exist.
# MOD TAH 200223: See if there is ../tables/ version of nbody.  If there
#  is use this version (same a svnav.dat below).  Use local copy if present.
   if( ! -e nbody ) then
      if( -e ../tables/nbody ) then
         ln -s ../tables/nbody .
#     Original code 
      else if ( ! -e ~/gg/tables/nbody ) then
         ln -s  ~/gg/tables/nutabl.${yr[1]}       nutabl. 
         ln -s  ~/gg/tables/soltab.${yr[1]}.$1    soltab. 
         ln -s  ~/gg/tables/luntab.${yr[1]}.$1    luntab.  
      else   
         ln -s  ~/gg/tables/nbody                 nbody 
      endif
   endif 

   ln -s  ~/gg/tables/gdetic.dat      gdetic.dat 
   ln -s  ~/gg/tables/guess_rcvant.dat      guess_rcvant.dat 
   ln -s  ~/gg/tables/hi.dat          hi.dat
   ln -s  ~/gg/tables/leap.sec        leap.sec
# MOD TAH 190627: Only link svnav.dat if a local copy does
#  not exist.  Allows use of igs_metadata.snx file. 
   if( ! -e svnav.dat ) then
      ln -s  ~/gg/tables/svnav.dat       svnav.dat
   endif
   ln -s  ~/gg/tables/svs_exclude.dat svs_exclude.dat    
   ln -s  ~/gg/tables/tform.dat       tform.dat    
   ln -s  ~/gg/tables/antmod.dat      antmod.dat 
   ln -s  ~/gg/tables/rcvant.dat      rcvant.dat 
   ln -s  ~/gg/tables/otl.grid        otl.grid 
   ln -s  ~/gg/tables/otl.list        otl.list 
   ln -s  ~/gg/tables/atl.grid        atl.grid 
   ln -s  ~/gg/tables/atl.list        atl.list
# MOD TAH 190606: Fixed typos in next two lines atlm. 
   ln -s  ~/gg/tables/atml.grid       atml.grid 
   ln -s  ~/gg/tables/atml.list       atml.list 
   ln -s  ~/gg/tables/met.grid        met.grid 
   ln -s  ~/gg/tables/met.list        met.list 
   ln -s  ~/gg/tables/map.grid        map.grid
   ln -s  ~/gg/tables/map.list        map.list         
   ln -s  ~/gg/tables/gpt.grid        gpt.grid 
   ln -s  ~/gg/tables/dcb.dat         dcb.dat 
   ln -s  ~/gg/tables/otlcmc.dat      otlcmc.dat   
   ln -s  ~/gg/tables/core            core

# MOD TAH 200503: Checked to see if local valid copy.  If so than link to gg/tables
#  not made.
   if( ! -e otides.dat ) then
      ln -sf ~/gg/tables/otides.dat            otides.dat
   endif

     
   if ($#argv < 6) then
     ln -s  ~/gg/tables/ut1.bull_b   ut1.    
     ln -s  ~/gg/tables/pole.bull_b  pole.   
   else 
      ln -s  ~/gg/tables/ut1.$6    ut1.    
      ln -s  ~/gg/tables/pole.$6   pole.   
   endif

  if ($#argv > 3) then
    /bin/rm -f    l$4${yr[3]}.$3                                              
     if ($#argv < 5) then
        ln -s ../tables/lfile.  l$4${yr[3]}.$3                       
     else 
        ln -s $5   l$4${yr[3]}.$3                                              
     endif
  endif
#   
endif
#ls -l `find . -type l -print`
# MOD TAH 971215: Only do ls -l if there are links in the directory
if ( `find  . -type l -print | wc | awk '{print $1}'` > 0 ) ls -l `find  . -type l -print`
