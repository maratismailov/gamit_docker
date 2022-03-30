#!/bin/csh -f      

#  Script links.tables to link from an experiment tables directory to global gamit/tables

if ($#argv == 0) then
   echo ' '
   echo '    LINKS.ARC creates standard GAMIT links necessary to run arc in the current directory'
   echo ' ' 
   echo ' FORMAT' 
   echo ' '
   echo '    links.arc frame year [ut1/pole]'
   echo ' '
   echo ' REQUIRED'               
   echo ' '  
   echo '    frame    = B1950 or J2000  (for soltab, luntab is used)'
   echo ' '
   echo '    year     = 1986-2000   (for nutabl, soltab, luntab if used)'  
   echo ' '
   echo '    Default with 2 arguments:   ut1.usno   pole.usno'
   echo ' '
   echo ' OPTIONAL'
   echo ' '        
   echo '    ut1/pole = bull_b_      IERS Bulletin B      (1-day interval)   valid  1 Jan 1992 - current '
   echo '               bull_b_5day  IERS Bulletin B_5day (5-day interval)   valid  1 Sep 1987 - 25 Nov 1993 '
   echo '               usno         IERS Bulletin A      (rapid-service)    valid  1 Jan 1992 - current '  
   echo '               vlbi         IRIS                                    valid  5 Jan 1984 - 28 Jun 95'
# MOD TAH 200326: Added default to keep local copy is none exists.
   echo '               <not specified>  Local pole. and ut1. will be kept, usno used if not present'
   echo ' '
   echo ' EXAMPLES'                               
   echo '    links.arc J2000 1994 '
   echo '    links.arc J2000 1994 vlbi'
   echo '  '
   echo '  '             
   echo " NOTE: If the link to 'nbody' is empty, indicating the absence of a PEP or JPL ephemeris, "
   echo "       then arc, model, and ngstot will revert to the old luntab., soltab., and nutabl. "
   echo "       Otherwise, the nbody ephemeris will be read and the nutation computed by calling "
   echo "       gamit/lib/MHB_2000 (copy of routine of same name in kf/gen_util/precess.f). "
   echo " " 
#   echo ' CURRENT LINKS:'      (now skipped: see below)
#   echo '  '

else      
# create standard links in a working directory
#

# Variable yr[1] = 4 char yr, yr[2] = 2 char yr, yr[3] = 1 char yr
   set yr = `sh_year -year $2`

# MOD TAH 200326: Don't delete yet.  Use -f in ln later.
#  /bin/rm -f ut1.
#  /bin/rm -f pole.
   /bin/rm -f nutabl.
   /bin/rm -f soltab.
   /bin/rm -f luntab.    
# MOD TAH 200223: Don't remove local cooy or link                 
#  /bin/rm -f nbody
   /bin/rm -f leap.sec
   /bin/rm -f leap.sec
# MOD TAH 190627: Only link svnav.dat if a local copy does
#  not exist.  Allows use of igs_metadata.snx file. 
#   /bin/rm -f svnav.dat

# These needed for sh_sp3fit (but not arc)
#   /bin/rm -f svs_exclude.dat   # MOD MAF 191206: Not necessary to remove and relink svs_exclude.dat if user wants to use local copy
   /bin/rm -f otlcmc.dat

  
# MOD TAH 000801: Changed the year argument to 4-digit version (instead
#  of original 2-digit version.
# MOD XXX: Someone commented out these lines with no note.  Change
#  OK with nbody linked.  Propogated this change to other link.xxx scripts. 
#  Change probably made 20190118 11:12.  Converted to test as well.
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
   ln -s  ~/gg/tables/leap.sec              leap.sec

# MOD TAH 190627: Only link svnav.dat if a local copy does
#  not exist.  Allows use of igs_metadata.snx file. 
# MOD TAH 191202: Try to link to ../tables/svnav.dat since in normal
#  processing this will be one used in sh_gamit.
   if( ! -e svnav.dat ) then
#     See if ../tables/version exists
      if( -e ../tables/svnav.dat ) then
         ln -s ../tables/svnav.dat .
      else
         ln -s  ~/gg/tables/svnav.dat .
      endif
   else
# MOD TAH 191202: if ../tables/svnav.dat exists test consistency
      if( -e ../tables/svnav.dat ) then
          set ndiff = `diff svnav.dat ../tables/svnav.dat | wc -l`
          if( $ndiff != 0 ) then
             echo -n "*** WARNING *** local svnav.dat differs from ../tables version. "
             echo " This may cause problems in sh_gamit"
             head -5 svnav.dat ../tables/svnav.dat
          endif
      endif 
   endif
# MOD MAF 200212: Duplicated svnav.dat if-tests, above, for svs_exclude.dat
   if ( ! -e svs_exclude.dat ) then
     if ( -e ../tables/svs_exclude.dat ) then
       ln -s ../tables/svs_exclude.dat .
     else
       ln -s ~/gg/tables/svs_exclude.dat .
     endif
   else
     if ( -e ../tables/svs_exclude.dat ) then
       set ndiff = `diff svs_exclude.dat ../tables/svs_exclude.dat | wc -l`
       if ( $ndiff != 0 ) then
         echo -n "*** WARNING *** local svs_exclude.dat differs from ../tables version. "
       endif
     endif
   endif
   ln -s  ~/gg/tables/otlcmc.dat            otlcmc.dat

# MOD TAH 180620: Add otides.dat for arc integration: Implemented active 190606.
# MOD TAH 200503: Checked to see if local valid copy.  If so than link to gg/tables
#  not made.
   if( ! -e otides.dat ) then
      ln -sf ~/gg/tables/otides.dat            otides.dat
   endif

# MOD RWK 130326: Link the sestbl. only if there is no local file
   if( -e sestbl. ) then
      echo "sestbl. exists--do not craeate a link"
   else
      ln -s  ../tables/sestbl.               sestbl.
   endif 
     
   if ($#argv < 3) then
# MOD TAH 20026: If local copy exist; just keep these
     if( ! -e ut1. || ! -e pole. ) then 
        echo "No local pole. and ut1. files. Using ~/gg/tables/usno"
        ln -sf  ~/gg/tables/ut1.usno ut1.    
        ln -sf  ~/gg/tables/pole.usno pole.  
     endif 
   else 
#    MOD TAH 190625: If local copy exists for pole.$3 and ut1.$3 use these
#    else try tables version
#    MOD TAH 200326: add f to ln options to force link and replace local copy.
     if( -e ut1.$3 && -e pole.$3 ) then
        ln -sf  ut1.$3    ut1.    
        ln -sf  pole.$3   pole.  
     else 
        ln -sf  ~/gg/tables/ut1.$3    ut1.    
        ln -sf  ~/gg/tables/pole.$3   pole.   
        if( ! -e  ~/gg/tables/ut1.$3 || ! -e  ~/gg/tables/pole.$3 ) then
           echo "Requested $3 pole. and ut1. file don't exisit: Reverting to usno"
           ln -sf  ~/gg/tables/ut1.usno ut1.    
           ln -sf  ~/gg/tables/pole.usno pole.   
        endif
     endif
   endif      # No third argument

#   
endif
#ls -l `find . -type l -print`
# MOD TAH 971215: Only do ls -l if there are links in the directory
# MOD TAH/RWK 080423: Skip this because of the large number of obs-file links
# if ( `find  . -type l -print | wc | awk '{print $1}'` > 0 ) ls -l `find  . -type l -print`

