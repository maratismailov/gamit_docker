#!/bin/csh -f      

#  Script sh_links.arc to create GAMIT links needed to run arc and/or sh_sp3fit in the current directory (e.g. /igs)
#  rwk 151231 from link.arc, modified for more flexible calling arguments and to use 'topt', like sh_links.tables
if ($#argv == 0) then
   echo ' '
   echo '    sh_links.arc creates standard GAMIT links necessary to run arc and/or sh_sp3fit in the current directory'
   echo ' ' 
   echo ' FORMAT' 
   echo ' '
   echo '    sh_links.arc -year [year] -frame [frame] -eop [eop] -topt [topt]'
   echo ' '
   echo ' REQUIRED'               
   echo ' '  
   echo '    year     = 4-digit year for nutabl, soltab, luntab'  
   echo ' '
   echo ' OPTIONAL'                                               
   echo ' ' 
   echo '    frame    = B1950 or J2000  for soltab, luntab; default J2000'
   echo ' '        
   echo '    eop  = bull_b_      IERS Bulletin B      (1-day interval)   valid  1 Jan 1992 - current '
   echo '           bull_b_5day  IERS Bulletin B_5day (5-day interval)   valid  1 Sep 1987 - 25 Nov 1993 '
   echo '           usno         IERS Bulletin A      (rapid-service)    valid  1 Jan 1992 - current '  
   echo '           usnd         IERS Bulletin A 1-day with predictions  '
   echo '           vlbi         IRIS                                    valid  5 Jan 1984 - 28 Jun 95'
   echo '           default usno '
   echo ' '    
   echo '    -topt defines the list of table files to be relinked if they already exist'
   echo '       none all ut1. pole. nutabl. soltab. luntab. leap.sec svnav.dat svs_exclude.dat otlcmc.dat '
   echo '       [Default = none]' 
   echo ' '
   echo ' EXAMPLES'                               
   echo '    sh_links.arc -year 2014 '
   echo '    sh_links.arc -year 2014 -eop usnd'
   echo '    sh_links.arc -year 2014 -topt ut1. pole. '
   echo '  '
#   echo ' CURRENT LINKS:'      (now skipped: see below) '
#   echo '  '

else      
# create standard links in a working directory
#

# Set the defaults
set frame = 'J2000'
set eop   = 'usno'
set topts = ''

# Decode the command line
  while ($#argv > 0 )
    set input = ( $argv )
    switch($input[1])
      case -y*:
        set year = $input[2]
# variable yr[1] = 4 char yr, yr[2] = 2 char yr, yr[3] = 1 char yr
        set yr = `sh_year -year $2`
      breaksw       
      case -e*:
        set eop = $input[2]  
      breaksw 
      case -t*:
        set topts =  (`echo $argv | cut -d- -f2`); shift topts 
      breaksw 
    endsw
    if ( $#argv > 0 ) shift argv
  end
alldone:
         
# Remove broken links
   if ( ! -e nutabl.	       )    \rm -r   nutabl.	       >& /dev/null
   if ( ! -e soltab.	       )    \rm -r   soltab.	       >& /dev/null
   if ( ! -e luntab.	       )    \rm -r   luntab.	       >& /dev/null
   if ( ! -e leap.sec	       )    \rm -r   leap.sec	       >& /dev/null  
   if ( ! -e ut1.	           )    \rm -r   ut1.	           >& /dev/null
   if ( ! -e pole.	           )    \rm -r   pole.	           >& /dev/null
   if ( ! -e svnav.dat        )    \rm -r   svnav.dat	       >& /dev/null
   if ( ! -e svs_exclude.dat  )    \rm -r   svs_exclude.dat   >& /dev/null
   if ( ! -e otlcmc.dat	   )    \rm -r   otlcmc.dat        >& /dev/null
   if ( ! -e core             )    \rm -r   core              >& /dev/null
                   
# Remove links listed in topts.

  foreach i ( $topts)
    if ( $i == all ) then
      \rm -r   nutabl.	     >& /dev/null
      \rm -r   soltab.	     >& /dev/null
      \rm -r   luntab.	     >& /dev/null
      \rm -r   leap.sec	     >& /dev/null      
      \rm -r   ut1.	     >& /dev/null
      \rm -r   pole.	     >& /dev/null
      \rm -r   svnav.dat     >& /dev/null
      \rm -r   svs_exclude.dat   >& /dev/null
      \rm -r   otlcmc.dat    >& /dev/null
      \rm -r   eq_rename	     >& /dev/null     
      \rm -r   core 	     >& /dev/null
     else
      \rm -r   $i >& /dev/null
     endif
  end

# Remove links for year-long directories that don't match the request  
  if ( $yr[1] != `ls -l nutabl.  | awk '{print $NF}' | awk -F"." '{print $NF}'` )    \rm  nutabl.  >&! /dev/null   
  if ( $yr[1] != `ls -l soltab.  | awk '{print $NF}' | awk -F"." '{print $(NF-1)}'` )  \rm  soltab.  >&! /dev/null   
  if ( $yr[1] != `ls -l luntab.  | awk '{print $NF}' | awk -F"." '{print $(NF-1)}'` )  \rm  luntab.  >&! /dev/null   

# Make new links if necessary

  ln -s  ~/gg/tables/nutabl.${yr[1]}        nutabl.           >& /dev/null
  ln -s  ~/gg/tables/soltab.${yr[1]}.$frame soltab.           >& /dev/null
  ln -s  ~/gg/tables/luntab.${yr[1]}.$frame luntab.           >& /dev/null
  ln -s  ~/gg/tables/leap.sec               leap.sec          >& /dev/null 
  ln -s  ~/gg/tables/ut1.$eop                ut1.             >& /dev/null 
  ln -s  ~/gg/tables/pole.$eop               pole.            >& /dev/null 
  ln -s  ~/gg/tables/svnav.dat              svnav.dat         >& /dev/null
  ln -s  ~/gg/tables/svs_exclude.dat        svs_exclude.dat   >& /dev/null
  ln -s  ~/gg/tables/otlcmc.dat             otlcmc.dat        >& /dev/null
  ln -s  ~/gg/tables/core                   core              >& /dev/null

# MOD RWK 130326: Link the sestbl. only if there is no local file
   if( -e sestbl. ) then
      echo "sestbl. exists--do not create a link"
   else
      ln -s  ../tables/sestbl.                 sestbl.
   endif 
#   
endif
#ls -l `find . -type l -print`
# MOD TAH 971215: Only do ls -l if there are links in the directory
# MOD TAH/RWK 080423: Skip this because of the large number of obs-file links
# if ( `find  . -type l -print | wc | awk '{print $1}'` > 0 ) ls -l `find  . -type l -print`

