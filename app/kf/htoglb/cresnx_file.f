CTITLE CRESNX_FILE
 
      subroutine cresnx_file (unit, unitc)
 
      implicit none

*     Routine to decode the file  blocks from SINEX:
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
*   unitc       - Commnts file unit number
 
      integer*4 unit, unitc
 
 
* LOCAL VARIABLES
 
*     Copy over the file/REFERENCE lines
      write(unit,'(a)') '+FILE/REFERENCE'
      call cp_comments(unit, unitc, 'FILE/REFERENCE')
      write(unit,'(a)') '-FILE/REFERENCE'
*     Copy over the file/COMMENT lines
      write(unit,'(a)') '+FILE/COMMENT'
      call cp_comments(unit, unitc, 'FILE/COMMENT')
      call cre_file_com(unit)
      write(unit,'(a)') '-FILE/COMMENT'

****  Now do the FILE

****  Thats all
      return
      end

CTITLE CRE_FILE_COM

      subroutine cre_file_com(unit)

      implicit none

*     Routine to write comments based on content of binary files
*
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'

* PASSED VARIABLES
      integer*4 unit

* LOCAL VARIABLES
      logical kbit   ! Tests of bit is set

      character*16 used(2)  ! (Applied/Not Applied string)

* MOD TAH 190628: Added verision comment
      integer*4 unitggv  ! Unit number for GGVersion file
     .,         ierr     ! IOSTAT error
     .,         trimlen  ! Function for length of string.

      character*256 home   ! User's home directory (assumed to have gg link)
     .,         ggver_file ! Full file name for GGVersion file/
     .,         line       ! Line read from file
     .,         ggv        ! Line with last non-blank line (currect version)

      character*8 ut1libr  ! String with UT1 Libration model applied.


***** Write information about models as we know them
      write(unit,100)
 100  format('* MODEL INFORMATION AS AVAILABLE',/,
     .       '* Full detail will depend on origin of sinex file')

      call sub_null(cspeopmod) 
      call sub_null(cetidemod) 
      call sub_null(cotidemod) 
      call sub_null(coatmlmod) 
      call sub_null(catmtdmod) 
      call sub_null(chydromod)
      call sub_null(cgnut)
      call sub_null(cggrav)
      call sub_null(ceradmod)
      call sub_null(cantradmod)
      call sub_null(cdryzen)
      call sub_null(cwetzen)
      call sub_null(cdrymap)
      call sub_null(cwetmap)
      call sub_null(cionsrc)
      call sub_null(cmagfield)

*     For loads see if they are applied
      used(1) = 'Not Applied'   ! Atm load
      used(2) = 'Not Applied'   ! Hydrology load
      ut1libr = ' '
      if( kbit(cload_mod, 9) ) used(1) = 'Applied'  ! Atmload
      if( kbit(cload_mod,25) ) used(1) = 'Applied'  ! Hydrology laod
      if( kbit(cgamit_mod, 5)) ut1libr = 'UT1-LIBR'

      write(unit,150) cspeopmod, ut1libr  , cetidemod, cotidemod,  
     .                coatmlmod, used(1)  , catmtdmod, chydromod,  
     .                used(2),   cgnut,     cggrav,    ceradmod,
     .                cantradmod,cdryzen,   cwetzen,   cdrymap,
     .                cwetmap,   cionsrc,   cmagfield
 150  format('* Short Period EOP  : ',a8,1x,a8,/,
     .       '* Solid Earth Tide  : ',a8,/,
     .       '* Ocean Tide        : ',a8,/,
     .       '* Atmospheric Load  : ',a8,' State ',a,/,
     .       '* Atmos. Tide Load  : ',a8,/,
     .       '* Hydrographic Load : ',a8,' State ',a,/,
     .       '* Nutation model    : ',a8,/,
     .       '* Gravity Field     : ',a8,/,
     .       '* Earth Alebdo      : ',a8,/,
     .       '* Antenna Thrust    : ',a8,/,
     .       '* Atm delay models  : ',4(a4,1x),/,
     .       '* 2nd order Ion     : ',2(a8,1x)   )
*     Report on solid Earth pole tide
      if( kbit(cgamit_mod,19) ) then    ! Pole tide applied
         if( kbit(cgamit_mod,21) ) write(unit,160)
 160     format('* SE Ptide IERS96   : Applied')
         if( kbit(cgamit_mod,23) ) write(unit,165)
 165     format('* SE Ptide IERS2010 : Applied')
         if( kbit(cgamit_mod,26) ) write(unit,170)
 170     format('* SE Ptide IERS2020 : Applied')

*        Report ocean pole tide if known.
         if( kbit(cgamit_mod,24) ) write(unit,180)
 180     format('* OPTide MP IERS96  : Applied')
         if( kbit(cgamit_mod,25) ) write(unit,185)
 185     format('* OPTide MP IERS2010: Applied')
         if( kbit(cgamit_mod,27) ) write(unit,190)
 190     format('* OPTide MP IERS2020: Applied')
      endif

* MOD TAH 190628: Get the GGVerion
      unitggv = 106   ! Unit should be open (105 used for comments)
*     File is $HOME/gg/com/GGVersion
      call getenv('HOME',home)
      ggver_file = trim(home) // '/gg/com/GGVersion'
*     Open and read file to the end (no blank line at end)
      open(unitggv, file=ggver_file, iostat=ierr, status='old')
      ggv = 'No GGVersion available'
      if( ierr.ne.0 ) then
         call report_error('IOSTAT',ierr,'open',ggver_file,0,'glbtosnx')
      else
         do while ( ierr.eq.0 ) 
            read(unitggv,'(a)',iostat=ierr) line
            if( ierr.eq.0 .and. trimlen(line).gt.0 ) ggv = line
         enddo
         close(unitggv)
       endif
       write(unit,220) ggv(1:38)   ! Takes to end of date field.
C      10.71.002 Fri Jun 28 13:01:30 EDT 2019 Sinex generation from back solution, Glonass metadata updates

 220   format('* GAMIT/GLOBK Vers  : ',a)
*
*     Thats all for now
      return
      end




 
 
