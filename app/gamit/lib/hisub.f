      Subroutine hisub( ihi,raw_up,raw_north,raw_east,antcod_in
     .  , htcod,site,iyr,idoy,isessn,offarp,warnings )
c
c PURPOSE: Calculate the vertical distance ('Height-of-Instrument')
c          from a ground monument to the antenna reference point (ARP)
c          given a field measurement from the monument to some specified
c          part of the antenna structure.  Called by MAKEX, FIXDRV,
c          MODEL or XTORX after call to sb station.info.   MODEL later
c          calls get_antpcv.f to obtain the offset between the ARP and
c          the L1 and L2 phase centers.
c
c PARAMETERS :
c          IN:  ihi       = unit number of hi.dat                          I*4
c               raw_up    = Measured vertical offset of antenna            R*8
c               raw_north = Measured north offset of antenna               R*8
c               raw_east  = Measured east offset of antenna                R*8
c               antcod_in = Antenna type                                   CH*6
c               htcod    = Type of antenna-height measurement             CH*5
c               site      = Site code                                      CH*4
c               iyr       = Year of session                                I*4
c               idoy      = Day of year of session                         I*4
c               isessn    = Session number of session                      I*4   
c               warnings  = if true, print warnings (T for processing, 
c                           F for mstinf2 and hfupd                        L*4
c
c         OUT:  offarp (3) = offset of ARP from monument (U N E)           R*8
c
Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego. 2005.   All rights reserved.
c
c Subroutines called : REPORT_STAT
c                      
c     Written   Oct 2005 by R. King, based on lib/hisub (see history)

c     This is a new version of subroutine lib/hisub.f, revised to read the vertical
c     and horizontal offsets from a table (hi.dat) rather than from values imbedded
c     in the routine.  Diagrams for the antennas are kept in the red 'Antennas' notebook
c     in MIT 54-611 or 54-612.  A source for many of these is the NOAA Geodynamics
c     Research Division web page: http://www.grdl.noaa.gov/GRD/GPS

c     If an antenna is not found in hi.dat, a warning will be issued and HI measurement
c     is assumed to be to the ARP.

      implicit none

      include '../includes/dimpar.h'

      logical warnings,found,eof

      character*4  site,upperc
      character*5  htcod,htcod1
      character*6  antcod_in,antcod,antcod1 
      character*20 anttyp
      character*80 prog_name
      character*256 message,line
           
      integer*4  iyr,idoy,isessn,ihi,len,rcpar,ioerr
                        
      real*8  raw_up,raw_north,raw_east,horiz,vert,offarp(3)

c Get the calling module name for report_stat   

      len =  rcpar(0,prog_name)              
      

c*      print *,'DEBUG ihi raw_up antcod_in htcod site warnings ',
c*     .    ihi,raw_up,antcod_in,htcod,site,warnings
                               
c Trap bad entries before calling julday

      if( iyr.lt.1980.or.iyr.gt.2050 ) then   
           write(message,'(a,i6,a,a4)') 'Unreasonable year ',iyr
     .       ,' for site ',site   
           call report_stat('FATAL',prog_name,'lib/hisub',' '
     .                     ,message,0)  
      endif
      

c Get the standard name for antennas from the alias table

      call ant_alias(antcod_in,antcod) 

c For certain related antennas, hi.dat has just one of them (these
c are not necesarily electrically the same, and hence are not true aliases
                               
      if( antcod_in.eq.'TI_100' .or. antcod_in.eq.'TI2000' ) 
     .      antcod = 'TI4000'
      if( antcod_in.eq.'LC299I' .or. antcod_in.eq.'LC399I' .or.
     .    antcod_in.eq.'LC_201' .or. antcod_in.eq.'LC202N' .or.
     .    antcod_in.eq.'LC202G' .or. antcod_in.eq.'LC302G' .or.
     .    antcod_in.eq.'LC302N' ) antcod = 'LC_201'
      if( antcod_in.eq.'LC303R' ) antcod = 'LC_303' 
      if (antcod_in.eq.'ASHG3R' .or. antcod_in.eq.'ASHGDA' .or.
     .    antcod_in.eq.'ATGE3B' .or. antcod_in.eq.'ATGE33' .or.
     .    antcod_in.eq.'ATGE3A' .or. antcod_in.eq.'ATGEA1' ) 
     .       antcod = 'ASHGD3'
      if( antcod_in.eq.'ASHDMC' .or. antcod_in.eq.'ATDMCB' .or.
     .    antcod_in.eq.'ATDMCD' .or. antcod_in.eq.'ATDMCF' .or.
     .    antcod_in.eq.'ATDMRA' .or. antcod_in.eq.'ATDMRE' .or.
     .    antcod_in.eq.'ATDMRE' .or. antcod_in.eq.'ATDMRB' .or.
     .    antcod_in.eq.'ATDM1C' .or. antcod_in.eq.'ATDM1E' .or.
     .    antcod_in.eq.'ATDMGG' .or. antcod_in.eq.'ASHGG3' .or.
     .    antcod_in.eq.'ATGGD1' .or. antcod_in.eq.'ATGGDA' .or.
     .    antcod_in.eq.'ATDMG2' ) antcod = 'ASHDMG'  
      if( antcod_in.eq.'ATMR2B' .or. antcod_in.eq.'ATMR2C' ) 
     .      antcod = 'ATMR2A'

c*      print *,'antcod_in antcod ',antcod_in,antcod
                    
c Translate DHPAB to DHARP except for the AOAD/M_B choke ring, which 
c needs a true pre-amp base, above the ARP, which is the bottom of the choke rings

      if( antcod_in.ne.'ROGAOA' .and. antcod_in.ne.'TRBMTA' .and.
     .    htcod.eq.'DHPAB' ) htcod = 'DHARP'

c If the measurement type is direct height to the ARP or L1-phase center transfer 
c the values and exit

      if( htcod.eq.'DHARP'.or.htcod.eq.'L1PHC' ) then

        offarp(1) = raw_up

c Otherwise, read the HI file to get the ARP offset (file opened by each module)

      else
        
        rewind(ihi,iostat=ioerr)
        if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/hisub'
     .      ,' ','Error rewinding antenna HI file',ioerr)
c*        print *,'DEBUG rewound ihi'

        found = .false.
        eof = .false.
        do while( .not.found .and. .not.eof )
          read( ihi,'(a)',iostat=ioerr ) line
          if( ioerr.eq.-1 ) then
            eof = .true.
          elseif (ioerr.ne.0 ) then
              call report_stat('FATAL',prog_name,'lib/hisub',' '
     .          ,'Error reading line from antenna HI file',ioerr)
          elseif ( line(1:1).ne.' ' .or.line(1:7).eq.'       ') then
            continue
          else
            read(line,'(1x,a20,1x,a6,1x,a5,2f7.2)',iostat=ioerr)
     .          anttyp,antcod1,htcod1,vert,horiz 
            if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .             ,'lib/hisub',' '
     .             ,'Error reading values from antenna HI file',ioerr)
            if( antcod1.eq.antcod .and. htcod1.eq.htcod ) then
               found = .true. 
               if( horiz.ne.0. ) then
c                check the validity of the slant height
                 if( raw_up.le.horiz ) then
                   write(message,'(a,f5.3,a,f5.3,a,a4)') 
     .                'Slant height (',raw_up
     .                ,') less than antenna radius (',horiz
     .                ,') for ',site
                   call report_stat('FATAL',prog_name,'lib/hisub',' '
     .               ,message,0)
                 endif
               endif             
               offarp(1) = dsqrt((raw_up**2)-(horiz**2)) - vert
c              print *,' HISUB raw horiz vert offarp ',
c     .                raw_up,horiz,vert,offarp(1)    
            endif
          endif
        enddo
                       
c*** change this: too dangerous to allow a measurement that might be only centimeters off
c        if( .not.found .and.warnings ) then
c          write(message,'(1X,a,a6,1x,a5,a,a4,1x,i4,1X,i3,1x,i2,a)')
c     .        'Antenna code (',antcod,htcod,') for :',upperc(site)
c     .       ,iyr,idoy,isessn,' not in hi.dat; assume ARP'   
c           call report_stat('WARNING',prog_name,'lib/hisub',' ',message
c     .                              ,0)
        if( .not.found  ) then
          write(message,'(1X,a,a6,1x,a5,a,a4,1x,i4,1X,i3,1x,i2,a)')
     .        'Antenna code (',antcod,htcod,') for :',upperc(site)
     .       ,iyr,idoy,isessn,' not in hi.dat'   
           call report_stat('FATAL',prog_name,'lib/hisub',' ',message
     .                              ,0)
        endif
c*** temporary to avoid compiler warning of unused variable 'warnings'
        if( warnings ) then 
           continue
        endif    
c***  
      
      endif

      offarp(2) = raw_north
      offarp(3) = raw_east

      return
      end

                               
