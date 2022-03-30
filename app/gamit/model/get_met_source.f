      Subroutine get_met_source( metsrc4 )

c     Decode the b-file hierarchical list of meteorlogical model sources
c     R. King  2007 January 10    

                
      implicit none

      include '../includes/dimpar.h'         
      include '../includes/units.h'
      include '../includes/model.h'
      
     
c  Input

c       metsrc4(5) : Array of values read from b-file                     
c         In common /lunits/ and /filenames/ in model.h
c       metfiln    : Name of RINEX met file (for p-file print)
c       iuw        : Unit number for RINEX met file (if 0, not available)  
c       iprnt      : Unit number of p-file 
c         In common /ufcom/ in model.h:
c       lmet       : Logical indicating presence of u-file values from a
c                    'met' grid or list file
c       lmap       : Logical indicating presenece of u-file values from a
c                    'map' grid or list file 
c       ntmap      : Number of epochs for met values on RINEX file or u-file 
                                            
      character*4 metsrc4(5)  
      integer*4 ioerr

c  Read from RINEX or u-file and stored in common /ufcom/:

c       met_time(maxmet) : Times of met values (decimal day-of-yr)
c       met_val(maxmet,3): Pressure, temperature, and humidity
c       map_name(8)      : 2-character codes for values from 'map' list or grid file
c       map_val(maxmap,8): Mapping functions coefficients (not used here), 
c                          ZHD, and (optionally) pressure and temperature from
c                          u-file (coded only for VMF1 grid or list so far) 

c  Output      
c         In common /metcom/:
c       metsrc   : 3-character primary source (RNX UFL GPT STP PTH)
c       metdef   : 3- character default source if value missing (no longer used)
c       zhd0, pres0, temp0, wetvar0 : Default values of ZHD (m), prsssure (hPa),
c                                     temperature (C), and humidity (0-1) to
c                                     be used if primary source not available 
c       zsrc, prsc, tsrc wsrc  :  Primary source of ZHD, pressure, temperature
c                                 and water vapor for SETUP printout   
c       sealoc  :  1-character flag indicating that pressure and temperature
c                  are for sea-level (S) or local (L)

       
c   Local variables

      integer*4 i,j,imap,nblnk,ngood   
      real*8 humid
      logical blank   

c   External function

      integer*4 match_name


    
c  RINEX has first preference
               
      if(metsrc4(1).eq.'RNX ') then
        if(iuw.gt.0) then     
          call read_metrnx
          sealoc = 'L' 
c         if a RINEX met file requested and available, it is used
c          get nominal values from the first valid entries  
          call check_met(maxmet,ntmet,'P',met_val(1,1),ngood)
          if( ngood.le.0 ) then 
c rwk 160426: this changed to FATAL
            call report_stat('FATAL','MODEL','get_met_source',' '
     .                    ,'No valid pressure values on RINEX file ',0)
c            metsrc = metdef           
c            psrc =   metdef   
          else
            pres0 = met_val(ngood,1)
            metsrc = 'RNX'
            psrc = 'RNX'
          endif
          call check_met(maxmet,ntmet,'T',met_val(1,2),ngood)
          if( ngood.le.0 ) then 
            call report_stat('WARNING','MODEL','get_met_source',' '
     .                ,'No valid temperature values on RINEX file ',0)
c rwk 160426: This changed from 'metdef' to 'GPT'
            tsrc =   'GPT'
          else
            temp0 = met_val(ngood,2)
            tsrc = 'RNX'
          endif
          call check_met(maxmet,ntmet,'H',met_val(1,3),ngood)
          if( ngood.le.0 ) then 
            call report_stat('WARNING','MODEL','get_met_source',' '
     .                    ,'No valid humidty values on RINEX file ',0)
              wsrc = '   '
          else
            wetvar0 = met_val(ngood,3)/100.
            wsrc = 'RNX'
          endif
        endif 
      endif   

c  If not RINEX, look next for a U-file 
           
      if( metsrc.ne.'RNX' .and. 
     .      match_name(5,4,metsrc4,'UFL ').gt.0 ) then
c         u-file source can be a met file or a mapping function file  
        if( lmet ) then     
           call report_stat('FATAL','MODEL','get_met_source',' '
     .         ,'Use of u-file met (not map) sources not yet coded ',0)
        elseif( lmap ) then    
c         ** Assume data are VMF1 and represent values on the orthographic surface
c         ** Eventually set this from model name or value tokens 
c            VMF1 values not yet available on grid:   ZWD Tm WP 
          metsrc = 'UFL'     
          sealoc = 'L'  
          imap = match_name(nmap,2,map_name,'ZH')
          if( imap.gt.0 ) then
c           get the default (usually first value)
            call check_met(maxmap,ntmap,'Z',map_val(1,imap),ngood)
            if( ngood.le.0 ) then 
              call report_stat('WARNING','MODEL','get_met_source',' '
     .                          ,'No valid ZHD values on u-file ',0)
              metsrc = 'GPT'           
              zsrc =   'GPT'     
              write(iprnt,'(a)') '  No valid ZHD values on u-file '
            else
              zhd0 = map_val(ngood,imap)
              zsrc = metsrc
            endif
          endif                                       
          imap = match_name(nmap,2,map_name,'SP')   
          if( imap.gt.0 ) then
c           get the default (usually first value)
            call check_met(maxmap,ntmap,'P',map_val(1,imap),ngood)
            if( ngood.le.0 .and. zsrc.ne.'UFL' ) then 
              call report_stat('WARNING','MODEL','get_met_source',' '
     .                        ,'No valid pressure values on u-file ',0)
              psrc =  metdef
            else
              pres0 = map_val(ngood,imap)
              psrc = metsrc
            endif
          endif                            
          imap = match_name(nmap,2,map_name,'TP') 
          if( imap.gt.0 ) then
c           get the default (usually first value)
            call check_met(maxmap,ntmap,'T',map_val(1,imap),ngood)
            if( ngood.le.0 ) then 
              tsrc =  metdef
            else
              temp0 = map_val(ngood,imap)
              tsrc = metsrc
            endif
          endif      
          if( zsrc.eq.'   ' ) then
            if( psrc.eq.'UFL' ) then
               write(iprnt,'(a)') 
     .           'ZHD not available on u-file, using u-file pressure'
             else
               write(iprnt,'(a)') 
     .       'Neither ZHD nor pressure available on u-file, using GPT'
               psrc = 'GPT'
             endif
          endif
        endif                   
c       get T and WV from GPT since not available from VMF1 grid
        if( tsrc.eq.'   ' ) tsrc = 'GPT'
        if( wsrc.eq.'   ' ) wsrc = 'GPT'
      endif
cd      print *,'GET_MET_SOURCE lmap psrc tsrc sealoc wsrc '
cd     .       ,                lmap,psrc,tsrc,sealoc,wsrc

c  If neither RINEX nor u-file, source is the GPT model or standard values

c      (only one of these should have been specified in the sestbl.)
      if( metsrc.eq.'   ' ) then  
         metsrc = 'GPT'   
         psrc = 'GPT'
         tsrc = 'GPT'
         wsrc = 'GPT'       
         sealoc = 'L'    
      endif    

c The last b-file token is the humidity: use this value throughout the
c session unless it's negative and a wet variable is available from the
c GPT2, RINEX, or U-file
 
      nblnk = 0                                
      blank = .false.
      do while (.not.blank .and.nblnk.le.5 )  
        nblnk = nblnk + 1 
        if( metsrc4(nblnk).eq.'    ' ) blank = .true.
      enddo    
      read(metsrc4(nblnk-1),'(f4.0)',iostat=ioerr) humid
      if(ioerr.ne.0) call report_stat('FATAL','MODEL','get_met_source'
     .  ,' ','Error reading humidity value from b-file ',ioerr)
      if( humid.ge.0 ) then
        wetvar0 = humid/100.d0  
        wsrc = 'INP'
      else      
        if( wsrc.eq.'RNX'.or. wsrc.eq.'UFL' .or.wsrc.eq.'GPT' ) then  
c         continue: wetvar0 already set
          continue
        else
          call report_stat('WARNING','MODEL','get_met_source',' '
     .   ,'Sestbl humidity < 0 but no valid value from RINEX or U-file '
     .      ,0)  
          wetvar0 = abs(humid)/100.
         endif
       endif
                            
      return
      end
