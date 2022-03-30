      Subroutine get_antinfo(debug) 

c     Get the antenna offsets and PCVs and print them to the p-file; called by setup 
c     at the beginning of a session, then again by update_coords if there is a change 
c     mid-session. PCV values are not stored here but rather 'saved' in get_antpcv,
c     called for every observation.  

c     R. King 100205  from code originally in setup.f; last modified 100906
     

c     Input
c       jd0      Julian day for observations   (model.h)
c            In common /lunits/ in model.h :   
c       ipnrt    unit number for print (p-file)
c       ipcv     unit number for antmod.dat file 
c            In common /rcvantcom/ of model.h :
c       anttyp  20-character antenna name including radome (for c- and h-file headers) 
c           (radome may be changed on output if not available on ANTEX file)

c     Output     
c       pcvminelev   Minimum elevation of PCV model (used to omit data in processing)  
c             In common /rcvantcom/ of model.h:
c       kstarts(5) kstops(5) : yr doy hr min sec of start of current station.info entry
c       kstartr(5) kstopr(5) : yr doy hr min sec of start/stop of current eq/rename   
c       anttyp               : 20-character antenna name including radome (for c- and h-file headers) 
c       radome_in            : 5-character radome code for measurements
c       atnmod_in            : antenna PCV model requested (AZEL ELEV or NONE)
c       antmod               : antenna model available and used (AZEL ELEV or NONE)
c       antmod_snx           : 10-character SINEX code for antenna model
c       offarp(3)            : offset of antenna ARP from monument (U N E) (m)
c       offl1(3) offl2(3)    : mean offset of antenna phase center from monument (U N E) (m)   
c       pcvminelev           : minimum elevation in PCV model from ANTEX file (deg)
c       elvtabl1(maxel),elvtabl2(maxel) : Azimuth-averaged PCVs (90-0 deg) (mm)
c       tablel1(maxel,maxaz),tablel2(maxel,maxaz): PCVs by elevation (90-0 deg) and 
c                                                azimuth (0-360) (mm)
c MOD TAH 200205: Changed pcoffl1 and pcoffl2 from get_antpcv to reflect the 
c       effects of rotated antenna.  No change to this routine needed. antdaz
c       in model.h and read from station.info (0 if no AntDAZ column in station.info
       
      implicit none  
            
      include '../includes/dimpar.h'    
      include '../includes/units.h'
      include '../includes/model.h'

      integer*4 mchkey,indx,minelv,ivalue,ifreq
     .        , ioerr,i,j,jj
      
      real*8 sign(2),rvalue,az,el,pcoffl1(3),pcoffl2(3),corrl1,corrl2
                      
      character*4 cvalue  
      character*16 amodel 
      character*256 message,line   

      logical found, eof, fcheck, radome_sub, debug
          
c  Read the ANTEX file to see what model is available for this antenna
                       
      el = 90. 
      az = 0.    
      newant = .true.  
      if( debug ) 
     .  print *,'GET_ANTINFO newant anttyp ',newant,anttyp
      call get_antpcv( jd0,el,az,debug
     .               , pcoffl1,pcoffl2,corrl1,corrl2 )

            
c  Set and print models for what was actually available
           

      if( debug) print *,'GET_ANTINFO radome_in anttyp '
     .        , radome_in, anttyp
      if( radome_in(1:4).ne.anttyp(17:20) )  
     .   write(iprnt,'(a,a5,a)') 
     .' **WARNING: requested radome ',radome_in
     .    ,' not found in antmod.dat, using NONE'
c     Set radome code in anttyp for c-file and h-file headers
c**      if( radome_in(1:4).eq.'UNKN' .or. radome_in(1:4).eq.'NONE' ) then
c**        anttyp(16:16) = ' ' 
      if( radome_in(1:4).ne.anttyp(17:20) ) then
        anttyp(16:16) = '-'
      else
        anttyp(16:16) = '+'
      endif  
      write(iprnt,'(a,a20)') 
     .    'After reading ANTEX file anttyp now ',anttyp
* MOD TAH 200512: Modified format
      write(iprnt,'(a,3(1x,a4),a,a4,2x,a10)')
     .  'Phase center variation PCV model for ANTEX frequencies ',atxfrq
     .   ,' is ',antmod,antmod_snx
      write(iprnt,'(1x,a20,a,3(f6.1,1x))')  anttyp
     . ,' L1 phase center offsets from ARP from antmod.dat [UNE] (mm): '
     . ,pcoffl1(3),pcoffl1(1),pcoffl1(2)
      write(iprnt,'(1x,a20,a,3(f6.1,1x))')  anttyp
     . ,' L2 phase center offsets from ARP from antmod.dat [UNE] (mm): '
     . ,pcoffl2(3),pcoffl2(1),pcoffl2(2)
      sign(1) = 1.d0
      sign(2) = 1.d0  
      minelv = idint(pcvminelev)      
      if( antmod_in.eq.'ELEV' ) then 
        write(iprnt,'(a)') 'L1 PCVs (mm)' 
        write(iprnt,'(a,45(i5,1x))') '     Elev:',(i,i=minelv,90,5)  
        write(iprnt,'(10x,19(f5.1,1x))') 
     .     (elvtabl1(i)*sign(1),i=19,1,-1)
        write(iprnt,'(a)') 'L2 PCVs (mm)' 
        write(iprnt,'(a,45(i5,1x))') '     Elev:',(i,i=minelv,90,5)  
        write(iprnt,'(10x,19(f5.1,1x))') 
     .     (elvtabl2(i)*sign(1),i=19,1,-1)
      elseif( antmod_in(1:2).eq.'AZ') then  
        write(iprnt,'(a)') 'L1 PCVs (mm)' 
        write(iprnt,'(a,45(i5,1x))') 'Az   Elev:',(i,i=minelv,90,5) 
        do j=1,naz
          write(iprnt,'(i3,7x,19(f5.1,1x))') 
     .      idint((j-1)*dazi),(tablel1(i,j)*sign(1),i=19,1,-1)
        enddo
        write(iprnt,'(a)') 'L2 PCVs (mm)' 
        write(iprnt,'(a,45(i5,1x))') 'Az   Elev:',(i,i=minelv,90,5) 
        do j=1,naz
          write(iprnt,'(i3,7x,19(f5.1,1x))') 
     .      idint((j-1)*dazi),(tablel2(i,j)*sign(1),i=19,1,-1)
        enddo    
      endif

c  Calculate the phase-center to monument offset (arp2pc values returned in mm)
                      
      if( debug ) 
     .   print *,'GET_ANTINFO offarp pcoffl1 pcoffl1 '
     .                      , offarp,pcoffl1,pcoffl1
c     GAMIT expects U N E in meters, but ANTEX is N E U
      offl1(1) = offarp(1) + pcoffl1(3)/1.d3
      offl1(2) = offarp(2) + pcoffl1(1)/1.d3
      offl1(3) = offarp(3) + pcoffl1(2)/1.d3            
      offl2(1) = offarp(1) + pcoffl2(3)/1.d3
      offl2(2) = offarp(2) + pcoffl2(1)/1.d3
      offl2(3) = offarp(3) + pcoffl2(2)/1.d3
      write(iprnt,'(a,3(f8.4,1x))')  
     .  ' Mean L1 offset from monument [UNE] (m): ', offl1
      write(iprnt,'(a,3(f8.4,1x))') 
     .  ' Mean L2 offset from monument [UNE] (m): ', offl2
        write(iprnt,'(1x)')
 

c     Check minimum elevation angle used by AUTCLN and warn if below table cutoff   

      if( antmod(1:2).eq.'EL'.or.antmod_in(1:2).eq.'AZ' ) then   
        if( fcheck('autcln.cmd') ) then  
          open(unit=43,file='autcln.cmd',status='old') 
          found = .false.              
          eof = .false.
          do while ( .not.found .and. .not.eof ) 
            read(43,'(a)',iostat=ioerr) line  
            if( ioerr.eq.-1 ) then
              eof =.true.
            elseif( ioerr.lt.0 ) then
              call report_stat('FATAL','MODEL','get_antinfo',' '
     .        ,'Error reading autcln.cmd file',ioerr)
            endif                            
            if( line(1:1).eq.' '.or. line(1:4).eq.'POST') then
              if( mchkey(line,'site_param',256,10).gt.0 .or.
     .            mchkey(line,'SITE_PARAM',256,10).gt.0 ) then
                found = .true.    
                indx = 2
c               read the elevation cutff (4th argument)
                call read_line(line,indx,'CH',ioerr,rvalue,cvalue)
                call read_line(line,indx,'CH',ioerr,rvalue,cvalue)
                call read_line(line,indx,'I4',ioerr,ivalue,cvalue)
                call read_line(line,indx,'I4',ioerr,ivalue,cvalue)
                if( ioerr.ne.0 ) 
     .            call report_stat('FATAL','MODEL','get_antinfo',' '
     .          ,'Too few arguments in autcln site_param command',ioerr) 
                found = .true.
              endif                      
c             autcln default is 0 deg cutoff
              if( .not.found ) ivalue = 0 
            endif
          enddo                     
          if( (pcvminelev-.01d0).gt.dfloat(ivalue)) then
            write(message,'(a,f4.1,a,i2,a,f4.1,a)') 
     .        'PCV table min elv = '
     .        ,pcvminelev,' but autcln cutoff = ',ivalue,', skip obs < '
     .        ,pcvminelev,' deg'
            call report_stat('WARNING','MODEL','get_antinfo',' '
     .                      ,message,0)
            write(iprnt,'(2a)') ' **WARNING: ',message
          endif
          close(43)
        else
          call report_stat('WARNING','MODEL','get_antinfo',' '
     . ,'autcln.cmd file missing, cannot check min elev of PCV table',0)
        endif
      endif
c     convert the autcln minimum elevation to a multiple of 5 for printout
      minelv = dint(pcvminelev)
      minelv = minelv - mod(minelv,5)  
c      rwk 061130: I can't remember why this is here, but it needs to check ivalue, not pcvminelev
cd      print *,'DEBUG pcvminelev  minelv ',pcvminelev,minelv
cd      if (minelv.lt.0 .or. minelv.gt.40 ) then
cd        write(message,'(a,i3,a)') 'Minimum elevation from PCV table (='
cd     .       ,minelv,') is unreasonable in p-file summary'
cd        call report_stat('FATAL','MODEL','get_antinfo',' ',message,0)
cd      endif
                      
C      stop
      return
      end



