      Subroutine read_antex(luin,antex_vers,anttyp,antsn,gnss,atxfrq,jd
     .                    , found_ant,found_f1,found_f2
     .                    , sinex_code,dazi,zen1,zen2,dzen
     .                    , offsetl1,offsetl2,elvtabl1,elvtabl2
     .                    , tablel1,tablel2 )


c     Input:
c       luin        logical unit number for ANTEX file  
c       antex_vers  format version (real*4)
c       anttyp      20-character code of requested antenna; if blank read the next entry
c       antsn       20-character code of requested serial number (Gnn for GPS, where nn=PRN) 
c       atxfrq(2)   ANTEX codes for requested frequencies (e.g.,G01,G02; E01,E05))
c       gnss        GNSS system for requested calibratons (G R C E I)
c       jd          requested epoch (PEP Julian day) (=0, dummy for converting files)

c     Output 
c       found_ant   requested antenna found      
c       found_f1    requested higher frequency found
c       found_f2    requested lower frequency found 
c       sinex_code  10-character id of PCV model for SINEX file 
c       dazi        increment of aziuth (=0. for elevation-only)
c       zen1        minimum zenith angle (max elev), usually = 0.
c       zen2        maxium zenith angle (min elev), usually 90.
c       dzen        increment of zenith angle of table 
c       offsetl1    N E U or X Y Z in mm for upper frequency (e.g. L1, L9)
c       offsetl2    N E U or X Y Z in mm for lower freqeuncy (e.g. L2, L5) 
c       elvtabl1    values of a non-azimuth-dependent model in mm for L1
c       elvtabl2    values of a non-azimuth-dependent model in mm for L2
c       tablel1     azimuth and elevation dependent values for L1
c       tablel2     azimuth and elevation dependent values for L2

      implicit none
   
      include '../includes/dimpar.h'

* MOD TAH 200209: Removed ifrq(20) since is not used and is the only
*     thing dimensioned with 20 and associated with the fatal commented
*     out below.
      integer*4 luin,jd,nzen,naz,nfrq,istart(5),istop(5)
     .        , jdstart,jdstop,ioerr,julday,len,rcpar,i,j,k  
                               
      real*4 sstart,sstop,antex_vers

      real*8 dazi,zen1,zen2,dzen,offsetl1(3),offsetl2(3)
     .     , elvtabl1(maxel),elvtabl2(maxel)
     .     , tablel1(maxel,maxaz),tablel2(maxel,maxaz)
     .     , offset(3),elvtab(maxel),table(maxel,maxaz)

      character*1 gnss
* MOD TAH 200512: Increased atxfrq from 2 to 3 to allow as secondary
*     low frequeny choice.  Also increased to 4 characters so that
*     4th character can show which one found (* added to name). 
      character*(*) atxfrq(3) 
      character*3 atxfrqx   ! Frequency code read from ANTEX file
      character*5 noazi             
      character*10 sinex_code,antid3,antid4
      character*20 anttyp,antsn,anttypx,antsnx   
      character*80 prog_name
      character*256 line,message  

* MOD TAH 200512: Added feature to have primary and secondary frequency
*     selection of the low frquency e.g., G05 as primary,. G02 as backup
      logical prime_found ! Set true if L2 primary found.

      logical found_ant,found_f1,found_f2,start_freq,debug
c**   logical radome_warning

c**      data radome_warning/.false./
       data debug/.false./
         
c**   Get the calling module name for report_stat 

      len = rcpar(0,prog_name)
                             
c*    Check the version number
     
      if( antex_vers.gt.1.4 ) call report_stat('WARNING',prog_name
     .  ,'lib/read_antex',' ','ANTEX version > 1.4 ',0)
                              
c*    Initialize variables check in logic
c      anttyp = ' '
c      antsn = ' ' 
                     
      if( debug ) print *,'READ_ANTEX antex_vers anttyp antsn jd '
     .                             ,  antex_vers,anttyp,antsn,jd
                          
c**   Read the antenna records (assumes that read_antex_head has been 
c     previously called and that there are no comments between end-of-header
c     and the first antenna entry
                 
      found_ant = .false. 
      found_f1 = .false.
      found_f2 = .false.
      prime_found = .false. 

      sinex_code = ' '
* MOD TAH 041227: Initialize jdstart and jdstop so that they will not be used
*     unless read from antex file
      jdstart = 0
      jdstop  = 0     
c      print *,'READ_ANTEX begin loop anttyp ',anttyp
      do while ( .not.found_ant ) 

        read(luin,'(a)',iostat=ioerr) line   
c        print *,'LINE 1',line
        if( ioerr.eq.-1 ) then   
c           call report_stat('WARNING',prog_name,'lib/read_antex',' '
c     .      ,'EOF on ANTEX file',ioerr )
           found_ant = .false.  
           return 
        elseif( ioerr.ne.0 ) then
          call report_stat('FATAL',prog_name,'lib/read_antex',' '
     .       ,'Error reading START OF ANTENNA line in ANTEX file',ioerr)   
        elseif( line(61:72).ne.'START OF ANT') then 
c         skip comments and other unexpected records
          continue
        else
c         must have START OF ANTENNA: next line should be antenna TYPE  
          read(luin,'(a)',iostat=ioerr) line  
c          print *,'LINE 2',line
          if( ioerr.ne.0.or.line(61:64).ne.'TYPE') 
     .       call report_stat('FATAL',prog_name,'lib/read_antex',' '
     .        ,'Error reading TYPE line',ioerr)
          read(line,'(a20,a20,a10,a10)',iostat=ioerr) anttypx,antsnx
     .         ,antid3,antid4
c         for SV antenna, antid3=svn  antid4=cospar #; not used for rcvr antenna
c         if radome missing, assume NONE 
          if( anttypx(1:5).ne.'BLOCK'.and.anttypx(1:5).ne.'GLONA'.and.
     .        anttypx(1:5).ne.'QZSS '.and.anttypx(1:5).ne.'GALIL'.and.
     .        anttypx(1:5).ne.'BEIDO'.and.anttypx(1:5).ne.'IRNSS'.and.
     .        anttypx(17:20).eq.'    ') then
             anttypx(17:20) = 'NONE' 
          endif
          if( ioerr.ne.0 ) then    
            call report_stat('FATAL',prog_name,'lib/read_antex',' '
     .       ,'Error reading antenna from TYPE line',ioerr)
          else                
cd           if( anttypx(1:3).eq.'AOA'.or.anttypx(7:9).eq.'IIA' ) then
cd              if( debug ) then
cd               print *,'antsn anttype ',antsn,anttype
cd     .            ,' antsnx anttypx antid3 ',antsnx,anttypx,antid3
cd              endif
cd            endif        
            if(debug) print *,'  read anttypx anttyp ',anttypx,anttyp
            if( (anttypx.eq.anttyp .or. anttyp(1:1).eq.' ') .and.
     .          (antsn.eq.antsnx  .or.  antsnx(1:1).eq.' ') ) then       
              found_ant = .true.             
              if( debug ) then
                print *,'antsn anttyp antid3 antid4 '
     .                  ,antsn,anttyp,antid3,antid4
                print *,' found_ant = T'
              endif
c             need to assign input antenna type and SN if it's blank
              if( anttyp(1:1).eq.' ' ) anttyp = anttypx 
              if( antsn(1:1).eq.' ' ) antsn = antsnx
c             skip calibration method line
              read(luin,'(a)') line 
c              print *,'LINE 3 ',line
              read(luin,'(a)') line 
c              print *,'LINE 4 ',line
              if( line(61:64).ne.'DAZI') then 
                 call report_stat('FATAL',prog_name,'lib/read_antex',' '
     .            ,'Missing DAZI line in ANTEX file',0)
              else
                read(line,'(2x,f6.1)') dazi
              endif
              read(luin,'(a)',iostat=ioerr) line 
c              print *,'LINE 5 ',line
              if( ioerr.ne.0.or.line(61:64).ne.'ZEN1' ) then 
                call report_stat('FATAL',prog_name,'lib/read_antex',' '
     .                        ,'Missing ZEN1 line in ANTEX file',ioerr)
              else
                read(line,'(2x,3f6.1)',iostat=ioerr) zen1,zen2,dzen  
                if( ioerr.ne.0 ) then          
                 call report_stat('FATAL',prog_name,'lib/read_antex',' '
     .                    ,'Error decoding ZEN1 ZEN2 DZEN', ioerr)
                else
                  nzen = nint((zen2-zen1)/dzen) + 1
                  if( nzen.gt.maxel ) then
                    write(message,'(a,i3,a,i3,a)') 
     .                 'Number of zenith values (',nzen
     .                ,') exceed maxel (',maxel,')' 
                    call report_stat('FATAL',prog_name,'lib/read_antex'
     .                              ,' ',message,0)
                  endif  
                endif
              endif
              read(luin,'(a)',iostat=ioerr) line   
c              print *,'LINE 6 ',line
              if( ioerr.ne.0.or.line(61:69).ne.'# OF FREQ') then  
                 call report_stat('FATAL',prog_name,'lib/read_antex',' '
     .                   ,'Missing # OF FREQ line in ANTEX file',ioerr)
              else
                read(line,'(i6)') nfrq 
* MOD TAH 200209: Removed the fatal here since nothing seems to be 
*               dimensioned with max nfrq of 20.
C               if( nfrq.gt.20 ) 
C    .            call report_stat('FATAL',prog_name,'lib/read_antex'
C    .                 ,' ','# frequencies > 20',ioerr)
              endif    
              do i=1,nfrq
c               the frequency information may be preceded by validity range,
c               SINEX code, and comments
                start_freq = .false.
                do while(.not.start_freq)
                  read(luin,'(a)',iostat=ioerr) line 
c                  print *,'LINE 7 ',line
                  if( ioerr.ne.0 ) then  
                    call report_stat('FATAL',prog_name,'lib/read_antex' 
     .             ,' ','Error reading START OF FREQUENCY in ANTEX file'
     .                ,ioerr)
                  elseif( line(61:70).eq.'VALID FROM' ) then
                    read(line,'(5i6,f13.7)',iostat=ioerr) istart,sstart 
                    if( ioerr.ne.0 ) then  
                     call report_stat('FATAL',prog_name,'lib/read_antex'
     .              ,' ','Error reading VALID FROM in ANTEX file',ioerr)
                    else
                      jdstart= julday(istart(2),istart(3),istart(1))
* MOD RWK 110818: Since checking dates of validity is done only if start and stop are non-zero,
*                 need to set a finite stop time (1 Jan 2100) in case VALID UNTIL is not preseent
                      jdstop = 2488070              
                       if(debug) print *,'set default jdstop ',jdstop
                    endif
                  elseif( line(61:71).eq.'VALID UNTIL' ) then
                    read(line,'(5i6,f13.7)',iostat=ioerr) istop,sstop
                    if(debug) print *,'read VALID UNTIL istop ',istop
                    if( ioerr.ne.0 ) then 
                     call report_stat('FATAL',prog_name,'lib/read_antex'
     .             ,' ','Error reading VALID UNTIL in ANTEX file',ioerr)
                    else
                      jdstop= julday(istop(2),istop(3),istop(1))
                       if(debug) print *,'set jdstop ',jdstop 
                    endif 
                  elseif( line(61:70).eq.'SINEX CODE' ) then
                    read(line,'(a10)',iostat=ioerr) sinex_code
                  elseif( line(61:67).eq.'COMMENT') then
                     continue
                  elseif( line(61:78).eq.'START OF FREQUENCY' ) then
                    start_freq = .true.                          
                    read(line,'(3x,a3)') atxfrqx
                  else     
                     call report_stat('FATAL',prog_name,'lib/read_antex'
     .                ,' ','Unexpected record before freq start',ioerr)
                  endif 
c               end do on finding start-frequency line
                enddo                                     
c               read the N E U or X Y Z offsets
                read(luin,'(3f10.2)',iostat=ioerr) offset   
                if( debug ) print *,'i atxfrq atxfrqx offset '
     .                             ,i,atxfrq(i),atxfrqx,offset  
                if( ioerr.ne.0) then
                  write(message,'(a,a20,1x,a,a3)') 
     .              'Error reading offsets for ',anttypx
     .                 ,' freq ',atxfrqx 
                  call report_stat('FATAL',prog_name,'lib/read_antex'
     .                  ,' ',message,ioerr)
                endif
c                 read the non-azimuth-dependent pattern
c MOD SCM and MJM 120104: To allow 0.5 degree antex files, need to be
c able to read 181 fields
c               read(luin,'(3x,a5,91f8.2)',iostat=ioerr) 
                read(luin,'(3x,a5,181f8.2)',iostat=ioerr) 
     .                  noazi,(elvtab(j),j=1,nzen) 
c                 print *,'LINE 9 elvtab ',elvtab(1)
                if( ioerr.ne.0.or.noazi.ne.'NOAZI') 
     .             call report_stat('FATAL',prog_name,'lib/read_antex'
     .              ,' ','Error reading NOAZI line in ANTEX file',ioerr)
                if( ioerr.ne.0.or.noazi.ne.'NOAZI') 
     .              call report_stat('FATAL',prog_name,'lib/read_antex'
     .              ,' ','Error reading NOAZI line in ANTEX file',ioerr)
c                 print *,'LINE 9B elvtabl2',elvtabl2(2)
c               if azimuth dependent values, read them now
                if( dazi.ne.0.0d0 ) then
c MOD SCM and MJMi 120104: to allow 0.5 degree antex files, move int()
c                  naz = 360/int(dazi) + 1 
                  naz = int(360/dazi) + 1 
c                  print *,'naz',naz,' nzen',nzen
                  do k=1,naz
c MOD SCM and MJM 120104: To allow 0.5 degree antex files, need to be
c able to read 181 fields
c                    read(luin,'(8x,91f8.2)',iostat=ioerr) 
                    read(luin,'(8x,181f8.2)',iostat=ioerr) 
     .                  (table(j,k),j=1,nzen)  
c                     print *,'table ',(table(j,k),j=1,nzen)
c                     print *,'LINE x1 ',line
                  enddo              
                endif
                read(luin,'(a)',iostat=ioerr) line  
c                print *,'LINE 10 ',line
                if( ioerr.ne.0.or.line(61:71).ne.'END OF FREQ')
     .           call report_stat('FATAL',prog_name,'lib/read_antex',' '
     .            ,'Error reading END OF FREQUENCY in ANTEX file',ioerr)
* MOD TAH: Check just first 3 characters
                if ( atxfrqx.eq.atxfrq(1)(1:3) ) then
                  do j=1,3                               
                    offsetl1(j) = offset(j)        
                  enddo  
                  do j=1,nzen
                    elvtabl1(j) = elvtab(j)
                    if( dazi.ne.0.d0 ) then
                      do k=1,naz
                        tablel1(j,k) = table(j,k)
                      enddo   
                    endif 
                  enddo  
                  found_f1 = .true.
                  if(debug) print *,'found_f1 atxfrqx offsets: '
     .                   ,found_f1,atxfrqx,offsetl1
* MOD TAH 200512: OK, save values if (a) this is primary frequency
*               or (b) primary not found year and this is secondary
*               elseif(atxfrqx.eq.atxfrq(2)) then
                elseif( atxfrqx.eq.atxfrq(2)(1:3) .or. 
     .            ( .not. prime_found .and. 
     .               atxfrqx.eq.atxfrq(3)(1:3) ) ) then
*                 If this is prime frequency, mark as found
                  if( atxfrqx.eq.atxfrq(2)(1:3) ) prime_found = .true.
                  do j=1,3                               
                    offsetl2(j) = offset(j)        
                  enddo                     
                  do j=1,nzen
                    elvtabl2(j) = elvtab(j)
                    if( dazi.ne.0.d0 ) then
                      do k=1,naz
                        tablel2(j,k) = table(j,k)
                      enddo   
                    endif 
                  enddo  
                  found_f2 = .true.                 
                  if(debug) 
     .            print *,'F2 AntexOff ',prime_found, atxfrq, anttyp, 
     .                          nint(offsetL2)
                endif 
c             end do on frequencies
              enddo                    
c             if the requested day is outside the range of validity, look for another entry     
              if( debug ) print *,'jd jdstart jdstop found_ant '
     .                              ,jd,jdstart,jdstop,found_ant
              if( jd.ne.0 .and. jdstart.ne.0 .and.jdstop.ne.0 ) then
                if( jd.lt.jdstart.or.jd.gt.jdstop ) found_ant = .false.
                if( debug ) print *,'found_ant reset ',found_ant
              endif    
c           end if on target antenna    
            endif  
c         end if on valid TYPE 
          endif 
c       end outermost if on read START OF ANTENNA
        endif    

c     end loop on antennas       
        if( found_ant .and. debug ) print *,'anttyp antsn sinex_code '
     .                                     , anttyp,antsn,sinex_code
      enddo 

*MOD TAH 201512: Mark which frequency was found
      if( found_f2 ) then
         if( prime_found ) then
             atxfrq(2)(4:4) = '*'
             atxfrq(3)(4:4) = ' '  
         else
             atxfrq(2)(4:4) = ' '
             atxfrq(3)(4:4) = '*'  
         endif
      endif 
                    
cd      if( anttyp(1:5).eq.'GLONA'.and.antsn(1:3).eq.'R03') then
cd        print *,'offsetl1 offsetl2 )',offsetl1,offsetl2
cd        stop
cd      endif
      if(debug) print *,'Exiting READ_ANTEX'         
      return
      end
      
    


