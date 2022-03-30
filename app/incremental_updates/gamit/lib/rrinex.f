      Subroutine RRINEX ( debug,iflag,rxver,gnss,nobtyp,rxobtyp,iobtypx
     .                  , nprn,isvid,rxtime,iwkn,sow,nepoch
     .                  , dofl1,dofl2,prgl1,prgl2,illi,issi
     .                  , anth,ante,antn
     .                  , fend,ferr )
                                   

c     Read one epoch from a RINEX format data file.

c       iflag = 0  : called by a scanning program, just read the file, don't return observables
c             = 1  : called by makex, return the observables (requires iobtypx to be set)
c
c     Assume that the RINEX file is already opened on unit urinex (in ../includes/makex.h),
c     Assume that the header has already been read.

c     MOD MAF 190404: Increase maximum line length to 1024 and associated maximum
c                     observables allowed (also "MAXOBT" in gamit/includes/makex.h)
c                     to 63, to accommodate long records in RINEX 3 files.
c     MOD RWK 150925: Allow reading of RINEX 3 as well as RINEX 2 values.  All SVs 
c     are read, but instead of returning all SVs and observables, the routine now 
c     returns only the SVs from the GNSS requested and only the 4 observables (2 phase 
c     and 2 PR) selected by MAKEX (and 2 backups for L1 and L2 pseudorange) to write 
c     on the x-file after reading the header record.

      implicit none

      include '../includes/dimpar.h' 
      include '../includes/makex.h'

      logical prior_l1,prior_p1,prior_p2

C       INPUT VALUES
                       
c       Flag for type of read
      integer*4 iflag
c       Set true to turn on debug print
      logical debug
c       GNSS values requested 
      character*1 gnss 
c       Indices in the rxobtyp array to be returned and used 
      integer*4 iobtypx(6)
c       Values from the RINEX header
c       RINEX version
      real*4 rxver
c       number of different types of observable quantities for the requested GNSS 
c       (required to be the same for all SVs for RINEX 2, can vary even within a GNSS for RINEX3
      integer nobtyp  
c       actual number of observables on a line
      integer nobtyp1
c       sat id letter designation and numbers (PRN)
      character*1 asvid(maxchn)
      integer*4 isvid(maxchn)     
c       Time system (usually GPS)
      character*3 rxtime       
c     labels for RINEX 2 and RINEX 3 observable types
      character*3  rxobtyp(maxobt)

c     OUTPUT VALUES
c       number of satellites at this epoch 
c       (for RINEX 2 initially all SVs, then reset to requested GNSS SVs only)
      integer*4 nprn
c       GPS week number for this epoch
      integer*4 iwkn
c       second of GPS week for this epoch
      real*8  sow  
c       number of this epoch
      integer     nepoch
C       L1, L2 doppler phase in cycles
      real*8  dofl1(maxchn), dofl2(maxchn)
C       L1, L2  pseudorange in meters
      real*8  prgl1(maxchn), prgl2(maxchn)
c       updated antenna offsets     
c       Phase and pseudorange signal strength and loss-of-lock flags
      integer*4 issi(maxchn,4),illi(maxchn,4)
      real*8 anth, ante, antn
c       value set to .true. at end of file
      logical fend
c       value set to .true. on error
      logical ferr

c       RINEX observations
      real*8 obs(maxchn,maxobt)

c       function to return day of year
      integer*4  idoy
c       function number of non-blank characters in a line
      integer*4 nblen

c       RINEX signal strength indicator- all observables for requested GNSS only
      integer*4 jssi(maxchn,maxobt)
C       RINEX loss-of-lock indicator -- all observables for  requested GNSS only
      integer*4 jlli(maxchn,maxobt)

c       error code (may be system dependent?)
      integer*4    ioerr
C       error count
      integer*4    nerr

      integer*4  nl  ! Number of lines that the RINEX 2 data
                     ! will occur on ((nobtyp-1)/5+1)
      integer    nel ! End entry for RINEX 2 line being read (increments
                     ! by 5 each line until last line)
     .,           k  ! Loop variable reading lines.

c       RINEX buffer
      character*1024 line

c       counter for requested SVs at each epoch (RINEX 2): used to set output nprn
      integer*4 nprnx

c       miscellaneous variables       
      character*1 recid
      real*8       sec,utcoff,ds,secdif,ltvel
      integer*4    ihr,min,iyr,j,itflag,month,idflag,i,iday,jdoy
     .        ,    jd,il,iend
      logical event                                      

c      functions
      integer*4 julday
      real*8 taiutc 

c       last good time tags
      integer*4 iwkn1
      real*8    sow1
      save iwkn1,sow1

c       variables for error reporting
      integer*4 len,rcpar
      character*80 prog_name
      character*256 message

      data iwkn1/0/
      data sow1 /0.0d0/    
      data ltvel/2.99792458d8/


c                DEFINITIONS OF RINEX DATA FORMATS
c                **********************************
c
c     Epoch/satellite record:

c     Epoch:   year  month  day  hour  minute  second
c     Flag :   0=OK,   1=power failure between previous and current epoch
c     Number of satellites in current epoch
c     Sys/PRN (e.g. G04) for this epoch (RINEX 2 only)
c

c     Observations:
                                          
c     RINEX 2:
c     Observation value
c     LLI:  Loss-of-lock indicator (important for phase only)
c           0 or blank = OK or not known,   1=loss of lock from previous epoch
c     Signal strength:  mapping of C/N ratio to range of 1-9 (to be determined)
c     5 observations per 80-character line, repeated on subsequent lines if > 5
                       
c     RINEX 3:
c     Sys/PRN (e.g. G04, R08) 
c     Observation values, LLI, signal strength in the order given by rxobtyp, no line-length limit

c     NOTE: Currently rxhed and rrinx return values only for the GNSS requested,
c           so values for other systems are skipped in the read

c     A note on conventions.
c     RINEX is defined to be in the PSEUDORANGE CONVENTION
c
c             d(phase)        1        d(range)
c               -----   = + ------     -------
c               d(t)        wavelength  d(t)
c
c     BUT X-files are in the DOPPLER CONVENTION
c
c             d(phase)        1        d(range)
c               -----   = - ------     -------
c               d(t)        wavelength  d(t)
c
c     dofl1 and dofl2 are returned in the DOPPLER convention.
         
c     get the calling program name
      len = rcpar(0,prog_name)


c     IMPLEMENTATION NOTES
c     ********************
                             
c     Case of Geotracer 2200 translated with GEOTRACER GPS Decoder Ver. 2.1
c     There are 9 observation types with L1 and D1 repeated: P1 P2 C1 L1 L2 L1 D1 D2 D1
c     lib/rrinex overwrites the first L1 with the second, but the second one has ISSI = 0
c     so the data are not kept.  Code added when the data are read to check for 
c     the prior existence of 'L1' and to ignore the second occurrence.
          
c     Initialization                                 

      nerr = 0
      ferr = .false.
      fend = .false.   
          
      if (debug  ) then          
        write (*,*) 'RRINEX: rxver iobtypx  nobtyp rxobtyp:'
     .        , rxver,(iobtypx(i),i=1,6),nobtyp,(rxobtyp(i),i=1,nobtyp)
      endif
      if (nepoch .eq. 0) then
         iwkn1 = 0
         sow1 = 0.0d0
      endif

      do i = 1,maxchn
         issi(i,1) = 0
         issi(i,2) = 0
         dofl1(i) = 0.0d0
         dofl2(i) = 0.0d0
         prgl1(i) = 0.0d0
         prgl2(i) = 0.0d0
         isvid(i) = 0
      enddo
c     non-zero values of antenna offsets used to key call to HISUB in MAKEX
      anth = 0.d0
      ante = 0.d0
      antn = 0.d0
         

c     Check dimensions (even though checked in rrxhed

      if( nobtyp.gt.maxobt ) then         
        write(message,'(a,i2,a,i2)') '# obs types =',nobtyp
     .       ,' > maxobt = ',maxobt  
       call report_stat('FATAL',prog_name,'lib/rrinex',' '
     .                 ,message,0)    
      endif
                         
cd      debug = .true.

c     Read a record into a buffer and check for file errors

 100  continue       
        read (unit   = urinex,
     .        iostat = ioerr,
     .        end    = 4000,
     .        fmt    = '(a80)') line
      if( debug ) write(*,*) 'RRINEX read ',line
                          
c     Five cases:
c        1. Error on read:  increment the error count error read the next line 
c        2. Error decoding epoch line:  increment the error count and read the next line
c        3. Good epoch line but number records (SVs) exceeds channels: 
c           Exit the program  
c        4. Good epoch line but with event flag set:  read the event info and then
c           the SV lines
c        5. Normal epoch line: decode and exit
c     After any read error, call rxerr (at the end of this file) to count the errors
c     and print out the line; after 100, stop the run.

c     Case 1: error on 'a' read
      if( ioerr.ne.0 ) then
        call rxerr (urinex,line,ioerr,nerr)
        if(debug) write(*,*) 'Error on read, ioerr ',ioerr
        goto 100
      endif
           
c     Decode the epoch line 

      if( rxver.lt.3.0 ) then   
        read (unit   = line,
     .        fmt    = '(5i3,f11.7,i3,i3)',
     .        iostat = ioerr) iyr,month,iday,ihr,min,sec,idflag,nprn 
        call fix_y2k(iyr)                                              
        if(debug) write(*,*) 'RINEX 2 decode epoch: '
     .     ,iyr,month,iday,ihr,min,sec,idflag,nprn
      else
        read (unit   = line,
     .        fmt    = '(a1,1x,i4,4(1x,i2.2),f11.7,2x,i1,i3)',
     .        iostat = ioerr) 
     .        recid,iyr,month,iday,ihr,min,sec,idflag,nprn
        if(debug) write(*,*) 'RINEX 3 decode epoch: '
     .        , recid,iyr,month,iday,ihr,min,sec,idflag,nprn
      endif

c     Case 2: error decoding the epoch line
      if( ioerr.ne.0 ) then 
        call rxerr(urinex,line,ioerr,nerr) 
        if(debug) write(*,*) 'Error on epoch decode, ioerr ',ioerr
        goto 100 
      endif

c     decoded ok

c     Case 3: too many SV records for GAMIT dimensions
      if( idflag.eq.0. and. nprn.gt.maxchn ) then 
        write(message,'(a,i2,a,i2,a,i4,4i3)')
     .    'Number of SVs on RINEX file (=',nprn
     .    ,') exceeds array dimensions (maxchn=',maxchn
     .    ,') ',iyr,month,iday,ihr,min
          call report_stat('FATAL',prog_name,'lib/rrinex',' '
     .                    ,message,0)    
      endif 
                           
c     Case 4:  RINEX 2 or 3 event information  
      if( idflag.gt.1 ) then
          if(debug) write(*,*) 'Event line nprn (nline):', nprn
          write(message,'(2a)')
     .      'Warning in RRINEX: obs file contains version 2 data--',
     .      'new antenna offset keys HISUB call; all else now ignored'
          call report_stat('WARNING',prog_name,'lib/rrinex',' '
     .                    ,message,0) 
          write(message,'(a,i4,1x,4i3)') 
     .        'Epoch: ',iyr,month,iday,ihr,min   
          call report_stat('WARNING',prog_name,'lib/rrinex',' '
     .                    ,message,0) 
          if( rxver.lt.3.0 ) then
c           RINEX 2 
            do i=1,nprn 
c             here nprn indicates number of event records, not SVs
              read( urinex,iostat=ioerr,end=4000,fmt='(a80)') line
              call report_stat('WARNING',prog_name,'lib/rrinex',' '
     .                        ,line,0)
              if (line(61:80) .eq. 'ANTENNA: DELTA H/E/N') then
                 read (line,'(3f14.4)',iostat=ioerr) anth,ante,antn
              endif
            enddo
            goto 100
          else      
c           RINEX 3
            event = .true.
            do while( event ) 
              read (unit =urinex,iostat=ioerr,end=4000,fmt='(a)') line
              if(debug) write(*,*) 'RINEX 3 event line: ',line
              if(line(1:1).eq.'>' ) then 
                event = .false.
              else
                if (line(61:80) .eq. 'ANTENNA: DELTA H/E/N') then
                  read (line,'(3f14.4)',iostat=ioerr) anth,ante,antn
                  call report_stat('WARNING',prog_name,'lib/rrinex',' '
     .                      ,'EVENT: new antenna offset',0)
                endif
              endif
            enddo
            backspace(urinex,iostat=ioerr ) 
          endif

      endif

c     Case 5: normal epoch record
      if( rxver.lt.3.0 ) call fix_y2k(iyr) 
        
c     Check for valid time tag
* MOD TAH 210101: Updated 2020 to 2100 (also 1985 to 1980 for consistency).
      if (iyr    .ge. 1980 .and. iyr    .le. 2100 .and.
     .    month  .ge.  1   .and. month  .le. 12 .and.
     .    iday   .ge.  1   .and. iday   .le. 31 .and.
     .    ihr    .ge.  0   .and. ihr    .le. 23 .and.
c**  .    min    .ge.  0   .and. min    .le. 59) then
c**    substitute this more lenient test to handle bad Leica converter (min =60)
     .    min    .ge.  0   .and. min    .le. 60) then
c        Assume that this is indeed an epoch record
c           (i.e., we are not lost due to file errors)
            
c        -- RINEX 2 ---             
        if( rxver.lt.3.0) then 
          if( nprn.gt.12 ) then
            iend = 12
          else
            iend = nprn
          endif   
c         decode the first line again, saving the svids
          read(line,'(5i3,f11.7,i3,i3,12(a,i2))',iostat=ioerr)
     .           iyr,month,iday,ihr,min,sec,idflag,nprn
     .           ,(asvid(i),isvid(i),i=1,iend) 
          call fix_y2k(iyr)
          if( debug ) then
             write(*,'(a)') 'RINEX 2 epoch line ',line
             write(*,'(a,i5,4i3,f11.7,i3,i3,12(a,i2))') 'Decoded: '
     .          ,iyr,month,iday,ihr,min,sec,idflag,nprn
     .          ,(asvid(i),isvid(i),i=1,iend)    
          endif
          if( ioerr.ne.0 ) then   
            call rxerr(urinex,line,ioerr,nerr)  
            if(debug) write(*,*) 'Error decoding epoch line',ioerr
            if( ioerr.eq.-1 ) then
              goto 4000
            else
              goto 100 
            endif 
          endif  
c         if there is more than 1 line, read and decode these now
          if(debug) write(*,*) 'nprn ',nprn,' > 12 read another SV line'
          if( nprn.gt.12) then
            nl = int((nprn-1)/12+1)
            do il = 2,nl
              if( nprn.gt. 2*il ) then
                iend = il*12
              else
                iend = nprn
              endif 
cd            print *,'nprn nl il iend ',nprn,nl,il,iend 
              read (unit=urinex,iostat=ioerr,fmt='(a80)') line 
              if(debug) write(*,*) 'Read another SV ID line ',line
              if (ioerr.ne.0 ) then
                call rxerr (urinex,line,ioerr,nerr)
                if(debug) write(*,*) 'Error on SV ID line ioerr ',ioerr
                if( ioerr.eq.-1 ) then
                  goto 4000
                else
                  goto 100 
                endif 
              endif  
              read(line,'(32x,12(a1,i2))',iostat=ioerr) 
     .            (asvid(i),isvid(i),i=(il-1)*12+1,iend)
              if( ioerr.ne.0 ) then 
                call rxerr(urinex,line,ioerr,nerr) 
                if( ioerr.eq.-1 ) then
                  goto 4000
                else
                  goto 100 
                endif 
              endif
            enddo
          endif

c         We have to hope that the observation records are in the order
c         and number expected (RINEX error recovery is difficult)

c         Read and decode the data lines
          nl = int((nobtyp-1)/5)+1  ! Number of lines for data
          do i=1,nprn
            do  k = 1, nl
*             Get end range of data to be read 
              nel = min0(5*k,nobtyp)  ! Need to use explicit name is 
                                      ! min is re-defined to be minutes.
*             Read line from file
              read(unit   = urinex,
     .             iostat = ioerr,
     .             end    = 4000,
     .             fmt    = '(a80)') line 
              if(debug) write(*,*) 'Data i k nl nel ioerr line '
     .                ,i,k,nl,nel,ioerr,line
*             Check for error on read
              if (ioerr .ne. 0) then   
                call rxerr (urinex,line,ioerr,nerr)
c               call this sat bad, because of file error
                isvid(i) = 0
              else    ! Read from file was OK, so extract
*                    ! values from line
                read(line,'(5(f14.3,i1,i1))', iostat  = ioerr)
     .                 (obs(i,j),jlli(i,j),jssi(i,j),j=(k-1)*5+1,nel)
                if(debug) write(*,*) 'Decoded RINEX 2 data line ioerr:'
     .           ,ioerr,(obs(i,j),jlli(i,j),jssi(i,j),j=(k-1)*5+1,nel)
                if(ioerr.ne.0) call rxerr (urinex,line,ioerr,nerr) 
cd                 print *,'obs 1-5 ',(obs(i,j),j=1,5)
              endif
            enddo
            if(debug) write(*,'(a,L,2i3,32(f14.3,2x,2i3))') 
     .          'ferr iprn nobtyp obs '
     .          ,ferr,i,nobtyp,(obs(i,j),jlli(i,j),jssi(i,j),j=1,nobtyp)
          enddo 
c         reduce the observation arrays to the requested GNSS only (reset nprn)
          call get_gnss(gnss,nprn,asvid,isvid,nobtyp,obs,jlli,jssi)


c      ---RINEX 3 --
        else                                 
cd          print *,'rxver 3 ',rxver
c         decode the epoch line
          read (unit   = line,
     .         fmt    = '(a1,1x,i4,4(1x,i2.2),f11.7,2x,i1,i3)',
     .        iostat = ioerr) 
     .        recid,iyr,month,iday,ihr,min,sec,idflag,nprn
          if(recid.ne.'>') then         
            call report_stat('WARNING',prog_name,'lib/rrinex',' ' 
     .            , 'RINEX 3 bogus epoch line, check code logic',0)
            message(1:40) = line(1:40)
            call report_stat('FATAL',prog_name,'lib/rrinex',' '
     .            ,message,0 )
          endif
c         read and decode the data lines - skip the lines that are not the requested GNSS
          nprnx = 0 
          do i=1,nprn
            read(urinex,'(a)',iostat=ioerr,end=4000) line
cd            print *,'read line  ioerr ',ioerr,line
            if(ioerr.ne.0 ) 
     .        call report_stat('FATAL',prog_name,'lib/rrinex',' '
     .                       ,'Error reading data line',ioerr )
            read(line,'(a1,i2.2)',iostat=ioerr) asvid(i),isvid(i)
cd            print *,'decoded ioerr i asvid isvid '
cd     .              ,ioerr,i,asvid(i),isvid(i)
            if( ioerr.ne.0 ) then
              call rxerr (urinex,line,ioerr,nerr)
c             call this sat bad, because of file error
              isvid(i) = 0              
            else   
              if( asvid(i).eq.gnss ) then
                nprnx = nprnx + 1
c               read from file was OK, so extract values from line  
c               hard-wired to a limit of 31 by the length of the line (512 chars) [superseded, below]
c               hard-wired to a limit of 63 by the length of the line (1024 chars)
                if( nobtyp.gt.63 ) call report_stat('FATAL',prog_name
     .              ,'lib/rrinex',' ','Only 63 observables allowed',0 )  
c               zero-out the full array and then read just the values present, which will 
c               vary by SV; since we don't know whether than last observable will have
c               jssi and jlli values, add a little bit to the line length in computing 
c               the number of obserables present
                do j=1,nobtyp
                  obs(nprnx,j) = 0.d0
                  jlli(nprnx,j) = 0
                  jssi(nprnx,j) = 0 
                enddo
                nobtyp1 = nblen(line)/16     
cd                print *,'nobtyp1 ' ,nobtyp1
cd                print *,'line:'
cd                print *,line
                read(line,'(a1,i2.2,31(f14.3,i1,i1))',iostat =ioerr) 
     .             asvid(nprnx),isvid(nprnx)
     .           ,(obs(nprnx,j),jlli(nprnx,j),jssi(nprnx,j),j=1,nobtyp1)
cd                print *,'decode line ioerr ',ioerr
                if(ioerr.ne.0) call rxerr (urinex,line,ioerr,nerr) 
                if( debug ) write(*,'(a1,i2.2,31(f14.3,i1,i1))')
     .             asvid(nprnx),isvid(nprnx)
     .           ,(obs(nprnx,j),jlli(nprnx,j),jssi(nprnx,j),j=1,nobtyp1)
              endif
            endif
          enddo      
          nprn = nprnx         
                         
c       endif on RINEX version
        endif
                
c       Select the observables to return and convert from pseudorange
c       convention to doppler convention. Also note that obs = 0 is no 
c       data, and should be flagged as bad.

        if (iflag.gt.0 .and. ioerr .eq. 0 ) then          
          do i=1,nprn
            dofl1(i) = -obs(i,iobtypx(1))  
            issi(i,1) = jssi(i,iobtypx(1))
            illi(i,1) = jlli(i,iobtypx(1)) 
            if(dabs(dofl1(i)) .lt. 0.001d0 ) issi(i,1) = 1
            dofl2(i) = -obs(i,iobtypx(2)) 
            issi(i,2) = jssi(i,iobtypx(2))
            illi(i,2) = jlli(i,iobtypx(2)) 
            if(dabs(dofl2(i)) .lt. 0.001d0 ) issi(i,2) = 1
            prgl1(i) =  obs(i,iobtypx(3))    
            issi(i,3) = jssi(i,iobtypx(3))
            illi(i,3) = jlli(i,iobtypx(3)) 
c           if P1 is missing, see if C1 is available
            if( prgl1(i).eq.0.d0 ) then
              if(iobtypx(5).ne.0) then
                prgl1(i) =  obs(i,iobtypx(5))    
                issi(i,3) = jssi(i,iobtypx(5))
                illi(i,3) = jlli(i,iobtypx(5)) 
              endif   
            endif
            prgl2(i) = obs(i,iobtypx(4)) 
            issi(i,4) = jssi(i,iobtypx(4))
            illi(i,4) = jlli(i,iobtypx(4)) 
c           if P2 is missing, see if C2 is available
            if( prgl2(i).eq.0.d0 ) then
              if(iobtypx(6).ne.0) then
                prgl2(i) =  obs(i,iobtypx(6))    
                issi(i,4) = jssi(i,iobtypx(6))
                illi(i,4) = jlli(i,iobtypx(6)) 
              endif   
            endif
          enddo
        endif  

c       Restore the leap second to Glonass pseudoranges-- no, do not 
c        if( rxtime.eq.'GPS'.and.gnss.eq.'R' ) then
c          jd = julday(month,iday,iyr)
c          utcoff = taiutc(jd) - 19.d0       
cd          print *,'RRINEX rxtime gnss utcoff ltvel '
cd     .           ,        rxtime,gnss,utcoff,ltvel
c          do i=1,nprn
c            prgl1(i) = prgl1(i) + utcoff*ltvel
c            prgl2(i) = prgl2(i) + utcoff*ltvel 
c          enddo
c        endif 


        if(debug) then  
          write(*,*) 'Observations selected L1 L2 P1 P2'
          do i=1,nprn
            write(*,*) isvid(i),dofl1(i),dofl2(i),prgl1(i),prgl2(i)
          enddo
        endif                                              
        if( debug) print *,'Before epoch-error checking ferr: ',ferr

c       Do some error checking here
c       Time since last good epoch, in seconds
        jdoy = idoy(iyr,month,iday)
        itflag= -4   
        call timcon(itflag,iwkn,sow,iyr,jdoy,ihr,min,sec,utcoff)
        ds = secdif(iwkn,sow,iwkn1,sow1)
        if (ds .gt. 0.0d0 ) then
           nepoch = nepoch + 1
           sow1 = sow
           iwkn1 = iwkn    
           if (debug) then 
             if( nprn.gt.0 ) then
               do i = 1,nprn
                  write (*,'(a,i5,1x,4(1pe15.4,1x))')
     .                 'RINEX: data: '
     .               ,isvid(i),dofl1(i),dofl2(i),prgl1(i),prgl2(i)
               enddo
             endif
           endif
        else if (ds .eq. 0.0d0) then   
          write(message,'(a,i5,3i3,f11.7,i3,i3,24i3)')
     .           'Time stands still: '
     .          , iyr,iday,ihr,min,sec,idflag,nprn,(isvid(i),i=1,nprn)
          call report_stat('WARNING',prog_name,'rrinex',' '
     .                           ,message,0)
          nerr = nerr+1
          ferr = .true.
        else if (ds .lt. 0.0d0) then  
          write(message,'(a,i5,4i3,f11.7,i3,i3,24i3)')
     .        'Time flows backwards: '
     .     , iyr,month,iday,ihr,min,sec,idflag,nprn,(isvid(i),i=1,nprn)  
           call report_stat('WARNING',prog_name,'lib/rrinex',' '
     .                        ,message,0)
          nerr = nerr+1
          ferr = .true.
c        else 
c          nerr = nerr + 1
c          ferr = .true.
        endif
      else 
         write(message,'(a,i5,4i3)') 'Bad time: ',iyr,month,iday,ihr
     .                             , min
         call report_stat('WARNING',prog_name,'lib/rrinex',' '
     .                   , message,0)
         write(message,'(a,a80)') 'Line read: ',line
         call report_stat('WARNING',prog_name,'lib/rrinex',' '
     .                   , message,0)
         nerr = nerr+1
         ferr = .true.
c     endif for valid time tag 
      endif
      if( debug) print *,'After error checking ferr: ',ferr

c     normal end         
      return

c     come here on end of file
4000  continue
      fend = .true.
      return
      end
       
c----------------------------------------------------------------------------------

      Subroutine rxerr (lunit,badline,ioerr,nerr)

c     Handle a RINEX file error
c     input:
c        lunit      logical unit for RINEX file
c        badline    offending line in file
c        ioerr      returned error message
c        nerr       total number of errors

      implicit none

      integer*4 lunit,ioerr,nerr,inqerr
      character*(*) badline
      integer*4 maxerr,rcpar,len

      character*80  prog_name,fname
      character*320  message

      data maxerr /100/
c     get calling program name and file name for report_stat
      len = rcpar(0,prog_name)
      inquire( unit=lunit, name=fname, iostat=inqerr )   
      if( inqerr.ne.0 )   
     .   call report_stat('FATAL',prog_name,'lib/rxerr',' '
     .                   ,'Cannot find RINEX filename',0)

      if (ioerr .gt. 0) then
                       
        write(message,'(a,a)') 'Bad line in RINEX file: ',badline(1:120)
        call report_stat('WARNING',prog_name,'lib/rxerr',fname
     .          ,message,ioerr)

         if (nerr .ge. maxerr) then
            write(message,'(a,i3)') 'Error count exceeds ',maxerr
            call report_stat('FATAL',prog_name,'lib/rxerr',fname
     .          ,message,ioerr)
         else    
            nerr = nerr + 1
            return
         endif
      else   
         return
      endif

      end
