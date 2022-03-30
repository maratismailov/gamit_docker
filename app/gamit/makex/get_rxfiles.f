Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego. 1997.   All rights reserved.

      Subroutine get_rxfiles( debug,site,iwknstart,sowstart,iwknstop
     .                      , sowstop,rx_doy_minus,rx_doy_plus,rxfiles )
           
c     Search a set of directories and find all the RINEX files for a given 
c     station that have data within a specified time span.  Coded to handle
c     up to 5 (local maxdir) directories and MAXFIL RINEX files (includes/makex.h maxfils)

c     R. King  7 Feb 97; last modified by rwk 20 May 98

c     Input
c                             
c       site        c*4     station name
c       iwknstart   i*4     GPS week number start of span
c       sowstart    r*8     GPS seconds-of-seek start of span
c       iwknstop    i*4     GPS week number end of span
c       sowstop     i*4     GPS seconsd-of-week end of span

c     Output
c
c       rxfiles(maxfil)  c*80     full path name for all RINEX files to be searched
c                                   and maxfil is set in ../includes/makex.h 


      implicit none
      
      include '../includes/makex.h'

      logical found_data(maxfil),fcheck,debug

      integer*4 iwknstart,iwknstop,iystart,iystop,idoystart,idoystop 
     .        , yrdays(2,maxfil),yr2,doy,nbuf,nrxs,nrx,ioerr
     .        , nblen,ic,count_arg,ib,lift_arg,lenmsg
     .        , isort(maxfil),rx_doy_minus,rx_doy_plus,nday
     .        , itflag,ihr,min,i,j,k
      integer*4 maxdir
      parameter (maxdir=5)              
      integer*4 trimlen

      character*4 site,sitelc
      character*80 wild,wildsave
      character*80 rxdirs(maxdir),rxsearch(maxfil),rxfiles(maxfil)
     .           , rxfbuf(maxfil)
      character*256 message

c     for debug only
                                               
      real*8 sowstart,sowstop,sec,utcoff
     .     , jdrx_search(maxfil),jdrx(maxfil)

      save rxdirs

      logical first_call

      data first_call/.true./
           
                      

      if( first_call ) then
         first_call = .false.   
c        check the array limits for the RINEX file list
         nrx = rx_doy_minus + rx_doy_plus + 1   
         if( nrx.gt.maxfil ) then     
           write(message,'(a,i2,a,i2,a)') 'Searching ',rx_doy_minus
     .       ,' before and ',rx_doy_plus
     .       ,' days after obs day for RINEX files'
           call report_stat('WARNING','MAKEX'
     .         ,'get_rxfiles',' ',message,0 )
           write(message,'(a,i3,a)') 'Number of RINEX files ',nrx
     .                 ,' exceeds MAXFIL'
           call report_stat('FATAL','MAKEX','get_rxfiles',' ',message,0)
         endif
c        decode the RINEX file line to get the directories to be searched    
         do i=1,maxdir
           rxdirs(i) = ' '
         enddo
c        the 'frinex' variable read by rbatch may contain up to maxdir directories, 
c        separated by blanks and not exceeding 79 characters.  frinex is
c        later overwritten as the full path of the current RINEX file
         ic = count_arg(frinex) 
         do i=1,ic
            if( ic.gt.maxdir) call report_stat('FATAL','MAKEX'
     .               ,'get_rxfiles',' ','Too many RINEX directories',0 )
            ib = lift_arg(frinex,rxdirs(i),i)
         enddo        
      endif             

c     determine the days to be searched: from start - rx_doy_minus to 
c     start + rx_doy_plus
                                            

      do i=1,2       
        do j=1,maxfil
          yrdays(i,j) = 0
        enddo 
      enddo  
      itflag = 4
      call timcon( itflag,iwknstart,sowstart,iystart,idoystart
     .           , ihr,min,sec,utcoff ) 
      call timcon( itflag,iwknstop,sowstop,iystop,idoystop
     .           , ihr,min,sec,utcoff ) 
c     fill the arrays of days to be searched  
      do j= 1, rx_doy_minus+rx_doy_plus+1 
        nday = j - (rx_doy_minus+1)
        yrdays(1,j) = iystart 
        yrdays(2,j) = idoystart + nday  
        call check_year( yrdays(1,j),yrdays(2,j))
c        print *,'j nday yrdays ',j,nday,(yrdays(i,j),i=1,2)
      enddo

c     first search the input directories to accumulate a list of all
c     RINEX files to be searched for data

c     initialize the number of RINEX files to be searched for data and
c     the file names
      nrxs = 0                     
      do i=1,maxfil
        rxsearch(i) = ' '
        rxfiles(i) =  ' '
        rxfbuf(i) = ' '
      enddo
     
c     set the site name to lowercase for directory search
      sitelc = site
      call lowers(sitelc)                                           

c     loop over directories
        
      i = 1     
      do while ( rxdirs(i)(1:1).ne.' ')
                                     
c        limit the search to the days specified 
         j=1       
         do while ( j.le.maxfil .and. yrdays(1,j).ne. 0 ) 
           yr2 =  mod(yrdays(1,j),100)
           doy = yrdays(2,j)   
           write(wild,'(a,a4,i3.3,a2,i2.2,a1)') 
     .       rxdirs(i)(1:nblen(rxdirs(i))),sitelc,doy,'?.',yr2,'o' 
c           write(wild,'(a4,i3.3,a2,i2.2,a1)') sitelc,doy,'?.',yr2,'o' 
c            print *,'wild ',wild   
            write(wildsave,'(a80)') wild
           call getdir(wild,maxfil,rxfbuf,nbuf) 
c           print *,'nbuf rxfbuf ',nbuf,(rxfbuf(k),k=1,3)
           do k=1,nbuf    
             rxfbuf(k) = rxfbuf(k)(1:trimlen(rxfbuf(k)))
             nrxs = nrxs + 1       
             if( nrxs.gt.maxfil ) goto 990   
             rxsearch(nrxs) = rxdirs(i)(1:nblen(rxdirs(i)))//rxfbuf(k)
c             print *,'rxsearch(nrxs) ',nrxs,rxsearch(nrxs)
           enddo  
           j = j+1 
         enddo
         i = i+1  
         
      enddo          
      if( nrxs.eq.0 ) then  
        write(message,'(a,a4,a,a80)') 'No RINEX files found for ',sitelc
     .      ,' from ',wildsave(1:nblen(wildsave))
        call report_stat('WARNING','MAKEX','get_rxfiles',' ',message,0) 
        write(uinfor,'(a)') message  
        return
      else         
c       don't allow more than 5 full-path names per message to avoid > 256 chars
        lenmsg = nrxs*(nblen(rxsearch(i))+1)  + 22
        write(uinfor,'(2a)') 'Searching for data in '
     .     ,(rxsearch(i)(1:nblen(rxsearch(i))+1),i=1,nrxs)  
        if( lenmsg.le.256 ) then 
           write(message,'(20a)') 'Searching for data in '  
     .        ,(rxsearch(i)(1:nblen(rxsearch(i))+1),i=1,nrxs)
        else
            write(message,'(a,i3,a)') 
     .          'Search for data in',nrxs,' RINEX files'
        endif
        call report_stat('STATUS','MAKEX','get_rxfiles',' ',message,0)  
      endif
        
c     now we have a list (full pathnames) of files to be searched--
c     keep only those that have data within the range specified

      do i = 1,nrxs
        frinex = rxsearch(i)
c       something screwy if this fails since we've just 'ls'ed the files
        if( .not.fcheck(frinex) ) goto 991  
        open(unit=urinex,file=frinex,status='OLD',iostat=ioerr)  
        if (ioerr.ne.0 ) goto 992
        call scan_rinex( debug,iwknstart,sowstart,iwknstop,sowstop
     .                 , found_data(i),jdrx_search(i) )   
c        print *,'frinex ',frinex
c        print *,' i jdrx_search ',i,jdrx_search(i)
        close(unit=urinex) 
      enddo


c     construct the output list to include only those files with requested data
c     and sort by time
             
      nrx = 0
      do i = 1,nrxs
        if( found_data(i) ) then
           nrx = nrx + 1
           rxfiles(nrx) = rxsearch(i)  
           jdrx(nrx) = jdrx_search(i)
        endif
      enddo    
      if( nrx.gt.0 ) then                  
        call indexx( nrx,jdrx,isort ) 
        do i=1,nrx
          rxfbuf(i) = rxfiles(i)
        enddo
        do i=1,nrx
          rxfiles(i) = rxfbuf(isort(i))
        enddo    
c       don't allow more than 5 full-path names per message to avoid > 256 chars
        lenmsg = nrxs*(nblen(rxsearch(i))+1)  + 22
        write(uinfor,'(2a)') 'Found data in '
     .     ,(rxfiles(i)(1:nblen(rxfiles(i))+1),i=1,nrx)  
        if( lenmsg.le.256 ) then 
           write(message,'(20a)') 'Found  for data in '  
     .        ,(rxfiles(i)(1:nblen(rxfiles(i))+1),i=1,nrx)
        else
            write(message,'(a,i3,a)') 
     .          'Found data in ',nrx,' RINEX files'
        endif
      else
        write(message,'(a,a4)') 'No data found on RINEX files for '
     .    ,sitelc
        call report_stat('WARNING','MAKEX','get_rxfiles',' ',message,0) 
        write(uinfor,'(a)') message  
        return
      endif
            
      return            

c     problem encountered

  990 call report_stat('FATAL','MAKEX','get_rxfiles',' ',
     .     'Number of RINEX files to be searched exceeds 10',0)                  
                          
  991 call report_stat('FATAL','MAKEX','get_rxfiles',rxfiles(i)
     .     ,'Cannot find selected RINEX file: ',0)

  992 call report_stat('FATAL','MAKEX','get_rxfiles',rxfiles(i)
     .     ,'Cannot open RINEX files: ',ioerr)
      return
      end                   
         

      Subroutine scan_rinex( debug,iwknstart,sowstart,iwknstop,sowstop
     .                     , found_data, jdrx )

c     R. W. King  10 February 1997

c     Scan a RINEX file to determine if there are data between two
c     input epochs.  Use the library header and data read routines,
c     even though they pass extraneous variables, because they have 
c     they have lots of error checking and elaborate logic to keep
c     synched on epoch records.

c     Input

c        iwkn(2)   i*4  start, stop week number of data to be found
c        sow(2)    r*8  start, stop seconds of week of data to be found
         
c     Output

c         found_data logical*4   true if data on this RINEX file
                                                               

c     RINEX file name 'fname' and logical unit number 'urinex' stored in 
c     common /fileinfo/ in ../includes/makex.h


      implicit none
              
      include '../includes/makex.h'

      integer*4 iwknstart,iwknstop,iwkno1,irxdoy,itflag
     .         ,julday,i

      real*8  sowstart,sowstop,sow1,utcoff,secdif,jdrx
c**      these used only for commented-out header code
c**      real*8 sow0 
c**      integer*4 idoy,iwkno0

      logical found_data
                                  
c     iflag = 0 to scan for epochs, do not return observables
      integer*4 iflag

c     RINEX defined items
      real*4 rxver
      character*20 rxpgm,rxusr,rxdat
c     comment
      integer irxcom
      character*60 rxcom(maxlin)
c     mark name
      character*60 rxmrk
c     observer
      character*20 rxobs
c     agency
      character*40 rxagy
c     receiever serial number, type and SW version
      character*20 rcvnum
      character*20 rctype
      character*20 rcvers
c     antenna serial number and type
      character*20 antnum
      character*20 anttyp
c     aproximate coordinates
      real*8 apx,apy,apz
c     antenna offsets
      real*8 anth,ante,antn
c     wavelength factors
      integer nwave1,nwave2
c     data interval in seconds
      real*8 rxint

c     start, stop times returned from rrxhed
      integer irxyr0,irxmo0,irxdy0,irxhr0,irxmn0
      real*8 rxsec0
c     data stop time
      integer irxyr1,irxmo1,irxdy1,irxhr1,irxmn1
      real*8 rxsec1                         
c     time type ('GPS' with mixed files) 
      character*3 rxtime

c     variables passed to/from rrinex
      logical debug
C     Value set to .true. at end of file
      logical fend
C     Value set to .true. on error
      logical ferr
c     GNSS system (not used in scan, so set artificially to 'G'
      character*1 gnss 
c     number of different types of observable quantities
      integer nobtyp     
c     observation types for selected gnss and indices in the array 
      character*3 rxobtyp(maxobt)                                  
      integer*4 iobtypx(6)
c     number of satellites
      integer*4 nprn
C     SV PRN numbers at each epoch
      integer*4 isvid(maxchn)
c     These observable quantities not used here
C     L1, L2 doppler phase in cycles
      real*8  dofl1(maxchn), dofl2(maxchn)
C     L1, L2  pseudorange in meters
      real*8  prgl1(maxchn), prgl2(maxchn)
c     number of this epoch
      integer     nepoch
c     loss of lock and signal strength indicators 
      integer*4 illi(maxchn,4),issi(maxchn,4)
               

c     initialize some variables
      found_data = .false.
      jdrx = 0.d0    
      iwkno1 = 0 
      sow1 = 0. 

c     set RINEX version to 1 until read, for initial check of end-of-header
      rxver = 1.0
                
c     read the header 
c        gnss not used in the scan, so set artificailly to GPS
      gnss = 'G'  
      if(debug) print *,'GET_RXFILES calling RRXHED '
      call rrxhed (debug,gnss,
     .   rxver,rxpgm,rxusr,rxdat,rxcom,irxcom,rxmrk,rxobs,rxagy,
     .   rcvnum,rctype,rcvers,antnum,anttyp,apx,apy,apz,
     .   anth,ante,antn,nwave1,nwave2,nobtyp,rxobtyp,rxint,rxtime,
     .   irxyr0,irxmo0,irxdy0,irxhr0,irxmn0,rxsec0,
     .   irxyr1,irxmo1,irxdy1,irxhr1,irxmn1,rxsec1)  
          
cd     print *,'GET_RXFILES after rrxhed rxver nobtyp ',rxver,nobtyp
c     scan the RINEX file to find any epoch within the requested window
      fend = .false. 
      nepoch = 0           
c     observales not returned in scanning 
      do i=1,maxobt
       rxobtyp(i) = ' '
      enddo
      iflag = 0 
      do i=1,6
        iobtypx(i) = 0 
      enddo
      if(debug) print *,'GET_RXFILES calling RRINEX '
      do while ( .not.fend )    
cd        print *,'before rrinex iflag gnss nobtyp ',iflag,gnss,nobtyp
        call rrinex( debug,iflag,rxver,gnss,nobtyp,rxobtyp,iobtypx
     .             , nprn,isvid,rxtime,iwkno1,sow1,nepoch
     .             , dofl1,dofl2,prgl1,prgl2,illi,issi
     .             , anth,ante,antn,fend,ferr )
cd           print *,'iwkno1 sow1 ',iwkno1, sow1
        if( secdif(iwkno1,sow1,iwknstart,sowstart).gt.0.d0 .and.
     .      secdif(iwkno1,sow1,iwknstop,sowstop).lt.0.d0 ) then
c              print *,'found data in window ',frinex
              found_data = .true.  
c             save the start time as a Julian date for get_rxfiles sort 
              itflag = 4    
              call timcon( itflag,iwkno1,sow1,irxyr0
     .                   , irxdoy,irxhr0,irxmn0,rxsec0,utcoff)
 
              call monday( irxdoy,irxmo0,irxdy0,irxyr0 )
              jdrx = dfloat( julday(irxmo0,irxdy0,irxyr0)) 
     .             + ( irxhr0 + ( irxmn0 + rxsec0/60.d0 )/60.d0 )/24.d0    
c               print *,'jdrx ',jdrx
              go to 999  
        endif 
      enddo  
        

  999 return
      end
           
      Subroutine add_day(year,day)

      implicit none 
     
      integer*4 year,day,year1,day1,ndyr
      logical leapyr

      ndyr = 365
      if( leapyr(year)) ndyr = 366
      day1 = day + 1
      year1 = year
      if( day1.gt.ndyr ) then
          day1 = day1 - ndyr
          year1 = year + 1
      endif
      day = day1
      year = year1
      return
      end
    
      Subroutine sub_day(year,day)
 
      implicit none 
     
      integer*4 year,day,year1,day1,ndyr
      logical leapyr

      ndyr = 365
      if( leapyr(year)) ndyr = 366
      day1 = day - 1
      year1 = year
      if( day1.lt.1 ) then
          day1 = day1 + ndyr
          year1 = year - 1
      endif
      day = day1
      year = year1
      return
      end


      Subroutine check_year(year,day)

      implicit none  

      integer*4 year,day,ndyr
      logical leapyr
      
      ndyr = 365
      if( leapyr(year)) ndyr = 366   
      if( day.lt.1 ) then   
        if( leapyr(year-1)) ndyr = 366   
        day = ndyr + day 
        year = year -1
      elseif( day.gt.ndyr ) then 
        if( leapyr(year)) ndyr = 366   
        day = day - ndyr
        year = year + 1
      endif
      return
      end

