      Program SCAN_RINEX

*     Scan a RINEX observation file and list the number and type of 
*     observataions for each SV.  R. King April 2020
      
*     Command line:   scan_rinex [file]

      Implicit none
                                     
      include '../includes/dimpar.h'
      include 'rinex.h'
                                                   
* In rinex.h:
*
*   Header
*  rx_version               : RINEX version number (r*4)
*  rx_rcvtyp                : 20-character receiver name
*  rx_rcvsw                 : 20-character firmware version
*  rx_anttyp                : 20-character antenna name including radome 
*  rx_ngss,rx_gnss(maxgnss) : GNSS systems found of those requested (c*1)
*  rx2_nobs                 : Number of observable types for RINEX 2
*  rx_obscod(maxgnss)       : Number of observable types for each system
*  rx2_pcncod               : Identifies cross-correlated and misreported L1 range observable
*  rx2_obscod(maxob11)       : RINEX 2 observable types mapped to 3-character RINEX 3 
*  rx_obscod(maxob11,maxgnss): 3-character observable types 
*   Data records (each epoch)
*  rx_nprn                - number of SV data returned
*  rx_prn(maxsch)         - PRNs selected (c*3)
*  rx_obs(maxsch,maxob11)  - observations return
                            
* Requested GNSS
      integer*4 ngnss
      character*1 gnss_list(maxgnss)

* Input/output names and units
      character*80 rnxfname,sumfile
      integer*4 lurnx/1/,lusum/2/
                   
* Arrays for used PRNs and observables
      
      integer*4 nprn,nobscod,nobs(maxob11,maxsch),iptr(maxgnss,maxob11)
      character*3 prn(maxsch),obscod(maxob11)
              
* Funtions
      integer*4 iclarg,nblen,isvarray,gindex 

* Other local       
      integer*4 iarg,ioerr,ignss,iepoch,ip,ig,io,i,j       
      character*128 message         
      logical fend,ferr,debug/.false./
      integer*4 ndatsat   ! Number of satellite with data
                          ! for a particular Observable
      integer*4 ndattot   ! Total number of measurements for
                          ! each observable

* Variables needed for sorting
      integer*4   sl(maxsch)     ! index to sorted PRN list
      integer*4   tmp_sl         ! Temporary storage for swap
      character*3 tmp_prn        ! Temporary stoage for PRN.
     .,           prncp(maxsch)  ! List to PRNs to be sorted    


* Decode the command line and open the files   
        
      iarg = iclarg(1,rnxfname)
      if(iarg.eq.0) then 
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) 'Program scan_rinex'
        write(*,*) ' '
        write(*,*) ' Count of observations on a RINEX file'
        write(*,*) ' '
        write(*,*) ' scan_rinex <rinex file> <out-file>/6 <gnss list>'
        write(*,*) ' '
        write(*,*) 'where <rinex file>  is the name of the RINEX file' 
        write(*,*) ' <out-file> is the output summary, and  <gnss list>'
        write(*,*) ' is a case-insensitive list of the systems to be '
        write(*,*) ' included (G R E  etc ) '
        write(*,*) ' '
        stop
      else 
        open(lurnx,file =rnxfname,status='old',iostat=ioerr)
        if(ioerr.ne.0) call report_stat('FATAL','SCAN_RINEX'
     .   ,'utils/scan_rinex',rnxfname,'Error opening input file',ioerr)
      endif  
      iarg = iclarg(2,sumfile)
* MOD TAH 200511: See if 6 or * passed so that piping can be used.
      if( sumfile(1:1).eq.'6' .or. sumfile(1:1).eq.'*' ) then
*        Set lusum to 6.
         lusum = 6
      elseif( iarg.eq.0.or.sumfile(2:2).eq.' ') then
        call report_stat('FATAL','SCAN_RINEX','utils/scan_rinex',' ' 
     .                  ,'Missing name of output file',0)
      else
        open(lusum,file = sumfile,status='unknown',iostat=ioerr)
        if(ioerr.ne.0) call report_stat('FATAL','SCAN_RINEX'
     .  ,'utils/scan_rinex',rnxfname,'Error opening output file',ioerr)
      endif
      ngnss = 0    
      iarg = iclarg(3,gnss_list(1))
      if(iarg.eq.0 ) then
        call report_stat('WARNING','SCAN_RINEX','utils/scan_rinex',' ' 
     .                  ,'Missing GNSS list, assume GPS-only',0)
        ngnss = 1
        gnss_list(1) = 'G'    
      else
        do i=1,maxgnss
          iarg = iclarg(i+2,gnss_list(i))
          call uppers(gnss_list(i))
          if(iarg.gt.0) then
            ngnss = ngnss + 1
          endif
        enddo
      endif               
      if(debug) print *,'ngnss gnss_list',ngnss,(gnss_list(i),i=1,ngnss)


* Read the RINEX header and print the information

      call read_rinex_header( lurnx,ngnss,gnss_list )
      if(debug) then               
         if( rx_version.lt.3. ) then 
           write(*,'(a,i3,50(1x,a2))') 'SCAN_RINEX rx2_nobs rx2_obscod '
     .           ,rx2_nobs,(rx2_obscod(i),i=1,rx2_nobs)
         else 
           write(*,'(a,i2,6i3)') 'SCAN_RINEX rx_ngnss rx_nobs '
     .     , rx_ngnss,(rx_nobs(i),i=1,rx_ngnss)         
           do ignss=1,rx_ngnss
             write(*,'(a,2i2,30(1x,a3))') 'ignss rx_nobs rx_obscod '
     .       ,ignss,rx_nobs(ignss),(rx_obscod(i,ignss),i=1,rx_nobs(i))
           enddo
        endif 
      endif
      write(lusum,'(3a,f4.2)') 'SCAN_RINEX summary of RINEX '
     .                              ,rnxfname,'  Version ',rx_version
      write(lusum,'(a,6(a1,1x))') 'Requested GNSS : '
     .                             ,(gnss_list(i),i=1,ngnss)
      write(lusum,'(/,a,a20)')  'Receiver : ',rx_rcvtyp
      write(lusum,'(a,a20)')    'Firmware : ',rx_rcvsw
      write(lusum,'(a,a,a20)')  'Antenna  : ',rx_anttyp
      if(rx_version.lt.3.) then 
        write(lusum,'(/,a,50(1x,a3))') 
     .     'RINEX 2 observables: ',(rx2_obscod(i),i=1,rx2_nobs)
      elseif( rx_version.ge.3.) then
        write(lusum,'(a)') 'Observables : '
        do j=1,rx_ngnss 
          if(debug) print *,'DEBUG rx_ngnss j rx_gnss rx_nobs '
     .          ,   rx_ngnss,j,rx_gnss,rx_nobs(j) 
          write(lusum,'(2x,a1,50(1x,a3))') 
     .        rx_gnss(j),(rx_obscod(i,j),i=1,rx_nobs(j))
        enddo  
      else
        call report_stat('FATAL','SCAN_RINEX','utils/scan_rinex',' '
     .                  , 'Unknown RINEX version',0)
      endif
        
* Fill the global observation code array from the RINEX 3 list
* for each GNSS, then set pointers from the global array to the 
* RINEX array for each GNSS 

      nobscod = 0     
      do j=1,ngnss 
        do i=1,rx_nobs(j) 
*         isvarray is written for PRNs but works as well for any c*3
          if( isvarray( rx_obscod(i,j),obscod,nobscod).eq.0 ) then
             nobscod = nobscod + 1
            obscod(nobscod) = rx_obscod(i,j)
          endif 
        enddo    
      enddo
      if(debug) then
        print *,'rx_version ',rx_version
        print *,'nobscod obscod ',nobscod,(obscod(i),i=1,nobscod)
      endif                              
      do i=1,nobscod
        do ig=1,ngnss
          do io=1,rx_nobs(ig)
            if(rx_obscod(io,ig).eq.obscod(i)) iptr(ig,io) = i
          enddo
        enddo
      enddo                                   
      if(debug) then      
        do ig = 1,ngnss
          do io=1,rx_nobs(ig) 
            print *,'ig rx_nobs io iptr ',ig,rx_nobs(ig),io,iptr(ig,io)
          enddo
        enddo 
      endif 

* Read through all the epochs and accumulate the statistics

      fend = .false.  
      iepoch  = 0 
      do while ( .not.fend )                 
                         
*       reinitialize the observables at each epoch
        do i=1,maxsch
          do j=1,maxob11
            rx_obs(i,j) = 0.d0
          enddo
        enddo 
        iepoch = iepoch  + 1
        call read_rinex_data(lurnx,ngnss,gnss_list,ferr,fend) 

       if (ferr) then 
        write (message,'(a,i4,f7.0)') 'Error in RINEX file at epoch '
     .                            ,rx_jd,rx_t
        call report_stat('WARNING','SCAN_RINEX','utils/scan_rinex',' '
     .          ,message,0)
        endif  
        if(debug.and.iepoch.gt.3) then 
          print *,'DEBUG set fend=.true.'
          fend = .true.
        endif 
cc        if(fend) then     
cc           write(message,'(a,i4,f8.1)') 'EOF on RINEX at '
cc     .        ,rx_jd,rx_t
cc           call report_stat('WARNING','SCAN_RINEX','utils/scan_rinex'
cc     .                       ,rnxfname,message,0)
cc        endif         
*       see if the global SV array needs to be incremented
        do i=1,rx_nprn    
          if( isvarray( rx_prn(i),prn,nprn).eq.0 ) then 
            nprn = nprn + 1
            prn(nprn) = rx_prn(i)
          endif
        enddo                                             
*       now increment the numbers for the current list of SVs
        if(debug) then 
          write(*,'(a,i3,30(1x,a3))') 'scan list: '
     .            ,nprn,(prn(i),i=1,nprn)
          write(*,'(a,i3,30(1x,a3))') 'rx   list: '
     .            ,rx_nprn,(rx_prn(i),i=1,rx_nprn) 
        endif 
        do i=1,nprn
          do ip=1,rx_nprn  
            if(rx_prn(ip).eq.prn(i)) then 
              if(debug) print *,'Found i prn ip prn '
     .           ,i,prn(i),ip,rx_prn(ip)   
              ig=gindex(rx_prn(ip)(1:1),gnss_list,ngnss)
              if(debug) print *,'ip prn1 ig ngnss gnss_list  '
     .           , ip,rx_prn(ip)(1:1),ig,ngnss,(gnss_list(j),j=1,ngnss)
              if(debug) print *,'Check PRN i ip prn ig gnss '
     .             ,i,ip,prn(i),ig,gnss_list(ig)
              do j=1,nobscod
                do io=1,rx_nobs(ig)
                  if(j.eq.iptr(ig,io).and. rx_obs(ip,io).ne.0.d0) then
                    nobs(j,i) = nobs(j,i) + 1   
                   if(debug) print *
     .               ,'Match match j io ig iptr rx_obs obscon nobs '
     .               ,j,io,ig,iptr(ig,io),rx_obs(ip,io)
     .               , obscod(iptr(ig,io)),nobs(j,i)
                  endif
                enddo 
              enddo
            endif
          enddo
        enddo 
        if(debug) then 
          print *,'SCAN_RINEX At epoch ',iepoch,' nprn ',nprn
          write(*,'(100(1x,a3))') (prn(i),i=1,nprn)
          write(*,'(100i4)') (nobs(1,i),i=1,nprn) 
        endif 
           
*     end loop through RINEX data
      enddo      
                                 
* Write out the observables found for RINEX 2

      if(rx_version.lt.3.) then 
        write(lusum,'(/,a,50(1x,a3))') 
     .     'Translation to RINEX 3 observables: '
        if(debug) then 
          write(*,'(a,i3,6(1x,a1,i3))') 
     .       ' SCAN_RINEX ngnss rx_gnss rx_nobs '
     .       ,ngnss,(rx_gnss(i),rx_nobs(i),i=1,ngnss)
        endif   
        do j=1,ngnss                       
          write(lusum,'(2x,a1,50(1x,a3))') 
     .        rx_gnss(j),(rx_obscod(i,j),i=1,rx_nobs(j))
        enddo  
      endif

* Write out the summary
                                                    
      write(lusum,'(/,a,/)') 'Summary: By Satellite'
      write(lusum,'(4x,200(3x,a3))') (obscod(i),i=1,nobscod)
* MOD TAH 20509: Replace space in name with 0 so strings can be
*           awk'd by token number.
      do i = 1,nprn
        if( prn(i)(2:2).eq.' ' ) prn(i)(2:2) = '0' 
*       Make copy for sorting
        prncp(i) = prn(i) 
        sl(i)    = i       ! Index that will get sorted and
                           ! used for output. 
      end do
* MOD TAH 200509: Sort the PRN strings.  Create index for 
*     sorted list.
      do i = 1,nprn-1
         do j = 1, nprn-1
            if( prncp(j).gt.prncp(j+1) ) then
*              Swap the prn and list
               tmp_prn    = prncp(j+1)
               tmp_sl     = sl(j+1)
               sl(j+1)    = sl(j)
               prncp(j+1) = prncp(j)
               sl(j)      = tmp_sl
               prncp(j)     = tmp_prn
            endif
         end do
      end do

*     Output the sorted list
      do i=1,nprn
         write(lusum,'(1x,a3,200(i6))') prn(sl(i)),
     .                (nobs(j,sl(i)),j=1,nobscod)   
      enddo

* MOD TAH 200520: Now summarize by frequency (allows test to see what is 
*       present.
      write(lusum,'(/,a,/)') 'Summary: By Observable'
      write(lusum,'(a)') ' OBS NPRN UPRN       %  Total Obs'
      do j = 1, nobscod
         ndatsat = 0
         ndattot = 0
         do i = 1,nprn
            if( nobs(j,i).gt.0 ) then
               ndatsat = ndatsat + 1
               ndattot = ndattot + nobs(j,i)
             endif
         end do
*        Add small number to nprn so we don't divide by zero
         write(lusum,'(1x,a3,1x,I4,1x,I4,1x,F7.1,1x,i8)') obscod(j),
     .         nprn, ndatsat, ndatsat/(nprn+0.0001)*100, ndattot
      end do
      write(lusum,'(1x)')

      write(*,'(/,a)') 'Scan complete'

      stop
      end

   

*************************************************************************

      Subroutine READ_RINEX_HEADER ( lun,ngnss_req,gnss_req ) 

*     Read a RINEX header and save the information needed to process 
*     the observations.  Robust for RINEX 3, some assumptions made for
*     RINEX 2.   Based on the original version for RINEX 1 by Kurt Feigl
*     (1991), modified by R. King for RINEX 2 GPS (1991-2007), then for
*     single-GNSS and RINEX 3 (2016), then for multi-GNSS 2018.
*     
*     The header inforamtion regarding receiver and antenna is not 
*     saved since it might be unreliable and GAMIT defers to station.info

*     File is assumed to be assigned a unit number and opened in the 
*     calling module

      implicit none

      include '../includes/dimpar.h' 
      include 'rinex.h'
 
* Input 

* lun - logical unit number for the file
      integer*4 lun                                 
    
* gnss_req          - number of GNSS requested
* gnss_req(maxgnss) - List of GNSS requested to be processed (c*1)
      integer*4 ngnss_req 
      character*(*) gnss_req(maxgnss)
  
* Output stored in includes/rinex.h
         
*  rx_version                : RINEX version number (r*4)
*  irxcom,rx_comment(maxlin) : Number and c*60 content of comment lines
*  rx_rcvtyp                 : 20-character receiver name
*  rx_rcvsn                  : 20-character receiver serial number    
*  rx_rcvsw                  : 20-character firmware version
*  rx_anttyp                 : 20-character antenna name including radome 
*  rx_antsn                  : 20-character antenna serial number 
*  rx_ngnss,rx_gnss(maxgnss) : GNSS systems found of those requested (c*1)
*  rx2_nobs                  : Number of observable types for RINEX 2
*  rx_obscod(maxgnss)        : Number of observable types for each system
*  rx2_pcncod                : Identifies cross-correlated and misreported L1 range observable
*  rx2_obscod(maxob11)        : RINEX 2 observable types mapped to 3-character RINEX 3 
*  rx_obscod(maxob11,maxgnss) : 3-character observable types 
*  rx_inter                  : Data interval in seconds (r*8)
*  rx_time                   : Time type ('GPS' or 'UTC') (c*3
*  rx_nphshft                : Number of phase shift entries
*  rx_phshtcod(maxob11)       : 4-character code for phase shifts (e.g. 'GL2W')
*  rx_phshft(maxob11)         : Phase shift to be applied (r*8) 
*    Note: RINEX 3 allows SV-dependent phase shifts, but we have not coded this
*
*  See read_rinex_data for the epoch-dependent variables.

* RINEX variables kept locally
                       
*  File type (filtyp should always be 'O'; rxsvs can be 'M' for mixed, 
*  or any one GNSS)
      character*1  filtyp,rxsvs
      
*  File origin identifiers
      character*20 rxprog,user,date
        
*  Mark name, mark number,observer, and agency
      character*60 mark
      character*20 marknum,observer
      character*40 rxagy

*  Approximate coordinates
      real*8 apx,apy,apz

*  Antenna offsets ( ARP - mark )
      real*8 anth,ante,antn

*  Wavelength factors  (RINEX 2)
      integer*4 nwave1,nwave2
*     --different wavelength factors for satellites 
      integer nsvwave,isvwave(maxsch),nw,nw1,nw2
     .       ,iwave1(maxsch),iwave2(maxsch) 

*  Flag for receiver clock offsets applied
      integer*4 ircvclk

*  SV system and observation types
      integer*4 nobsx 
      character*1 gnssrx 
      character*3 obscodx(maxob11)
  
*  Data start and stop times
      integer irxyr0,irxmo0,irxdy0,irxhr0,irxmn0 
     .      , irxyr1,irxmo1,irxdy1,irxhr1,irxmn1
      real*8 rxsec0,rxsec1

*  Number and names of satellites if in header
      integer irxnsv
      character*1 rxstyp(maxsch)
      integer irxsvn(maxsch)

*  Number of observations for each satellite
      integer irxnob(maxsch,maxob11)

* Other local variables

*     variables for status reporting
      integer*4 ioerr,nerr,len,rcpar
      character*80 fname,prog_name
      character*256 message

*     strings for buffering data
      character*20 buff20
      character*80 buff80,blnk80

*     L1-only flag
      logical l1_only,ok_obscod  

      integer irxsat,n1,n2,itype,i,j

* Functions
      
      integer*4 ignss,gindex     

* Debug
      logical debug/.false./                     

* Get calling program name and RINEX file name for report_stat

       len = rcpar(0,prog_name)
       inquire( unit=lun, name=fname, iostat=ioerr )
       if(ioerr.ne.0) call report_stat('WARNING',prog_name
     .    ,'lib/read_rinex_header'
     .   ,' ','Cannot get name of RINEX file for error reporting',ioerr)
     
* Initialize buffers to avoid using info from previous file if line missing
      write (blnk80,'(80x)')
      write (rxprog ,'(20x)')
      write (user ,'(20x)')
      write (mark ,'(20x)')
      write (observer ,'(20x)')
      write (rxagy ,'(40x)')
      write (rx_rcvsn,'(20x)')
      write (rx_rcvtyp,'(20x)')
      write (rx_rcvsw,'(20x)')
      write (rx_anttyp,'(20x)')
      write (rx_antsn,'(20x)')
      apx = 0.d0
      apy = 0.d0
      apz = 0.d0
      anth = 0.d0
      ante = 0.d0
      antn = 0.d0
      nwave1 = 0
      nwave2 = 0
      irxyr0 = 0
      irxmo0 = 0
      irxdy0 = 0
      irxhr0 = 0
      irxmn0 = 0
      rxsec0 = 0.d0
      irxyr1 = 0
      irxmo1 = 0
      irxdy1 = 0
      irxhr1 = 0
      irxmn1 = 0
      rxsec1 = 0.d0
      nsvwave = 0 
      do i=1,maxob11
       obscodx(i) = '   '
      enddo    
      rx_time = '   '
      do i=1,maxlin
        rx_comment(i) = ' ' 
      enddo

*  Set RINEX version to 1 until read, for initial check of end-of-header
      rx_version = 1.0
      irxcom = 0
      irxsat = 1 
      rx_inter = 0.
      nerr = 0              
      

*    --rwk 060707: Read first line of RINEX header separately to guard against a bogus blank first line
      read( unit=lun,iostat=ioerr,fmt='(a)') buff80
      if( buff80(61:73).eq.'rinex version' ) then
         call report_stat('WARNING',prog_name,'lib/read_rinex_header'
     .      ,fname,'First line of RINEX lower case',ioerr) 
      endif
      call uppers(buff80)       
      if( buff80(61:73).ne.'RINEX VERSION' )  then 
           call report_stat('FATAL',prog_name,'lib/read_rinex_header'
     .      ,fname,'Bogus first line of RINEX file',ioerr) 
      else
*       If the version number is an integer written under version 1 format
*       ( I6 ) it will not be read correctly by the version 2 format ( f9.2)
*       so read the version number free-format before reading the rest of
*       the line.  rwk 060526
        read( buff80(1:9),*,iostat=ioerr) rx_version  
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .      ,'lib/read_rinex_header'
     .      ,fname,'Error reading RINEX version number ',ioerr) 
          read (buff80,'(20x,a1,19x,a1,19x)',iostat=ioerr) filtyp,rxsvs
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .        ,'lib/read_rinex_header',fname
     .        ,'Error reading RINEX file type ',ioerr) 
      endif                         
      if(debug) print *,'READ_RINEX_HEADER rx_version rxsvs '
     .                                    ,rx_version,rxsvs
                
* Initialize the counter for the number of GNSS in the file
*  --for RINEX 3, we'll count these from the header; for RINEX 2
*    we have to count them as we go in read_rinex_data.
      rx_ngnss = 0 
                 
* Come here to read another line 
 10   continue                
        buff80 = ' '
        read (unit   = lun,
     .        iostat = ioerr,
     .        fmt    = '(a)') buff80
        BUFF20 = buff80(61:80)       
        call uppers(BUFF20)
        if (ioerr .ne. 0) then
            call report_stat('FATAL',prog_name,'lib/read_rinex_header'
     .            ,fname,'Error reading header for RINEX file',ioerr)
        elseif (
     .    (rx_version.le.1.1.and. buff20.eq.'                    ') .or.
*          set to accept END OF HEADER with RINEX 1 (non-standard)
*     .     (irxver.eq.2 .and. buff20(:13).eq.'END OF HEADER') .or.  
     .     ( buff20(:13).eq.'END OF HEADER') .or.
     .      ioerr.eq.-1 ) then
          goto 20

        elseif (BUFF20(1:19) .eq. 'PGM / RUN BY / DATE') then
           read (buff80,'(3a20)',iostat=ioerr)  rxprog,user,date

        elseif (BUFF20(1:7) .eq. 'COMMENT') then
          read (buff80,'(a60)',iostat=ioerr) rx_comment(irxcom+1)
*           MOD TAH 980115: +1 added to irxcom in if to ensure no bound
*                              violation.
          if( irxcom+1.lt.maxlin ) irxcom = irxcom + 1

        elseif (BUFF20(1:11) .eq. 'MARKER NAME') then
          read (buff80,'(a60)',iostat=ioerr) mark

        elseif (BUFF20(1:13) .eq. 'MARKER NUMBER') then
          read (buff80,'(a20)',iostat=ioerr) marknum

        elseif (BUFF20(1:17) .eq. 'OBSERVER / AGENCY') then
          read (buff80,'(a20,a40)',iostat=ioerr)  observer,rxagy

        elseif (BUFF20(1:19) .eq. 'REC # / TYPE / VERS') then
          read (buff80,'(3a20)',iostat=ioerr)  
     .      rx_rcvsn,rx_rcvtyp,rx_rcvsw

        elseif (BUFF20(1:12) .eq. 'ANT # / TYPE') then
          read (buff80,'(2a20)',iostat=ioerr) rx_antsn,rx_anttyp

        elseif (BUFF20(1:19) .eq. 'APPROX POSITION XYZ') then
          read (buff80,'(3f14.4)',iostat=ioerr) apx,apy,apz

        elseif (BUFF20 .eq. 'ANTENNA: DELTA H/E/N') then
          read (buff80,'(3f14.4)',iostat=ioerr) anth,ante,antn

        elseif (BUFF20 .eq. 'WAVELENGTH FACT L1/2') then
*         check to see if different factors for different SVs  
          read(buff80(18:18),'(i1)') nw
          if( nw.ne.0 ) then  
            nw1 = nsvwave + 1
            nw2 = nw1 + nw - 1
            if( nw2.gt.maxsch ) call report_stat('FATAL',prog_name
     .          ,'lib/read_rinex_header',fname,'Bad WAVELENGTH FACT',0)
            read (buff80,'(2i6,i6,7(3x,a1,i2))',iostat=ioerr)  
     .              n1,n2,nw,(isvwave(i),i=nw1,nw2) 
            do i=nw1,nw2 
              iwave1(i) = n1
              iwave2(i) = n2
            enddo
            nsvwave = nw2
          else
            read (buff80,'(2i6,i6)',iostat=ioerr) n1,n2
            nwave1 = n1
            nwave2 = n2
          endif
                      
        elseif (buff20(1:19) .eq. '# / TYPES OF OBSERV') then 
*         RINEX 1/2 will have 'TYPES OF OBSERV', RINEX 3 will have 'OBS TYPES'
*         so we don't have to check which version here
          read (buff80,'(i6)',iostat=ioerr) nobsx
          if( nobsx.gt.maxob11 ) then   
            write(message,'(a,i2,a,i2)') '# obs types =',nobsx
     .             ,' > maxob11 = ',maxob11  
            call report_stat('FATAL',prog_name,'lib/read_rinex_header'
     .               ,fname,message,0)
          elseif( nobsx.le.9 ) then      
            read (buff80,'(i6,9(4x,a2))',iostat=ioerr)
     .              nobsx,(obscodx(j)(1:2),j=1,nobsx)   
          else       
            read (buff80,'(i6,9(4x,a2))',iostat=ioerr) nobsx
     .               ,(obscodx(j)(1:2),j=1,9)     
*           nobsx > 9, continue to next line    
            if( nobsx.gt.9 ) then
              read (unit=lun,iostat=ioerr,fmt='(a)') buff80
              if( ioerr.ne.0.or.buff80(65:72).ne.'TYPES OF' ) 
     .           call report_stat('FATAL',prog_name
     .              ,'lib/read_rinex_header'
     .              ,fname,'Unexpected TYPES OF OBSERV line 2',ioerr)
*              MOD TAH 120125: Read up to 18 first and then the remainder
               read (buff80,'(6x,9(4x,a2))',iostat=ioerr)
     .               (obscodx(j)(1:2),j=10,min(18,nobsx))       
               if( ioerr.ne.0 ) 
     .            call report_stat('FATAL',prog_name
     .             ,'lib/read_rinex_header',fname
     .             ,'Error reading TYPES OF OBSERV line 2',ioerr)
            endif  
*           MOD TAH 120125   See if there are still more (19-nobsx)
            if( nobsx.gt.18 ) then  ! Read next line
              read (unit=lun,iostat=ioerr,fmt='(a)') buff80
              if( ioerr.ne.0.or.buff80(65:72).ne.'TYPES OF' ) 
     .           call report_stat('FATAL',prog_name
     .              ,'lib/read_rinex_header'
     .              ,fname,'Unexpected TYPES OF OBSERV line 2',ioerr)
              read (buff80,'(6x,9(4x,a2))',iostat=ioerr)
     .                      (obscodx(j)(1:2),j=19,min(27,nobsx))       
              if( ioerr.ne.0 ) 
     .           call report_stat('FATAL',prog_name
     .              ,'lib/read_rinex_header',fname
     .              ,'Error reading TYPES OF OBSERV line 3',ioerr)
            endif
          endif 
* *       Observables are required to be uppercase, but at least one old version (1.7)
* *       of TRRINEXO produced lowercase ids, so put this in for safety (Weiping Jiang/RWK 031016)
          do i=1,nobsx
             call uppers(obscodx(i))
          enddo                                         
          if(debug) write(*,'(a,f4.2,3i3,30(1x,a3))')
     .        'READ_RINEX_HEADER rx_version nwave1 nwave2 nobsx obscodx'
     .           , rx_version,nwave1,nwave2,nobsx,(obscodx(i),i=1,nobsx)
*         Check for cross-correlated or misreported L1 range observable
          call check_pcncod(rx_rcvtyp,rx2_pcncod)    
          if(debug) print *,'READ_RINEX_HEADER rx2_pcncod ',rx2_pcncod
          rx2_nobs = nobsx
          do i=1,rx2_nobs
            rx2_obscod(i) = obscodx(i)
          enddo    
          if(debug) write(*,'(a,i3,30(1x,a3))') 
     .        'READ_RINEX_HEADER rx2_nobs rx2_obscod '
     .        , rx2_nobs,(rx2_obscod(i),i=1,rx2_nobs)
*         Map the RINEX 2 two-character global observables to 
*         RINEX 3 three-character GNSS-dependent observables
          do i = 1,ngnss_req 
            if(debug) write(*,'(a,a1,2i4)') ' GNSS i rx2_nobs '
     .                ,  gnss_req(i),i,rx2_nobs
            call map_rinex2(gnss_req(i),i)
          enddo 
          if(debug) then 
            write(*,'(a)') 'READ_RINEX_HEADER after MAP_RINEX' 
            do i=1,ngnss_req                       
              write(*,'(a1,20(1x,a3))') 
     .           gnss_req(i),(rx_obscod(i,j),j=1,rx_nobs(i))
            enddo
          endif 
            
        elseif ( buff20(1:19).eq.'SYS / # / OBS TYPES' ) then
*         RINEX 3 
          read (buff80,'(a1,i5)',iostat=ioerr) gnssrx,nobsx   
          if(debug) then
            print *,'READ_RINEX_HEADER buff80 ',buff80
            print *,'  gnssrx  nobsx ',gnssrx,nobsx
            print *,' decoding '
          endif 
          if( nobsx.gt.maxob11 ) then   
            write(message,'(a,i2,a,i2)') '# obs types =',nobsx
     .            ,' > maxob11 = ',maxob11  
            call report_stat('FATAL',prog_name,'lib/read_rinex_header'
     .                      ,fname,message,0)
          elseif( nobsx.le.13 ) then
            read (buff80,'(1x,i5,13(1x,a3))',iostat=ioerr)
     .             nobsx,(obscodx(j),j=1,nobsx)   
          else 
            read (buff80,'(1x,i5,13(1x,a3))',iostat=ioerr)
     .            nobsx,(obscodx(j),j=1,13)     
*           nobsx > 13, continue to next line    
            if( nobsx.gt.13 ) then
              read (unit=lun,iostat=ioerr,fmt='(a)') buff80
              if( ioerr.ne.0.or.buff80(71:79).ne.'OBS TYPES' ) 
     .          call report_stat('FATAL',prog_name
     .                   ,'lib/read_rinex_header',fname
     .                   ,'Unexpected OBS TYPES line 2',ioerr)
*                   Read up to 26 first and then the remainder
              read (buff80,'(6x,13(1x,a3))',iostat=ioerr)
     .                (obscodx(j),j=14,min(26,nobsx))       
              if( ioerr.ne.0 ) 
     .           call report_stat('FATAL',prog_name
     .                 ,'lib/read_rinex_header'
     .                 ,fname,'Error reading OBS TYPES line 2',ioerr)
              endif
*             see if there are still more (27-nobsx)
              if( nobsx.gt.26 ) then  ! Read next line
                 read (unit=lun,iostat=ioerr,fmt='(a)') buff80
                 if( ioerr.ne.0.or.buff80(71:79).ne.'OBS TYPES' ) 
     .              call report_stat('FATAL',prog_name
     .                   ,'lib/read_rinex_header'
     .                   ,fname,'Unexpected OBS TYPES line 2',ioerr)
                 read (buff80,'(6x,13(1x,a3))',iostat=ioerr)
     .                  (obscodx(j),j=27,min(27,nobsx))       
                 if( ioerr.ne.0 ) 
     .             call report_stat('FATAL',prog_name
     .               ,'lib/read_rinex_header',fname
     .               ,'Error reading OBS TYPES line 3',ioerr)
                endif  
              endif  
              if(debug) write(*,'(a,f4.2,i3,30(1x,a3))') 
     .        'READ_RINEX_HEADER rx_version nwave1 nwave2 nobsx obscodx'
     .           , rx_version,nwave1,nwave2,nobsx,(obscodx(i),i=1,nobsx)
              ignss = gindex(gnssrx,gnss_req,ngnss_req) 
              if(debug) print *,'gnssrx ignss ',gnssrx,ignss 
              if(ignss.gt.0) then
                rx_ngnss = rx_ngnss + 1 
* RWK 205011: Update next 2 lines   
C               rx_gnss(rx_ngnss) = gnssrx
C               rx_nobs(rx_ngnss) = nobsx 
                rx_gnss(ignss) = gnssrx
                rx_nobs(ignss) = nobsx 
                do i=1,nobsx
                  rx_obscod(i,ignss) = obscodx(i)
                enddo      
                if(debug) write(*,'(a,2i2,6(1x,a3))') 
     .            'READ_RINEX_HEADER rx_ngnss ignss rx_obscod ',rx_ngnss
     .             ,ignss,(rx_obscod(i,ignss),i=1,rx_nobs(rx_ngnss))
              endif  

        elseif (BUFF20(1:8) .eq. 'INTERVAL') then  
               read (buff80,'(f10.0)',iostat=ioerr)  rx_inter  

        elseif (BUFF20(1:17) .eq. 'TIME OF FIRST OBS') then
              read (buff80,'(5i6,f13.7,5x,a3)',iostat=ioerr)
     .          irxyr0,irxmo0,irxdy0,irxhr0,irxmn0,rxsec0,rx_time
              call fix_y2k(irxyr0)
                
        elseif (BUFF20(1:16) .eq. 'TIME OF LAST OBS') then
              read (buff80,'(5i6,f12.6)',iostat=ioerr)
     .           irxyr1,irxmo1,irxdy1,irxhr1,irxmn1,rxsec1
             call fix_y2k(irxyr1)
  
        elseif (BUFF20(1:16) .eq. '# OF SATELLITES') then
               read (buff80,'(i6)',iostat=ioerr) irxnsv

        elseif (BUFF20(1:14) .eq. 'PRN / # OF OBS') then
              read (buff80,'(3x,a1,i2,9i6)',iostat=ioerr)
     .           rxstyp(irxsat),irxsvn(irxsat)
     .         , (irxnob(irxsat,i),i=1,nobsx)
              if( irxsat.lt.maxsch ) irxsat = irxsat + 1
            
        elseif (BUFF20(1:19) .eq. 'RCV CLOCK OFFS APPL' ) then
c             warn but do nothing more for this entry
              read (buff80,'(2x,i4)') ircvclk
              if( ircvclk.ne.0 ) 
     .           call report_stat('WARNING',prog_name
     .              ,'lib/read_rinex_header',fname
     .              ,'Receiver clock offsets applied',0)

*       endif on type of line read
        endif
                    
*           go read another line        
            goto 10                 

*         come here on end of header 
  20      continue
           
      if(debug) then 
        do i=1,nobsx
           print *,i,obscodx(i)
        enddo
      endif 

*  Check that the RINEX version is supported with expected inputs

      if (rx_version .gt. 3.05 ) then
         write (message,'(a,f9.2,a)') 'RINEX version = ',rx_version
     .                    ,' greater than 3.05'
         call report_stat('WARNING',prog_name,'lib/read_rinex_header'
     .                   ,fname,message,0)
      endif                                  
      if (filtyp .ne. 'O') then
         write (message,'(a,a1)')  'File type not O ',filtyp
         call report_stat('WARNING',prog_name,'lib/read_rinex_header'
     .                   ,fname,message,0)
      endif
      if ( anth.eq. 0.d0 ) then
         call report_stat('WARNING',prog_name,'lib/read_rinex_header'
     .                   ,fname,'Antenna height is zero',0)
      endif
 
*   Missing time time

      if( rx_time.eq.'   ' ) rx_time = 'GPS'

       
* Trap a number of special cases indicating invalid header information or observations


*    RINEX 1 or 2, change special-case observation types
      if( rx_version.lt.3.0 ) then
        l1_only = .true.        
        do itype=1,nobsx
*         P2 called C2 for serial P-code Trimble  
*         rwk mod 060111: Do this only for old observations since there is now a true C2
*         if( obscodx(itype).eq.'C2 ' ) obscodx(itype) = 'P2 '
          if( obscodx(itype).eq.'C2 ' .and.irxyr0.lt.2004 )  
     .        obscodx(itype) = 'P2 '
          if( obscodx(itype).eq.'L2 ' ) l1_only = .false.
        enddo                          
        if( l1_only ) call report_stat('WARNING',prog_name
     .     ,'lib/read_rinex_header',' '
     .     ,'Single frequency observations',0)   
      endif   

*     if multiple wavelength factors input, set single values  (RINEX 2 only?) 
      if( nsvwave.gt.0 ) then  
         nwave1 = 1
         nwave2 = 1
         do i=1,nsvwave
           if( iwave1(i).ne.nwave1 ) then  
*             L1 should always be 1 (except Macrometer II) 
             call report_stat('FATAL',prog_name
     .          ,'lib/read_rinex_header',fname
     .          ,'L1 wavelength factor not 1 ',0)
           endif
           if( iwave2(i).ne.nwave2 ) then
             call report_stat('WARNING',prog_name
     .         ,'lib/read_rinex_header',fname
     .         ,'At least one SV L2 wavelength factor =2, set all =2',0)
             iwave2(i) = 2
           endif
         enddo
         if( iwave2(1).eq.2 ) nwave2 = 2
      endif

*    Fix for column misalignment 

      if ( rx_version.lt.3.0 .and. (nwave1.lt.1.or.nwave1.gt.2) )  then
         write (message,'(a,i2)') 'Invalid nwave1 ',nwave1
         call report_stat('WARNING',prog_name,'lib/read_rinex_header'
     .                   ,fname,message,0)   
         if( nwave1.gt.0.and.nwave1.ne.2 ) then 
           call report_stat('WARNING',prog_name,'lib/read_rinex_header'
     .                     ,fname,'  Setting nwave1=1 ',0)  
           nwave1 = 1 
         endif
      endif
      if (.not.l1_only .and. rx_version.lt.3.0 .and. 
     .    (nwave2.lt.1.or.nwave2.gt.2) ) then
         write (message,'(a,i2)') 'Invalid nwave2 ',nwave2
         call report_stat('WARNING',prog_name,'lib/read_rinex_header'
     .                   ,fname,message,0)  
         if( rx_version.lt.3.0 .and. nwave1.gt.0.and.nwave1.ne.2 ) then
           call report_stat('WARNING',prog_name,'lib/read_rinex_header'
     .                    ,fname,'  Setting nwave1=1 ',0)  
           nwave1 = 1 
         endif
      endif

*    Unreasonable antenna offsets 

      if (anth .gt. 2.d0 ) then
         write (message,'(a,d12.3)') 'Antenna height > 2 meters: ',anth
         call report_stat('WARNING',prog_name,'lib/read_rinex_header'
     .                   ,fname,message,0)
      endif
      if (dabs(ante) .gt. 0.002) then
         write (message,'(a,f10.1)') 'East antenna offset > 2 mm: '
     .                               ,ante*1000
        call report_stat('WARNING',prog_name,'lib/read_rinex_header'
     .                  ,fname,message,0)

      endif
      if (dabs(antn) .gt. 0.002) then
         write (message,'(a,f10.1)') 'North antenna offset > 2 mm: '
     .                                ,antn*1000
        call report_stat('WARNING',prog_name,'lib/read_rinex_header'
     .                  ,fname,message,0)
      endif      
               
*    Bad observable in header with an L1-only Ashtech receiver

      if(rx_version.lt.3.0.and.rx_rcvtyp(1:15).eq.'ASHTECH GG-XXIV') 
     .      then
        ok_obscod = .true.
        do itype=1,nobsx
          if(obscodx(itype).eq.'P1 ') ok_obscod = .false.
          if(obscodx(itype).eq.'P2 ') ok_obscod = .false.
          if(obscodx(itype).eq.'L2 ') ok_obscod = .false.
          if(obscodx(itype).eq.'D2 ') ok_obscod = .false.
        enddo
        if( .not.ok_obscod) call report_stat('WARNING',prog_name
     .    ,'lib/read_rinex_header',fname
     .  ,'Rcvr type indicates L1-only but obs type P or L2 in header',0)
      endif                                              
 
*    Check for Ashtech codeless receiver uncorrected for the time-tag error

      if( rx_rcvtyp(1:9).eq.'ASHTECH' ) call check_fixash(lun,fname)
             
*  
      return
      end

*------------------------------------------------------------------------------------------------------*

      Subroutine MAP_RINEX2 ( gnss,index )

*     Map RINEX 2 and mislabeled RINEX 3 observables into correct RINEX 3
*     R. King September 2019

      implicit none
                                   
      include '../includes/dimpar.h'
      include 'rinex.h'                   
                          
*  index  : index in the rx_gnss array of the requested GNSS
      character*1 gnss
      integer*4 index   

*  Local
      integer*4 iblnk,i 
      character*256 message 
      logical blanks,debug/.false./


*     We use the following mappings (some are ambiguous, all to be rechecked!!):

*  RINEX2      RINEX3
*        
*          G       R      C       E       I      J
*   L1    L1C     L1C    L1X     L1X            L1C
*   L2    L2C     L2C    L2I                    L2X       
*   L5    L5X            L5X     L5X     L5A    L5X
*   L6                   L6I     L6X            L6X
*   L7                   L7I     L7X
*   L8                   L8X     L8X 
*   L9                                   L9A   
*   C1    C1C     C1C    C1X     C1X            C1C
*   C2    C2X     C2C    C2I                    C2X
*   P1    C1W     C1P            
*   P2    C2W     C2P                           
*   C5    C5X            C5X     C5X     C5A    C5X 
*   C6                   C6I     C6X            C6X 
*   C7                   C7I     C7X
*   C8                   C8X     C8X 
*   C9                                   C9A

* Warn the user that mappings may be ambiguous

      if( rx_version.lt.3.0 ) then
        write(message,'(a,a1,a)') 
     .     'RINEX 2 observable translations for GNSS ',gnss,' ambiguous'
        call report_stat('WARNING','lib/read_rinex_header','map_rinex2'
     .                  ,' ',message,0)
      endif

* Loop over all the RINEX 2 observables, translating for the input GNSS
* the 2-character codes to a 3-character code.  If a frequency is not
* applicable to this GNSS, set the code to '   ', then use these at
* the end to reset the number of observables for this GNSS  

      if(debug) write(*,'(a,2i4,a1)') 'MAP_RINEX2 index rx2_nobs gnss '
     .     ,index,rx2_nobs,gnss 
      do i=1,rx2_nobs
                    
        if( gnss.eq.'G' ) then
          if (rx2_obscod(i).eq.'L1 ') then
            rx_obscod(i,index) = 'L1C'
          elseif (rx2_obscod(i).eq.'C1') then
            rx_obscod(i,index) = 'C1C' 
          elseif (rx2_obscod(i).eq.'P1') then 
            rx_obscod(i,index) = 'C1W' 
          elseif (rx2_obscod(i).eq.'D1 ') then
            rx_obscod(i,index) = 'D1C'
          elseif (rx2_obscod(i).eq.'S1 ') then
            rx_obscod(i,index) = 'S1C'        
          elseif (rx2_obscod(i).eq.'L2 ') then 
            rx_obscod(i,index) = 'L2W'
          elseif (rx2_obscod(i).eq.'C2 ') then 
            rx_obscod(i,index) = 'C2C'
          elseif (rx2_obscod(i).eq.'P2 ') then
              rx_obscod(i,index) = 'C2W'
          elseif (rx2_obscod(i).eq.'D2 ') then
             rx_obscod(i,index) = 'D2C'
          elseif (rx2_obscod(i).eq.'S2 ') then
             rx_obscod(i,index) = 'S2C'
          elseif (rx2_obscod(i).eq.'L5 ') then 
            rx_obscod(i,index) = 'L5Q'
          elseif (rx2_obscod(i).eq.'C5 ') then
            rx_obscod(i,index) = 'C5Q'     
          elseif (rx2_obscod(i).eq.'D5 ') then
            rx_obscod(i,index) = 'D5Q'
          elseif (rx2_obscod(i).eq.'S5 ') then
            rx_obscod(i,index) = 'S5W'
          else
            rx_obscod(i,index) = '   '
          endif   

        elseif( gnss.eq.'R' ) then
          if (rx2_obscod(i).eq.'L1 ') then 
            rx_obscod(i,index) = 'L1C' 
          elseif (rx2_obscod(i).eq.'P1 ') then
            rx_obscod(i,index) = 'C1P'
          elseif (rx2_obscod(i).eq.'C1 ') then
             rx_obscod(i,index) = 'C1C'
          elseif (rx2_obscod(i).eq.'D1 ') then
            rx_obscod(i,index) = 'D1C'
          elseif (rx2_obscod(i).eq.'S1 ') then
            rx_obscod(i,index) = 'S1C'
          elseif (rx2_obscod(i).eq.'L2 ') then
            rx_obscod(i,index) = 'L2C'
          elseif (rx2_obscod(i).eq.'C2 ') then
            rx_obscod(i,index) = 'C2C'
          elseif (rx2_obscod(i).eq.'P2 ') then
            rx_obscod(i,index) = 'C2P'  
          elseif (rx2_obscod(i).eq.'D2 ') then 
            rx_obscod(i,index) = 'D2C'
          elseif (rx2_obscod(i).eq.'S2 ') then
            rx_obscod(i,index) = 'S2C'
          elseif (rx2_obscod(i).eq.'L3 ') then
            rx_obscod(i,index) = 'L3X'
          elseif (rx2_obscod(i).eq.'C3 ') then 
            rx_obscod(i,index) = 'C3X'   
          elseif (rx2_obscod(i).eq.'D3 ') then 
            rx_obscod(i,index) = 'D3C'
          elseif (rx2_obscod(i).eq.'S3 ') then
            rx_obscod(i,index) = 'S3C'
          else
            rx_obscod(i,index) = '   '
          endif 

        elseif( gnss.eq.'C' ) then             
*           RINEX 2 or 3 mislabeled Beidou second frequency (B2 but RINEX C7/L7/S7) 
*           No: assume '2' is correct and codes for Beidou B-12 
*           if( rx2_obscod(i)(2:2).eq.'2' )  rx2_obscod(i)(2:2) = '7'
          if (rx2_obscod(i).eq.'L1 ') then
            rx_obscod(i,index) = 'L1X'             
          elseif (rx2_obscod(i).eq.'C1 ') then 
            rx_obscod(i,index) = 'C1X'       
          elseif (rx2_obscod(i).eq.'D1 ') then
            rx_obscod(i,index) = 'D1X'             
          elseif (rx2_obscod(i).eq.'S1 ') then 
            rx_obscod(i,index) = 'S1X'             
          elseif (rx2_obscod(i).eq.'L2 ') then
            rx_obscod(i,index) = 'L2I'             
          elseif (rx2_obscod(i).eq.'C2 ') then
            rx_obscod(i,index) = 'C2I'  
          elseif (rx2_obscod(i).eq.'D2 ') then 
            rx_obscod(i,index) = 'D2X'             
          elseif (rx2_obscod(i).eq.'S2 ') then 
            rx_obscod(i,index) = 'S2X'             
          elseif (rx2_obscod(i).eq.'L5 ') then 
            rx_obscod(i,index) = 'L5X'             
          elseif (rx2_obscod(i).eq.'C5 ') then 
            rx_obscod(i,index) = 'C5X' 
          elseif (rx2_obscod(i).eq.'D5 ') then 
            rx_obscod(i,index) = 'D5X'             
          elseif (rx2_obscod(i).eq.'S5 ') then 
            rx_obscod(i,index) = 'S5X'              
          elseif (rx2_obscod(i).eq.'L6 ') then 
            rx_obscod(i,index) = 'L6I'             
          elseif (rx2_obscod(i).eq.'C6 ') then
            rx_obscod(i,index) = 'C6I' 
          elseif (rx2_obscod(i).eq.'D6 ') then 
            rx_obscod(i,index) = 'D6X'             
          elseif (rx2_obscod(i).eq.'S6 ') then 
            rx_obscod(i,index) = 'S6X'              
          elseif (rx2_obscod(i).eq.'L7 ') then 
            rx_obscod(i,index) = 'L7I' 
          elseif (rx2_obscod(i).eq.'C7 ') then 
            rx_obscod(i,index) = 'C7I' 
          elseif (rx2_obscod(i).eq.'D7 ') then 
            rx_obscod(i,index) = 'D7I'             
          elseif (rx2_obscod(i).eq.'S7 ') then
            rx_obscod(i,index) = 'S7I'  
          elseif (rx2_obscod(i).eq.'L8 ') then 
            rx_obscod(i,index) = 'L8X'   
          elseif (rx2_obscod(i).eq.'C8 ') then 
            rx_obscod(i,index) = 'C8X'
          elseif (rx2_obscod(i).eq.'D8 ') then 
            rx_obscod(i,index) = 'D8X'             
          elseif (rx2_obscod(i).eq.'S8 ') then
            rx_obscod(i,index) = 'S8X'    
          else
            rx_obscod(i,index) = '   '
          endif

        elseif( gnss.eq.'E' ) then
          if (rx2_obscod(i).eq.'L1 ') then 
            rx_obscod(i,index) = 'L1X'  
          elseif (rx2_obscod(i).eq.'C1 ') then 
            rx_obscod(i,index) = 'C1X'  
          elseif (rx2_obscod(i).eq.'D1 ') then
            rx_obscod(i,index) = 'D1X'  
          elseif (rx2_obscod(i).eq.'S1 ') then
            rx_obscod(i,index) = 'S1X'   
          elseif (rx2_obscod(i).eq.'L5 ') then 
            rx_obscod(i,index) = 'L5X' 
          elseif (rx2_obscod(i).eq.'C5 ') then
            rx_obscod(i,index) = 'C5X'        
          elseif (rx2_obscod(i).eq.'D5 ') then
            rx_obscod(i,index) = 'D5X' 
          elseif (rx2_obscod(i).eq.'S5 ') then 
            rx_obscod(i,index) = 'S5X'            
          elseif (rx2_obscod(i).eq.'L6 ') then 
            rx_obscod(i,index) = 'L6X' 
          elseif (rx2_obscod(i).eq.'C6 ') then
            rx_obscod(i,index) = 'C6X'        
          elseif (rx2_obscod(i).eq.'D6 ') then 
            rx_obscod(i,index) = 'D6X' 
          elseif (rx2_obscod(i).eq.'S6 ') then 
            rx_obscod(i,index) = 'S6X' 
          elseif (rx2_obscod(i).eq.'L7 ') then
            rx_obscod(i,index) = 'L7X' 
          elseif (rx2_obscod(i).eq.'C7 ') then
            rx_obscod(i,index) = 'C7X'        
          elseif (rx2_obscod(i).eq.'D7 ') then
            rx_obscod(i,index) = 'D7X' 
          elseif (rx2_obscod(i).eq.'S7 ') then
            rx_obscod(i,index) = 'S7X' 
          elseif (rx2_obscod(i).eq.'L8 ') then
            rx_obscod(i,index) = 'L8X' 
          elseif (rx2_obscod(i).eq.'C8 ') then
            rx_obscod(i,index) = 'C8X'        
          elseif (rx2_obscod(i).eq.'D8 ') then
            rx_obscod(i,index) = 'D8X' 
          elseif (rx2_obscod(i).eq.'S8 ') then
            rx_obscod(i,index) = 'S8X' 
          else
            rx_obscod(i,index) = '   ' 
          endif

        elseif( gnss.eq.'I' ) then 
          if (rx2_obscod(i).eq.'L5 ') then
            rx_obscod(i,index) = 'L5A'   
          elseif (rx2_obscod(i).eq.'C5 ') then
            rx_obscod(i,index) = 'C5A'     
          elseif (rx2_obscod(i).eq.'D5 ') then
            rx_obscod(i,index) = 'D5X' 
          elseif (rx2_obscod(i).eq.'S5 ') then
            rx_obscod(i,index) = 'S5X' 
          elseif (rx2_obscod(i).eq.'L9 ') then
            rx_obscod(i,index) = 'L9A'
          elseif (rx2_obscod(i).eq.'C9 ') then
            rx_obscod(i,index) = 'C9A'   
          elseif (rx2_obscod(i).eq.'D9 ') then
            rx_obscod(i,index) = 'D9X' 
          elseif (rx2_obscod(i).eq.'S9 ') then
            rx_obscod(i,index) = 'S9X'   
          else
            rx_obscod(i,index) = '   '
          endif

        elseif( gnss.eq.'J' ) then 
          if (rx2_obscod(i).eq.'L1 ') then
            rx_obscod(i,index) = 'L1C'   
          elseif (rx2_obscod(i).eq.'C1 ') then
             rx_obscod(i,index) = 'C1C'     
          elseif (rx2_obscod(i).eq.'D1 ') then
            rx_obscod(i,index) = 'D1C' 
          elseif (rx2_obscod(i).eq.'S1 ') then
            rx_obscod(i,index) = 'S1C'
          elseif (rx2_obscod(i).eq.'L2 ') then
            rx_obscod(i,index) = 'L2X'   
          elseif (rx2_obscod(i).eq.'C2 ') then
            rx_obscod(i,index) = 'C2X'     
          elseif (rx2_obscod(i).eq.'D2 ') then
            rx_obscod(i,index) = 'D2X' 
          elseif (rx2_obscod(i).eq.'S2 ') then
            rx_obscod(i,index) = 'S2X'
          elseif (rx2_obscod(i).eq.'L5 ') then
            rx_obscod(i,index) = 'L5X'   
          elseif (rx2_obscod(i).eq.'C5 ') then
            rx_obscod(i,index) = 'C5X'     
          elseif (rx2_obscod(i).eq.'D5 ') then
            rx_obscod(i,index) = 'D5X' 
          elseif (rx2_obscod(i).eq.'S5 ') then
            rx_obscod(i,index) = 'S5X'           
          elseif (rx2_obscod(i).eq.'L6 ') then
            rx_obscod(i,index) = 'L6X'   
          elseif (rx2_obscod(i).eq.'C6 ') then
            rx_obscod(i,index) = 'C6X'     
          elseif (rx2_obscod(i).eq.'D6 ') then
            rx_obscod(i,index) = 'D6X' 
          elseif (rx2_obscod(i).eq.'S6 ') then
            rx_obscod(i,index) = 'S6X'
          else
            rx_obscod(i,index) = '   '
          endif 
       
        endif 
        if(debug) write(*,'(a,2(i3,1x,a3))') 
     .    '  i rx2_obscod index rx_obscod '
     .       , i,rx2_obscod(i),index,rx_obscod(i,index) 
      enddo    

* Now reset the RINEX 3 number and observables for this GNSS 
                 
      rx_nobs(index) = rx2_nobs     
      blanks = .true.               

      do while( blanks )
       
   10   iblnk = 0 
        do i=1,rx_nobs(index)
          if( rx_obscod(i,index).eq.'   ' ) then 
            iblnk=i
          endif
        enddo  
        if( iblnk.ne.0 ) then 
          rx_nobs(index) = rx_nobs(index) - 1
          do i=iblnk,rx_nobs(index)
            rx_obscod(i,index) = rx_obscod(i+1,index)
          enddo 
          goto 10 
        else
          blanks = .false.
        endif

      enddo                                       
      if(debug) write(*,'(a,i3,50(1x,a3))') 
     .    '  rx_nobs rx_obscod '
     .    , rx_nobs(index),(rx_obscod(i,index),i=1,rx_nobs(index))
      return
      end

*--------------------------------------------------------------------- 
        
      Subroutine READ_RINEX_DATA( lun,ngnss_req,gnss_req,ferr,fend )
                             
*     Read the epoch and SV lines of a RINEX 2 or RINEX 3 file
*     R. King 6 December 2018 from old lib/rrinex.f

      implicit none           

      include '../includes/dimpar.h'
      include 'rinex.h'

*  Input  

*   lun                - logical unit number, assumed open 
*   ngnss_req          - number of GNSS requested
*   gnss_req(maxgnss)  - list of GNSS requested ( c*1 ) 
      integer*4 lun,ngnss_req
      character*1 gnss_req(*)

*  Output

*   ferr - true if any read errors detected        
*   fend - true if end-of-file
      logical ferr,fend   
                   
*  Values stored in includes/rinex.h; all only for SVs included in the gnss_req list

*    i*4 rx_ngnss               - Of those requested, the number of systems on the RINEX file 
*    c*1 rx_gnss                - List of systems on the RINEX file                       
*    i*4 rx_jd                  - PEP JD of epoch read
*    r*8 rx_t                   - Second-of-day of epoch read
*    i*4 rx_nprn                - number of SV data returned
*    c*3 rx_prn(maxsch)         - PRNs selected 
*    r*8 rx_obs(maxsch,maxob11)  - observations return
*    i*4 rx_issi(maxsch,maxob11) - RINEX signal strengh indicator
*    i*4 rx_illi(maxsch,maxob11) - loss-of-lock indicator
*    log rx_data_ok(maxsch)     - False if a read error 

* Local

*   iyr,ihr,mon,iday,ihr,min, sec, jd - Time on the epoch line  
      integer*4 iyr,month,iday,ihr,min,jd             
      real*8    sec    
*   anth, antn, ante - Antenna offsets from an event record (ignored)
      real*8 anth,antn,ante 
*   idflag               - 0 if epoch ok, 1 if power failure > 1 special event
      integer*4 idflag
*   nprnrx - number of SVs in the record 
      integer*4 nprnrx 
*   prnrx = SVs in record
      character*3 prnrx(maxsch)  
*   observation values each PRN (RINEX 2) 
      real*8 obsrx(maxob11)
      integer*4 issirx(maxob11),illirx(maxob11)
*   iprn - counter for selected SVs  (becomes rx_nprn in rinex.h)
      integer*4 iprn
*   recid - 1st column in a RINEX 3 file (e.g. '>' for an epoch line)
      character*1 recid 
*   nl - number of lines for RINEX 2 data each SV
      integer*4 nl
*   nel, iend  - end indices for RINEX line (increments by 5 until last line)
      integer*4 nel,iend 
*   loop variables
      integer*4 i,j,k,ii,il 
*   string length from rcpar 
      integer*4 len  
*   line buffer
      character*512 line
*   report_stat message
      character*256 message 
*   error flag, and number
      integer*4 ioerr,nerr
*   name of calling program
      character*80 prog_name   
*   time variable
      real*8 t               
*   event - true if an event has occurred
      logical event         
      logical debug/.false./ 

* Function
      integer*4 rcpar,julday

      if(debug) write(*,'(a,i3,50(1x,a2))') 
     .    'READ_RINEX_DATA ngnss_req gnss_req '
     .   ,ngnss_req,(gnss_req(i),i=1,ngnss_req)
                
* Get the calling program
      len = rcpar(0,prog_name)

* Read the next epoch line into a buffer and check

      fend = .false.                     
  100 read (lun,'(a)',iostat=ioerr,end=999) line                       
      event = .false.
      if( debug ) write(*,*) 'READ_RINEX_DATA read ',line
                          
*     Five cases:
*        1. Error on read:  increment the error count error read the next line 
*        2. Error decoding epoch line:  increment the error count and read the next line
*        3. Good epoch line but number records (SVs) exceeds channels: 
*           Exit the program  
*        4. Good epoch line but with event flag set:  read the event info and then
*           the SV lines
*        5. Normal epoch line: decode and exit
*     After any read error, call rxerr (at the end of this file) to count the errors
*     and print out the line; after 100, stop the run.

* Case 1: error on 'a' read
      if( ioerr.ne.0 ) then
        call rxerr (lun,line,ioerr,nerr)
        if(debug) write(*,*) 'READ_RINEX_DATA Error on read, ioerr '
     .    ,ioerr
        goto 100
      endif
           
* Decode the epoch line 

      if( rx_version.lt.3.0 ) then   
        read (unit   = line
     .       ,fmt    = '(5i3,f11.7,i3,i3)'
     .       ,iostat = ioerr
     .       ,end = 999) iyr,month,iday,ihr,min,sec,idflag,nprnrx
        call fix_y2k(iyr)                                              
        if(debug) write(*,*) 'RINEX 2 decode epoch: '
     .     ,iyr,month,iday,ihr,min,sec,idflag,nprnrx
      else
        read (unit   = line,
     .        fmt    = '(a1,1x,i4,4(1x,i2.2),f11.7,2x,i1,i3)',
     .        iostat = ioerr, end = 999 ) 
     .        recid,iyr,month,iday,ihr,min,sec,idflag,nprnrx
        if(debug) write(*,*) 'RINEX 3 decode epoch: '
     .        , recid,iyr,month,iday,ihr,min,sec,idflag,nprnrx
      endif

* Case 2: error decoding the epoch line
      if( ioerr.ne.0 ) then 
        call rxerr(lun,line,ioerr,nerr) 
        if(debug) write(*,*) 
     .    'READ_RINEX_DATA Error on epoch decode, ioerr ',ioerr
        goto 100 
      endif

*     decoded ok 

* Case 3: too many SV records for GAMIT dimensions
      if( idflag.eq.0. and. nprnrx.gt.maxsch ) then 
        write(message,'(a,i2,a,i2,a,i4,4i3)')
     .    'Number of SVs on RINEX file (=',nprnrx
     .    ,') exceeds array dimensions (maxsch=',maxsch
     .    ,') ',iyr,month,iday,ihr,min
          call report_stat('FATAL',prog_name,'lib/read_rinex_data',' '
     .                    ,message,0)    
      endif 
                           
* Case 4:  RINEX 2 or 3 event information  
      if( idflag.gt.1 ) then
        if(debug) write(*,*) 
     .    'READ_RINEX_DATA Event line nprnrx (nline):', nprnrx
        write(message,'(2a)') 'RINEX file contains version 2 data-'
     .    ,'new antenna offset keys HISUB call; all else now ignored'
        call report_stat('WARNING',prog_name,'lib/read_rinex_data',' '
     .                  ,message,0) 
        write(message,'(a,i4,1x,4i3)') 
     .      'Epoch: ',iyr,month,iday,ihr,min   
        call report_stat('WARNING',prog_name,'lib/read_rinex_data',' '
     .                  ,message,0) 
        if( rx_version.lt.3.0 ) then
c         RINEX 2                                                                          
          do i=1,nprnrx 
c           here nprn indicates number of event records, not SVs
            read(lun,'(a80)',iostat=ioerr,end=999) line
            call report_stat('WARNING',prog_name,'lib/read_rinex_data'
     .                      ,' ',line,0)
            if (line(61:80) .eq. 'ANTENNA: DELTA H/E/N') then
              read (line,'(3f14.4)',iostat=ioerr,end=999) anth,ante,antn
            endif
          enddo
* MOD TAH 201113: By adding the condition below that case 5 should 
*         onlu execute when idflag==0, this goto is not needed and
*         should never by used. 
CCCCCC    goto 100
        else      
c         RINEX 3
* MOD TAH 201113: Not clear why this code operates this way since
*         nprnrx should contain number of records in event??
          event = .true.
          do while( event ) 
            read (lun,'(a)',iostat=ioerr,end=999) line
            if(ioerr.eq.-1 ) then
              fend = .true.
              return
            elseif(ioerr.ne.0) then
              call report_stat('FATAL',prog_name,'read_rinex_epoch',' '
     .                      ,'Error reading RINEX 3 event lines',ioerr)
            else
              if(debug) write(*,*) 'RINEX 3 event line: ',line
              if(line(1:1).eq.'>' ) then 
                event = .false.
              else
                if (line(61:80) .eq. 'ANTENNA: DELTA H/E/N') then
                  read (line,'(3f14.4)',iostat=ioerr,end=999) 
     .                  anth,ante,antn
                  call report_stat('WARNING',prog_name
     .         ,'lib/read_rinex_data',' ','EVENT: new antenna offset',0)
                endif
              endif
            endif
          enddo
          backspace(lun,iostat=ioerr ) 
c         endif on RINEX 2 or 3
        endif
c     endif on event line
      endif

* Case 5: normal epoch record
 
* Check for valid time tag and convert the times   
* MOD TAH 201113: The rest of this code should only run if the
*     idflag == 0 meaning this is data record.  Adding this
*     condition to the if statement, removes need for goto
*     above.    
      if( rx_version.lt.3.0 ) call fix_y2k(iyr)
      if (iyr    .ge. 1985 .and. iyr    .le. 2030 .and.
     .    month  .ge.  1   .and. month  .le. 12 .and.
     .    iday   .ge.  1   .and. iday   .le. 31 .and.
     .    ihr    .ge.  0   .and. ihr    .le. 23 .and.
     .    min    .ge.  0   .and. min    .le. 60 .and. idflag.le.1 ) then
        jd = julday(month,iday,iyr)
        t = 3600.d0*ihr + 60.d0*min + sec

* RINEX 2:   Read and decode the data lines - skip the lines that are not the requested GNSS 
         
        if( rx_version.lt.3.0) then 

          if( nprnrx.gt.12 ) then
            iend = 12
          else
            iend = nprnrx
          endif   
*         Decode the 1st epoch line again, then read the others, saving the svids
          read(line,'(5i3,f11.7,i3,i3,12(a3))',iostat=ioerr,end=999)
     .           iyr,month,iday,ihr,min,sec,idflag,nprnrx
     .           ,(prnrx(i),i=1,iend)   
          call fix_y2k(iyr)  
          if( debug ) then
             write(*,'(a)') 'RINEX 2 epoch line ',line
             write(*,'(a,i5,4i3,f11.7,i3,i3,12(a3))') 'Decoded: '
     .          ,iyr,month,iday,ihr,min,sec,idflag,nprnrx
     .          ,(prnrx(i),i=1,iend)    
          endif
          if( ioerr.ne.0 ) then   
            call rxerr(lun,line,ioerr,nerr)  
            if(debug) write(*,*) 
     .         'READ_RINEX_DATA Error decoding epoch line',ioerr
            if( ioerr.eq.-1 ) then
              fend = .true.
              return
            else
              goto 100 
            endif 
          endif  
*         If there is more than 1 line, read and decode these now
          if(debug) write(*,*) 'nprnrx ',nprnrx
     .         ,' > 12 read another SV line'
          if( nprnrx.gt.12) then
            nl = int((nprnrx-1)/12+1)
            do il = 2,nl
              if( nprnrx.gt. 2*il ) then
                iend = il*12
              else
                iend = nprnrx
              endif 
              if(debug)  print *,'READ_RINEX_DATA nprnrx nl il iend '
     .            ,nprnrx,nl,il,iend 
              read (unit=lun,iostat=ioerr,fmt='(a80)',end=999) line 
              if(debug) write(*,*) 'Read another SV ID line ',line
              if (ioerr.ne.0 ) then
                call rxerr (lun,line,ioerr,nerr)
                if(debug) write(*,*) 
     .             'READ_RINEX_DATA Error on SV ID line ioerr ',ioerr
                if( ioerr.eq.-1 ) then
                  fend = .true.
                  return
                else
                  goto 100 
                endif 
              endif  
              read(line,'(32x,12(a3))',iostat=ioerr,end=999) 
     .            (prnrx(i),i=(il-1)*12+1,iend)
              if( ioerr.ne.0 ) then 
                call rxerr(lun,line,ioerr,nerr) 
                if( ioerr.eq.-1 ) then
                  fend = .true.
                  return
                else
                  goto 100 
                endif 
              endif
            enddo
*         endif on nprnrx > 12 for RINEX 2
          endif            
          do i=1,nprnrx
*           add the GNSS if an old-style RINEX file
            if(prnrx(i)(1:1).eq.' ') prnrx(i)(1:1) = 'G'
          enddo    
                                                                                 
*         We'd like to know which of the GNSS requested are actually
*         available on the RINEX file (rx_ngnss and rx_gnss(i)) in 
*         rinex.h).  For RINEX 3, we have this information in the 
*         header, but for RINEX 2, we have to gather it from the 
*         data records.  So call this routine, appended to the bottom
*         of read_rinex_data, to increment the rinex.h variables.
          call count_gnss( ngnss_req,gnss_req,nprnrx, prnrx )

*         We have to hope that the observation records are in the order
*         and number expected (RINEX error recovery is difficult)

*         Read and decode the data lines 
          nl = int((rx2_nobs-1)/5)+1  ! Number of lines for observables, each SV
          iprn = 0      ! counter for prns for requested GNSS SVs 
          do i=1,nprnrx
            rx_data_ok(i) = .true.
            do  k = 1, nl
*             Get end range of data to be read 
              nel = min0(5*k,rx2_nobs)  ! Need to use explicit name if
                                      ! min is re-defined to be minutes.
*             Read line from file
              read(lun,'(a80)',iostat=ioerr,end=999) line  
              if(debug) write(*,*) 
     .             'READ_RINEX_DATA Data i k nl nel ioerr line '
     .                ,i,k,nl,nel,ioerr,line
              if(ioerr.eq.-1) then
                fend = .true.
                return
              elseif(ioerr .ne. 0) then   
                call rxerr (lun,line,ioerr,nerr)
*               call this sat bad, because of file error
                rx_data_ok(i) = .false. 

              else  ! Read OK  
                read(line,'(5(f14.3,i1,i1))',iostat=ioerr,end=999)
     .            (obsrx(j),illirx(j),issirx(j)
     .                        ,j=(k-1)*5+1,nel)
                if(debug) write(*,'(a,i5,1x,a3,5(f14.3,i1,i1))')
     .             'RINEX 2 data line ioerr prnrx, decoded:'
     .                  , ioerr,prnrx(i)
     .                ,(obsrx(j),illirx(j),issirx(j)
     .                       ,j=(k-1)*5+1,nel)
                if(ioerr.ne.0) then 
                   call rxerr (lun,line,ioerr,nerr) 
                endif
              endif     
*           enddo over lines read for one SV - see if requested
            enddo 
            do ii=1,ngnss_req
              if( prnrx(i)(1:1).eq.gnss_req(ii) ) then          
                iprn = iprn + 1
                rx_prn(iprn) = prnrx(i)
                do j=1,rx2_nobs
                  rx_obs(iprn,j) = obsrx(j)
                  rx_issi(iprn,j) = issirx(j) 
                  rx_illi(iprn,j) = illirx(j)
                enddo 
c                if(debug) write(*,'(a,5i3,30(f14.3))') 
c     .           'i iprn rx2_nobs nl k rx_obs '
c     .           ,i,iprn,rx2_nobs,nl,k,(rx_obs(iprn,j),j=1,rx2_nobs)
                if(debug) write(*,'(a,a3,32(f14.3,2x,2i3))') 
     .          'Selected ',rx_prn(iprn),(rx_obs(iprn,j)
     .                ,rx_illi(iprn,j),rx_issi(iprn,j),j=1,rx2_nobs)
              endif
            enddo
       
*         enddo over SVs on the observation lines
          enddo   
          rx_nprn = iprn                               
          if(debug) then    
            write(*,'(a,2i3)') 'READ_RINEX_DATA rx_nprn rx2_nobs '
     .                       ,  rx_nprn,rx2_nobs 
            write(*,'(50(1x,a3))') (rx2_obscod(j),j=1,rx2_nobs)
            do i=1,rx_nprn
            write(*,'(a3,32f14.3)') rx_prn(i),(rx_obs(i,j),j=1,rx2_nobs)
            enddo
          endif 
       
* RINEX3 :  Read and decode the data lines - skip the lines that are not the requested GNSS
    
        else

          iprn = 0 
          do i=1,nprnrx
            read(lun,'(a)',iostat=ioerr,end=999) line
            if(debug) print *,'data line ',line 
            if (ioerr.eq.-1 ) then
              fend = .true.
              if(debug) print *,'READ_RINEX_DATA read line  ioerr '
     .            ,ioerr,line
            elseif (ioerr.ne.0 ) then
              call report_stat('FATAL',prog_name,'lib/read_rinex_data'
     .                   ,' ','Error reading data line',ioerr )
            else
              read(line,'(a3)',iostat=ioerr,end=999) prnrx(i)
              if(debug)  print *,'ioerr read prnrx ok ',ioerr,prnrx(i)
              if( ioerr.ne.0 ) then     
                call rxerr (lun,line,ioerr,nerr)
              else  
                do j=1,ngnss_req
                  if(prnrx(i)(1:1).eq.gnss_req(j)) then          
                    iprn = iprn + 1    
                    read(line,'(a3,31(f14.3,i1,i1))'
     .                 ,iostat=ioerr,end=999) 
     .                   rx_prn(iprn)
     .                 , (rx_obs(iprn,k),rx_illi(iprn,k),rx_issi(iprn,k)
     .                     ,k=1,rx_nobs(j))
                    if(ioerr.ne.0) then
                       call rxerr (lun,line,ioerr,nerr) 
                    endif
                    if( debug ) write(*,'(a,a3,31(f14.3,i1,i1))') 
     .               'Selected ',rx_prn(iprn),(rx_obs(iprn,k)
     .                ,rx_illi(iprn,k),rx_issi(iprn,k),k=1,rx_nobs(j))
                  endif
                enddo
              endif
            endif
          enddo

*       endif on RINEX version
        endif   

*       endif on valid epoch, with observations read into the rinex.h arrays
        endif     

*       Reset the number of SVs to only those selected                      
        rx_nprn = iprn

        return                         

* Come here on end of file
  999  fend = .true.
       return 
        end 

*-----------------------------------------------------------------------

      Subroutine COUNT_GNSS ( ngnss_req,gnss_req,nprnrx,prnrx )

*     For RINEX 2, the GNSS available is not recorded on the header, so 
*     increment the name and codes as they are encountered for the first time

      implicit none              

      include '../includes/dimpar.h'
      include 'rinex.h'

* Input     
*   ngnss_req  - number of GNSS requested
*   gnss_req   - list of GNSS requested                           
*   nprnrx     - number of PRNs in the data record
*   prnrx      - PRNs in the data record (c*3)
      integer*4 ngnss_req,nprnrx
      character*3  gnss_req(*),prnrx(maxsch)

* Output - in rinex.h
*   rx_ngnss          - Current count of the GNSS on the file
*   rx_gnss(*)  - Current GNSS codes (c*1)

* Local
      integer*4 newgnss,i,j    
      logical debug/.false./

* Function
      logical use_gnss

* Need to initialize the first GNSS to make the logic work
      if( rx_ngnss.eq.0 ) then
        rx_ngnss = 1
        rx_gnss(1) = prnrx(1)(1:1)
      endif  
      if(debug) print *,'Entering COUNT_GNSS nprnrx prnrx '
     .  ,nprnrx,(prnrx(i),i=1,nprnrx)    
      if(debug) print *,'  ngnss_req gnss_req '
     .   ,ngnss_req,(gnss_req(i),i=1,ngnss_req)
      do i=1,nprnrx                 
        if(debug) print *,' prnrx(1:1) ',i,prnrx(i)(1:1)           
*  Is the PRN GNSS in the request list?
        if( use_gnss(prnrx(i)(1:1),ngnss_req,gnss_req) ) then
*  If so, is it already in the list found on the RINEX file?
          if( use_gnss(prnrx(i)(1:1),rx_ngnss,rx_gnss) ) then
            if(debug) print *,'already in list rx_ngnss ',rx_ngnss
            continue
          else                                                                
            rx_ngnss = rx_ngnss + 1  
            rx_gnss(rx_ngnss) = prnrx(i)(1:1)
            if(debug) print *,'add to list rx_ngnss rx_gnss '
     .                ,rx_ngnss,rx_gnss(rx_ngnss)
          endif
        endif
      enddo

      return
      end


*-------------------------------------------------------------------------

      Subroutine RXERR (lunit,badline,ioerr,nerr)

*     Handle a RINEX file error

*        Input:
*          lunit      logical unit for RINEX file
*          badline    offending line in file
*          ioerr      returned error message
*          nerr       total number of errors

      implicit none

      integer*4 lunit,ioerr,nerr,inqerr
      character*(*) badline
      integer*4 maxerr,rcpar,len

      character*80  prog_name,fname
      character*320  message

      data maxerr /100/

* Get calling program name and file name for report_stat

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
*--------------------------------------------------------------

      Function GINDEX(gnss,gnss_array,n)

*     Get the index in gnss_array of the input gnss
*     R. King 19 January 2019

      implicit none

      integer*4 gindex,n,i 
                                          
      character*1 gnss,gnss_array(*)
                 
      gindex = 0 
      do i=1,n
       if(gnss.eq.gnss_array(i)) gindex = i
      enddo     
      return
      end

*----------------------------------------------------------

      Function ISVARRAY( prn,prn_array,n )
      
* Determine the position in the PRN array of a satellite
* R. King August 2018
                          
      integer isvarray,n,i
      character*3  prn,prn_array(*)

cd      print *,'ISVARRAY  prn n ',prn,n    
cd      write(*,'(a,128(1x,a3))') 'PRN array   ',(prn_array(i),i=1,n)
      isvarray = 0
      do i=1,n
        if(prn.eq.prn_array(i)) then
          isvarray = i
        endif
      enddo 
cd      print *,'  isvarray ',i 
      return
      end 

*--------------------------------------------------------------------

      Subroutine CHECK_FIXASH( lurnx,fname )

*   Check the header of the RINEX file to see if the receiver is 
*   a codeless Ashtech that has not had its time-tag corrected.

*   R. King June 2019 

*   In 10.7 and earlier, the corrections could be applied in 
*   MODEL, but now we're forcing the user to correct the RINEX
*   file with script sh_fixash before running GAMIT.

*   The receivers with a bad time tag, ~45 microsecond offset,
*   are usually identified by receiver types L-XII, LM-XII, LM-XIIC, 
*   or M-XII with firmware 6A or 7A.  Receivers LM-XII3 and LM-XII
*   firmware 7AP6, 6C, 6G, or 6M have P2 available and do not have a 
*   timing error.  

                                      
      implicit none             
                                  
      include '../includes/dimpar.h'
      include 'rinex.h'

      integer*4 lurnx,mchkey,i
  
      character*(*) fname
      character*7 codeless_rcvrs(4)/'L-XII  ','LM-XII ','LM-XIIC'
     .                             ,'M-XII  '/
      character*18 ashfixline
      character*20 rctype,rcvers,buff20 
      character*80 prog_name,buff80,message 
                
      logical rcmatch,versmatch,corrected
     

* See if this is a candidate receiver and firmware

      rcmatch = .false.                         
      versmatch = .false.
      do i=1,4                                   
        if( mchkey(rctype,codeless_rcvrs(i),20,7).gt.0 ) 
     .     rcmatch = .true.
      enddo
      if( rcmatch ) then                                        
        if(rcvers(1:2).eq.'6A'.or.
     .    (rcvers(1:2).eq.'7A'.and.rcvers(1:4).ne.'7AP6') ) then 
          versmatch = .true.
        endif
      endif                            

* Check for the message 'ASHTECH L-XII clock offset fixed by  45   usec;
* written by sh_fixash as a header COMMENT.  

      if( rcmatch.and.versmatch ) then
        corrected = .false.          
        ashfixline = 'clock offset fixed'
        rewind(lurnx)
        do i=1,irxcom
          if(mchkey(rx_comment,ashfixline,60,18).gt.0) corrected=.true.
        enddo
        if(corrected) then
          write(message,'(a)') 
     .      'Ashtech codeless with corrected time tag'  
          call report_stat('WARNING',prog_name,'lib/check_fixash'
     .                    ,fname,message,0)
        else  
          write(message,'(a)') 
     .      'Ashtech codeless with uncorrected time tag'
          call report_stat('FATAL',prog_name,'lib/check_fixash'
     .                    ,fname,message,0)
        endif   
      endif
 
      return
      end                          

*----------------------------------------------------------------------

      Subroutine CHECK_PCNCOD( rcvtyp,pcncod )

*     For older receivers and RINEX 2, determine if the L1 range observable is
*     cross-correlating (C1C) or reported as C1 (C1C) but is really P1 (C1W).
*     This routine replaces the entries in the now-defunct rcvant.dat file.

*     R. King  16 Novebmer 2019

      implicit none

* Input:  rcvtyp   c*20   RINEX-standard receiver name
      character*20 rcvtyp

* Output: pcncod   c*1    Code indicating type of pseudorange observable
      character*1  pcncod 
*    
*  'P' : receiver is cross-correlating and requires correction of P2' and C1 
*          Rogue SNR, Trimble 4000, etc.
*  'C' : receiver is non-cross-correlating but reports C1 instead of P1
*         Trimble 4700, 5700, Leica RS500, CRS1000, SR9600, etc. unless AS is off
*  'N' : receiver is non-cross-correlating and reports true P1, P2 
           
* Cross-correlating receivers

      if( rcvtyp(1:9).eq.'ROGUE SNR' .or. 
     .    rcvtyp(1:17).eq.'AOA ICS-4000Z   ' .or.
     .    rcvtyp(1:12).eq.'TRIMBLE 4000' .or. 
     .    rcvtyp(1:12).eq.'TRIMBLE 4400' .or.
     .    rcvtyp(1:12).eq.'TRIMBLE 4600' .or.
     .    rcvtyp(1:12).eq.'TRIMBLE 4800' .or. 
     .    rcvtyp(1:12).eq.'TRIMBLE 7400' .or.
     .    rcvtyp(1:7) .eq.'MINIMAC' .or.
     .    rcvtyp(1:9) .eq.'TOPCON TT' .or.
     .    rcvtyp(1:12).eq.'TOPCON GP-DX' .or.
     .    rcvtyp(1:11).eq.'LEICA SR299' .or.
     .    rcvtyp(1:11).eq.'LEICA SR399' .or.
     .    rcvtyp(1:16).eq.'SPP GEOTRACER100' .or. 
     .    rcvtyp(1:14).eq.'GeoTracer 2000' ) then
        pcncod = 'P'
          
* Non-cross-correlating but reports C1 

      elseif( rcvtyp(1:7).eq.'TRIMBLE' .or.
     .        rcvtyp.eq.'ASHTECH S-XII       ' .or.
     .        rcvtyp.eq.'ASHTECH L-XII       ' .or.
     .        rcvtyp.eq.'ASHTECH M-XII       ' .or.
     .        rcvtyp.eq.'ASHTECH LM-XII3     ' .or.
     .        rcvtyp.eq.'ASHTECH P-XII3      ' .or.
     .        rcvtyp.eq.'ASHTECH 3DF-XXIV    ' .or.
     .        rcvtyp.eq.'ASHTECH D-XII       ' .or.
     .        rcvtyp.eq.'ASHTECH MS-XII      ' .or.
     .        rcvtyp.eq.'ASHTECH SUPER-CA    ' .or.
     .        rcvtyp.eq.'ASHTECH 3DF-XXIV    ' .or.
     .        rcvtyp.eq.'TOPCON GP-R1        ' .or.
     .        rcvtyp.eq.'TOPCON GP-R1D       ' .or.
     .        rcvtyp.eq.'TOPCON GP-R1DP      ' .or.
     .        rcvtyp.eq.'TOPCON GP-SX1       ' .or.
     .        rcvtyp(1:6).eq.'CHC N7' .or.
     .        rcvtyp.eq.'HEM P320 ECLIPSE II ' .or.
     .        rcvtyp(1:5).eq.'LEICA' .or.
     .        rcvtyp(1:3).eq.'NOV' .or.       
     .        rcvtyp(1:3).eq.'SOK' .or.
     .        rcvtyp(1:3).eq.'SPP' .or. rcvtyp(1:7).eq.'SPECTRA' .or.
     .        rcvtyp(1:7).eq.'RAYMAND' .or.
     .        rcvtyp(1:6).eq.'RVMETS' .or.
     .        rcvtyp(1:8).eq.'CHAMPION' .or.         
     .        rcvtyp(1:3).eq.'CMC' ) then
        pcncod = 'C'

      else
        pcncod = 'N'

      endif   

      return
      end
                     
*-------------------------------------------------------------------------

      Logical function USE_GNSS(gnss,ngnss_sel,gnss_sel)

*     Determine if the input GNSS code is included in the selection array
*     R. King 30 August 2018

      implicit none

      character*1 gnss,gnss_sel(*)

      integer*4 ngnss_sel,i
       
      use_gnss = .false.
      do i=1,ngnss_sel
        if( gnss.eq.gnss_sel(i) ) then
          use_gnss = .true.
        endif
      enddo
      return
      end
        
     

 


      


