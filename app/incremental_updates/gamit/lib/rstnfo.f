      Subroutine RSTNFO(lun, sitcod, yr_in, doy_in, sod_in
     .                 , span, sname, anth, antn, ante, antdaz 
     .                 , rcvcod, antcod
     .                 , htcod, radome, swver, rcvers, rcvrsn, antsn
     .                 , istart, istop )   

c     Subroutine to read the new-style station.info table.

c     R. King  27 April 2001; last modified R. King 8 September 2010

c MOD TAH 200203: Added AntDAZ token for antenna  Alignment from True N
c MOD TAH 200524: Updated to allow a gap between end and start times of
c                 up 6hrs.  This happens with firmware updates where the
c                 gap can be 1-minute to a few hours.  There is probably no
c                 data in gap but GAMIT is still generating epochs during
c                 the gap (makex and maybe model fatals). Cleaned up 
c                 indenting around modified code.

c  Input/Output Variables
c                 
c   Input:    
c          
c    lun      I*4   logical unit number for station.info
c    sitcod   C*4   4-character station code for requested information
c    yr_in    I*4   year of epoch requested (4-digits)
c    doy_in   I*4   day-of-year of epoch requested
c    sod_in   I*4   seconds-of-day of epoch requested  
c    span     I*4   length of session (sec) if call for header value;
c                     =0 for single-epoch call
c
c   Output:
c
c    sname      C*16  16-character station name
c    anth       R*8   antenna height from file (type specified by antcod)
c    antn       R*8   north offset of antenna
c    ante       R*8   east offset of antenna
c    antdaz     R*8   Antenna deviation from true North (Alignment from True N in log)
c                     Added TAH 200203 for repro3
c    rcvcod     C*6   6-character receiver type
c    antcod     C*6   6-character antenna type
c    htcod      C*5   5-character type of antenna-height measurement
c    radome     C*5   5-character radome type
c    swver      R*4   receiver software version
c    istart(5)  I*4   start time (yr doy hr min sec) for values (yr 9999 if not found)
c    istop(5)   I*4   stop time  (yr doy hr min sec) for values 

c The number and type of data items in station.info file will vary according
c to applications and are specified by labels in a header line identified by 
c the keyword *SITE and including at least the labels for the other four 
c required entries:
c 
c   Session Start     Session Stop        RcvCod           AntCod
c                                (or    Receiver Type     Antenna Type)
c     
c Observations will usually include also an antenna height (anth) and
c possibly north and east offsets (antn ante) and receiver firmware (RINEX 
c 20-characters or GAMIT decimal code). 
c
c Following is a complete list of the allowed labels, the token (variable
c name) used in the code to represent them, and a description of the 
c data items they identify:
c
c Label         Token  Format   Default      Description
c ------        -----  -----    ----------  ----------------------------------------------
c *SITE         sitcod  a4                  4-char station code 
c Station Name  sname   a16                 16-char station name
c Session Start  start  2i4,3i3  0 0 0 0 0  start time (yr doy hr min sec)
c Session Stop  stop    2i4,3i3  0 0 0 0 0  stop time (yr doy hr min sec)
c Sessn         sessn   i1       0          session number
c Ant Ht        anth    f7.4     0.0        height of antenna ARP above the monument
c Ant N         antn    f7.4     0.0        N offset of ant ARP from the center of the monument
c Ant E         ante    f7.4     0.0        E offset of ant ARP from the center of the monument
c RcvCod        rcvcod  a6                  6-char/acter GAMIT receiver code
c Receiver Type rctype  a20                 IGS (RINEX) receiver name    
c Receiver SN   rcvrsn  a20                 receiver serial number    
c SwVer         swver   f5.2     0.0        GAMIT firmware code (real)
c Vers          rcvers  a20                 receiver firmware version from RINEX
c AntCod        antcod   a6                 6-charr GAMIT antenna code     
c Antenna SN    antsn   a20                 20-char antenna serial number
c Antenna Type  anttyp  a15                 IGS antenna name, 1st 15 chars of RINEX antenna field
c Dome          radome   a5                 IGS radome name, last 5 chars of RINEX antenna field
c HtCod         htcod    a5     DHARP       5-char GAMIT code for type of height measurement
c AntDAZ	antdaz  f5.0     0.0        Alignment from True N (deg).  TAH 2020203. 

c Any line with a non-blank character in the first c column is read as a comment.  
c Any item specified in the list of labels must be non-blank ('-----' for dummy value).
c For readability, items are separated by two blank spaces.  An example is shown below:
c                                                    
c *SITE  Session Start      Session Stop       Ant Ht  RcvCod  Antcod  HtCod
c  BLHL  1998 253 12 00 00  1998 254 12 00 00  1.4116  TRMSSE  TRMSST  SLTGP
c
c Because of the different logic associated with having an explicit stop time,
c the current version of rstnfo (and all of GAMIT) is not compatible with old-style
c station.info files.  If an old-style file is read, it will be detected by 
c subroutine read_stnfohead (below) and the run will be terminated.


      implicit none
         
c       station.info variables
                                
      character*4 sitcod
      character*5 htcod,radome  
      character*6 rcvcod,antcod   
      character*7 span_match
      character*15 anttyp
      character*16 sname
      character*20 rctype,rcvers,rcvrsn,antsn
* MOD TAH 200317: Increased size of comment to 36 from 30 to allow
*     logfile name and source to be saved.
      character*36 comment  
      real*4 swver
      real*8 anth,antn,ante 
      real*8 antdaz  ! Antenna aligment from True N (deg).
      integer*4 istart(5),istop(5),span
        
c       arrays for values and pointers
c MOD TAH 200203: Still slots avaiable so no change in dimenrion.
      integer*4 nlist
      character*6 item_list(20)
      character*20 values(20)

c       other variables

      logical first_call,found,found_comment,eof,warnings
                            
      integer*4 lun,yr_in,doy_in,sod_in
     .        , ih,im,is,len,rcpar,ioerr,i,j

      real*8 xsod,sec  

      character*4 site_test,upperc 
      character*80 prog_name
      character*256 line,requested_line,message

c MOD TAH 200524: Added to allow gaps when span is zero (ie., 
c     specific epoch not coverd by station.info entry.  Also
c     modified so that whole station.info file is not re-read
c     with each call.
      integer*4 num_in_gap  ! Number of seconds processed inside
                            ! gap (fatal if more than 6-hrs).  This
                            ! inteval covers firmware update gaps (I think).
      character*256 line_next    ! Next line in station.info when a data point
                            ! falls in a gap between entries.
      character*20 val_next(20)  ! Values for the items in the next station.info
                            ! line.
      integer*4 istr_next(5) ! Start time of next station.info line 
      integer*4 itime(5)     ! Time of current epoch being checked.

c     function
      integer*8 itimdif
                  
      data first_call /.true./    
      data nlist/0/
      data num_in_gap  / 0 /   ! Initialize 

      save first_call,prog_name,nlist,item_list, num_in_gap
      save line_next, istr_next 
 
c   Get program name for report_stat calls
      len = rcpar(0,prog_name)         

c   Set warnings = true for call to decode_values (set false by
c   mstinf2 for merging station.info files
       warnings =.true.
             
c   Initialize entries that might be missing
      rcvrsn = ' ' 
      antsn = ' ' 
      rcvers = ' ' 
                    
c   Rewind the file and read the headers

* MOD TAH 200526: Only re-wind the file and start again if
*     we have not already the file.  When scan=0, called
*     for single epoch and if num_in_gap is non-zero then
*     are waiting for the time to use the next entry.
      if( num_in_gap.eq. 0 ) then              
         rewind(lun,iostat=ioerr) 
         call read_stnfohead( lun,line )
      endif

c     if first call get the token list; 
c     otherwise, just position the file to read data
      if( first_call ) then
        do i=1,20
          item_list(i)=' '
        enddo  
        call get_tokens( line,nlist,item_list )
        first_call = .false.
      endif  

c     make sure the site code is uppercase
      call uppers(sitcod)
         
c___________________________________________________________________________
       
c  Get the requested information from the file

c     Kinematic: read the next line and return what's there
c     Read the file for the matching values for the input station

* MOD TAH 200525: If we are in a gap, see if new observation is after
*     start of next line
      if( num_in_gap.gt.0 ) then
*        See we have passed the start of next station.info line
         itime(1) = yr_in 
         itime(2) = doy_in   
         call ds2hms( itime(1),itime(2),dble(sod_in),itime(3),
     .                itime(4),sec) 
         itime(5) = int(sec)  
         if( itimdif(itime,istr_next).gt.0 ) then  ! We have found next line
*           Copy ove the line so that we don't need to read file.
            found = .true.
            requested_line = line_next
            write(message,'(a,a4,i5,i4,3i3,1x,I3,1x,a)') 
     .           'Station.info entry ',sitcod,itime,num_in_gap,
     .           'epochs in meta data gap'
            call report_stat('WARNING',prog_name,'/lib/rstnfo',' '
     .                     , message,0) 

            num_in_gap = 0  ! Reset in case there is another gap later in
         else  ! We are still are still in the gap so use current meta data
            RETURN
         endif 
      else      ! Normal original code   
c       initialize the 'found' variable : note that with the new start/stop logic allowing 
c       entries to start late or stop early in a session for a particular station, we can not
c       longer conveniently check for duplicates or overlapping entries
        found = .false.
        
c       begin loop on all records of station.info
                             
        eof = .false.  

        do while( .not.eof .and. .not.found ) 
                                            
c        read until we find a non-comment line
         found_comment = .true.
         do while ( found_comment )            
           read( lun,'(a)',iostat=ioerr,end=90 ) line 
           if( ioerr.eq.-1 ) then
              eof = .true.        
           elseif( ioerr.ne.0 ) then 
              call report_stat('FATAL',prog_name,'lib/rstnfo',' '
     .                 ,'Error reading station.info ',ioerr)   
           endif
           if( line(1:1).eq.' ')  found_comment=.false. 
         enddo
         
c        get the values and check for a station match  
         call read_values( line,nlist,item_list,values, comment )
         site_test = upperc(line(2:5))
         if( site_test.eq.sitcod ) then    
               
c          check for requested time within the bounds of the entry time 
           do i=1,nlist 
             if( item_list(i).eq.'start' ) then
              read(values(i),'(2i4,3i3)',iostat=ioerr) (istart(j),j=1,5) 
              if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .             ,'lib/rstnfo',' '
     .             ,'Error reading start time from station.info ',ioerr)
             elseif( item_list(i).eq.'stop' ) then   
               read(values(i),'(2i4,3i3)',iostat=ioerr) (istop(j),j=1,5)
               if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .             ,'lib/rstnfo',' '
     .             ,'Error reading stop time from station.info ',ioerr)
c              change year = 9999 to something reasonable for time conversions
               if( istop(1).eq.9999 ) then
                  istop(1) = 2100 
                  istop(2) = 1
                  do j=3,5
                    istop(j) = 0
                  enddo  
               endif 
             endif         
           enddo   
c          check for bogus start/stop times
           if( itimdif(istop,istart).lt.0 ) then  
             write(message,'(a,i5,i4,3i3,a,i5,i4,3i3,a,a4,i4,a)')
     .          'Stop ',istop,' earlier than start ',istart,' for '
     .          ,sitcod,doy_in,' in station.info'
             call report_stat('FATAL',prog_name,'lib/rstnfo',' '
     .                      ,message,0)            
           endif                    

           call check_span( yr_in, doy_in, sod_in, span, istart, istop
     .                    , span_match ) 
     
           if( span_match.eq.'span_in' ) then
             found = .true.
             requested_line = line  
c          if this is a header value, check also for partial overlap: keep but issue warning
           elseif ( span.gt.0 .and. span_match.eq.'early__' ) then
             found = .true.
             requested_line = line    
             write(message,'(a,a4,i5,i4,3i3,1x,i5,i4,3i3,a)') 
     .             'Station.info entry ',sitcod,istart,istop
     .             ,' starts late for session but may be ok for station'
             call report_stat('WARNING',prog_name,'/lib/rstnfo',' '
     .                     , message,0)    
           elseif ( span.gt.0 .and. span_match.eq.'late___' ) then
             found = .true.
             requested_line = line    
             write(message,'(a,a4,i5,i4,3i3,1x,i5,i4,3i3,a)') 
     .             'Station.info entry ',sitcod,istart,istop
     .             ,' ends early for session but may be ok for station'
             call report_stat('WARNING',prog_name,'/lib/rstnfo',' '
     .                     , message,0) 
* MOD TAH 200524: Allow some gap for early end even when span equals zero. 
           elseif ( span.eq.0 .and. span_match.eq.'near_en' ) then
             found = .true.   ! Allow line to be used.
             requested_line = line    
             if( num_in_gap.eq.0 ) then ! print warning
                write(message,'(a,a4,i5,i4,3i3,1x,i5,i4,3i3,1x,a)') 
     .             'Station.info entry ',sitcod,istart,istop,
     .             'ends early: gap allowed'
                call report_stat('WARNING',prog_name,'/lib/rstnfo',' '
     .                     , message,0)
                num_in_gap = 1
* MOD TAH 200524: Read the next line of station.info.  We will use 
*               current meta data until the new time is passed the
*               start of the next entry
                read( lun,'(a)',iostat=ioerr,end=90 ) line_next
*               Now get the start time form list (assume no comment)
                call read_values( line_next,nlist,item_list,val_next,
     .                            comment )
                do i=1,nlist 
                    if( item_list(i).eq.'start' ) then
                       read(val_next(i),'(2i4,3i3)',iostat=ioerr) 
     .                    (istr_next(j),j=1,5) 
                       if( ioerr.ne.0 ) call report_stat('FATAL',
     .                      prog_name,'lib/rstnfo',' ',
     .                     'Error reading start time line_next ',
     .                      ioerr)
                    endif
                enddo
*               Now we wait for data epoch to pass this the istr_next time.
 
             elseif( num_in_gap.gt. 240 ) then  ! 2hrs hours@30sec: Too long
                write(message,'(a,a4,i5,i4,3i3,1x,i5,i4,3i3,a)') 
     .             'Station.info entry ',sitcod,istart,istop
     .             ,' Too long in gaps'
                call report_stat('FATAL',prog_name,'/lib/rstnfo',' '
     .                     , message,num_in_gap)
             else
                num_in_gap = num_in_gap + 1
             endif  
           endif

         endif 

c     --end of loop over data records
        enddo

      endif      ! num_in_gap zero

   90 if( found ) then   
                    
c       read the values from the line into an array  
        call read_values( requested_line,nlist,item_list,values,comment) 
        call decode_values( nlist,item_list,values
     .                    , doy_in, sitcod, sname
     .                    , anth, antn, ante, antdaz, rcvcod, antcod  
     .                    , rctype, rcvrsn, rcvers, anttyp, antsn
     .                    , htcod, radome, swver, istart, istop
     .                    , warnings ) 

c** rwk 050812 remove this from here; use subsitution only in reading antmod.dat (antex) file       
c              blank or ----- are set UNKN in decode_values
c        if( radome.eq.'UNKN ') then   
c          write(message,'(a,a4,1x,2i4,3i3,a)')
c     .   'Radome ',sitcod,yr_in,doy_in,ih,im,is,' unknwn, assume NONE'
c          call report_stat('WARNING',prog_name,'lib/rstnfo',' '
c     .                    ,message,0)   
c        endif 

      else
             
        xsod = dble(sod_in)               
        call ds2hms(yr_in,doy_in,xsod,ih,im,sec)
        is = int(sec) 
        write(message,'(a,a4,1x,2i4,3i3,a)')   
     .   'No match for ',sitcod,yr_in,doy_in,ih,im,is,' in station.info'    
c       if called by makexp, may have linked files with data outside the span, so
c       issue only a warning.    
        if( prog_name(1:6).eq.'makexp' ) then 
          call report_stat('WARNING',prog_name,'lib/rstnfo',' '
     .                    ,message,0)
        else
          call report_stat('FATAL',prog_name,'lib/rstnfo',' '
     .                    ,message,0)
        endif
      endif
       
      return 
      end

c***********************************************************************

      Subroutine read_stnfohead ( lun,line )

c     read station.info through the header records
                                     
c     Input:  lun   -  logical unit number
c     Output: line  -  line containing data tokens 

      implicit none
                    
      logical first_call,label_line

      integer*4 lun,ioerr,len,rcpar,i
                 
      character*5 sitlab  
      character*80 prog_name
      character*256 line,message
      
      data first_call/.true./

      save first_call
  
c     read first three lines to see if oldstyle format
 
      if( first_call ) then 
c       get program name for report_stat calls
        len = rcpar(0,prog_name)   
c         oldstyle will have 'TRCK' in line 3
c         newstyle will have '*SITE' in header, which may have only 1 line
c            or an indetermine number of comment lines
        do i=1,3
          read( lun,'(a)',iostat=ioerr) line
          if( ioerr.eq.-1 ) then 
            continue
          elseif( ioerr.ne.0 ) then
             call report_stat('FATAL',prog_name,'lib/rstnfo',' '
     .        ,'Error reading header lines of station.info   ',ioerr)
          endif
        enddo
        if( line(1:4).eq.'TRCK' .or. line(2:5).eq.'TRCK' ) then  
c         old-style -- stop the program  
          write(message,'(2a)') 'Old-style station.info incompatible'
     .         ,' with GAMIT 10.34, convert with conv_stnfo'
         call report_stat('FATAL',prog_name,'lib/rstnfo',' ',message,0)
        endif 
        first_call = .false.
      endif
                 
c     'TRCK' not found, must be new format; look for '*SITE' line
        
      rewind ( lun )    
      label_line = .false.
      do while (.not.label_line)
        line = ' '      
        sitlab = '*SITE'
        read( lun,'(a)',iostat=ioerr ) line 
        if( ioerr.ne.0 )  
     .     call report_stat('FATAL',prog_name,'lib/rstnfo',' '
     .         ,'Error reading token lines of station.info',ioerr)
        if( line(1:5).eq.sitlab ) then 
          label_line = .true.
        endif
      enddo

      return
      end       

c********************************************************************************

      Subroutine get_tokens( line,nlist,item_list )

c     Decode header labels and assign items and columns
c     R. King 9 Jan 2001

      implicit none
          
c       arrays for values and pointers

      integer*4 nlist
      character*6 item_list(20)
         
c        other variables
 
      logical found_rcv,found_ant
              
      integer*4 labcol(20),len,rcpar,mchcol,i
                        
      character*256 line,prog_name  
                      
      
c   Get program name for report_stat calls
      len = rcpar(0,prog_name)
              
c   Read the labels from the header line and set the token list 
c     (save the column number to sort the list at the end)

c     get first the five required entries
      nlist = 1 
c     sitcod must be first and start in column 2  or we wouldn't have gotten this far
      labcol(nlist) = 2     
      item_list(nlist) = 'sitcod'
      i = mchcol(line,'Start',256,5)
      if( i.gt.0 ) then
        nlist = nlist + 1
        labcol(nlist) = i 
        item_list(nlist) = 'start ' 
      else
       call report_stat('FATAL',prog_name,'lib/rstnfo',' '
     .          ,'Session Start missing from station.info',0)
      endif
      i = mchcol(line,'Stop',256,4)
      if( i.gt.0 ) then
        nlist = nlist + 1
        labcol(nlist) = i  
        item_list(nlist) = 'stop  '
      else
       call report_stat('FATAL',prog_name,'lib/rstnfo',' '
     .          ,'Session Stop missing from station.info',0)
      endif                                              
      i = mchcol(line,'Ant Ht',256,6)
      if( i.gt.0 ) then
        nlist = nlist + 1
        labcol(nlist) = i  
        item_list(nlist) =  'anth  '  
      else
      call report_stat('FATAL',prog_name,'lib/rstnfo',' '
     .          ,'Antenna height missing from station.info',0)
      endif                       
      found_rcv = .false.
      i = mchcol(line,'RcvCod',256,6)
      if( i.gt.0 ) then
        nlist = nlist + 1
        labcol(nlist) = i 
        item_list(nlist) = 'rcvcod'
        found_rcv = .true.
      endif
      i = mchcol(line,'Receiver Type',256,13)
      if( i.gt.0 ) then
        nlist = nlist + 1
        labcol(nlist) = i  
        item_list(nlist) =  'rctype'
        found_rcv = .true.
      endif
      if( .not.found_rcv )
     .    call report_stat('FATAL',prog_name,'lib/rstnfo',' '
     .    ,'Receiver name or code missing from station.info',0)  
      found_ant = .false.
      i = mchcol(line,'AntCod',256,6)
      if( i.gt.0 ) then
        nlist = nlist + 1
        labcol(nlist) = i 
        item_list(nlist) = 'antcod'
        found_rcv = .true.
      endif
      i = mchcol(line,'Antenna Type',256,12)
      if( i.gt.0 ) then
        nlist = nlist + 1
        labcol(nlist) = i  
        item_list(nlist) =  'anttyp'
        found_rcv = .true.
      endif
      if( .not.found_rcv )
     .    call report_stat('FATAL',prog_name,'lib/rstnfo',' '
     .    ,'Antenna name or code missing from station.info',0)
c     now the optional entries
      i = mchcol(line,'Station',256,7)
      if( i.gt.0 ) then
        nlist = nlist + 1
        labcol(nlist) = i  
        item_list(nlist) =  'sname '
      endif                       
      i = mchcol(line,'Sessn',256,5)
      if( i.gt.0 ) then
        nlist = nlist + 1
        labcol(nlist) = i  
        item_list(nlist) =  'sessn '
      endif                       
      i = mchcol(line,'Ant N',256,5)
      if( i.gt.0 ) then
        nlist = nlist + 1
        labcol(nlist) = i  
        item_list(nlist) =  'antn  '
      endif                       
      i = mchcol(line,'Ant E',256,5)
      if( i.gt.0 ) then
        nlist = nlist + 1
        labcol(nlist) = i  
        item_list(nlist) =  'ante  '
      endif                       
      i = mchcol(line,'Receiver SN',256,11)
      if( i.gt.0 ) then
        nlist = nlist + 1
        labcol(nlist) = i  
        item_list(nlist) =  'rcvrsn'
      endif                       
      i = mchcol(line,'Antenna SN',256,10)
      if( i.gt.0 ) then
        nlist = nlist + 1
        labcol(nlist) = i  
        item_list(nlist) =  'antsn '
      endif   
      i = mchcol(line,'SwVer',256,5)
      if( i.gt.0 ) then
        nlist = nlist + 1
        labcol(nlist) = i  
        item_list(nlist) =  'swver '
      endif                       
      i = mchcol(line,'Vers',256,4)
      if( i.gt.0 ) then
        nlist = nlist + 1
        labcol(nlist) = i  
        item_list(nlist) =  'rcvers'
      endif                       
      i = mchcol(line,'Dome',256,4)
      if( i.gt.0 ) then
        nlist = nlist + 1
        labcol(nlist) = i  
        item_list(nlist) =  'radome'
      endif                          
      i = mchcol(line,'HtCod',256,5)
      if( i.gt.0 ) then
        nlist = nlist + 1
        labcol(nlist) = i  
        item_list(nlist) =  'htcod '
      endif

* MOD TAH 200203: Added AntDAZ (antdaz) 
      i = mchcol(line,'AntDAZ',256,6)
      if( i.gt.0 ) then
        nlist = nlist + 1
        labcol(nlist) = i  
        item_list(nlist) =  'antdaz'
      endif

     
c     sort the list by column number to get the order right
      call sort_intch (nlist,labcol,item_list ) 

      return
      end          
   
c****************************************************************************************
 
       
      Subroutine read_values( line,nlist,item_list,values,comment )

c     read the items of a station.info line into a character array
c     under format control
c     R. King 28 April 2001

      implicit none

c       arrays for values and pointers

      integer*4 nlist,ncol
      character*6 item_list(20)
      character*20 values(20)
                   
c       other variables
                            
      logical first_call

      character*7 fmti(20) 
* MOD TAH 200317: Make strings use passed length.
      character*(*) comment
      character*80 prog_name
      character*128 format 
      character*256 line
              
      integer*4 rcpar,len,ioerr,i  
      integer*4 lencom  ! Length of string used to hold comments at
                        ! end of stinfo lines (TAH 200317)  
      integer*4 lenpn   ! Length of program name (used to be len) 
                    
      data first_call/.true./
      save first_call,format
                 
c     get the program name for report_stat call
* MOD TAH 200317: LEN is instrinsic routine can should not be used
*     as variable name.  Replaced len with lenpn
          
      lenpn = rcpar(0,prog_name)  

c     MSTINF can call rstnfo with differently formatted station.info files, so
c     must regenerate the format 
      if( prog_name(1:6).eq.'mstinf' ) first_call = .true.
         
c     construct the format for reading the lines      
      if( first_call ) then      
        first_call = .false.
	 
        if( item_list(1).eq.'sitcod' ) then
          fmti(1) = '(1x,a4' 
c         ncolumn counts the values columns (not the format) and set to the last
c         column of the field (i.e, the next field begins at ncol + 3); it is used 
c         only to determine the column for reading in-line comments 
          ncol = 5
        else
          call report_stat('FATAL',prog_name,'/lib/rstnfo',' '
     .               ,'First token in label line not *SITE',0)
        endif
        do i= 2,nlist
          if( item_list(i).eq.'sname'  ) then
            fmti(i)=  ',2x,a16'    
            ncol = ncol + 18           
          elseif( item_list(i).eq.'sessn' ) then
            fmti(i) = ',2x,a1' 
            ncol = ncol + 3
          elseif( item_list(i).eq.'start'  ) then
            fmti(i)= ',2x,a17' 
            ncol = ncol + 19
          elseif( item_list(i).eq.'stop'  ) then
            fmti(i)= ',2x,a17' 
            ncol = ncol + 19 
          elseif( item_list(i).eq.'anth' ) then 
            fmti(i) = ',2x,a7' 
            ncol = ncol + 9    
          elseif( item_list(i).eq.'htcod') then
            fmti(i) = ',2x,a5'
            ncol = ncol + 7
          elseif( item_list(i).eq.'antn' ) then
            fmti(i) = ',2x,a7'  
            ncol = ncol + 9
          elseif( item_list(i).eq.'ante' ) then
            fmti(i) = ',2x,a7'
            ncol = ncol + 9
          elseif( item_list(i).eq.'rcvcod' ) then
            fmti(i) = ',2x,a6' 
            ncol = ncol + 8
          elseif( item_list(i).eq.'rctype' ) then
            fmti(i) = ',2x,a20'
            ncol = ncol + 22   
          elseif( item_list(i).eq.'rcvrsn' ) then
            fmti(i) = ',2x,a20'       
            ncol = ncol + 22              
          elseif( item_list(i).eq.'swver' ) then
            fmti(i) = ',2x,a5' 
            ncol = ncol + 7     
          elseif( item_list(i).eq.'rcvers' ) then
            fmti(i) = ',2x,a20'
            ncol = ncol + 22      
          elseif( item_list(i).eq.'antcod' ) then
            fmti(i) = ',2x,a6' 
            ncol = ncol + 8
          elseif( item_list(i).eq.'anttyp' ) then
            fmti(i) = ',2x,a15' 
            ncol = ncol + 17    
          elseif( item_list(i).eq.'antsn' ) then
            fmti(i) = ',2x,a20'
            ncol = ncol + 22   
          elseif( item_list(i).eq.'radome' ) then
            fmti(i) = ',2x,a5' 
            ncol = ncol + 7   
          elseif( item_list(i).eq.'antdaz' ) then
            fmti(i) = ',2x,a5' 
            ncol = ncol + 8   
          else
             call report_stat('FATAL',prog_name,'/lib/rstnfo',
     .            item_list(i),
     .           'Station.info data token not in standard list',0) 
          endif   
        enddo  
        write(format,'(21a7)',iostat=ioerr) 
     .     (fmti(i),i=1,nlist),')'   
      endif  

c     read the values
      read( line,format,iostat=ioerr ) (values(i),i=1,nlist)
      if( ioerr.ne.0 )  call report_stat('FATAL',prog_name,'lib/rstnfo'
     .   ,' ','Error reading data values from station.info',ioerr)   
                     

c     read any comments appended to the line (30 characters max)
* MOD TAH 200317: Use the length of the comment variable passed.
*     Use lencom instead of 30
      lencom = len(comment)

      comment = ' ' 
* MOD TAH 200317: Replaced 256 with actual line length.  Removed test
*     on length of line and had read to end of line.
C     if( (ncol+lencom).le.len(line) ) then
      read( line(ncol+1:),'(a)',iostat=ioerr ) comment
      call trimlead(comment)   ! Remove any leading blanks from comments 
C      endif   

      return
      end

c***********************************************************************      

      Subroutine decode_values( nlist,item_list,values
     .                        , doy_in, sitcod, sname
     .                        , anth, antn, ante, antdaz, rcvcod, antcod
     .                        , rctype, rcvrsn, rcvers, anttyp, antsn
     .                        , htcod, radome, swver, istart, istop
     .                        , warnings )
                                
c     Transfer the items in the values array into the expected variable names
c     R. King 28 April 2001   

c     If warnings = T, call report_stat to report missing entries

c MOD TAH 200203: Added AntDAZ token for antenna  Alignment from True N
      
      implicit none

c       station.info variables
                                
      character*4 sitcod
      character*5 htcod,radome,aswver
      character*6 rcvcod,antcod
      character*20 rcvrsn
      character*15 anttyp
      character*16 sname
      character*20 rctype,rcvers,antsn,anttyp20
      real*4 swver
      real*8 anth,antn,ante 
      real*8 antdaz  ! Antenna aligment from True N (deg).
      integer*4 isessn,istart(5),istop(5)
        
c       arrays for values and pointers

      integer*4 nlist
      character*6 item_list(20)
      character*20 values(20)
       
c       other variables
         
      logical warnings                                                
      character*1 pcncod
      character*5 char5
      character*6 char6
      character*20 char20       
      character*80 prog_name
      character*256 message    
      integer*4 doy_in,rcpar,len,ioerr,i,j  

c      get the program name for report_stat call
          
      len = rcpar(0,prog_name)

c       initialize values 

      sname = ' '  
      isessn = 1
      anth = 0.0d0
      antn = 0.0d0
      ante = 0.0d0 
      antdaz = 0.d0 
      rcvcod = ' ' 
      rctype = ' ' 
      rcvrsn = ' ' 
      swver = 0.0
      rcvers =  ' '
      antcod = ' ' 
      anttyp = ' '   
      anttyp20 = ' ' 
      htcod = ' ' 
      antsn = ' '  
      radome = 'UNKN '
      htcod = 'DHARP'
      do i=1,5 
        istart(i) = 0
        istop(i) = 0
      enddo

c       decode the values
                     
      do i = 1, nlist

        if( item_list(i).eq.'sitcod' ) then
           read(values(i),'(a4)',iostat=ioerr) sitcod  
           if( ioerr.ne.0 ) goto 10
        elseif( item_list(i).eq.'sname ' ) then
           read(values(i),'(a16)',iostat=ioerr) sname 
           if( ioerr.ne.0 ) goto 10
        elseif( item_list(i).eq.'sessn ' ) then
           read(values(i),'(i1)',iostat=ioerr ) isessn
           if( ioerr.ne.0 ) goto 10
        elseif( item_list(i).eq.'anth  ' ) then
           read(values(i),'(f7.4)',iostat=ioerr) anth 
           if( ioerr.ne.0 ) goto 10
        elseif( item_list(i).eq.'antn  ' ) then
           read(values(i),'(f7.4)',iostat=ioerr) antn 
           if( ioerr.ne.0 ) goto 10
        elseif( item_list(i).eq.'ante  ' ) then
           read(values(i),'(f7.4)',iostat=ioerr) ante 
           if( ioerr.ne.0 ) goto 10 
        elseif( item_list(i).eq.'antdaz' ) then
           read(values(i),'(f5.0)',iostat=ioerr) antdaz 
           if( ioerr.ne.0 ) goto 10 
        elseif( item_list(i).eq.'rcvcod' ) then
           read(values(i),'(a6)',iostat=ioerr) rcvcod   
           if( ioerr.ne.0 ) goto 10
        elseif( item_list(i).eq.'rctype' ) then
           read(values(i),'(a20)',iostat=ioerr) rctype 
           if( ioerr.ne.0 ) goto 10
        elseif( item_list(i).eq.'rcvrsn' ) then
           read(values(i),'(a20)',iostat=ioerr) rcvrsn 
           if( ioerr.ne.0 ) goto 10
        elseif( item_list(i).eq.'swver ' ) then   
           read(values(i),'(a5)',iostat=ioerr) aswver  
cd           write(*,'(a,a5,i5)'),'aswver ',aswver,ioerr 
           if( ioerr.ne.0 ) goto 10
           if( aswver(3:4).eq.'--' ) then 
             swver = 0.0                  
           else
             read(aswver,'(f5.2)',iostat=ioerr) swver 
           endif
        elseif( item_list(i).eq.'rcvers' ) then
           read(values(i),'(a20)',iostat=ioerr) rcvers
           if( ioerr.ne.0 ) goto 10
        elseif( item_list(i).eq.'antcod' ) then
           read(values(i),'(a6)',iostat=ioerr) antcod  
           if( ioerr.ne.0 ) goto 10
        elseif( item_list(i).eq.'anttyp' ) then
           read(values(i),'(a15)',iostat=ioerr) anttyp
           if( ioerr.ne.0 ) goto 10
        elseif( item_list(i).eq.'antsn' ) then
           read(values(i),'(a20)',iostat=ioerr) antsn 
           if( ioerr.ne.0 ) goto 10 
        elseif( item_list(i).eq.'htcod ' ) then
           read(values(i),'(a5)',iostat=ioerr) htcod 
           if( htcod.eq.'DHPAB' ) htcod = 'DHARP'
           if( ioerr.ne.0 ) goto 10
        elseif( item_list(i).eq.'radome') then
           read(values(i),'(a5)',iostat=ioerr) radome    
           if( ioerr.ne.0 ) goto 10       
           if( radome.eq.'     '.or.radome.eq.'-----') then
             radome = 'UNKN '      
           else      
c rwk 101111: Make sure 'NONE' is uppercase (not always true in SOPAC file)
             call uppers(radome)   
           endif 
        elseif( item_list(i).eq.'start ' ) then
           read(values(i),'(2i4,3i3)',iostat=ioerr) (istart(j),j=1,5) 
           if( ioerr.ne.0 ) goto 10
c          check values for reasonableness to trap common typos
           if( istart(1).lt.1980.or.istart(1).gt.9999 ) then
             write(message,'(a,i6,a,a4)') 'Unreasonable start year ('
     .              ,istart(1),') for site ',sitcod
             call report_stat('FATAL',prog_name,'lib/rstnfo',' '
     .                     ,message,ioerr)
           endif
        elseif( item_list(i).eq.'stop  ' ) then
           read(values(i),'(2i4,3i3)',iostat=ioerr) (istop(j),j=1,5)  
           if( ioerr.ne.0 ) goto 10   
           if( istop(1).lt.1980.or.istop(1).gt.9999 ) then
             write(message,'(a,i6,a,a4)') 'Unreasonable stop year ('
     .              ,istop(1),') for site ',sitcod
             call report_stat('FATAL',prog_name,'lib/rstnfo',' '
     .                    ,message,ioerr)
           endif
        endif   
10      if( ioerr.ne.0 ) then               
           write(message,'(a,a6,a,a4,a)') 'Error decoding '
     .        ,item_list(i),' for ',sitcod,' from station.info'
           call report_stat('FATAL',prog_name,'lib/rstnfo',' '
     .                     ,message,ioerr)
        endif
      enddo
        
c     assign corresponding RINEX/GAMIT receiver, antenna, and firmware tokens   
      if((rcvcod.eq.'      '.or.rcvcod.eq.'------') .and.
     . (rctype(1:5).eq.'     '.or.rctype(1:5).eq.'-----') ) then
        if( warnings ) then
          write(message,'(a,a4,i5,i4,3i3,a)') 
     .      'Neither RcvCod nor Receiver Type found for ',sitcod
     .      ,(istart(i),i=1,5),' in station.info'
          call report_stat('WARNING',prog_name,'lib/rstnfo',' '
     .                    ,message,0)  
        endif
      elseif (rctype(1:5).eq.'     '.or.rctype(1:5).eq.'-----') then
        call read_rcvant(1,2,char6,char20,char5,rcvcod,rctype,pcncod)
      elseif (rcvcod.eq.'      '.or.rcvcod.eq.'------') then   
        call read_rcvant(2,2,char6,char20,char5,rcvcod,rctype,pcncod) 
      else
        call report_stat('FATAL',prog_name,'lib/rstnfo',' '
     .    ,'Receiver entry inconsistency in station.info',0)
      endif  
                  
      if((antcod.eq.'      ' .or.antcod.eq.'------') .and.
     .    (anttyp(1:5).eq.'     '.or.
     .    anttyp(1:5).eq.'-----')) then
        if( warnings ) then 
          write(message,'(a,a4,i5,i4,3i3,a)') 
     .      'Neither AntCod nor Antenna Type found for ',sitcod
     .      ,(istart(i),i=1,5),' in station.info'
          call report_stat('WARNING',prog_name,'lib/rstnfo',' '
     .      ,message,0)    
        endif
      elseif (anttyp(1:5).eq.'     '.or.anttyp.eq.'-----') then   
        call read_rcvant(1,1,antcod,anttyp20,radome,char6,char20
     .    ,pcncod)
        anttyp = anttyp20(1:15)
      elseif (antcod.eq.'      '.or.antcod.eq.'------') then 
        anttyp20(1:15) = anttyp    
        call read_rcvant(2,1,antcod,anttyp20,radome,char6,char20
     .                  ,pcncod) 
      else
         call report_stat('FATAL',prog_name,'lib/rstnfo',' '
     .      ,'Antenna entry inconsistency in station.info',0)
      endif
   
c     Check for blunder to avoid negative sqrt in hisub 
      if( anth.eq.0.d0 .and. htcod(1:2).eq.'SL' ) then
         write(message,'(a,a4,i4,a)') 'Slant ht = 0. for '
     .        ,sitcod,doy_in,' in station.info'
         call report_stat('FATAL',prog_name,'lib/rstnfo',' '
     .                   ,message,0)
      endif 

      return
      end  

*******************************************************************************

      subroutine sort_intch (ndat,iarray,carray )

c     Sort a character array using the values of an associated
c     integer array, smallest to largest.  This version for up to 
c     20 6-character station.info header tokens.    R, King 021212

      implicit none
                       
      integer*4 ndat,iarray(20),iswtch,i,k  
      character*6 carray(20),cswtch

      if(ndat.le.1) go to 100
      do 50 k=1,32000
         iswtch=0
         do 30  i=2,ndat
            if(iarray(i).ge.iarray(i-1)) go to  30
            iswtch=iarray(i)   
            cswtch=carray(i)
            iarray(i)=iarray(i-1) 
            carray(i)=carray(i-1)
            iarray(i-1)=iswtch 
            carray(i-1)=cswtch
 30      continue
         if(iswtch.eq.0) go to 100
 50   continue
c
 100  continue
c
      return
      end
c*********************************************************************************

      integer function mchcol(strinx,striny,ix,iy)

c     Find a string within another string and return the starting column number   
c     R. King 12 Dec 2002, based on D. Dong routine mchkey.f (in gamit/lib)
c
c     Input:
c        stringx:  string to be searched
c        stringy:  candidate string to be found
c        ix     :  length of stringx    
c        iy     :  length of stringy
c           
c     Output:  Value of mchcol (-1 if not found, +col if found)

      implicit none

      character*(*) strinx,striny
      integer ix,iy,i,i1
c
      mchcol = -1
      i1 = 0
      do 10 i = iy,ix
         i1 = i1+1
         if (striny(1:iy).eq.strinx(i1:i)) goto 20
 10   continue
      return
c
 20   mchcol = i1
      return
      end

c*******************************************************************************

      Subroutine check_span( year, doy, sod, span, istart, istop
     .                     , span_match )
                       
      implicit none
               
      character*7 span_match  

      integer*4 year, doy, sod, span, istart(5), istop(5)
     .        , itime(5)
      integer*8 itimdif  

      real*8 sec

c     Routine to determine if a requested time is within the range of a
c     station.info entry.  The 'span' is the length of the range to be
c     tested, in seconds.  For header values, span is the session length
c     and the check is whether the station.info entry falls anywhere within
c     the session.  If the data span begins too early or extends beyond
c     the end of the station.info span, set a return flag to issue a warning
c     (though there may be no problem if the actual data span for a station
c     is shorter than the x-file session length.  For epoch values, span is 
c     zero, and the epoch must be within the station.info entry start and stop times.  
c     
c     There are four possibilities for span_match:

c      not_in_  No overlap between the epoch or span and the station.info entry
c      span_in  The epoch or span falls entirely with the station.info entry
c      early__  The epoch is earlier than the station.info entry but the
c                  the session span overlaps with the station.info entry
c      late___  The (start) epoch is within the the station.info entry but the 
c                  session span extends beyond the station.info stop time. 
c MOD TAH 200524: Added case to handle gaps in station.info entries
c      past_en  The span is 0 and the data time is <6hrs past the end of span.
c               This case is kept track of to see when next record starts  

c      The tolerance for these checks is 5 seconds (originally 60s, so be careful
c      that we haven't introduced a bug.  rwk 100331    

c
      itime(1) = year
      itime(2) = doy                  
      call ds2hms( itime(1),itime(2),dble(sod),itime(3),itime(4),sec) 
      itime(5) = int(sec)  

      span_match = 'not_in_'

      if( span.gt.0 ) then  

        if( itimdif(itime,istart).gt.(-5) .and.
     .      itimdif(itime,istop).lt.(5-span) ) then  
          span_match = 'span_in' 
        elseif( itimdif(itime,istart).le.(-5).and. 
     .          itimdif(itime,istop).lt.0 .and.
     .          itimdif(itime,istart).gt.(-span) ) then
          span_match = 'early__'     
        elseif( itimdif(itime,istart).gt.(-5) .and.  
     .          itimdif(itime,istop).lt.0 .and.
     .          itimdif(itime,istop).gt.(-span) ) then
          span_match = 'late___'  
        else    
c         no overlap  
        endif

      elseif( span.eq.0 ) then
        if( itimdif(itime,istart).gt.(-5) .and.
     .      itimdif(itime,istop).lt.5 ) then
          span_match = 'span_in'          
* MOD TAH 200524: See if end time is just past the end of the
*     span (i.e., allow for gap due to say firmware switch).
*     Allo up to 6-hour gap (7200 seconds)
        elseif( itimdif(itime,istart).gt.(-5) .and.
     .      itimdif(itime,istop).lt.7200 ) then
          span_match = 'near_en' 
        endif      
      endif   

      return
      end





