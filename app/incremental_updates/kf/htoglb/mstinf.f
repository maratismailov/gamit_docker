      Program MSTINF

      implicit none

*     This program has two purposes:   
*     (a) Merge station.info files into one file with a reference file
*         that takes precedence
*     (b) Add new information into station.info based on rinex file headers
*         (after checking to see if information is already present).
*
*     Constructed by rwk March 2002 for new-format station.info, from original
*     program MSTINF written by tah; ability to add information from the RINEX
*     header, dummy in original MSTINF, added.
           
*     Calling arguments:
*           
*     Required:
*     -f [file]              Reference station.info file 
*     Optional:  
*     -w [file]              Output station.info file  (if omitted, simply compare)
*     -s [file1, file2, ..]  One or more station.info files to be merged   
*     -r [file]              RINEX file to be read for station.info information (optional) 
*     -i [file]              IGS log file to be read for station.info information (optional) 
*     -x [file]              SINEX file to be read for station.info information (optional)  
*     -st [file]             GIPSY sta_rcvr and sta_svec files (must have both: rcv, ant) to be read for station.info information (optional)
*     -xs [site]             Site to be extracted from SINEX file (if omitted, extract all)
*     -o                     Allow overwriting of input station.info filename    
*     -rep [value]           Replace option (default diff)
*                               'all' add a new entry even if values are the same (recommended with IGS log files)
*                               'none' check values but do not add new entries
*                               'diff' add a new entry only if the values are different 
*     -rxo                   Treat all RINEX entries as open-ended (may produce erroneous results if later entries)
*     -no                    Do not write out duplicate or overridden entries (default is to write as comment '-'
*     -ns                    Do not sort the added sites (use only when confident of time order and no duplicates)
*     -c                     Copy all comments from input station.info files    
*     -t                     Tolerance in seconds for deciding whether two start times are the same (default 120)
*     -u [site1, site2, ..]  List of stations to be used from reference station.info
*     -l [file]              File containing list of stations to be used from reference station.info
*     -h [value]             Height (m) above which RINEX values will be treated as slant heights  
*     -apr                   Write out an apr file from the coordinates in the IGS or SINEX file (name: 'mstinf.apr')            
*     -nogaps                Don't allow gaps in the entries from IGS logs; used to prevent spans with no
*                            metadata where data recording has occurred (error in logs). 
*     -debug                 Writes IGS log file as read to find bad date lines.

      include 'mstinf.h' 
      include '../../gamit/includes/makex.h'


* LOCAL VARIABLES     
* ierr -- IOSTAT error
* i    -- Loop counter

      integer*4 uigs,usnx,ustar,ustaa,ierr, i   
                  
      character*80 line
      
         
* Function
      integer trimlen

*     Write the program name
c*      write(*,'(a)') 'MSTINF: Merge or update station.info files' 
      call report_stat('CLEAR','MSTINF',' ',' ' ,' ',0)
      call report_stat('STATUS','MSTINF','htoglb/mstinf',' '
     .  ,'Merge or update station.info files: Version 2.3 2015-11-08',0)

*     Decode the runstring
      call decode_stinf_run     
               
*     Open the output apr file if requested
      if( write_apr ) then
        uapr = 50
        open( uapr, file='mstinf.apr', status='unknown', iostat=ierr)
        if( ierr.ne.0 )  call report_stat('FATAL','MSTINF'
     .   ,'htoglb/mstinf',' '
     .    ,'Error opening output apr file', ierr)  
       if( num_rxf.gt.0 ) write(uapr,'(a,i3,a)') 
     .    '* Coordinates from ',num_rxf,' RINEX files'
       if( num_igf.gt.0 ) write(uapr,'(a,i3,a)')
     .      '* Coordinates from ',num_igf,' IGS log files'
       if( num_sxf.eq.1 ) write( uapr,'(2a)') 
     .      '* Coordinates from ',sx_file(1)
       call report_stat('STATUS','MSTINF','htoglb/mstinf',' '
     .      ,'Apriori coordinates will be written to mstinf.apr',0)
      endif
         
*     Initialize the array indicating old or new entries
      do i=1,max_str
        newflag(i) = .false.
      enddo

*     Now read the reference station.info file
      open(100, file=ref_stinf, status='old', iostat=ierr)
      if( ierr.ne.0 )  call report_stat('FATAL','MSTINF'
     .   ,'htoglb/mstinf',' '
     .    ,'Error opening reference station.info file', ierr)  
      call read_stinf(100, 'REF', 0)
      close(100)       
      
*     Now loop over the other station.info files
      do i = 1, num_stf
          open(100, file=upd_stinf(i), status='old', iostat=ierr) 
          if(ierr.ne.0) call report_stat('FATAL','MSTINF'
     .      ,'htoglb/mstinf'
     .      ,upd_stinf(i),'Error opening input station.info',ierr)
          if( ierr.eq.0 ) then
              call read_stinf(100, 'UPD', i )
          endif
          close(100)
      end do 
  
*     Assign the site_indx array for later insertion and sorting
      do i=1,num_str
        site_indx(i) = i
      enddo
  
*     Now update global array with data from any RINEX files requested    
      do i = 1, num_rxf     
         urinex = 100
         open(100, file=rx_file(i), status='old', iostat=ierr)  
         if(ierr.ne.0) call report_stat('FATAL','MSTINF'
     .     ,'htoglb/mstinf',rx_file(i),'Error opening RINEX file',ierr)
         if( ierr.eq.0 ) then 
             call upd_from_rinex(i) 
         endif
         close(100)
      end do      

*     Now update the global array with data from any IGS log files requested
      do i = 1, num_igf     
         uigs = 100
         open(uigs, file=ig_file(i), status='old', iostat=ierr)  
         if(ierr.ne.0) then 
           call report_stat('FATAL','MSTINF'
     .     , 'htoglb/mstinf',ig_file(i),'Error opening IGS log file'
     .     , ierr)
         else  
           call report_stat('STATUS','MSTINF'
     .     , 'htoglb/mstinf',ig_file(i),'Opened IGS log file',0)
           read(uigs,'(a)') line 
           rewind(uigs)
           call upd_from_igslog(uigs,i)    
         endif
         close(100)
      end do      

*     Now update the global array with data from any SINEX files requested
      do i = 1, num_sxf     
         usnx = 100
         open(usnx, file=sx_file(i), status='old', iostat=ierr)  
         if(ierr.ne.0) call report_stat('FATAL','MSTINF'
     .     , 'htoglb/mstinf',ig_file(i),'Error opening SINEX file'
     .     , ierr)
         if( ierr.eq.0 ) then 
             call upd_from_sinex(usnx,i )    
         endif
         close(100)
      end do            

*     Now update the global array with data from any GIPSY sta files requested
      if( sr_file(1:2).ne.'  '.or. sa_file(1:2).ne.'  ' ) then
         if( sr_file(1:2).ne.'  ') then
           ustar = 100
           open(ustar, file=sr_file, status='old', iostat=ierr)  
           if(ierr.ne.0) call report_stat('FATAL','MSTINF'
     .     , 'htoglb/mstinf',sr_file,'Error opening GIPSY receiver file'
     .     , ierr)
         endif
         if( sr_file(1:2).ne.'  ') then
           ustaa = 101
           open(ustaa, file=sa_file, status='old', iostat=ierr)  
           if(ierr.ne.0) call report_stat('FATAL','MSTINF'
     .     , 'htoglb/mstinf',sa_file,'Error opening GIPSY antenna file'
     .     , ierr)
         endif
         if( ierr.eq.0 ) then 
             call upd_from_gipsy(ustar,ustaa)    
         endif
         if( sr_file(1:2).ne.'  ') close(100)
         if( sa_file(1:2).ne.'  ') close(101)
      endif 
          
*     Now sort and merge the entries we have in station.info
      call merge_stinf  
      
*     Now we are complete:  Write out the station.info to a file
      if( trimlen(out_stinf).ne.0 ) then
        if( allow_overwrite ) then
            open(200, file=out_stinf, status='unknown', iostat=ierr) 
            if(ierr.ne.0) call report_stat('FATAL','MSTINF'
     .         ,'htoglb/mstinf',out_stinf
     .       ,'Error opening output station.info (overwrite ok)',ierr)
        else
            open(200, file=out_stinf, status='new', iostat=ierr) 
            if(ierr.ne.0) call report_stat('FATAL','MSTINF'
     .         ,'htoglb/mstinf',out_stinf,
     .     'Error opening output station.info (overwrite not specified)'
     .       ,ierr)
        endif     
        call write_stinf(200)
        close(200)    
      endif

****  Thats all
      end
      
CTITLE DECODE_STINF_RUN 

      subroutine decode_stinf_run

      implicit none

*     Routine to decode the runstring for mstinf.  Note: Has some
*     code to support merging rinex headers but this is not used.

      include 'mstinf.h'

* LOCAL VARIABLES
* len_run  -- Length of runsting
* nr       -- Number of runstring entry
* rcpar    -- Function to return runstring
* trimlen  -- Length of string
* ierr     -- IOSTAT error number

      integer*4 len_run, nr, rcpar, trimlen, ierr, i
     .        ,  num_srf, num_saf

* decode_stf, decode_use -- Set true when station.info file or 
*     use site names are to be read

      logical decode_stf, decode_use

* runstring  -- Runstring returned by rcpar
* word       -- First few characters of runstring  
      
      character*256 message
      character*128 runstring  
      character*4   word

****  Loop over runstring getting all the elements.  The decoding
*     is driven by the option string passed 

*     Initialize variables for the program:
      allow_overwrite = .false.
      copy_comments   = .false.   
      nosort = .false. 
      nowrite = .false.
      num_stf = 0
      num_use = 0  
      num_srf = 0 
      num_saf = 0 
      decode_use = .false.
      decode_stf = .false.
      ref_stinf = ' '
      out_stinf = ' '
      max_slant = 0    
      dup_tol = 120.d0/86400.d0/365.25d0   
      snx_site = ' '  
      replace_entry = 'diff' 
      rx_open = .false. 
      write_apr = .false.
      do i=1,max_rxf
        rx_file(i) = ' '
      enddo
      do i=1,max_igf
        ig_file(i) = ' '
      enddo
      do i=1,max_sxf
        sx_file(i) = ' ' 
      enddo  
      sr_file = ' ' 
      sa_file = ' '
      dump_log = .false.

*     Loop decoding the runstring
      nr = 0
      len_run = 1
      do while ( len_run.gt.0 )
         nr = nr + 1
         len_run = rcpar(nr, runstring )
         if( len_run.gt.0 ) then

*            See if we have new option
             word = runstring           
             call casefold(word)     
*            See if file overwrite will be allowed
             if( word(1:2).eq.'-O' ) then
                 allow_overwrite = .true.

             else if( word(1:2).eq.'-C' ) then
                 copy_comments = .true.

*            See if station.info file names will follow.  Set 
*            the decode option to say get station.info names
             else if( word(1:3).eq.'-S ' ) then
                 decode_stf = .true.
                 decode_use = .false.

*            See if RINEX file name given next 
             else if( word(1:3).eq.'-R ' ) then  
                 nr = nr + 1
                 len_run = rcpar(nr,rx_file(1))
                 num_rxf = 1

*            See if IGS log file name given next 
             else if( word(1:2).eq.'-I' ) then  
                 nr = nr + 1
                 len_run = rcpar(nr,ig_file(1))
                 num_igf = 1

*            See if SINEX file name given next 
             else if( word(1:3).eq.'-X ' ) then  
                 nr = nr + 1
                 len_run = rcpar(nr,sx_file(1))
                 num_sxf = 1                   
       
*            See if GIPSY receiver file name given next
             else if( word(1:3).eq.'-SR' ) then  
                 nr = nr + 1
                 len_run = rcpar(nr,sr_file)
                 num_srf = 1

*            See if GIPSY antenna file name given next
             else if( word(1:3).eq.'-SA' ) then  
                 nr = nr + 1
                 len_run = rcpar(nr,sa_file)
                 num_saf = 1
  
*            See if site name from SINEX file given next
             else if( word(1:3).eq.'-XS' ) then
                 nr = nr + 1
                 len_run = rcpar(nr,snx_site)   
                 call uppers(snx_site)                      

*            See if reference station.info file name given next
             else if( word(1:2).eq.'-F' ) then
                 nr = nr + 1
                 len_run = rcpar(nr,ref_stinf) 

*            See if output station.info file name given next
             else if( word(1:2).eq.'-W' ) then
                 nr = nr + 1
                 len_run = rcpar(nr,out_stinf)    
                        
*            See if list of stations to be used will follow.  Set the
*            decode option to say get station-include names
             else if ( word(1:2).eq.'-U') then 
                 decode_stf = .false.
                 decode_use = .true. 
                 
*            See if station-include list file given next
             else if( word(1:2).eq.'-L' ) then
                 nr = nr + 1
                 len_run = rcpar(nr,runstring)
*                Now open and read the use list
                 call read_use( runstring )
                                        
*            See if max slant height argument given next  
             else if ( word(1:2).eq.'-H' ) then
                  nr = nr + 1
                  len_run = rcpar(nr,runstring)
                  read(runstring,*) max_slant

*            See if no-sort option given next
             else if( word(1:3).eq.'-NS' ) then
                  nosort =.true.

*            Set if replace option given next
             else if( word(1:3).eq.'-RE' ) then
                nr = nr + 1
                len_run = rcpar(nr,replace_entry)
                                 
*            Set if option to make RINEX entries open-ended given next
             else if( word(1:4).eq.'-RXO' ) then
                rx_open = .true.

*            Set if nowrite option given next
             else if( word(1:4).eq.'-NOW' ) then
                 nowrite = .true.
                                                       
*            See if write_apr option is given next
             else if( word(1:4).eq.'-APR' ) then
                 write_apr = .true.
       
*            See if no_gaps option is given next
             else if( word(1:4).eq.'-NOG' ) then
                 no_gaps = .true. 

*            See if tolerance for detecting duplicates is given next
             else if( word(1:2).eq.'-T' ) then
                  nr = nr + 1
                  len_run = rcpar(nr,runstring)
                  read(runstring,*) dup_tol     
                  dup_tol =  dup_tol/86400.d0/365.25d0  

*            Depending on decode option, read the names of station.info
*            files or station-include list.
             else if ( decode_stf ) then
                 num_stf = num_stf + 1
                 call check_max(num_stf, max_stf, 
     .                   'Number of input station.info files')
                 upd_stinf(num_stf) = runstring 

             else if ( decode_use ) then
                 num_use = num_use + 1
                 call casefold(runstring)
                 call check_max(num_use, max_site, 
     .                   'Number of use sites')
                 use_names(num_use) = runstring(1:4)
* MOD TAH 180802: Added debug option to find bad dates
             else  if ( word(1:4) .eq.'-DEB' ) then
                 dump_log = .true.
             else   ! Unknown command
                 call report_stat('WARNING','MSTINF','htoglb/mstinf'
     .              ,runstring(1:len_run),'Unknown runstring entry ',0)
             end if
         end if
      end do

****  Make sure we have enough information for run
      if( trimlen(ref_stinf).eq.0  ) then 
           call report_stat('WARNING','MSTINF','htoglb/mstinf',' '
     .         ,'Reference station.info missing from command line',0)
           call proper_runstring('mstinf.hlp','mstinf',1)
      endif

****  Let user know what is happening  
      write(message,'(a,i4,a,i4,a)')
     .    'There will be ',num_stf, ' station.info files and '
     .    ,num_rxf,' RINEX headers merged'  
      call report_stat('STATUS','MSTINF','htoglb/mstinf',' ',message
     .             ,0)
      write(message,'(a,i4,a,i4,a)')
     .    'There will be ',num_igf, ' IGS log files and '
     .    ,num_sxf,' SINEX headers merged'  
      call report_stat('STATUS','MSTINF','htoglb/mstinf',' ',message
     .             ,0)
      write(message,'(a,i4,a,i4,a)')
     .    'There will be ',num_srf,' GIPSY receiver files and '
     .    ,num_saf,' GIPSY antenna files merged'  
      call report_stat('STATUS','MSTINF','htoglb/mstinf',' ',message
     .             ,0)


      if( max_slant.gt.0 ) then    
          write(message,'(a,f7.3,a)') 'Heights greater than ',max_slant
     .         ,' m in rinex files will be taken to be slant heights'   
           call report_stat('WARNING','MSTINF','htoglb/mstinf',' '
     .                      ,message,0)
      endif

****  Thats all     
      return
      end

CTITLE READ_STINF

      subroutine read_stinf(unit, option, fn)

      implicit none

*     Routine to read station information file

      include 'mstinf.h'

* PASSED VARIABLES
* unit   -- Unit number of read file with 
* fn     -- File number (0 if reference file read)

      integer*4 unit, fn

* option -- Type read either initial reference (REF) or update with
*           addition files (UPD)
      character*(*) option

* LOCAL VARIABLES
* nl  -- Local counter of record number
* ierr -- IOSTAT errors
* i,j  -- Loop counters
* trimlen -- Length of string
* js   -- Station number for checking against list
* indx -- Pointer in string  
* doy_in -- day-of-year input to DECODE_VALUES, used for REPORT_STAT

      integer*4 nl, ierr, i, j, trimlen, js, indx, doy_in

* match   -- Logical set true when station.info entry matches one 
*     that we have previously encountered in the reference file

      logical match  

* warnings -- flag for call to decode_values in gamit/lib/rstnfo; set false here 
*             to turn off report_stat calls
      
       logical warnings

* line  -- Line read from file
* message -- Buffer for report_stat messages

      character*256 line,message 

* nlist -- number of items (columns) in station.info
* item_list -- tokens of items  
* values -- values of items      

      integer*4 nlist
      character*6 item_list(20) 
      character*20 values(20)     
         

* anth, antn, ante -- Antenna offsets up, north, east
      
      real*8 anth, antn,ante      
      real*8 antdaz  ! Antenna aligment from True N (deg).

* function in gamit/lib
 
      integer*4 itimdif   

***** If this is the reference file then initialize the record
*     and comment counters
      if( option(1:1).eq.'R' ) then
         num_str = 0
         num_rsi = 0
         num_comments = 0
         nl = 0
      else
*        Save number of lines
         nl = num_str
      end if
           
***** Read the header and decode the tokens (these routines in gamit/lib/rstnfo.f)
                                  
      call read_stnfohead( unit,line ) 
      call get_tokens( line, nlist, item_list )   

*     If this is reference file then save the token list for output file
      if( option(1:1).eq.'R' ) then
        ref_nlist = nlist
        do i = 1, nlist
          ref_item_list(i) = item_list(i)
        enddo
      endif   

****  Loop over the file until we hit an eof   
      ierr = 0
      do while ( ierr.eq.0 )
         read(unit,'(a)',iostat=ierr) line
         if( ierr.eq.0 .and. trimlen(line).gt.0 .and.
     .       line(1:1).eq.' ' ) then

*            OK, valid line: decode
             nl = nl + 1
             call check_max(nl, max_str,
     .                      'Number of station.info records')  
             call read_values(line,nlist,item_list,values,comment(nl))  
             doy_in = istart(1,nl)  
* MOD RWK 050719: Add warnings variable to call and set false to avoid report_stat calls
*                 reading all of station.info                         
             warnings = .false.  
cd             print *,'calling DECODE_VALUES nl rcvcod ',nl,rcvcod(nl)
             call decode_values( nlist, item_list, values, doy_in  
     .                         , sitcod(nl), sname(nl)
     .                         , anth, antn, ante, antdaz, rcvcod(nl)
     .                         , antcod(nl),  rctype(nl), rcvrsn(nl) 
     .                         , rcvers(nl), anttyp(nl), antsn(nl)
     .                         , htcod(nl), radome(nl), swver(nl)
     .                         , istart(1,nl), istop(1,nl), warnings )   
             dUNE(1,nl) = anth
             dUNE(2,nl) = antn
             dUNE(3,nl) = ante  
* MOD TAH 200203: Save the antenna azimuth
             dAntAZ(nl) = antdaz

             call casefold(sitcod(nl))  
c??             call trimlead(comment(nl)) 
             if( itimdif(istart(1,nl),istop(1,nl)).gt.0 ) then
               write(message,'(a,a4,1x,i4,1x,i3,a,i2)') 
     .            'Start later than stop for ',sitcod(nl),istart(1,nl)
     .                ,istart(2,nl),' in file ',nl
               call report_stat('WARNING','MSTINF','htoglb/mstinf',' '
     .              ,message,0) 
             endif
           
* rwk 090410: No longer ignore matching site/yr/doy entries, but rather 
*            flag as new to be compared later
             if( option(1:3).ne.'REF' ) then
               newflag(nl) = .true.
*              See if we should get this station
               if( num_use.gt.0 ) then 
                 indx = 0
                 call get_cmd(sitcod(nl), rsi_names, num_rsi,js,indx)
                 if( js.gt.0 ) then
                    nl = nl - 1 
                 end if
               end if
             else
*              See if we have a list of sites to use from the reference
*              station.info and if this is one of them
               if( num_use.gt.0 ) then 
                 indx = 0
                 call get_cmd(sitcod(nl), use_names, num_use,js,indx)
                 if( js.le.0 ) then
*                  See if we have this station name already
                   indx = 0
                   call get_cmd(sitcod(nl), rsi_names, num_rsi,js,indx)
                   if( js.le.0 ) then
                      num_rsi = num_rsi + 1
                      call check_max(num_rsi, max_site,
     .                        'number of excluded reference sites')
                      rsi_names(num_rsi) = sitcod(nl)
                   endif
                   nl = nl - 1
                 endif
               endif
             endif    

* rwk/tah 090731: Remove any entries that have no receiver or no antenna
             if((rcvcod(nl)(1:2).eq.' '.and.rctype(nl)(1:2).eq.'--').or.
     .          (antcod(nl)(1:2).eq.' '.and.anttyp(nl)(1:2).eq.'--'))
     .         then
                nl = nl - 1
             endif

         else if( ierr.eq.0 .and. line(1:1).ne.' ' .and.
     .       num_comments.lt.max_comments ) then
             num_comments = num_comments + 1
             stinf_com(num_comments) = line(1:132)
         end if
      end do

* Temporary(?): rename deprecated ASHTECH MICROZ
      if( rctype(nl)(1:14).eq.'ASHTECH MICROZ' ) then
         call report_stat('WARNING','MSTINF','htoglb/mstinf',' '
     ., 'Obsolete rcvr type ASHTECH MICROZ replaced by ASHTECH UZ-12',0)
        rctype(nl) = 'ASHTECH UZ-12       ' 
      endif

*     Tell user what is happening 
      if( option(1:1).eq.'R' ) then  
                       
         write(message,'(a,a,i7,a)') 'Reference file '
     .          ,ref_stinf(1:trimlen(ref_stinf)), nl,' entries found'
         call report_stat('STATUS','MSTINF','htoglb/mstinf',' '
     .                   ,message,0)
         num_ref = nl
      else      
         write(message,'(a,a,1x,i6,a)') 'Update file '
     .         ,upd_stinf(fn)(1:trimlen(upd_stinf(fn))),nl-num_str
     .         ,' entries found'
         call report_stat('STATUS','MSTINF','htoglb/mstinf',' '
     .                   ,message,0)
      end if
      num_str = nl 

****  Thats all
      return 
      end

CTITLE CHECK_MAX

      subroutine check_max( n, max, type )

      implicit none

*     Routine to check if we exceed the maxiumum number in some
*     quanity

* PASSED VARIABLES
* n, max  -- Current number and max for quanity

      integer*4 n, max

* type    -- String with description of quanity

      character*(*) type  

* LOCAL VARIABLE

      character*256 message


****  See if we have exceeded limit
      if( n.gt.max ) then                                
         write(message,'(2a,i6)') type,' exceeded; limit ',max
         call report_stat('FATAL','MSTIN2F','htoglb/mstinf',' '
     .       ,message,0)
          stop 'Limits exceeded'
      end if

      return 
      end


CTITLE MERGE_STINF

      subroutine merge_stinf

      implicit none

*     Routine to merge the station.info entries already read and
*     to eliminate any duplicates.
*               
      include 'mstinf.h' 
             
* FUNCTION IN gamit/lib
* itimdif --integer function to difference time as y,doy,h,m,s

* LOCAL VARIABLES   
* i,j,k  -- Loop counters
* js     -- Site number, returns -1 if new site
* indx   -- Pointer in string
* ns     -- Generic station number
* num_ent -- Number of entries at specific site 
* max_ent -- Maximum number of entries for any one site
* ns_with_max -- Site with maximum
* ik, ikp     -- indices for k and k+1 entry at a site 
* jk jkm jkp  -- indices for j, j-1, j+1 entries
* iswap       -- Swap value for exchanging k amd k+1 is needed.
* last        -- Last entry for each station    
* istype      -- type of span overlap
* message -- Buffer for report_stat call  

      integer*4 i,j,k, js, indx, ns, max_ent, ns_with_max, ik, ikp
     .        , jk, jkm, jkp, istype, iswap, num_ent, last, itimdif  

* tk, tkp     -- Decimal years corresponding to these value

      real*8 tk, tkp    

      character*256 message  

* MOD TAH 200205: Added list of sites that are in the new entries
*     read (not reference station.info) so that they can be
*     removed from reference entries if -rep replace option 
*     is used.
* MOD TAH 200417: Changed dimensioning of new_sites to max_site, added
*     fatal call rather then just stop with message.
      integer*4 new_sites(max_site)    ! Site numbers in new entries 
     .,         num_new               ! number of sites in new entries
      logical   match                 ! No match in list if new sites

      logical different,debug
      data debug/.false./
                      
****  First scan and get all the station names
      num_site = 0
      num_new  = 0
      do i = 1, num_str
         indx = 0
         call get_cmd(sitcod(i), site_names, num_site, js, indx)
         if( js.eq.-1 ) then
             num_site = num_site + 1
             call check_max(num_site, max_site,'Number of stations')
             site_names(num_site) = sitcod(i)
             js = num_site
         end if
         if( i.gt. num_ref ) then  ! We are past the referece station.info
                                   ! lines, so this is a site in updated entries.
C            See if this site is already in new list
             match = .false.
             do j = 1,num_new
                if( js .eq. new_sites(j) ) then
                   match = .true.
                   exit
                endif
             enddo
             if( .not.match ) then  ! No matchng site number so add
                 num_new = num_new + 1
                 if( num_new.gt.max_site ) 
     .              call report_stat('FATAL','MSTINF','htoglb/mstinf',
     .                    ' ','Too many new sites',num_new)
                 new_sites(num_new) = js
             end if
          endif

      end do
                    
      write(message,'(a,i6,a)') 'Station.info files have '
     .      ,num_site,' sites' 
      call report_stat('STATUS','MSTINF','htoglb/mstinf',' '
     .                ,message,0)

*     Now sort the site names  
**** rwk 090223: allsite_indx appears never to be used in current code
c**      call sort_snames(site_names, allsite_indx, num_site) 

*     Now sort the names in the station.info records 
      call sort_snames(sitcod, site_indx, num_str )
           

*     Now find the start of each station
      ns = 1
      site_start(1) = 1
      max_ent = 0
      do i = 1, num_str-1
         if ( sitcod(i+1).ne.sitcod(i) ) then
             ns = ns + 1
             site_start(ns) = i+1
         end if
         num_ent = site_start(ns) - site_start(ns-1)
         if( num_ent.gt.max_ent ) then
             max_ent = num_ent
             ns_with_max = ns - 1
         end if
      end do   
      write(message,'(3a,i5)') 'Site with most entries '
     .     , site_names(ns_with_max),' Maximum ',max_ent
      call report_stat('STATUS','MSTINF','htoglb/mstinf',' '
     .                 ,message,0)
           
CD DEBUG   
      if( debug ) then
         print *,'List before sort (i site_indx(i) )'
         do i = 1, num_str
          ik = abs(site_indx(i))
          write(*,'(2i4,1x,l1,1x,a4,2(1x,i4,1x,i3))') 
     .      i,site_indx(i),newflag(ik),sitcod(i)
     .      ,istart(1,ik),istart(2,ik),istop(1,ik),istop(2,ik)
        enddo    
      endif

****  Now within each site, sort list in ascending time order by start 
      do i = 1, num_site

*        Loop over entries for this sites, rearranging if need be
*        Get end index (differs if last site)
         if( i.eq.num_site ) then
            last = num_str
         else
            last = site_start(i+1) - 1
         end if  
         if( debug )  print *
     .      ,'sorting ',site_names(i),' 1s last ',site_start(i),last
         do j = site_start(i), last-1
           do k = site_start(i), last-(j-site_start(i)+1)
              ik = abs(site_indx(k))
              tk = istart(1,ik)+ (istart(2,ik) + (istart(3,ik) + 
     .           (istart(4,ik) + 
     .            istart(5,ik) /60.d0) /60.d0) /24.d0 )/365.25d0 
              ikp = abs(site_indx(k+1))   
              tkp = istart(1,ikp)+ (istart(2,ikp) + (istart(3,ikp) + 
     .           (istart(4,ikp) + 
     .            istart(5,ikp) /60.d0) /60.d0) /24.d0 )/365.25d0 
              if( tk.gt.(tkp+dup_tol) ) then
*                 Time at k is greater than time at k+1, so switch
*                 indices
                  iswap = site_indx(k)
                  site_indx(k) = site_indx(k+1)
                  site_indx(k+1) = iswap    
              endif
            enddo
         enddo
      enddo

*     List should now be sorted.
CD DEBUG   
      if( debug ) then
        print *,'List after sort (i site_indx(i) )'
        do i = 1, num_str
          ik = abs(site_indx(i))
          write(*,'(2i4,1x,l1,1x,a4,2(1x,i4,1x,i3))') 
     .      i,site_indx(i),newflag(ik),sitcod(i)
     .      ,istart(1,ik),istart(2,ik),istop(1,ik),istop(2,ik)
        enddo 
      endif
CD END DEBUG  
   

*     Check if new entries should replace old ones 
      do i = 1, num_site
        if( i.eq.num_site ) then
          last = num_str
        else
          last = site_start(i+1) - 1
        endif                        
        if( debug )  print *,'check loop start site i site_start last '
     .      ,i,site_start(i),last   
        if(last.gt.num_str)   stop
        do j = site_start(i), last    
          jk = abs(site_indx(j))
          if( (site_start(i)-last+1).eq.1 ) then
*           only 1 entry, no action required
          elseif( newflag(jk) ) then 
*           first compare with the entry before this one if it is old and
*           hasn't been previously flagged 
            jkm = abs(site_indx(j-1)) 
            if( j.ne.site_start(i).and..not.newflag(jkm).and.
     .                                site_indx(j-1).gt.0 ) then   
               if( debug ) then
                  print *,'j jk start stop ',j,jk,istart(1,jk)
     .               ,istart(2,jk),istop(1,jk),istop(2,jk)
                  print *,'j jkm start stop ',j,jkm,istart(1,jkm)
     .               ,istart(2,jkm),istop(1,jkm),istop(2,jkm)
               endif
               call span_type( istart(1,jk),istop(1,jk)
     .                       ,istart(1,jkm),istop(1,jkm),istype)    
               if( debug ) print *,'istype ',istype        
               if(istype.eq.1) then
*               new entry matches old: if different keep new, warn, and flag old
                call check_change(sitcod(j),jk,jkm,different)
                if( different ) then                                          
                  if( replace_entry.ne.'none' ) then
                    call print_warning(sitcod(j),istype
     .                 ,istart(1,jk),istop(1,jk),'replace')
                     site_indx(j-1) = -abs(site_indx(j-1))
                  else          
                    call print_warning(sitcod(j),istype
     .                 ,istart(1,jk),istop(1,jk),'ignore ')
                    site_indx(j) = -abs(site_indx(j))
                  endif
                else 
                  if( replace_entry.eq.'all ' ) then  
                    site_indx(j-1) = -abs(site_indx(j-1)) 
                  else
                    site_indx(j) = -abs(site_indx(j)) 
                  endif
                endif          
              elseif( istype.eq.2 ) then
*               new entry starts before and ends within span: 
*               if different adjust old and keep both; if same adjust new span and flag old
                call check_change(sitcod(j),jk,jkm,different)
                if( different ) then  
                  if( replace_entry.eq.'none' ) then
                     call print_warning(sitcod(j),istype
     .                       ,istart(1,jk),istop(1,jk),'ignore ')
                  else   
                    call print_warning(sitcod(j),istype
     .                       ,istart(1,jk),istop(1,jk),'update ')
                    do k=1,5
                      istart(k,jkm) = istop(k,jk)
                    enddo     
                  endif
                else   
                  do k=1,5
                    istop(k,jk) = istop(k,jkm)
                  enddo   
                  call print_warning(sitcod(j),istype
     .                      ,istart(1,jk),istop(1,jk),'update ')
                  site_indx(j-1) = -abs(site_indx(j-1))
                endif
              elseif( istype.eq.4 ) then
*               new entry starts at or within and ends before station.info span: 
*              if different replace old after new start; if same flag new 
               call check_change(sitcod(j),jk,jkm,different)
               if( different ) then  
                 if( replace_entry.eq.'none' ) then    
                    call print_warning(sitcod(j),istype
     .                     ,istart(1,jk),istop(1,jk),'ignore ')
                    site_indx(j) = -abs(site_indx(j))   
                 else  
                   call print_warning(sitcod(j),istype
     .                      ,istart(1,jk),istop(1,jk),'update ') 
                   do k=1,5
                     istop(k,jkm) = istart(k,jk)
                   enddo
                 endif
               else     
                 if( replace_entry.eq.'all ' ) then
                    site_indx(j-1) = -abs(site_indx(j-1))
                 else
                   site_indx(j) = -abs(site_indx(j))
                 endif
               endif 
             elseif( istype.eq.5 ) then
*              new entry starts at and ends after station.info span: 
*              replace old with new 
                call check_change(sitcod(j),jk,jkm,different)
                if( replace_entry.eq.'none' ) then    
                   call print_warning(sitcod(j),istype
     .                      ,istart(1,jk),istop(1,jk),'ignore ')
                    site_indx(j) = -abs(site_indx(j))
                else  
                   call print_warning(sitcod(j),istype
     .                       ,istart(1,jk),istop(1,jk),'update ')
                   site_indx(j-1) = -abs(site_indx(j-1))       
                endif    
              elseif( istype.eq.6 ) then  
*               new entry starts within or at end of station.info span and ends afterward: 
*               if different, keep both and adjust spans; if same, keep new flag old
                call check_change(sitcod(j),jk,jkm,different)
                if( different ) then  
                  if( replace_entry.eq.'none' ) then    
                     call print_warning(sitcod(j),istype
     .                      ,istart(1,jk),istop(1,jk),'ignore ')
                     site_indx(j) = -abs(site_indx(j))
                  else  
                    call print_warning(sitcod(j),istype
     .                       ,istart(1,jk),istop(1,jk),'update ')
                    do k=1,5
                       istop(k,jkm) = istart(k,jk)
                    enddo          
                  endif
                else
                  do k=1,5
                    istart(k,jk) = istart(k,jkm)
                  enddo  
                  site_indx(j-1) = -abs(site_indx(j-1))
                endif   
*             endif for span type 
              endif
*           endif for checking new vs previous entry
            endif                 
*           now compare with following entry if it is old and has not
*           been previously flagged
            jkp = abs(site_indx(j+1))
            if( j.ne.last .and. .not.newflag(jkp) .and.
     .                            site_indx(j+1).gt.0 ) then 
               call span_type( istart(1,jk),istop(1,jk)
     .                       ,istart(1,jkp),istop(1,jkp),istype)    
               if(debug) print *,'j jk jkp istype ',j,jk,jkp,istype
              if(istype.eq.1) then
*               new entry matches old: if different keep new, warn, and flag old
                call check_change(sitcod(j),jk,jkp,different)
                if( different ) then 
                  if( replace_entry.eq.'none' ) then
                    site_indx(j) = -abs(site_indx(j))   
                    call print_warning(sitcod(j),istype
     .                       ,istart(1,jk),istop(1,jk),'ignore ')
                  else
                     call print_warning(sitcod(j),istype
     .                 ,istart(1,jk),istop(1,jk),'replace')
                     site_indx(j+1) = -abs(site_indx(j+1))
                  endif
                else
                  if( replace_entry.eq.'all ') then   
                    site_indx(j-1) = -abs(site_indx(j-1))  
                  else
                    site_indx(j) = -abs(site_indx(j))   
                  endif
                endif
              elseif( istype.eq.2 ) then
*               new entry starts before and ends within span: 
*               if different adjust old and keep both; if same adjust new span and flag old
                call check_change(sitcod(j),jk,jkp,different)
                if( different ) then  
                  if( replace_entry.eq.'none' ) then
                     call print_warning(sitcod(j),istype
     .                       ,istart(1,jk),istop(1,jk),'ignore ')
                  else   
                    call print_warning(sitcod(j),istype
     .                       ,istart(1,jk),istop(1,jk),'update ')
                    do k=1,5
                      istart(k,jkp) = istop(k,jk)
                    enddo     
                  endif
                else   
                  do k=1,5
                    istop(k,jk) = istop(k,jkp)
                  enddo   
                  call print_warning(sitcod(j),istype
     .                      ,istart(1,jk),istop(1,jk),'update ')
                  site_indx(j+1) = -abs(site_indx(j+1))
                endif
              elseif( istype.eq.3 ) then
*               new entry starts before and ends at end of station.info span 
*               if different replace old with new; if same adjust new span, flag old
                call check_change(sitcod(j),jk,jkp,different) 
                if( replace_entry.eq.'none' ) then
                  call print_warning(sitcod(j),istype
     .                       ,istart(1,jk),istop(1,jk),'ignore ')
                  site_indx(j) = -abs(site_indx(j))
                else  
                  call print_warning(sitcod(j),istype
     .                      ,istart(1,jk),istop(1,jk),'update ')
                  site_indx(j+1) = -abs(site_indx(j+1))
                endif
              elseif( istype.eq.4 ) then
*              new entry starts at and ends before station.info span: 
*              if different replace old after new start; if same flag new 
                call check_change(sitcod(j),jk,jkp,different)
                if( different ) then  
                  if( replace_entry.eq.'none' ) then    
                     call print_warning(sitcod(j),istype
     .                      ,istart(1,jk),istop(1,jk),'ignore ')
                     site_indx(j) = -abs(site_indx(j))
                  else  
                    call print_warning(sitcod(j),istype
     .                       ,istart(1,jk),istop(1,jk),'update ')
                    if( istop(1,jkp).lt.2100 ) then   
                      do k=1,5
                       istart(k,jkp) = istop(k,jk)
                      enddo    
                    else
                      do k=1,5
                        istop(k,jk) = istop(k,jkp)
                      enddo  
                    endif   
                  endif
                else  
                  if(replace_entry.eq.'all ') then
                    site_indx(j+1) = -abs(site_indx(j+1))
                  else  
                    site_indx(j) = -abs(site_indx(j))
                  endif
                endif                               
              elseif( istype.eq.5 ) then
*              new entry starts at and ends after station.info span: 
*              replace old with new 
                call check_change(sitcod(j),jk,jkp,different)
                if( replace_entry.eq.'none' ) then    
                   call print_warning(sitcod(j),istype
     .                      ,istart(1,jk),istop(1,jk),'ignore ')
                    site_indx(j) = -abs(site_indx(j))
                else  
                   call print_warning(sitcod(j),istype
     .                       ,istart(1,jk),istop(1,jk),'update ')
                   site_indx(j+1) = -abs(site_indx(j+1))
                endif    
              elseif( istype.eq.7 ) then
*               new entry is wider than station.info span: replace old with new    
                if( replace_entry.eq.'none') then  
                  call print_warning(sitcod(j),istype
     .                      ,istart(1,jk),istop(1,jk),'ignore ')
                  site_indx(j) = -abs(site_indx(j))      
                else
                  call print_warning(sitcod(j),istype
     .                   ,istart(1,jk),istop(1,jk),'update ')
                  site_indx(j+1) = -abs(site_indx(j+1))
                endif
*             endif for span type
              endif
*           endif of checking new vs following entry
            endif                                   
*         endif on newflag
          endif
*       enddo on loop through entries for this site
        enddo
*     enddo on loop through sites
      enddo
                                  

CD DEBUG 
      if( debug ) then  
        print *,'List after new/old merge (i site_indx(i) )'
        do i = 1, num_str
          ik = abs(site_indx(i))
          write(*,'(2i4,1x,l1,1x,a4,2(1x,i4,1x,i3))') 
     .      i,site_indx(i),newflag(ik),sitcod(i)
     .      ,istart(1,ik),istart(2,ik),istop(1,ik),istop(2,ik)
        enddo    
      endif


*     Finally loop through all the entries again, adjusting start and 
*     stop times to avoid overlaps     
*   rwk this doesn't work with a complicated file, and may not be necessary in
*        any case, with the current sort scheme
      
c     do i =1, num_str-1 
c       if( site_indx(i).gt.0 ) then     
c          print *,'i sitcod sitcod1 ',i,sitcod(i),sitcod(i+1)
c         if( sitcod(i).eq.sitcod(i+1) ) then 
c           ik = abs(site_indx(i))    
c            print *,'i,ik istart istop '
c     .          ,i,ik,(istart(j,ik),j=1,5),(istop(j,ik),j=1,5)     
c           ikp = abs(site_indx(i+1))
c            print *,'i, ikp istart istop '
c     .          ,i,ikp,(istart(j,ikp),j=1,5),(istop(j,ikp),j=1,5)     
c           if( itimdif(istop(1,ik),istart(1,ikp)).gt.0 ) then  
c             do j=1,5
c               istop(j,ik) = istart(j,ikp)   
c             enddo  
c              print *,'ikp istop ',ikp,(istop(j,ikp),j=1,5)
c           endif
c         endif 
c       endif
c     enddo
                                       
CD DEBUG   
      if( debug ) then
        print *,'List after final sort (i site_indx(i) )'
        do i = 1, num_str
          ik = abs(site_indx(i))
          write(*,'(2i4,1x,l1,1x,a4,2(1x,i4,1x,i3))') 
     .      i,site_indx(i),newflag(ik),sitcod(i)
     .      ,istart(1,ik),istart(2,ik),istop(1,ik),istop(2,ik)
        enddo    
      endif

****  Thats all
      return 
      end       
  
CPRINT_WARNING

      Subroutine print_warning(site,istype,istart,istop,action)

      implicit none

      character*4 site      
      character*7  action
      character*256 message
                 
      integer*4 istart(5),istop(5),istype,i
                               
      write(message,'(a,i1,a,a4,2(1x,2i4,3i3),1x,a7)') 
     .   'New entry type ',istype,' for ',site
     .   ,(istart(i),i=1,5),(istop(i),i=1,5)
     .   ,action
      call report_stat('WARNING','MSTINF','htoglb/mstinf'
     .   ,' ',message,0)
      return
      end


CTITLE DYEAR

      real*8 function dyear(itime)
             
*     Convert YYYY DDD HH MM SS to decimal year for  for sorting
      integer*4 itime(5)
      
      dyear = itime(1) + 
     .         ( itime(2) + 
     .          (itime(3) + (itime(4) + 
     .           itime(5) /60.d0) /60.d0) /24.d0 )/365.25d0 
      return
      end

CTITLE CHECK_CHANGE
                             
*     Compare two sets of station.info entries to see if they are different.

      Subroutine check_change(sitcod1,i1,i2,different)

      implicit none
  
* Input    
      character*4 sitcod1
*       sitcod1   name of site (for warning printout)
      integer*4 is,i1,i2  
*       i1, i2   indices in the station.info arrays of the two entries to be compared
       
* Output
      logical different
*        set true if there are any differences

* Arrays are all in common /stinf_recs/ in mstinf.h
      include 'mstinf.h'

* Function
      real*8 dyear
     
* Local    
      integer*4 i   
      real*8 t1s,t1e,t2s,t2e
      character*256 message
      logical debug/.false./
               
      if( debug ) then
        print *,'CHECK_CHANGE i1 i2 ',i1,i2 
        print *,sitcod1,rcvcod(i1),istart(1,i1),istart(2,i1)
        print *,sitcod1,rcvcod(i2),istart(1,i2),istart(2,i2)
      endif
  
   

* First see if there are differences  

      different = .false. 
      if( htcod(i1).ne.htcod(i2) ) then
         write(message,'(a4,1x,2i4,1x,a,a5,4x,a5)')       
     .      sitcod1,istart(1,i1),istart(2,i1),
     .        'HtCod mismatch (old, new): '
     .        ,htcod(i2),htcod(i1)
          call report_stat('WARNING','MSTINF'
     .                       ,'htoglb/mstinf',' ',message,0)
        different = .true.
      endif
      if( dUNE(1,i1).ne.dUNE(1,i2)  .or.
     .    dUNE(2,i1).ne.dUNE(2,i2)  .or.
     .    dUNE(3,i1).ne.dUNE(3,i2) ) then
        write(message,'(a4,1x,2i4,1x,a,3f7.4,2x,3f7.4)')            
     .      sitcod1,istart(1,i1),istart(2,i1),
     .        'Ant Offset mismatch (old, new): '
     .        ,(dUNE(i,i2),i=1,3),(dUNE(i,i1),i=1,3)
          call report_stat('WARNING','MSTINF'
     .                       ,'htoglb/mstinf',' ',message,0)
        different = .true.
      endif
* MOD TAH 200203; Check Azimuth value
      if( dAntAZ(i1).ne.dAntAZ(i2) ) then
        write(message,'(a4,1x,2i4,1x,a,f5.1,2x,F5.1)')            
     .      sitcod1,istart(1,i1),istart(2,i1),
     .        'Ant Azimuth mismatch (old, new): '
     .        ,dAntAZ(i2), dAntAZ(i1)
          call report_stat('WARNING','MSTINF'
     .                       ,'htoglb/mstinf',' ',message,0)
        different = .true.
      endif

      if( rcvcod(i1).ne.rcvcod(i2) ) then
         write(message,'(a4,1x,2i4,1x,a,a6,4x,a6)')    
     .      sitcod1,istart(1,i1),istart(2,i1),
     .        'RcvCod mismatch (old, new): '
     .        , rcvcod(i2),rcvcod(i1)
          call report_stat('WARNING','MSTINF'
     .                       ,'htoglb/mstinf',' ',message,0)
        different = .true.
      endif      
      if( rcvers(i1).ne.rcvers(i2))  then
        write(message,'(a4,1x,2i4,1x,a,a20,4x,a20)')       
     .      sitcod1,istart(1,i1),istart(2,i1),
     .        'RcVers mismatch (old, new): '
     .        , rcvers(i2),rcvers(i1)
          call report_stat('WARNING','MSTINF'
     .                       ,'htoglb/mstinf',' ',message,0)
        different = .true.
      endif
      if(  rcvrsn(i1).ne.rcvrsn(i2) ) then
        write(message,'(a4,1x,2i4,1x,a,a20,4x,a20)')  
     .      sitcod1,istart(1,i1),istart(2,i1),
 
     .        'RcvrSN mismatch (old, new): '
     .        , rcvrsn(i2),rcvrsn(i1)
          call report_stat('WARNING','MSTINF'
     .                       ,'htoglb/mstinf',' ',message,0)
        different = .true.
      endif                    
      if( antcod(i1).ne.antcod(i2) ) then     
         write(message,'(a4,1x,2i4,a,a6,4x,a6)')  
     .      sitcod1,istart(1,i1),istart(2,i1),
     .        'AntCod mismatch (old, new): '
     .        , antcod(i2),antcod(i1)
          call report_stat('WARNING','MSTINF'
     .                       ,'htoglb/mstinf',' ',message,0)
        different = .true.
      endif
      if(  antsn(i1).ne.antsn(i2) ) then
        write(message,'(a4,1x,2i4,1x,a,a20,4x,a20)')   
     .      sitcod1,istart(1,i1),istart(2,i1),
     .        'AntSN mismatch (old, new): '
     .        , antsn(i2),antsn(i1)
          call report_stat('WARNING','MSTINF'
     .                       ,'htoglb/mstinf',' ',message,0)
        different = .true.
      endif                    
      if( radome(i1).ne.radome(i2) ) then    
        write(message,'(a4,1x,2i4,1x,a,a20,4x,a20)')   
     .      sitcod1,istart(1,i1),istart(2,i1),
     .        'Radome mismatch (old, new): '
     .        , radome(i2),radome(i1)
          call report_stat('WARNING','MSTINF'
     .                       ,'htoglb/mstinf',' ',message,0)
        different = .true.
      endif                    
      if( different) then
       if( debug ) then
       print *,'after found difference '
       print *,' 1 start stop ',(istart(i,i1),i=1,5),(istop(i,i1),i=1,5)
       print *,'  ',rcvcod(i1),rcvers(i1),rcvrsn(i1),antcod(i1)
     .                ,antsn(i1),radome(i1),(dUNE(i,i1),i=1,3)
       print *,' 2 start stop ',(istart(i,i2),i=1,5),(istop(i,i2),i=1,5)
       print *,'  ',rcvcod(i2),rcvers(i2),rcvrsn(i2),antcod(i2)
     .           ,antsn(i2),radome(i2),(dUNE(i,i2),i=1,3)                   
      endif
      endif

      return
      end 

CTITLE SPAN_TYPE
    
      Subroutine span_type( istart1,istop1,istart2,istop2,istype )

      implicit none

*   Determine which of six cases applies to times spans being compared:
 
*   There are four 6 cases to consider (assume entry 1 is new, entry 2 the existing station.info):   
*     0. No overlap in spans     
*     1. New entry span matches old
*          If different, keep new, warn, and flag old                 
*     2. New entry starts before and ends within station.info span 
*          If same, adjust new span, flag old
*          If different, keep both, adjust old times, and warn
*     3. New entry starts before and ends at end of station.info span
*          Replace old with new
*     4. New entry starts at or within and ends before station.info span 
*          If same, flag new
*          If different, replace old after new start; if finite old, warn 
*     5. New entry starts at and ends after station.info span
*          Replace old with new             
*     6. New entry starts at end or within and ends after station.info span
*          If same, adjust new span and flag old
*          If different, keep both and adjust spans 
*     7. New entry wider than station.info span            
*          Replace with new entry
* Open ended entries in both station.info and IGS log files have stop time 9999 999 0 0 0
                 
* Input
      integer*4 istart1(5),istop1(5),istart2(5),istop2(5)

* Output
      integer*4 istype
           
* Local       
      real*8 t1s,t1e,t2s,t2e
      character*256 message
* Tolerance for determine coincident times is 10 minutes
      real*8 tol/1.9d-5/
      logical debug/.false./
           
* Function
      real*8 dyear

*   Determine the type of overlap (use decimal years for convenience)
                            
      t1s = dyear(istart1)
      t1e = dyear(istop1)
      t2s = dyear(istart2)
      t2e = dyear(istop2)                 
CDEBUG   
      if( debug ) then
        print *,'SPAN_TYPE '
        print *,' 1 ',istart1,istop1,t1s,t1e
        print *,' 2 ',istart2,istop2,t2s,t2e
      endif
                 
c*  I don't understand why I changed to the next line
c*  from the following one for the no-overlap case
c*      if( t1e.gt.t2e .and. (t1s-tol).gt.t2e .or.
c*     .    t1s.lt.t2s .and. (t1e+tol).lt.t2s ) then
c*  This appears to have been close but I've modified slightly
c*      if( (t1s+tol).gt.t2e .or. (t1e-tol).lt.t2s ) then  
c*  Try this version
      if( ((t2s+tol).ge.t1e) .or. ((t1s+tol).ge.t2e) ) then
        istype = 0                    
      elseif( dabs(t1s-t2s).lt.tol .and. dabs(t1e-t2e).lt.tol ) then
        istype = 1 
      elseif( (t1s+tol).lt.t2s .and. 
     ,        ( (t1e-tol).gt.t2s. and.(t1e-tol).lt.t2e) ) then
        istype = 2
      elseif( (t1s+tol).lt.t2s .and. dabs(t1e-t2e).lt.tol ) then
        istype = 3
      elseif( (t1s+tol).gt.t2s .and. (t1e+tol).lt.t2e ) then   
        istype = 4     
      elseif( dabs(t1s-t2s).lt.tol .and. t1e.gt.t2e ) then
        istype = 5
      elseif( (t1s-tol).lt.t2e .and. (t1e+tol).gt.t2e ) then
        istype = 6
      elseif( (t1s+tol).lt.t2s .and. (t1e-tol).gt.t2e  ) then   
        istype = 7
      else
         write(message,'(a,2f12.6,a,2f12.6)') 
     .     'Logic problem, new span =',t1s,t1e,' old span=',t2s,t2e 
         call report_stat('FATAL','MSTINF','htoglb/mstinf',' '
     .                   ,message,0)
      endif                      
CDEBUG 
      if( debug )  print *,'SPAN_TYPE istype ',istype
      return
      end


CTITLE SORT_SNAMES
 
      subroutine sort_snames( names, list, num)

      implicit none
 
*     This routine uses an exchange sort algormithm to sort
*     the list names into ascending alphabetical order.  
*     There are num values in list and names
 
*   num     - Number of values to be sorted
*   list(num)  - List to be sorted in to ascending order.
 
      integer*4 num, list(num)

*   names(num) -- Names to be sorted
      character*4 names(num)
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
*   smallest_one    - pointer to least value

      integer*4 i,j, smallest_one

*   swap   - swap string
      character*4  swap
 


****  Generate index list     
* This now done in mstinf
c      do i = 1, num
c         list(i) = i
c      end do
 
****  Start loop using exchange sort
      do i = 1, num
          smallest_one = i
          do j = i+1, num
              if( names(j).lt. names(smallest_one) ) then
                  smallest_one = j
              end if
          end do 
 
*****     See if we should swap
          if( smallest_one.gt. i ) then
              swap = names(smallest_one)
              names(smallest_one) = names(i)
              names(i) = swap
*             Now swap the index
              j = list(smallest_one)
              list(smallest_one) = list(i)
              list(i) = j

          end if
      end do
 
***** Thats all.  Now sorted in ascending order
      return
      end

CTITLE WRITE_STINF

      subroutine write_stinf(unit) 

      implicit none

*     Routine to write out the updated station.info file

      include 'mstinf.h'

* PASSED VARIABLES
* unit   -- Unit number for output file

      integer*4 unit

* LOCAL VARIABLES
* i, j  -- Loop counter  
* ierr -- IOSTAT error
* ik, il   -- Point to record to be printed
* trimlen -- Length of string   
* ioerr -- i/o return error    

      integer*4 i,j, ierr, ik, il, trimlen, ioerr

* out_OK  -- Set true when line should be written out (this for duplicates
*            that differ)  
* first_call -- Set false after first call to WSTNFO

      logical out_OK, first_call
         
* ccom    -- Comment character for a line
     
      character*(1) ccom   
           
* anth antn ante -- antenna offsets up, north, east

      real*8 anth,antn,ante
      real*8 antdaz  ! Antenna aligment from True N (deg).

* comments for header

      integer*4 maxcmt,ncmt
      parameter (maxcmt=1)
      character*132 comments(maxcmt)  

*  user name and date
        
      character*16 uname
      integer*4 irunt(6),ihnsec      

*  buffer for report_stat call

      character*256 message

      logical debug/.false./

      
**** Report what we are doing
      
      write(message,'(3a,i5,a)') 'Output file ',out_stinf(1:20)
     .   ,' has ', num_str,' entries' 
      call report_stat('STATUS','MSTINF','htoglb/mstinf',' '
     .                 ,message,0)
 

**** Get the user name and date 
      call getdat(irunt(1),irunt(2),irunt(3))
      call gettim(irunt(4),irunt(5),irunt(6),ihnsec )
      call getusr(uname)
        
        
**** Construct the header for the new file
      
      write(unit,'(a,a16,a,i4,a,i2.2,a,i2.2,2x,i2.2,a,i2.2)'
     .     ,iostat=ioerr) 
     . '# Station.info written by MSTINF user '
     .   ,uname,' on ',irunt(1),"-",irunt(2),"-",irunt(3),irunt(4)
     ,  ,":",irunt(5)   
      write(unit,'(2a)') '* Reference file : '
     .     ,  ref_stinf(1:trimlen(ref_stinf))
      do i = 1, num_stf 
        write(unit,'(2a)') '* Merged station.info file : '
     .       , upd_stinf(i)(1:trimlen(upd_stinf(i))) 
      enddo 
      do i = 1, num_igf
        write(unit,'(2a)') '* IGS log file : '
     .       , ig_file(i)(1:trimlen(ig_file(i))) 
      enddo
      do i = 1, num_sxf
        write(unit,'(2a)') '* SINEX file : '
     .       , sx_file(i)(1:trimlen(sx_file(i)))                                                
      enddo         
      if( sr_file(1:1).ne.' ') then 
        write(unit,'(a,2(1x,a))') '* GIPSY sta_ files : '
     .       , sr_file(1:trimlen(sr_file))
     .       , sa_file(1:trimlen(sa_file))
      endif

*     Now write the comments that we already had
      if( copy_comments ) then
         do i = 1, num_comments
            if( trimlen(stinf_com(i)).gt.1 .and. 
     .          stinf_com(i)(1:1).ne.'.' )
     .      write(unit,'(a)') stinf_com(i)(1:trimlen(stinf_com(i))) 
         end do
      endif  
      ncmt = 1
      comments(1) = '*'
       
CD DEBUG   
      if( debug ) then
        print *,'List in write stinf (i site_indx(i) )'
        do i = 1, num_str
          ik = abs(site_indx(i))
          write(*,'(2i4,1x,l1,1x,a4,2(1x,i4,1x,i3),1x,a20,1x,a15)')
     .      i,site_indx(i),newflag(ik),sitcod(i)
     .      ,istart(1,ik),istart(2,ik),istop(1,ik),istop(2,ik)
     .      ,rctype(ik),anttyp(ik)
        enddo    
      endif


**** Now write out the lines of the merged file
       
      first_call = .true.  
       
      do i = 1, num_str
         ik = abs(site_indx(i))             
*        See if we should comment this entry because its a duplicate
         ccom = ' '
         out_OK = .true.      
         if( site_indx(i).lt.0 ) then

*            See if should print the line (if it is the same then don't print)  
*  rwk 090224: by-pass this code, checks done in merge_stinf
c            out_OK = .false.
c            il = abs(site_indx(i+1)) 
c            do j = 1, 3
c               if( dUNE(j,ik).ne.dUNE(j,il) ) out_OK = .true.
c            end do
c            if( rcvcod(ik).ne.rcvcod(il) ) out_OK = .true.
c            if( antcod(ik).ne. antcod(il) ) out_OK = .true.
c            if( htcod(ik).ne. htcod(il) ) out_OK = .true.
c            if( sname(ik).ne. sname(ik) ) out_OK = .true.

             ccom = '-'   
         end if

*        Now write out the line  
         
         if( ccom.eq.' ' .or. .not.nowrite ) then
           anth = dUNE(1,ik)   
           antn = dUNE(2,ik)
           ante = dUNE(3,ik) 
           antdaz = dAntAZ(ik)        
cd           print *,'calling wstnfo ref_nlist ik anttyp '
cd     .          ,ref_nlist,ik,anttyp(ik)
cd         print *,'before wstnfo swver radome :',swver(ik),radome(ik)
           call wstnfo( unit, ref_nlist, ref_item_list
     .       , first_call, maxcmt, ncmt, comments
     .       , ccom, sitcod(i), sname(ik), anth, antn, ante, antdaz
     .       , rcvcod(ik), antcod(ik), htcod(ik), radome(ik), swver(ik)
     .       , rctype(ik), rcvrsn(ik), rcvers(ik), anttyp(ik), antsn(ik)  
     .       , sn(ik), istart(1,ik), istop(1,ik), comment(ik) )   
           first_call = .false.
         endif

*        See if first character changes.  Leave blank if it does
         if( i.lt. num_str ) then
             if( sitcod(i)(1:1).ne.sitcod(i+1)(1:1) ) then
                 write(unit,'(a)',iostat=ierr) '.'
             end if
         end if
      end do

****  Thats all 
      return
      end

CTITLE READ_USE

      subroutine read_use( file_name )

      implicit none

*     Routine to read a list of 4-char codes to be used from the reference
*     station.info and not used from the merged files

      include 'mstinf.h'

* PASSED VARIABLES
* file_name  -- Name of file with list of sites

      character*(*) file_name

* LOCAL VARIABLES
* ierr  -- IOSTAT error
* indx  -- Pointer in string
* trimlen -- Length of string
* js      -- Station number in case of duplicates

      integer*4 ierr, indx, trimlen, js  ,i

* line  -- Line read from file
* word  -- First word extacted from line

      character*80 line
      character*5  word


****  Open the input file to get list
      open(101,file=file_name, status='old', iostat=ierr)  
      if(ierr.ne.0) call report_stat('FATAL','MSTINF','htoglb/mstinf'
     .     ,file_name,'Error opening station-use file',ierr)
      
*     Now read list from file     
      do while ( ierr.eq.0 )
         read(101,'(a)',iostat=ierr) line  
         if( ierr.eq.0 .and. line(1:1).eq.' ' .and. 
     .       trimlen(line).gt.0 ) then

*            Get the name from the string
             indx = 0
             call GetWord(line, word, indx)
             call casefold(word)

*            See if we already have this site name 
             indx = 0
             call get_cmd(word, use_names, num_use, js, indx)
             if( js.lt.0 ) then
                 num_use = num_use + 1
                 call check_max(num_use, max_site, 
     .                   'Number of use sites')
                 use_names(num_use) = word(1:4)
             end if
         end if
      end do    

****  Thats all
      close(101)
      return
      end

CTITLE UPD_FROM_RINEX

      Subroutine upd_from_rinex ( irxf )

      implicit none
     
*     Routine to get station.info entries from a RINEX header. 
*     rwk March 2002 from program upd_stnfo in gamit/utils. 
*     Revised significantly with addition of reading RINEX start/stop time.  R. King August 2003 

*     Revised significantly again to simplify, removing the logic to determine
*     if a new value is to be used and rather relying on merge_stinf to replace
*     the old value if necessary, remove duplicates, and adjust start 
*     and stop times.   R. King  February 2008
           
      include 'mstinf.h' 
      include '../../gamit/includes/makex.h'
              
* LOCAL VARIABLES
                                
      logical span_match,station_match,mismatch
      integer iyr, iyr2, iday, kcol, nblen, irxf, ioerr 
     .      , idoy, irxdoy0, irxdoy1,itimdif  
     .      , irxstart(5),irxstop(5),jstart(5),jstop(5)
     .      , ispan_str,i,j
      real*8 apr_epoch,coord(6)
      character*1 pcncod
      character*4 rxsite
      character*5 default_slcod 
      character*6 char6                      
      character*12 rxfname
      character*20 char20
      character*256 message
                   
* Function
      real*8 decyrs

* Variables passed to/from RINEX that are not in common
      real*4 rxver
      character*20 rxpgm,rxusr,rxdat
c     comment
      integer irxcom
      character*60 rxcom(max_comments)
c     mark name
      character*60 rxmrk
c     observer
      character*20 rxobs
c     agency
      character*40 rxagy
c     aproximate coordinates
      real*8 apx,apy,apz   
c     full antenna type (20-characters, vs 15 for station.info)
      character*20 rxanttyp 
c     number of different types of observable quantities
      integer nobtyp
c     labels for RINEX observable types
      character*3  rxobtyp(maxobt)

c     antenna offsets
      real*8 anth,ante,antn
      real*8 antdaz  ! Antenna aligment from True N (deg).
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

C      data blank20/'                    '/
    
    
       logical debug/.true./

* DUMMY HERE
      character*1 gnss/' '/
      character*3 rxtime/'   '/
             
***** Report what we're doing
                  
      write(message,'(2a)') 
     .   'Checking station.info entries for RINEX file ',rx_file(irxf) 
      call report_stat('STATUS','MSTINF','htoglb/mstinf',' '
     .                ,message,0)  
      call report_stat('WARNING','MSTINF','htoglb/mstinf',' '
     .                ,message,0)
 
*     Extract the 11-character RINEX file name, site, year, and day-of-year 
*     from the full-path RINEX file name
             
      kcol = nblen( rx_file(irxf ) )        
      rxfname = rx_file(irxf)(kcol-11:kcol-1)
      rxsite = rx_file(irxf)(kcol-11:kcol-8)
      call uppers(rxsite)

      read(rx_file(irxf)(kcol-7:kcol-5),'(i3)',iostat=ioerr) iday 
      if( ioerr.ne.0 ) call report_stat('FATAL','MSTINF'
     .     ,'htoglb/mstinf',rx_file(irxf)
     .     ,'Error reading day-of-yr from RINEX filename',ioerr)
      read(rx_file(irxf)(kcol-2:kcol-1),'(i2)',iostat=ioerr) iyr2
      if( ioerr.ne.0) call report_stat('FATAL','MSTINF','htoglb/mtinf2'
     .    ,rx_file(irxf),'Error reading year from RINEX filename',ioerr)
      iyr = iyr2
      call fix_y2k(iyr) 

*     Increment the station.info entry index (decrement later if RINEX values not used) 

      num_str = num_str + 1   
      newflag(num_str) = .true.  

*     Set the default entries    
                             
*     track code will need to be input for this to be used for kinematic  
      htcod(num_str) = 'DHARP'
      sn(num_str) = 0  
      comment(num_str)=
     .       ' mstinf: '//rx_file(irxf)(1:nblen(rx_file(irxf)))   
     

*     Read the RINEX header
          
      call rrxhed ( debug,gnss,
     .   rxver,rxpgm,rxusr,rxdat,rxcom,irxcom,rxmrk,rxobs,rxagy,
     .   rcvrsn(num_str),rctype(num_str),rcvers(num_str),
     .   antsn(num_str),rxanttyp,apx,apy,apz,
     .   anth,ante,antn,nwave1,nwave2,nobtyp,rxobtyp,rxint,rxtime,
     .   irxyr0,irxmo0,irxdy0,irxhr0,irxmn0,rxsec0,
     .   irxyr1,irxmo1,irxdy1,irxhr1,irxmn1,rxsec1) 
      antdaz = 0.0   ! Not defined in RINEX format
                              
*     Get the start and stop times for the RINEX file
   
*     If the start and stop times are on the header, use them     
      if( irxyr0.ne.0 .and. irxyr1.ne.0 ) then 
        irxdoy0 = idoy(irxyr0,irxmo0,irxdy0) 
        irxstart(1) = irxyr0
        irxstart(2) = irxdoy0
        irxstart(3) = irxhr0
        irxstart(4) = irxmn0
        irxstart(5) = int(rxsec0)  
        irxdoy1 = idoy(irxyr1,irxmo1,irxdy1) 
        irxstop(1) = irxyr1
        irxstop(2) = irxdoy1
        irxstop(3) = irxhr1
        irxstop(4) = irxmn1
        irxstop(5) = int(rxsec1)   
*     otherwise, scan the full RINEX file for the times
      else
       call get_rinex_times(rx_file(irxf),rxver,nobtyp,irxstart,irxstop)    
      endif  

*     Check for a bad start date
      if( irxstart(1).ne.iyr .or. irxstart(2).ne.iday ) then
        write(message,'(a,i4,1x,i3,a,i4,1x,i3)') 'RINEX start date ('
     .       ,irxstart(1),irxstart(2),') does not match file name date '
     .        ,iyr,iday
        call report_stat('WARNING','MSTINF','htoglb/mstinf'
     .                   ,rx_file(irxf),message,0)
      endif   
            
*     If the RINEX times are valid, use them for station.info; otherwise use the 
*     beginning and end of the day     

      if( irxstart(1).gt.1980 .and. irxstart(1).lt.2100 ) then
        do i=1,5
         istart(i,num_str) = irxstart(i)
        enddo  
      else
        call report_stat('WARNING','MSTINF','htoglb/mstinf'
     .        ,rx_file(irxf)
     .        ,'RINEX start time bogus, set = start of day',0) 
        istart(1,num_str) = iyr
        istart(2,num_str) = iday
        do i=3,5
          istart(i,num_str) = 0
        enddo
      endif
      if( irxstop(1).gt.1980 .and. irxstart(1).lt.2100 ) then
        do i=1,5
         istop(i,num_str) = irxstop(i)
        enddo       
      else
        call report_stat('WARNING','MSTINF','htoglb/mstinf'
     .        ,rx_file(irxf)
     .        ,'RINEX stop time bogus, set = end of day',0) 
        istop(1,num_str) = iyr
        istop(2,num_str) = iday 
        istop(3,num_str) = 24
        do i=4,5
          istop(i,num_str) = 0
        enddo
      endif
      
*     If specified in the command file, set the stop time open-ended
          
      if( rx_open ) then
        istop(1,num_str) = 9999
        istop(2,num_str) = 999 
        do i=3,5
          istop(i,num_str) = 0
        enddo  
      endif
   
*     Guess the receiver type and firmware version and antenna type

*       --use full 20 characters of antenna type for maximum flexibility  
      call guess_rcvant( 1,rxfname,rctype(num_str), rxanttyp
     .                 , radome(num_str), rcvers(num_str)
     .                 , rcvcod(num_str), swver(num_str)
     .                 , antcod(num_str), comment(num_str)  ) 
*       --now reassign according to new IGS standards 
c      print *,'DEBUG UPDATE_STINF rctype:',rctype(num_str)
c      print *,'                   rcvcod:',rcvcod(num_str) 
c      print *,'                   rcvers:',rcvers(num_str)
c      print *,'                    swver:',swver(num_str) 
*       --get standard IGS values from GAMIT codes
      call read_rcvant(1,2,char6,char20,radome(num_str),rcvcod(num_str)
     .                ,rctype(num_str),pcncod ) 
      call read_rcvant(1,1,antcod(num_str),rxanttyp,radome(num_str)
     .                ,char6,char20,pcncod)
c      print *,'after read_rcvant rctype:',rctype(num_str)
c      print *,'after read_rcvant rxanttyp:',rxanttyp
      anttyp(num_str) = rxanttyp(1:15)  
      radome(num_str) = rxanttyp(17:20) 

      if( radome(num_str)(2:5).eq.'    ' ) radome(num_str)='UNKN '
*     remove ^ added by guess_rcvant for firmware version
      call sub_char ( rcvers(num_str),'^',' ')   
c      print *,'after sub_char rcvers:',rcvers(num_str)
                       
*    Assign the other values

      sitcod(num_str) = rxsite 
      sname(num_str) = rxmrk(1:16)
      dUNE(1,num_str) = anth
      dUNE(2,num_str) = antn
      dUNE(3,num_str) = ante 
      dAntAZ( num_str) = antdaz
      htcod(num_str) = 'DHARP'
*     see if we need to modify the height type code
      if( max_slant.gt.0 .and. dUNE(1,num_str).gt.max_slant ) then
         htcod(num_str) = 'SLBGP'
         if(index(antcod(num_str),'DM').gt.0) htcod(num_str)='SLBCR'
*        If default has been passed then use this.
         if( default_slcod(1:1).ne.' ' ) 
     .       htcod(num_str) = default_slcod
         comment(num_str) = '?'
      endif 
       
*     Write an entry to 'mstinf.apr' if requested
      if( write_apr ) then
         apr_epoch = decyrs(iyr,iday,0.d0)
         write(uapr,'(2x,2a4,3(2x,f13.4),3(2x,f8.4),2x,f8.3)')
     .         rxsite,'_GPS', apx,apy,apz,0.,0.,0.,apr_epoch
      endif
          
*     Update the site_indx array used for sorting
      site_indx(num_str) = num_str 
                    
      return
      end
   
 
CTITLE GET_RINEX_TIMES

      Subroutine get_rinex_times(rx_file,rxver,nobtyp,irxstart,irxstop)    

      implicit none

*     Original code extracted mostly from gamit/lib/rrinex.f, but with 
*     data records skipped and different error handling; changed by RWK
*     090608 to call rrinex since the reading got more complicated.
      
      include '../../gamit/includes/makex.h'

      integer*4 irxstart(5),irxstop(5),iwkn,nobtyp,nprn,nrxerr,itflag,i
      real*4 rxver
      real*8 sow,sec    
      real*8 anth, ante, antn                
      character*128 rx_file
      character*256 message              

      logical first_epoch,fend,ferr,debug    
      data debug/.false./
               
*     RINEX quantities returned from rrinex but not used here
      character*1 asvid(maxchn)
      integer*4 isvid(maxchn),illi(maxchn,maxobt),issi(maxchn,maxobt)
     .        , nepoch,iflag
      character*3  rxobtyp(maxobt)
      real*8  dofl1(maxchn), dofl2(maxchn)
      real*8  prgl1(maxchn), prgl2(maxchn)
      real*8 utcoff
      
*     DUMMY
      character*1 gnss/' '/
      character*3 rxtime/'   '/
      integer*4 iobtypx(6)/6*0/
      

       
*     Initialization

      nrxerr = 0          
      first_epoch = .true.  
      do i=1,5
       irxstart(i) = 0
       irxstop(i) = 0
      enddo

*     Read through the file, saving the first and last epochs

      do while ( .not.fend )   
                     
c*        call rrinex
c*    .      ( debug,iflag,fend,ferr,nprn,asvid,isvid,iwkn,sow,rxobtp
c*     .      ,nobtyp,dofl1,dofl2,prgl1,prgl2,illi,issi,anth,ante,antn
c*     .      , nepoch) 
        call rrinex( debug,iflag,rxver,gnss,nobtyp,rxobtyp,iobtypx
     .              , nprn,isvid,rxtime,iwkn,sow,nepoch
     .              , dofl1,dofl2,prgl1,prgl2,illi,issi
     .              , anth,ante,antn
     .              , fend,ferr )  
        if( fend ) then
          continue
        elseif( ferr ) then
          nrxerr = nrxerr + 1          
        elseif( first_epoch )  then  
          itflag = 4    
          call timcon( itflag,iwkn,sow,irxstart(1),irxstart(2)
     .               , irxstart(3),irxstart(4),sec,utcoff)
          irxstart(5) = int(sec)
          first_epoch = .false.  
        elseif( anth.ne.0.d0 .or. ante.ne.0.d0 .or. antn.ne.0d0 ) then
          call report_stat('WARNING','MSTINF','htoglb/mstinf',rx_file
     .      ,'Antenna offset change, fix station.info manually',0 ) 
        endif
                   
      enddo

*     Save last epoch read
                 
      itflag = 4    
      call timcon( itflag,iwkn,sow,irxstop(1),irxstop(2),irxstop(3)
     .            , irxstop(4),sec,utcoff)
      irxstop(5) = int(sec)
               
*     Issue a warning if errors encountered

      if( nrxerr.gt.0 )  then
        write(message,'(i6,a)') nrxerr,' read errors on RINEX file'
        call report_stat('WARNING','MSTINF','htoglb/mstinf',rx_file
     .                  ,message,0 ) 
      endif


      return
      end


CTITLE GUESS_RCVANT

      Subroutine guess_rcvant( icall,fname,rctype,anttyp,radome,rcvers
     .                       , rcvcod,swver,antcod,comment )    

      implicit none
       
*     Get the correct GAMIT receiver, antenna, and firmware codes for an input 
*     string extracted from a RINEX, IGS, or SINEX file, using user-provided
*     entries in file guess_rcvant.dat.  Then (outlide of this routine) get the 
*     corresponding IGS codes.

*     Two modes:
 
*        icall = 1 :  Get rcvr, antenna, and firware codes (RINEX)
*              = 2 :  Get only the swver (IGS log or SINEX)

*      fname     :  Name of original file (for report_stat calls)

*      Input strings from origin file:  rctype, anttyp, radome, rcvers

*      Output GAMIT codes:  rcvcod,swver,antcod

*      Comments to add to station.info:  comment

                                         
      integer*2 idum2(18)
      integer*4 ltab,ioerr,mchkey,nrcv,nant,nrad,nswv,indx,len
     .        , lenf,icall,i
             
      integer maxcod           
      parameter(maxcod=1024)
                     
*     Function
      integer*4 nblen

      real*4 swver
                   
      character*1 pcncod 
      character*3 keywrd,lowerc  
      character*5 swvcod,swvcods(maxcod),radome,radomes(maxcod)
     .          , radstrs(maxcod),default_radome,default_slcod      
     .          , char5
      character*6 rcvcod, antcod, rcvcods(maxcod),antcods(maxcod)
     .          , code, default_rccod, default_ancod,char6
      character*20 rctype,  anttyp, rcvers, string,char20 
     .           , rcvstrs(maxcod),antstrs(maxcod),swvstrs(maxcod)   
      character*132 comment 
c     name of original file can be 11 characters (RINEX) or fewer (IGS log or SINEX)            
      character*(*) fname
      character*256 line,message

      logical eof,found

*     Initialize the defaults

      default_rccod = ' ' 
      default_ancod = ' '
      default_radome = 'NONE '  
      default_slcod = ' ' 
      swver = 0.0 
   
*     Set the length of the file name for warnings     
      lenf = nblen(fname)

*     Open the file of correspondences
                   
      ltab = 91
      open (unit   =  ltab,
     .      file   =  'guess_rcvant.dat',
     .      status =  'old',
     .      iostat =  ioerr)
      if (ioerr .ne. 0) call report_stat('FATAL','MSTINF'
     .   ,'htoglb/mstinf',fname(1:lenf),'Error opening guess_rcvant.dat'
     .   ,ioerr) 
       
*     Load the guess_rcvant.dat entries into arrays

      nrcv = 0
      nant = 0
      nswv = 0  
      nrad = 0 
      eof = .false.              
      do while (.not.eof ) 
        read( ltab,'(a)',iostat=ioerr,end=10 ) line
        if( ioerr.ne.0 ) 
     .     call report_stat('FATAL','MSTINF'
     .                    ,'htoglb/mstinf',fname(1:lenf)
     .                ,'Error reading guess_rcvant.dat for file ',ioerr)
        if( line(1:1).eq.' ' ) then
          indx = 2
          call read_line(line,indx,'CH',ioerr,idum2,keywrd)
          call read_line(line,indx,'CH',ioerr,idum2,string)
          call read_line(line,indx,'CH',ioerr,idum2,code)
          if( ioerr.gt.0 ) call report_stat('FATAL','MSTINF'
     .                        ,'htoglb/mstinf',fname(1:lenf)
     .                ,'Error reading line from _rcvant.dat',ioerr)
          if( lowerc(keywrd).eq.'rcv' ) then
            nrcv = nrcv + 1
            if( nrcv.gt.maxcod ) call report_stat('FATAL','MSTINF'
     .     ,'htoglb/mstinf','rcvant.dat','Too many rcv codes ',0)
            rcvstrs(nrcv) = string
            rcvcods(nrcv) = code
*           See if default given
            if( string(1:7).eq.'default' ) default_rccod = code

          elseif( lowerc(keywrd).eq.'ant' ) then
            nant = nant + 1 
            if( nant.gt.maxcod ) call report_stat('FATAL','MSTINF'
     .      ,'htoglb/mstinf','rcvant.dat','Too many ant codes ',0)
            antstrs(nant) = string
            antcods(nant) = code
*           See if default given
            if( string(1:7).eq.'default' ) default_ancod = code
                           
          elseif( lowerc(keywrd).eq.'rad' ) then
            nrad = nrad + 1 
* MOD TAH 210218: Removed extra ' ' from calling string.
            if( nrad.gt.maxcod ) call report_stat('FATAL','MSTINF'
     .          ,'htoglb/mstinf','rcvant.dat'
     .          ,'Too many ant codes ',0)
            radstrs(nrad) = string
            radomes(nrad) = code
*           See if default given
            if( string(1:7).eq.'default' ) default_radome = code
                                 
          else if( lowerc(keywrd).eq.'slt' ) then
             default_slcod = code

          elseif( lowerc(keywrd).eq.'swv' ) then
            nswv = nswv + 1     
            if( nrcv.gt.maxcod ) call report_stat('FATAL','MSTINF'
     .     ,'htoglb/mstinf','rcvant.dat','Too many swv codes ',0)
            swvstrs(nswv) = string   
            swvcods(nswv) = code   
          endif                               
        endif
      enddo                                                             
  10  if( nrcv.eq.0 .or. nant.eq.0 .or.nswv.eq.0 ) 
     .     call report_stat('FATAL','MSTINF','htoglb/mstinf',' ' 
     .    ,'Zero entries in _rcvant.dat table',0)     
      close(ltab)      

      
*     Get the receiver code   (skip if IGS log or SINEX file)
                            
      if( icall.eq.1 ) then   
        rcvcod = 'unknwn'   
*       see if RINEX entry is a correct IGS name      
        call read_rcvant(2,2,char6,char20,char5,rcvcod,rctype,pcncod)
        if( rcvcod.eq.'      ') then
          write(message,'(a,a6,a)') 'RINEX rcvr name ',rctype,
     .        ' non-standard, read guess_rcvant.dat '
          call report_stat('WARNING','MSTINF','htoglb/mstinf'
     .      ,fname(1:lenf),message,0 ) 
*         must substitute ^ for blank in RINEX string (and use ^ in table) to get uniqueness
          call sub_char ( rctype,' ','^')   
cd        print *,'rctype ',rctype
          do i=1,nrcv
            len = nblen(rcvstrs(i))
            if ( mchkey(rctype,rcvstrs(i),20,len).gt.0 ) 
     .         rcvcod = rcvcods(i)    
cd           print *,'len rcvstrs rcvcod',len,rcvstrs(i),rcvcod
          enddo  
          if( rcvcod.eq.'unknwn' ) then
*           see if a default type has been given by user
            if( default_rccod(1:1).ne.' ' ) then
              rcvcod = default_rccod
              write(message,'(a,a6,a,a20)')  'Using default rcvcod ( ',
     .                               rcvcod,' ) for RINEX entry ',rctype
              call report_stat('WARNING','MSTINF','htoglb/mstinf'
     .                       ,fname(1:lenf),message,0)
              comment(nblen(comment)+1:) = '?'        
            else
             write(message,'(a,a20)') 
     .               'Cannot get rcvcod from RINEX entry: ',rctype
             call report_stat('FATAL','MSTINF','htoglb/mstinf'
     .                    ,fname(1:lenf),message,0)
            endif
          endif
        endif
      endif  

*     Get the antenna code
         
      if( icall.eq.1 ) then
        antcod = 'unknwn'   
        if( anttyp(17:20).eq.'    ' ) then
          radome = 'UNKN ' 
        else
          radome(1:4) = anttyp(17:20)
        endif
*       see if RINEX entry is a correct IGS name 
        call read_rcvant(2,1,antcod,anttyp,radome,char6,char20,pcncod) 
        if( antcod.eq.'      ') then
          write(message,'(a,a6,a)') 'RINEX antenna name ',anttyp,
     .        ' non-standard, read guess_rcvant.dat '
          call report_stat('WARNING','MSTINF','htoglb/mstinf'
     .      ,fname(1:lenf),message,0 ) 
*         must substitute ^ for blank in RINEX string (and use ^ in table) to get uniqueness
          call sub_char ( anttyp,' ','^')
          do i=1,nant
            len = nblen(antstrs(i))
            if(mchkey(anttyp,antstrs(i),20,len).gt.0) antcod=antcods(i)
          enddo   
          if( antcod.eq.'unknwn' ) then
*           see if a default type has been given
            if( default_ancod(1:1).ne.' ' ) then
              antcod = default_ancod           
              write(message,'(a,a5,a,a20)') 'Using default antcod ( '
     .                            ,antcod,' ) for RINEX entry ',anttyp
              call report_stat('WARNING','MSTINF','htoglb/mstinf'
     .                    ,fname(1:lenf),message,0)        
              comment(nblen(comment)+1:) = '?'
            else
              write(message,'(a,a20)') 
     .           'Cannot get antcod from RINEX entry: ',anttyp
              call report_stat('FATAL','MSTINF','htoglb/mstinf'
     .                      ,fname(1:lenf),message,0)
            endif
          endif 
*         Radomes 
          do i=1,nrad
            len = nblen(radstrs(i))
             if ( mchkey(anttyp(16:20),radstrs(i),20,len).gt.0 )
     .            radome = radomes(i)
           enddo 
        endif   
      endif

*     The firmware version needs to be gotten from guess_rcvant.dat in all cases
*     since the correspondence between the long version and a unique GAMIT
*     decimal version is not recorded any  (this may change with a new 
*     rcvant.dat format in the future)

*     check to see if the firmware version is a decimal number before trying to read it
      call check_num(rcvers,ioerr) 
      if( ioerr.ne.0 ) then  
        call sub_char ( rcvers,' ','^')    
        found = .false.
        do i=1,nswv
         len = nblen(swvstrs(i))                                         
         if ( mchkey(rcvers,swvstrs(i),20,len).gt.0 ) then
           swvcod=swvcods(i) 
           found = .true.   
         endif
        enddo   
        if( found ) read(swvcod,'(f5.0)',iostat=ioerr) swver 
      else                                        
*       but if the decimal number is too large for the station.info field width
*       (f5.2), set it to 0.0
        read(rcvers,'(f20.0)',iostat=ioerr ) swver
        if( swver.gt. 99.) swver = 0.0
        write(message,'(a,f9.0,a)') 
     .     'swver (',swver,') from rcvers > 99., set to zero'
        call report_stat('WARNING','MSTINF','htoglb/mstinf'
     .                  ,fname(1:lenf),message,0)
      endif
      if ( swver.eq.0. ) then
        write(message,'(a,a20,a)') 
     .    'Firmware version ',rcvers
     .      , 'not in guess_rcvant.dat; set to 0.0'
        call report_stat('WARNING','MSTINF','htoglb/mstinf'
     .                  ,fname(1:lenf),message,0)
        comment(nblen(comment)+1:) = '?' 
      endif
      return
      end


CTITLE UPD_FROM_IGSLOG

      Subroutine upd_from_igslog ( uigs,ilogf ) 

      implicit none         
     
*     Routine to get station.info entries from an IGS log file. 
*     rwk Sept 2007, Feb 2007 from surboutine update_stinf (now upd_from_rinex). 
                                                              
*       iuigs is the unit number, set in the main program
*       ilogf  is the index in the list of IGS log files requested

*     Routine adds all log entries to the global array, relying on 
*     merge_stinf to replace old entries if requested
   
* COMMON VARIABLES

      include 'mstinf.h' 
      include '../../gamit/includes/makex.h'  
              
* LOCAL VARIABLES
              
      integer*4 uigs,ilogf,len
      integer rcvr_start(5,max_log),rcvr_stop(5,max_log)
     .      , ant_start(5,max_log),ant_stop(5,max_log)
     .      , mstart(5,max_log),mstop(5,max_log)
     .      , irpntr(max_log),iapntr(max_log),ioerr
     .      , num_save,nr,na,ir,ia,ik,ikp,i,j

      integer*4 trimlen ! Length of non-blank string
                              
      real*4 section,last_section
      real*8 anth1(max_log),antn1(max_log),ante1(max_log),sod,tk,tkp
     .     , decyrs,dUNE1(3,max_log),coord(6),coord_epoch
* MOD TAH 200203: Added AntDAZ to list of values from station.info
      real*8 dAntAZ1(max_log)  ! Antenna aligment from True N (deg).
     .,      antdaz1(max_log)  ! Copy (for consistency between dUNE and ant[hne]1.
                  
      character*1 pcncod
      character*4 sitcod1,sitcodlc,lowerc                   
      character*6 char6
      character*5 radome1(max_log),char5   
* MOD TAH 200213: Changed from  character*10 to character*20 for consistency
      character*20 rcvrsn1(max_log)
      character*15 anttyp1(max_log)
      character*16 sname1,snametmp
      character*17 rstart,rstop,astart,astop
      character*20 rctype1(max_log),rcvers1(max_log)
     .            ,antsn1(max_log),anttyp20,char20
      character*24 char24
      character*80 line
      character*256 message
     
      logical found,finished,debug,new_section
                          
      integer*4 mchkey 

      data debug/.false./
    
              
* FUNCTION

      integer*4 itimdif,nblen

*      Report what we're doing
                  
      write(message,'(2a)') 
     . 'Updating station.info entries from IGS log file ',ig_file(ilogf)
      call report_stat('STATUS','MSTINF','htoglb/mstinf',' '
     .                ,message,0)  
      call report_stat('WARNING','MSTINF','htoglb/mstinf',' '
     .                ,message,0)

*     Initialize the times for facilate checking
      do i=1,max_log
        rcvr_start(1,i) = 9999
        rcvr_start(2,i) = 999
        rcvr_start(3,i) = 0
        rcvr_start(4,i) = 0
        rcvr_start(5,i) = 0
        ant_start(1,i) = 9999
        ant_start(2,i) = 999
        ant_start(3,i) = 0
        ant_start(4,i) = 0
        ant_start(5,i) = 0
      enddo                
      do i=1,6
        coord(i) = 0.d0 
      enddo  
      coord_epoch = 2008.0     
                                                 
      
*     Read the site description 
                 
      found = .false.
      do while (.not.found ) 
        read(uigs,'(a)',iostat=ioerr) line   
* MOD TAH 180802: Added dump Line
        if( dump_log ) write(*,'(a," | ",a)') 'DUMP DES ',trim(line)   
        if(debug) print *,'read site description line ',line
        if( ioerr.ne.0 ) call log_error(ig_file(ilogf)
     .     ,'site description',' ',0,ioerr,'F') 
        if( line(6:14).eq.'Site Name' ) then
           read(line(33:48),'(a16)',iostat=ioerr) snametmp 
          sname1 = snametmp
          if( ioerr.ne.0 ) then
             call report_stat('FATAL','MSTINF','htoglb/mstinf'
     .         , ig_file(ilogf)
     .         , 'Error reading full site name from IGS log file',ioerr) 
          else
            found = .true.
          endif
        endif             
      enddo   
       
*     Find the 4-character site code and check it against the file name

      found = .false.
      do while (.not.found ) 
        read(uigs,'(a)',iostat=ioerr) line                  
* MOD TAH 180802: Added dump Line
        if( dump_log ) write(*,'(a," | ",a)') 'DUMP CID ',trim(line)   
        if(debug) print *,'read character ID line ',line
        if( ioerr.ne.0 ) call log_error(ig_file(ilogf)
     .     ,'Four Character ID line'
     .      ,' ',0,ioerr,'F') 
        if( line(6:14).eq.'Four Char' ) then
           read(line(33:36),'(a4)',iostat=ioerr) sitcod1
          if( ioerr.ne.0 ) then
             call log_error(ig_file(ilogf),'site code ',' ',0,ioerr,'F')
           else
             found = .true. 
             call uppers(sitcod1)   
           endif
        endif             
      enddo   
*     The file name may be a full path, so search for the 4 characters anywhere in the string 
      sitcodlc = lowerc(sitcod1)
      if( mchkey(ig_file(ilogf),sitcod1,128,4).eq.0 .and.
     .    mchkey(ig_file(ilogf),sitcodlc,128,4).eq.0 ) then
        write(message,'(a,a4,a)') '4-character code in log file ('
     .    ,sitcod1,') does not match file name)'
        call report_stat('FATAL','MSTINF','htoglb/mstinf'
     .      ,ig_file(ilogf),message,0)
      endif

*     Read all of the receiver and antenna info to see how many entries 
*     of each we have (unlike RINEX files, receiver and antenna info 
*     is not paired)
                   
*     Skip sections 1 
      finished = .false.
      do while (.not.finished)
        read(uigs,'(a)',iostat=ioerr) line   
* MOD TAH 180802: Added dump Line
        if( dump_log ) write(*,'(a," | ",a)') 'DUMP CRD ',trim(line)   
        if(debug) print *,'reading finding site coordr section ',line
        if( ioerr.ne.0 )  call log_error(ig_file(ilogf)
     .     ,'finding rcvr section',' ',0,ioerr,'F') 
        if( line(1:1).eq.'2') finished = .true.
      enddo

*     Site coordinate entries (Section 2)     
      finished = .false.                       
      last_section = 2.0
      nr=0
      do while(.not.finished)
        found = .false.
        do while (.not.found)     
          read(uigs,'(a)',iostat=ioerr) line  
* MOD TAH 180802: Added dump Line
        if( dump_log ) write(*,'(a," | ",a)') 'DUMP CRS ',trim(line)   
          if(debug) print *,'reading coords lines ',line
          if( ioerr.eq.-1) then
            finished = .true.   
            found = .true.
          elseif( ioerr.ne.0 ) then
            call log_error(ig_file(ilogf),'line for coordinates'
     .       ,' ',0,ioerr,'F')    
          elseif(  new_section(line,last_section) ) then
            finished = .true.  
            found = .true.
          elseif(line(8:20).eq.'X coordinate') then
            found =.true. 
            nr = nr + 1       
            read(line,'(32x,f13.5)',iostat=ioerr ) coord(1)       
            if( ioerr.ne.0) call log_error(ig_file(ilogf),'X coord'
     .                                  ,' ',nr,ioerr,'W')                                       
            read(uigs,'(a)',iostat=ioerr) line  
            read(line,'(32x,f13.5)',iostat=ioerr ) coord(2)       
            if( ioerr.ne.0) call log_error(ig_file(ilogf),'Y coord'
     .                                  ,' ',nr,ioerr,'W')   
            read(uigs,'(a)',iostat=ioerr) line              
            read(line,'(32x,f13.5)',iostat=ioerr ) coord(3)       
            if( ioerr.ne.0) call log_error(ig_file(ilogf),'Z coord'
     .                                  ,' ',nr,ioerr,'W')   
            if( write_apr) 
     .         write(uapr,'(2x,2a4,3(2x,f13.4),3(2x,f8.4),2x,f8.3)')
     .         sitcod1,'_GPS', (coord(i),i=1,6),2008.00
           endif
         enddo
      enddo
                        
*     Receiver entries (Section 3)                             
      finished = .false.                       
      last_section = 3.0
      nr=0
      do while(.not.finished)
        found = .false.
        do while (.not.found)     
          read(uigs,'(a)',iostat=ioerr) line     
* MOD TAH 180802: Added dump Line
        if( dump_log ) write(*,'(a," | ",a)') 'DUMP REC ',trim(line)   
          if(debug) print *,'reading rcvr line ',line
          if( ioerr.eq.-1) then
            finished = .true.   
            found = .true.
          elseif( ioerr.ne.0 ) then
            call log_error(ig_file(ilogf),'line for rcvr info'
     .       ,' ',0,ioerr,'F')    
          elseif(  new_section(line,last_section) ) then
            finished = .true.  
            found = .true.
          elseif(line(1:1).eq.'3'.and.line(3:3).ne.'x'.and.
     .          line(3:3).ne.' '.and.line(3:3).ne.'0' .and.
     .          line(33:36).ne.'(A20' ) then
            found =.true. 
            nr = nr + 1                                 
            if( nr.gt.max_log ) call report_stat('FATAL','MSTINF'
     .       ,'htoglb/mstinf',ig_file(ilogf)      
     .       ,'Number of receiver entries exceeds MAXLOG',0) 
            backspace( uigs )  
            rctype1(nr) = ' ' 
            call read_igs_rcvr( uigs,nr,ig_file(ilogf)
     .                        , rctype1(nr),rcvrsn1(nr),rcvers1(nr)
     .                        , rcvr_start(1,nr),rcvr_stop(1,nr) )   
           endif
         enddo
      enddo    
      if( nr.le.0 ) call report_stat('FATAL','MSTINF','htoglb/mstinf'
     .   ,ig_file(ilogf),'Number of receiver entries = 0',0) 
                 
*     Antenna entries (Section 4)
      finished = .false. 
      last_section = 4.0
      na=0
      do while(.not.finished)
        found = .false.
        do while (.not.found)     
          read(uigs,'(a)',iostat=ioerr) line   
* MOD TAH 180802: Added dump Line
        if( dump_log ) write(*,'(a," | ",a)') 'DUMP ANT ',trim(line)   
          if( debug ) write(*,'(2a,i5)') 'Read antenna line: '
     .                   ,line ,ioerr
          if( ioerr.eq.-1) then
            finished = .true.  
            found = .true.
          elseif( ioerr.ne.0 ) then  
            call log_error(ig_file(ilogf),'line for antenna info'
     .      ,' ',0,ioerr,'F') 
cxx          elseif(line(1:1).eq.'5'.or.line(3:3).eq.'x') then
cxx            elseif( section_number.gt.4 .or. line(3:3).eq.'x') then 
          elseif( new_section(line,last_section) ) then
            finished = .true.
            found = .true.          
            if(debug) write(*,'(2a,f6.2,l1)') 'line last_section finish'
     .                             ,  line,last_section,finished
          elseif(line(1:1).eq.'4'.and.line(3:3).ne.'x'.and.
     .           line(3:3).ne.' '.and.line(3:3).ne.'0'.and.
     .           line(33:36).ne.'(A20' ) then
            found =.true.
            na = na + 1
            backspace (uigs)   
            anttyp1(na) = ' '
* MOD TAH 200203: Added AntDAZ to list of values from station.info
            call read_igs_ant( uigs,na,ig_file(ilogf) 
     .                       , anttyp1(na),antsn1(na),radome1(na)
     .                       , anth1(na),antn1(na),ante1(na)
     .                       , antdaz1(na)
     .                       , ant_start(1,na),ant_stop(1,na) ) 
            dUNE1(1,na) = anth1(na)
            dUNE1(2,na) = antn1(na)
            dUNE1(3,na) = ante1(na)
            dAntAZ1(na)  = antdaz1(na)
          endif
        enddo
      enddo     
      if( na.le.0 ) call report_stat('FATAL','MSTINF','htoglb/mstinf'
     .   ,ig_file(ilogf),'Number of antenna entries = 0',0) 

 
      if( debug .or. dump_log ) then
        write(*,'(a)') 'Receiver entries: '
        do i=1,nr          
          write(*,'(i3,2(2i5,3i4),1x,3a20)') i,
     .       (rcvr_start(j,i),j=1,5),(rcvr_stop(j,i),j=1,5)
     .       ,rctype1(i),rcvrsn1(i),rcvers1(i)
        enddo
        write(*,'(a)') 'Antenna entries '
        do i=1,na
          write(*,'(i3,2(2i5,3i4),1x,a15,a20,a5,3f8.3,1x,f5.0)') i
     .     ,(ant_start(j,i),j=1,5),(ant_stop(j,i),j=1,5)
     .    ,anttyp1(i),antsn1(i),radome1(i),anth1(i),antn1(i),ante1(i),
     .     antdaz1(i)
        enddo
      endif                                  

*     Now pair up the receiver and antenna entries and insert them at the 
*     end of the global arrays

* MOD TAH 161118: Make sure we have a long name
      if( trimlen(sname1).eq.0 ) sname1 = sitcod1 // ' GPS Site'
           
* MOD TAH 200203: Added AntDAZ to list of values from station.info
      call merge_rcv_ant( ig_file(ilogf),nr,na,sitcod1,sname1
     .      , rcvr_start,rcvr_stop,rctype1,rcvrsn1,rcvers1
     .      , ant_start,ant_stop,anttyp1,antsn1,radome1,dUNE1,
     .        dAntAZ1 )

      return
      end                        

                    
CTITLE  FUNCTION NEW_SECTION

      logical function new_section( line,last_section ) 

      implicit none

      character*(*) line
      integer*4 ioerr
      real*4 last_section,section
      logical debug/.false./

*     Function returns true if the integer part of the section number is not
*     the same as the last section or if the decimal part is 'x'; either is
*     a sign to UPD_FROM_IGSLOG that this section is finished.  In this case
*     the routine also updates last_section.  Reading 0.0 in this field simply
*     means it is blank and further tests should be skipped.

                            
      if( debug ) print *,'NEW_SECTION ',last_section,line
      new_section = .false.
      if( line(1:1).eq.' ' ) then
        continue
      elseif( line(3:3).eq.'x'.or.line(4:4).eq.'x') then
        new_section = .true.
      else
        read(line(1:4),'(f4.0)',iostat=ioerr ) section
        if( ioerr.ne.0 ) then           
*         a non-integer value is read but not caught by the first two tests 
          new_section = .true.   
          last_section = section     
        elseif ( int(section).ne.int(last_section) ) then
          new_section = .true.
          last_section = section    
        endif     
      endif
      if( debug )  print *,'ioerr new_section ',ioerr,new_section
      return 
      end


CTITLE READ_IGS_RCVR

      Subroutine read_igs_rcvr( uigs,nr,igfile,rctype,rcvrsn,rcvers
     .                        , rcvr_start,rcvr_stop ) 

      implicit none
              
      integer*4 uigs,nr
          
*     VALUES FOR THIS RECEIVER

      integer*4 rcvr_start(5),rcvr_stop(5)       
      character*(*) igfile,rctype,rcvrsn,rcvers
      

*     LOCAL
                           
      integer*4 ioerr   
      integer*4 trimlen   ! Length of string. 
      character*17 rstart,rstop 
      character*80 line

      logical found,endr,debug
  
      data debug/.false./
                                                         
      read(uigs,'(32x,a20)',iostat=ioerr) rctype
      if( ioerr.ne.0) call log_error(igfile,'rcvr type'
     .                 ,rctype,nr,ioerr,'F') 
      if( debug ) write(*,'(a,i3,1x,a20)') 'read receiver type '
     .                                      ,nr,rctype       
      endr = .false.
      do while( .not.endr ) 
        line = ' '
        read(uigs,'(a)',iostat=ioerr) line
        if( ioerr.ne.0 ) call report_stat('FATAL','MSTINF'
     .             ,'htoglb/mstinf',igfile
     .             ,'Unexpected error reading receiver line',ioerr)
        if( debug ) write(*,'(2a)') 'READ_IGS_RCVR line: ',line
        if(line(6:11).eq.'Serial') then
          read(line,'(32x,a10)',iostat=ioerr) rcvrsn 
          if( ioerr.ne.0) call log_error(igfile,'rcvr sn'
     .                                      ,rcvrsn,nr,ioerr,'F') 
c          if( debug ) write(*,'(a,i3,1x,a10)') 'read serial number'
c     .                                   nr,rcvrsn   
           if( debug ) print *,'read serial number ',nr,rcvrsn
        elseif(line(6:13).eq.'Firmware') then
          read(line,'(32x,a20)',iostat=ioerr ) rcvers       
          if( ioerr.ne.0) call log_error(igfile,'rcvers'
     .                                  ,rcvers,nr,ioerr,'F') 
c          if( debug ) write(*,'(a,i3,1x,a20)') 'read firmware '
c     .                               ,nr,rcvers)             
           if( debug ) print *,'read firmware ',nr,rcvers
        elseif(line(12:19).eq.'nstalled') then
          read(line,'(32x,a16)',iostat=ioerr) rstart    
          if(ioerr.ne.0) call log_error(igfile,'rcvr start'
     .             ,rstart,nr,ioerr,'F') 
          if(debug) write(*,'(a,i3,1x,a16)') 'read start',nr,rstart 
          call convert_igst( 'rvr','S',rstart,rcvr_start)
        elseif(line(12:17).eq.'emoved') then
          read(line,'(32x,a16)',iostat=ioerr) rstop  
         if( ioerr.ne.0) call log_error(igfile,'rcvr stop'
     .                                      ,rstop,nr,ioerr,'F') 
          if(debug)  write(*,'(a,i3,1x,a16)') 'read stop ',nr,rstop
          call convert_igst( 'rvr','E',rstop,rcvr_stop) 

* MOD TAH 161222: Test for blank line (to avoid case where
*       and entry have been continued to next line e.g.,
* Radome Serial Number     : DOME, Outer Fiber Reinforced Plastic Dome,
*                            Half Sphere, 145cm Diameter, date: 2002-04-25
*       elseif(line(11:12).eq.' '.or.line(1:1).ne.' ' ) then
        elseif( trimlen(line).eq.0 ) then
*         endof subsection should have a blank line, but might not
          endr = .true.
        endif                 
      enddo                                                
      return
      end
           


CTITLE READ_IGS_ANT

      Subroutine read_igs_ant( uigs,na,igfile,anttyp,antsn,radome
     .         , anth,antn,ante,antdaz, ant_start,ant_stop ) 

      implicit none
              
      integer*4 uigs,na
      character*(*) igfile
       
*     VALUES FOR THIS ANTENNA
                        
      integer*4 ant_start(5),ant_stop(5)        
      real*8 anth,antn,ante 
      real*8 antdaz  ! Antenna aligment from True N (deg).
      character*5 radome
      character*15 anttyp
      character*20 antsn,anttyp1

*     LOCAL
                    
      integer*2 idum2(18)       
      integer*4 ioerr,indx 
      integer*4 trimlen   ! Length of string. 
      character*17 astart,astop
      character*80 line

      logical found,enda,debug

      data debug/.false./
                           
c     read(uigs,'(32x,a15,a5)',iostat=ioerr) anttyp,radome  
c     read(uigs,'(t33,a15,t48,a5)',iostat=ioerr) anttyp,radome   
      read(uigs,'(32x,a20)',iostat=ioerr) anttyp1
      if( debug ) print *,'read anttyp1 ',anttyp1
c     replace (illegal) tab between the antenna and radome name
      call sub_char(anttyp1,char(9),' ')
      indx = 1
      call read_line(anttyp1,indx,'CH',ioerr,idum2,anttyp)
      if( debug ) print *,'read_line indx anttyp ',indx,anttyp
      call read_line(anttyp1,indx,'CH',ioerr,idum2,radome)
      if( debug) print *,'read_line indx radome ',indx,radome
c     skip the error check for the radome since usually read later
c     if(ioerr.ne.0) call log_error(igfile,'ant type'
c    .           ,anttyp,na,ioerr,'W')  
      enda = .false.
      do while( .not.enda ) 
        line = ' '
        read(uigs,'(a)',iostat=ioerr) line                
        if( debug ) write(*,'(2a)') 'READ_IGS_ANT line: ',line
        if( ioerr.ne.0 ) call report_stat('FATAL','MSTINF'
     .             ,'htoglb/mstinf',igfile
     .             ,'Unexpected error reading antenna line',ioerr)
        
        if(line(6:11).eq.'Serial') then
          read(line,'(32x,a20)',iostat=ioerr) antsn     
          if( ioerr.ne.0) call log_error(igfile,'ant sn '
     .             ,antsn,na,ioerr,'W') 
          if( debug ) write(*,'(a,i3,1x,a20)') 'read serial number'
     .                       ,na,antsn
        elseif( line(14:19).eq.'ARP Up' ) then 
          read(line,'(32x,f8.4)',iostat=ioerr) anth 
          if( ioerr.ne.0) then 
             call log_error(igfile,'ant ht',' ',na,ioerr,'W') 
             anth = 0.
          endif
          if( debug ) write(*,'(a,i3,1x,f8.4)') 'read ant ht'
     .                       ,na,anth
        elseif( line(18:26).eq.'North Ecc' ) then   
          if( line(33:33).eq.'('.or.line(33:33).eq.'F'.or.
     .        line(33:33).eq.' ') then
              antn = 0.d0   
          else
            read(line,'(32x,f8.4)',iostat=ioerr) antn
            if( ioerr.ne.0 ) call log_error(igfile,'ARP N line '
     .                                   ,'',na,ioerr,'F') 
          endif
          if( debug ) write(*,'(2a)') 'read ant north',line
        elseif( line(18:26).eq.'North Ecc' ) then   
          if( line(33:33).eq.'('.or.line(33:33).eq.'F'.or.
     .        line(33:33).eq.' ') then
              antn = 0.d0   
          else
            read(line,'(32x,f8.4)',iostat=ioerr) antn
            if( ioerr.ne.0 ) call log_error(igfile,'ARP N line '
     .                                   ,'',na,ioerr,'F') 
          endif
          if( debug ) write(*,'(2a)') 'read ant north',line
        elseif( line(18:25).eq.'East Ecc' ) then   
          if( line(33:33).eq.'('.or.line(33:33).eq.'F'.or.
     .        line(33:33).eq.' ') then
              ante = 0.d0   
          else
            read(line,'(32x,f8.4)',iostat=ioerr) ante
            if( ioerr.ne.0 ) call log_error(igfile,'ARP E line '
     .                                   ,'',na,ioerr,'F') 
          endif
          if( debug ) write(*,'(2a)') 'read ant east ',line
* MOD TAH 200203: Read antenna azimuth
C    Alignment from True N    : 0 deg

        elseif( line(21:26).eq.'True N' ) then   
          if( line(33:33).eq.'('.or.line(33:33).eq.'F') then
              antdaz = 0.d0   
          else
* MOD TAH 200309: Casefold line because Yes/yes are often in log
            call casefold(line)
            call sub_char(line,'YES',' 0.0')
            read(line(32:),*,iostat=ioerr) antdaz
* MOD TAH 200213: Forgive -1 errors due to bad encoding
            if( ioerr.eq.-1 ) then
                call log_error(igfile,'True N line '
     .                ,'',na,ioerr,'W') 
                antdaz = 0 
            elseif( ioerr.ne.0 ) then
* MOD TAH 200309: Changed to warning because Az is not yet needed.
                call log_error(igfile,'True N line '
     .                ,line,na,ioerr,'W') 
                antdaz = 0 
            endif
         endif
          if( debug ) write(*,'(2a)') 'read antdaz ',line
        elseif( line(14:24).eq.'Radome Type' ) then
          read(line,'(32x,a5)',iostat=ioerr) radome   
          if( ioerr.ne.0 ) then
            call log_error(igfile,'Radome',radome,na,ioerr,'W') 
            radome = 'UNKN '
          endif
          if( debug ) write(*,'(a,a5)') 'read radome type ',radome
        elseif( line(6:19).eq.'Date Installed'.or.
     .          line(6:19).eq.'Date installed' ) then 
          read(line,'(32x,a16)',iostat=ioerr) astart    
          if( ioerr.ne.0 ) call log_error(igfile,'Ant start'
     .            ,astart,na,ioerr,'F') 
          if( debug ) write(*,'(a,a16)') 'read ant date start ',astart   
          call convert_igst( 'Ant','S',astart,ant_start )
        elseif( line(6:17).eq.'Date Removed'.or.
     .          line(6:17).eq.'Date removed' ) then 
          read(line,'(32x,a16)',iostat=ioerr) astop
          if( ioerr.ne.0 ) call log_error(igfile,'Ant stop'
     .            ,astop,na,ioerr,'F') 
          if( debug ) write(*,'(a,a16)') 'read ant date stop ',astop   
          call convert_igst( 'Ant','E',astop,ant_stop )

* MOD TAH 161222: Test for blank line (to avoid case where
*       and entry have been continued to next line e.g.,
* Radome Serial Number     : DOME, Outer Fiber Reinforced Plastic Dome,
*                            Half Sphere, 145cm Diameter, date: 2002-04-25
*       elseif(line(11:12).eq.' '.or.line(1:1).ne.' ' ) then
        elseif( trimlen(line).eq.0 ) then
*         end of subsection should have a blank line but might not
          enda = .true.
        endif                 
      enddo            
      if( debug ) print *,'End of secction in read_igs_ant '
      return
      end
      

CTITLE  LOG_ERROR

      Subroutine log_error(fname,readmsg,item,nentry,ioerr,W_OR_F)

      implicit none

*       File name, message identifying location in the program and item being read
      character*(*) fname,readmsg,item
                        
*       Number of entry in group
      integer*4 nentry

*       IOSTAT error code
      integer*4 ioerr

*       WARNING or FATAL
      character*1 W_or_F    

      character*256 message

      write(message,'(5a,i3)') 'Error reading ',readmsg,' Item: ',item
     .     ,' #Entry:',nentry
      if( W_or_F.eq.'W' ) then
             call report_stat('WARNING','MSTINF','htoglb/mstinf'
     .                       ,fname,message,ioerr)       
      else
             call report_stat('FATAL','MSTINF','htoglb/mstinf'
     .                       ,fname,message,ioerr)       
      endif
      end
       
                 
CTITLE UPD_FROM_GIPSY

      Subroutine upd_from_gipsy ( ustar, ustaa )

      implicit none
     
*     Routine to get station.info entries from two GIPSY sta_ files,
*     one for receivers and one for antennas 
*     rwk January 2013 from surboutine update_from_snx 
                                                              
*       ustar is the receiver-file unit number, set in the main program
*       ustaa is the antenna-file unit number, set in the main program
                   
*     Note that this program uses the comments at the end of the lines, which
*     have official IGS names.  For pure sta_ files, we'll need to add a call to
*     read a file that translates GIPSY 9-character antenna names to IGS 20-
*     character names.
*     
*     Routine adds all log entries to the global array, relying on 
*     merge_stinf to replace old entries if requested
  
* COMMON VARIABLES

      include 'mstinf.h' 
      include '../../gamit/includes/makex.h'  
              
* LOCAL VARIABLES
              
      integer*4 ustar,ustaa,rcvr_start(5,max_log)
     .        , rcvr_stop(5,max_log),ant_start(5,max_log)
     .        , ant_stop(5,max_log),ir,ia,nr,na,nrs
     .        , iyr,imo,iday,ihr,imin,ioerr,num_save,i,j

      real*4 anth1,antn1,ante1,sec
      real*8 dUNE1(3,max_log)
c MOD TAH 200203: Added AntDAZ token for antenna  Alignment from True N
      real*8 dAntAZ1(max_log)  ! Antenna Azimuth; not in GIPSY file but
                               ! need as dummy
      character*4 sitcod1(max_snx_site),sitcodr,sitcoda,site,char4
      character*5 radome1(max_log),htcod1,char5,antsn5
      character*9 goacod      
      character*15 anttyp1(max_log)
      character*16 sname1(max_snx_site)
      character*20 rctype1(max_log),rcvrsn1(max_log)
     .            ,rcvers1(max_log),rctypex,antsn1(max_log)
     .            ,igscod
      character*24 char24
      character*128 fname
      character*256 message,line,linex
     .            , rlines(max_snx_site),alines(max_snx_site)
     
      logical eof,debug,finished_rcvf,same_site
      data debug/.false./
              
* FUNCTION               
      integer*4 nblen,itimdif,idoy
      real*8 decyrs

*      Report what we're doing
                  
      write(message,'(4a)') 
     . 'Updating station.info entries from GIPSY files '
     .    ,sr_file(1:nblen(sr_file)),' ',sa_file(1:nblen(sa_file))
      call report_stat('STATUS','MSTINF','htoglb/mstinf',' '
     .                ,message,0)  
      call report_stat('WARNING','MSTINF','htoglb/mstinf',' '
     .                ,message,0)
                         
*     Set a generic file name for comments
      fname =  ' '
      fname(1:10) = 'GIPSY sta_'   
              
*     Read the receiver and and files one sites at a time. merge 
*     the two sets of entries and insert them into the global array

*        read one line of the receiver file, then find all other
*        receiver lines and all antenna lines matching the site code
                                 
      finished_rcvf = .false.
      do while(.not.finished_rcvf )
                   
        nr = 1     
        rlines(nr) = ' '  
        read(ustar,'(a)',iostat=ioerr) rlines(nr)  
        if(debug) write(*,*) '1st ioerr nr rlines ',ioerr,nr,rlines(nr)
        site = rlines(nr)(2:5)                                     
        if( ioerr.eq.-1 ) then
          finished_rcvf =. true. 
          nr = nr - 1
        else                
          same_site=.true.
          do while(same_site)       
            nr = nr + 1
            read(ustar,'(a)',iostat=ioerr) rlines(nr)
            if(debug) write(*,*) '2nd ioerr nr rlines '
     .         ,ioerr,nr,rlines(nr)
            if( ioerr.eq.-1 ) then
              finished_rcvf = .true.
              same_site = .false.
              nr = nr -1
            elseif ( rlines(nr)(2:5).ne.site ) then
              same_site = .false. 
              nr = nr -1 
              backspace(ustar)
              if(debug) write(*,*) 'new site ',rlines(nr+1)(2:5)
     .               ,'nr same_site ',nr,same_site
            endif
          enddo   
          do i=1,nr 
           read(rlines(i),'(1x,a4,1x,i4,4i3,f6.0,14x,2a20)'
     .         ,iostat=ioerr)sitcodr,iyr,imo,iday,ihr,imin,sec
     .         , rctype1(i),rcvers1(i)     
           rcvrsn1(i) = ' '
c rwk trap a bad entry here since julday doesn't identify the site entry
           if( iyr.lt.1980.or.iyr.gt.2100 ) then
             write(message,'(a,a4,i5,4i3)') 
     .      'Bad receiver start date for ',sitcoda,iyr,imo,iday,ihr,imin
             call report_stat('FATAL','MSTINF','htoglb/mstinf',' ' 
     .                     ,message,ioerr)  
           endif
           if( ioerr.ne.0 ) then
               write(message,'(2a)') 
     .             'Error decoding sta_rcvr receiver line '
     .              ,rlines(i)(1:nblen(rlines(i)))
                   call report_stat('FATAL','MSTINF','htoglb/mstinf',' ' 
     .                              ,message,ioerr)  
           elseif( sitcodr.ne.site ) then
             write(message,'(a,a4,1x,a4)') 
     .         'Something wrong: sitcodr.ne.site ',sitcodr,site
                   call report_stat('FATAL','MSTINF','htoglb/mstinf',' ' 
     .                              ,message,ioerr)                    
           else         
             sitcod1(i) = sitcodr
             sname1(i) = sitcodr
             rcvr_start(1,i) = iyr
             rcvr_start(2,i) = idoy(iyr,imo,iday)
             rcvr_start(3,i) = ihr
             rcvr_start(4,i) = imin
             rcvr_start(5,i) = ifix(sec) 
             rcvr_stop(1,i) = 9999
             rcvr_stop(2,i) = 999 
             rcvr_stop(3,i) = 0
             rcvr_stop(4,i) = 0
             rcvr_stop(5,i) = 0 
           endif
          enddo         
          if( debug ) then  
            write(*,*) nr,' rcvr entries for ',sitcod1(1)
            do i=1,nr
              write(*,'(i5,1x,a4,1x,a16,2(5i5),3(1x,a20))') 
     .           i,sitcod1(1),sname1(1),(rcvr_start(j,i),j=1,5)
     .          ,(rcvr_stop(j,i),j=1,5),rctype1(i),rcvrsn1(i),rcvers1(i)
            enddo
          endif
          rewind(ustaa)        
          if(debug) write(*,*) 'rewound ustaa ',ustaa 
          eof = .false.
          if(debug)  print *,'DEBUG site ',site
          na = 0
          do while (.not.eof ) 
            read(ustaa,'(a)',iostat=ioerr) line
cd            print *,'alines ',line
            if( ioerr.eq.-1 ) then
               eof = .true.                  
            else                          
c             skip commented lines and also entries with no antenna
              if( line(2:5).eq.site.and.line(1:1).eq.' '.and.
     .            line(48:51).ne.'ZERO'.and.line(53:56).ne.'ZERO') then
                na = na + 1
                alines(na) = line  
              endif
            endif
          enddo     
          if(debug)  print *,'na ',na 
          do i=1,na              
            antsn1(i) = ' ' 
            read(alines(i),'(1x,a4,6x,i4,4i3,f6.0,14x,a9,1x,2f11.0,11x
     .         ,f11.0)',iostat=ioerr) sitcoda,iyr,imo,iday,ihr,imin,sec
     .          , goacod,ante1,antn1,anth1
            if( debug ) write(*,*) 'ioerr,i alines ',ioerr,i,alines(i)
c rwk trap a bad entry here since julday doesn't identify the site entry
           if( iyr.lt.1980.or.iyr.gt.2100 ) then
             write(message,'(a,a4,i5,4i3)') 
     .      'Bad antenna start date for ',sitcoda,iyr,imo,iday,ihr,imin
             call report_stat('FATAL','MSTINF','htoglb/mstinf',' ' 
     .                     ,message,ioerr)  
           endif
            call goa_to_igs(goacod,igscod)
            anttyp1(i) = igscod(1:15)
            radome1(i)(1:4) = igscod(17:20)
            radome1(i)(5:5) = ' '   
            if(debug) write(*,'(a,1x,a9,1x,a15,1x,a5)')
     .          'goa igs dome ',goacod,igscod(1:15),radome1(i)
c           serial numbers not reliably available
            antsn1(i) = ' ' 
            if( ioerr.ne.0 ) then
              write(message,'(a,i4)') 
     .          'Error decoding sta_svec antenna line '
              call report_stat('FATAL','MSTINF','htoglb/mstinf',' ' 
     .                     ,message,ioerr)  
            else                 
              sitcod1(i) = sitcoda 
              sname1(i) = sitcoda
              ant_start(1,i) = iyr
              ant_start(2,i) = idoy(iyr,imo,iday)
              ant_start(3,i) = ihr
              ant_start(4,i) = imin
              ant_start(5,i) = ifix(sec) 
              ant_stop(1,i) = 9999
              ant_stop(2,i) = 999 
              ant_stop(3,i) = 0
              ant_stop(4,i) = 0
              ant_stop(5,i) = 0  
              dUNE1(1,i) = anth1 
              dUNE1(2,i) = antn1
              dUNE1(3,i) = ante1 
              dAntAZ1(i) = 0.0     ! Not in GIPSY file.  
            endif  
          enddo      
          if( debug ) then  
             write(*,*) na,' antenna entries for ',sitcod1(1)
             do i=1,na       
               write(*,'(i5,1x,a4,1x,a16,2(5i5),1x,a15,1x,a20,1x
     .            ,a5,3f10.4)') i,sitcod1(i),sname1(i)
     .             ,(ant_start(j,i),j=1,5),(ant_stop(j,i),j=1,5)    
     .             , anttyp1(i),antsn1(i), radome1(i),(dUNE1(j,i),j=1,3)
             enddo
          endif

           
*         The merging routine expects the entries to be time ordered, but 
*         GIPSY uses reverse time order, so reverse the lists
              
          if( debug ) then 
            print *,'unsorted rcv entries  nr ',nr
            do i=1,nr
              write(*,'(i3,a4,1x,2i5,1x,a20)') 
     .          i,site,(rcvr_start(j,i),j=1,2),rctype1(i)
            enddo
          endif                   
cd*** Extreme DEBUG
          if( debug ) then
            read(ustar,'(a)',iostat=ioerr) linex 
            write(*,*) 'before rev_rcvlist ustar line ',ioerr,linex
            backspace(ustar)
          endif

          call rev_rcvlist(max_log,nr,rcvr_start,rcvr_stop,rctype1
     .                    ,rcvrsn1,rcvers1)       
   
cd*** Extreme DEBUG
          if( debug ) then
            read(ustar,'(a)',iostat=ioerr) linex 
            write(*,*) 'after rev_rcvlist ustar line ',ioerr,linex
            backspace(ustar)
          endif

         if( debug ) then 
            print *,'sorted rcv entries  nr ',nr
            do i=1,nr
              write(*,'(i3,a4,1x,2i5,1x,a20)') 
     .          i,site,(rcvr_start(j,i),j=1,2),rctype1(i)
            enddo
          endif 
          
          if( debug ) then 
            print *,'unsorted ant entries  na ',na
            do i=1,na
              write(*,'(i3,a4,1x,2i5,2(1x,a20))') 
     .          i,site,(ant_start(j,i),j=1,2),anttyp1(i),antsn1(i)
            enddo
          endif             

cd*** Extreme DEBUG
          if( debug ) then
            read(ustar,'(a)',iostat=ioerr) linex 
            write(*,*) 'before rev_antlist ustar line ',ioerr,linex
            backspace(ustar)
          endif


          call rev_antlist( max_log,na,ant_start,ant_stop
     .                    ,  anttyp1,antsn1,radome1,dUNE1 )   
cd          print *,'after rev_antlist radome1(2):',radome1(2)     

cd*** Extreme DEBUG
          if( debug ) then
            read(ustar,'(a)',iostat=ioerr) linex 
            write(*,*) 'after rev_antlist ustar line ',ioerr,linex
            backspace(ustar)
          endif



          if( debug ) then 
            print *,'sorted ant entries  na ',na
            do i=1,na
              write(*,'(i3,a4,1x,2i5,2(1x,a20))') 
     .          i,site,(ant_start(j,i),j=1,2),anttyp1(i),antsn1(i)
            enddo
          endif 
cd          print *,'ant_start(1,1-2) ',ant_start(1,1),ant_start(2,1)


*         Now pair up the receiver and antenna entries and insert them at the 
*         end of the global arrays (global counter is num_str, in common)
                                            
         if(debug ) then
            write(*,'(2a,2i3,a4)') 'Calling merge fname nr na sitcod '
     .                      ,fname(1:16),nr,na,sitcod1(1)
            do i=1,nr
              write(*,*) (rcvr_start(j,i),j=1,5),rctype1(i)
            enddo
            do i=1,na
              write(*,*) (ant_start(j,i),j=1,5),anttyp1(i),dUNE(1,i),
     .                    dAntAZ(i)
            enddo  
            print *,'site nr na num_str ',sitcod1(1),nr,na,num_str
          endif  
cd        print *,'calling merge anttyp1 antsn1 dUNE1 '
cd        do i=1,na   
cd           write(*,'(1x,a15,1x,a20)') anttyp1(i),antsn1(i)
cd     .      ,dUNE1(j,i),j=1,3)
cd        enddo                         
cd         print *,'calling ant_start(1-2) ',ant_start(1,1),ant_start(2,1)                      
cd*** Extreme DEBUG
          if( debug ) then
            read(ustar,'(a)',iostat=ioerr) linex 
            write(*,*) 'before merge ustar line ',ioerr,linex
            backspace(ustar)
          endif
          call merge_rcv_ant( fname,nr,na,sitcod1(1),sname1(1)
     .               , rcvr_start,rcvr_stop,rctype1,rcvrsn1,rcvers1
     .               , ant_start,ant_stop,anttyp1,antsn1,radome1,dUNE1,
     .                 dAntAZ1 )
        endif                          
          if( debug ) then
            read(ustar,'(a)',iostat=ioerr) linex 
            write(*,*) 'after merge ustar line ',ioerr,linex
            backspace(ustar)
          endif

cd        print *,'after merge_rcv_ant radome1(2):',radome1(2)
cd        print *,'  radome(2):',radome(2)
        
                    
*       completed all rcvr and antenna entries for this site, read the
*       next one in the receiver sta_rcvr file
        enddo  

      return
      end      
                 
CTITLE GOA_TO_ANT
c
c PURPOSE: This subroutine reads a tranalstion table (on first entry) and 
c          gets the IGS antenna and radome codes from the GIPSY-OASIS codes.
c 
c R. King        18 March 2014
c
      Subroutine goa_to_igs( goa_in, igs_out )
c
c VARIABLES:                         
c     goa_in   - GIPSY-OASIS 9-character code                      C*9
c     igs_out  - IGS 20--character antenna code                    C*20

      implicit none   
                   

      integer*4 max_goa
      parameter (max_goa=300) 

      integer*4 lu,igoa,ngoa,ioerr,i
                                    
      character*4 radome
      character*9 goa_in,goa_ant(max_goa)
      character*20 igs_out,igs_ant(max_goa)
      character*80 line
      character*256 message

      logical eof,found,first_call
      data first_call/.true./
                             
  
      save ngoa,goa_ant,igs_ant,first_call 

               
c   if first-call, read all the values into storage
      if( first_call ) then             
        open(103,file='ant_goa.dat',status='old',iostat=ioerr)
        if( ioerr.ne.0) call report_stat('FATAL','MSTINF','goa_to_igs'
     .     ,'ant_goa.dat','Cannot open GIPSY conversion file',ioerr)
c       find the antenna entries: in the current rendition, they are
c       the first non-blank, non-# lines
        eof = .false.    
        ngoa = 0 
        do while( .not.eof )  
         read(103,'(a)',iostat=ioerr) line 
         if( ioerr.eq.-1 ) then
            eof = .true.     
          elseif( ioerr.ne.0 ) then
            call report_stat('FATAL','MSTINF','goa_to_igs',' '
     .                  ,'Error reading GIPSY antenna code file',0)
          elseif( line(1:1).ne.' '.and.line(1:1).ne.'#')  then
            ngoa = ngoa + 1
            if( ngoa.gt.max_goa ) then
              write(message,'(a,i4)') 
     .        'Number of entries in GIPSY antenna-code file > ',max_goa
               call report_stat('FATAL','MSTINF','goa_to_igs',' '
     .                  ,message,0)
            endif
            read(line,'(a20,1x,a9 )',iostat=ioerr) 
     .             igs_ant(ngoa),goa_ant(ngoa)
cd            write(*,'(a,a,a9,1x,a20)')  'line goa_ant igs_ant '
cd     .           ,line,goa_ant(ngoa),igs_ant(ngoa)
            if( ioerr.ne.0 )              
     .        call report_stat('FATAL','MSTINF','goa_to_igs',' '
     .                ,'Error decoding GIPSY antenna code entries',0)
          endif
        enddo  
        first_call = .false.
        close(103)      
cd        print *,'close 103'
      endif       

  
c  find the match

      igoa = 1       
      found = .false.
      do while (.not.found .and.igoa.le.ngoa) 
cd        print *,'igoa goa_in goa_ant ',igoa,goa_in,goa_ant(igoa)
        if( goa_in.eq.goa_ant(igoa) ) then
          igs_out = igs_ant(igoa)
          found = .true.
        endif 
        igoa = igoa + 1
      enddo             
cd      print *,'GOA_TO_IGS ', igoa,igs_out 
      if( .not.found ) then
       write(message,'(a,a9,a)') 
     .   'GOA code ',goa_in,' not found in ant_goa.dat'
       call report_stat('FATAL','MSTINF','goa_to_igs',' ',message,0)
      endif 
      return
      end

CTITLE  MERGE_RCV_ANT
              
c     Merge receiver and antenna entries for a single station,
c     needed for IGS log files and GIPSY sta_rcvr and sta_svec files in 
c     which receiver and antenna entries are not necessarily paired by
c     dates. 
c
c     Use 'max_log' nominally set to 100, to dimension the number
c     of allowable receiver or antenna entries for a site from the 
c     IGS log or sta files. Insert the merged entries into the
c     global array (counter is num_str).   R. King 130124
* MOD TAH 200203: Added dAntAZ1 to list of values from station.info

      Subroutine merge_rcv_ant( fname,nr,na,sitcod1,sname1
     .      , rcvr_start,rcvr_stop,rctype1,rcvrsn1,rcvers1
     .      , ant_start,ant_stop,anttyp1,antsn1,radome1,dUNE1,
     .        dAntAZ1 )

      implicit none
        

*       Primary variables and index num_str in common
      include 'mstinf.h' 
      
      integer*4 nr,na,rcvr_start(5,max_log),rcvr_stop(5,max_log)
     .        , ant_start(5,max_log),ant_stop(5,max_log),num_save
     .        , ik,ikp,ia,ir,i,j
* MOD TAH 210218: Declated pcncor C*1 and not I*4.
      character*1 pcncod                    
   
      character*4 sitcod1
      character*5 radome1(max_log),char5
      character*6 char6  
* MOD TAH 200213: Changed from  character*10 to character*20 for consistency
      character*20 rcvrsn1(max_log)
      character*15 anttyp1(max_log)   
      character*16 sname1
      character*20 rctype1(max_log),rcvers1(max_log),antsn1(max_log)
     .           , anttyp20,char20
      character*24 char24                              
      character*128 fname 
      character*132 char132,message

      real*8 anth1(max_log),antn1(max_log),ante1(max_log)
     .     , dUNE1(3,max_log),sod,tk,tkp,decyrs
      real*8 dAntAZ1(max_log)
          
      logical finished, debug/.false./,too_late/.false./

c  FUNCTION
      integer*4 itimdif,nblen
        
      debug = dump_log

c     save the begining value of the global index to count entries and set the GAMIT codes
      num_save = num_str
      num_str = num_str + 1
      site_indx(num_str) = num_str
        
c     set a time in the future for the array position after the end to
c     facilitate time checking
      rcvr_start(1,nr+1) = 2100
      rcvr_start(2,nr+1) = 0
      ant_start(1,na+1) = 2100
      ant_start(2,na+1) = 0 
                           
c     start with the later of the two first entries (cannot track w/ rcvr or ant alone).
c     the logic here takes account of the fact that the separateness of the receivr and
c     antenna entries in IGS log files and especially the GIPSY rcv and svec files may 
c     have multiple entries of antennas prior to the first receiver, or vice versa, an 
c     illogical but all-too-common case
      ir = 1
      ia = 1   
      if(debug) write(*,'(a,1x,a4,2(2i5,3i3))') 
     .    'sitcod1 rcvr_start ant_start'
     .    ,sitcod1,(rcvr_start(i,1),i=1,5),(ant_start(i,1),i=1,5)
      if( itimdif(rcvr_start(1,1),ant_start(1,1)).eq.0) then    
        ir = 1
        ia = 1                                           
        if( debug )                 
     .    write(*,'(a,3i3,5i5,2x,5i5,1x,a20,1x,a20)') 
     .    'Common rcv/ant start times num_str ir ia'
     .   ,num_str,ir,ia,(rcvr_start(i,1),i=1,5),(rcvr_stop(i,1),i=1,5)
     .      ,rctype1(1),anttyp1(1)
        call insert_times( 'start',rcvr_start(1,1) )        
        sitcod(num_str) = sitcod1   
        call insert_rcvr_info( rctype1(1),rcvrsn1(1),rcvers1(1)) 
        call insert_ant_info( anttyp1(1),antsn1(1),radome1(1)
     .                      , dUNE1(1,1), dAntAz1(1) )                                                                          
      elseif( itimdif(rcvr_start(1,1),ant_start(1,1)).gt.0) then 
        if( debug )                 
     .    write(*,'(a,3i3,1x,2i4,3i3,2x,2i4,3i3,1x,a20,1x,a20)') 
     .    'Rcvr latest num_str ir ia'
     .   ,num_str,ir,ia,(rcvr_start(i,1),i=1,5),(rcvr_stop(i,1),i=1,5)
     .      ,rctype1(1),anttyp1(1)
        call insert_times( 'start',rcvr_start(1,1) )        
        sitcod(num_str) = sitcod1   
        call insert_rcvr_info( rctype1(1),rcvrsn1(1),rcvers1(1)) 
        print *,'INSERTING ir  num_str',ir,num_str
     .    ,rcvr_start(1,ir),rcvr_start(2,ir),rctype1(ir)
        too_late = .false.
        ia = 1  
        do while (.not.too_late.and.ia+1.le.na+1 )
          if( debug ) print *,'ia' ,ia
          if( itimdif(ant_start(1,ia+1),rcvr_start(1,1)).gt.0 ) then
            too_late = .true.                                     
            if( debug ) then
              write(*,'(a,i3,2i5,3i3)')  
     .            'ia+1 too late: ia ant_start ',ia
     .             ,(ant_start(i,ia),i=1,5)
              write(*,'(a,2i3,a15,a20,a5,f10.4,1x,f5.0)') 
     .        'INSERTING ia num_str anttyp antsn radome dUNE'
     .        ,ia,num_str,anttyp1(ia),antsn1(ia),radome1(ia),
     .         dUNE1(1,ia),dAntAZ1(ia)
            endif
            call insert_ant_info( anttyp1(ia),antsn1(ia),radome1(ia)
     .                        , dUNE1(1,ia), dAntAz1(ia) )     
          endif                         
          ia = ia + 1
        enddo 
        ia = ia - 1
  
      elseif( itimdif(ant_start(1,1),rcvr_start(1,1)).gt.0) then 
        if( debug )                 
     .    write(*,'(a,3i3,1x,2i4,3i3,2x,2i4,3i3,1x,a20,1x,a20)') 
     .    'Antenna latest num_str ir ia'
     .   ,num_str,ir,ia,(ant_start(i,1),i=1,5),(ant_stop(i,1),i=1,5)
     .      ,rctype1(1),anttyp1(1)
        call insert_times( 'start',ant_start(1,1) )        
        sitcod(num_str) = sitcod1    
        call insert_ant_info( anttyp1(1),antsn1(1),radome1(1)
     .                        , dUNE1(1,1), dAntAz1(1) )                             
        if( debug ) write(*,'(a,2i3,a15,a20,a5,f10.4)') 
     .      'INSERTING ia num_str anttyp antsn radome dUNE'
     .      ,ia,num_str,anttyp1(ia),antsn1(ia),radome1(ia),dUNE1(1,ia)
        too_late = .false.
        ir = 1  
        do while (.not.too_late.and.ir+1.le.nr+1 )
          if( itimdif(rcvr_start(1,ir+1),ant_start(1,1)).gt.0 ) then
            too_late = .true.   
            if( debug ) then
              write(*,'(a,i3,2i5,3i3)')  
     .            'ir+1 too late: ir rcvr_start ',ir
     .            ,(rcvr_start(i,ir),i=1,5)
              write(*,'(a,2i3,3a20)') 
     .          'INSERTING ir num_str rctype rcvrsn rcvers '
     .          ,ir,num_str,rctype(ir),rcvrsn1(ir),rcvers1(ir)     
            endif
            call insert_rcvr_info( rctype1(ir),rcvrsn1(ir),rcvers1(ir))
          endif
          ir = ir + 1
        enddo   
* MOD TAH 180802: decremented ir to make code same as ia
        ir = ir - 1
      endif

c     set stop time the earliest of the two entries
      if( itimdif(rcvr_stop(1,ir),ant_stop(1,ia)).le.0 ) then
c        times are the same or receiver stops first
         call insert_times( 'stop ',rcvr_stop(1,ir))
      else
c       antenna stops first
        call insert_times('stop ',ant_stop(1,ia))
      endif
      if( debug ) write(*,'(a,3i4,2(2i5,3i3),a20,a15)') 
     .   'Initial entries inserted ir ia num_str start stop '
     .    , ir,ia,num_str,(istart(i,num_str),i=1,5)
     .    ,(istop(i,num_str),i=1,5),rctype(num_str),anttyp(num_str)

      newflag(num_str) = .true.   
c     now add entries whenever the rcvr or ant change; don't worry about stop 
c     times, as they will be sorted out in the merge
      finished = .false.          
      do while(.not.finished)  
        if( debug ) write(*,'(a,2(i3,2i5,3i3))') 
     .              'Next ir start ia start '
     .             ,ir+1,(rcvr_start(i,ir+1),i=1,5)
     .             ,ia+1,(ant_start(i,ia+1),i=1,5)
        if( (ir+1).le.nr .and.
     .      itimdif(rcvr_start(1,ir+1),ant_start(1,ia+1)).le.0) then
c          next epoch is a receiver entry     
          if( debug ) print *,'Rcvr is next entry '
          ir = ir + 1
          num_str = num_str + 1  
          call check_max(num_str, max_str, 
     .                   'Number of IGS log or GIPSY entries')
          site_indx(num_str) = num_str 
          newflag(num_str) = .true.  
          sitcod(num_str) = sitcod1   
          if( debug)  print *,'INSERTING ir  num_str',ir,num_str
     .    ,rcvr_start(1,ir),rcvr_start(2,ir),rctype1(ir)
          call insert_times( 'start',rcvr_start(1,ir) ) 
          call insert_rcvr_info( rctype1(ir),rcvrsn1(ir),rcvers1(ir) )
c         see what antenna entry to use
          if( (ia+1).le.na.and.
     .        itimdif(rcvr_start(1,ir),ant_start(1,ia+1)).eq.0) then
c           times of next entries match
            ia = ia + 1 
            if( debug ) print *,'Next ant time matches rcvr ia = ',ia
            if( debug)  print *,'INSERTING ia  num_str',ia,num_str
     .         ,ant_start(1,ia),ant_start(2,ia),anttyp1(ia),radome1(ia)
            call insert_ant_info( anttyp1(ia),antsn1(ia),radome1(ia)
     .                       , dUNE1(1,ia), dAntAz1(ia) )
          else
c           use the last antenna entry
            if( debug ) print *,'Next ant time later than rcvr ia = ',ia 
            if( debug)  write(*,'(a,3i4,2i5,3i3,1x,a15,1x,a4)')
     .             'INSERTING ir ia  num_str',ir,ia,num_str
     .            ,(rcvr_start(i,ir),i=1,5),anttyp1(ia),radome1(ia)
            call insert_ant_info( anttyp1(ia),antsn1(ia),radome1(ia)
     .                       , dUNE1(1,ia), dAntAz1(ia) )
          endif
          if( debug ) print *,'Inserted num_str ir ia ',num_str,ir,ia
     .      ,(rcvr_start(i,ir),i=1,5),rctype1(1),anttyp1(1)
        elseif (ia.lt.na ) then
c         next epoch is an antenna entry
          if( debug ) print *,'Ant is next entry '
          ia = ia + 1
          num_str = num_str + 1   
          call check_max(num_str, max_str, 
     .                   'Number of IGS log or GIPSY entries')
          site_indx(num_str) = num_str   
          newflag(num_str) = .true.         
          sitcod(num_str) = sitcod1                            
         if( debug)  write(*,'(a,3i4,2i5,3i3,1x,a15,a4)')
     .        'INSERTING ir ia  num_str',ir,ia,num_str
     .         ,(ant_start(i,ia),i=1,5),anttyp1(ia),radome1(ia)
          call insert_times( 'start',ant_start(1,ia) ) 
          call insert_ant_info( anttyp1(ia),antsn1(ia),radome1(ia)
     .                       , dUNE1(1,ia),dAntAz1(ia) )
c         since we've checked the equal-times case, and used the rcvr entry
c         as primary, the receiver must not have changed; use the last receiver entry       
          if( debug)  write(*,'(a,2i4,2i5,3i3)')
     .         'INSERTING ir num_str',ir,num_str
     .         ,(rcvr_start(i,ir),i=1,5),rctype1(ir)
          call insert_rcvr_info( rctype1(ir),rcvrsn1(ir),rcvers1(ir) )   
        else
          finished = .true.         
          if( debug ) print *,'Finished ir nr ia na ',ir,nr,ia,na
        endif         
c       set stop time the earliest of the two entries
        if( itimdif(rcvr_stop(1,ir),ant_stop(1,ia)).le.0 ) then
c          times are the same or receiver stops first
           call insert_times( 'stop ',rcvr_stop(1,ir))
        else
c         antenna stops first
          call insert_times('stop ',ant_stop(1,ia))
        endif
      enddo     
             
*     Adjust the start and stop times to avoid overlaps or (optionally) gaps
                                    
      if( debug ) print *,'Adjusting times, no_gaps ',no_gaps
        
      do i = num_save+1,num_str-1      
        ik = abs(site_indx(i))    
        if( debug ) print *,'i,ik istart istop '
     .          ,i,ik,(istart(j,ik),j=1,5),(istop(j,ik),j=1,5)     
        ikp = abs(site_indx(i+1))
        if( debug ) print *,'i, ikp istart istop '
     .        ,i,ikp,(istart(j,ikp),j=1,5),(istop(j,ikp),j=1,5)     
        if( itimdif(istop(1,ik),istart(1,ikp)).gt.0 ) then  
          do j=1,5
            istop(j,ik) = istart(j,ikp)   
          enddo     
          if( debug ) print *,'overlap fix ikp istop '
     .         ,ikp,(istop(j,ikp),j=1,5)
        elseif ( no_gaps .and. 
     .           itimdif(istop(1,ik),istart(1,ikp)).lt.0 ) then 
          do j=1,5
            istop(j,ik) = istart(j,ikp)   
          enddo     
          if( debug ) print *,'gap fix ikp istop '
     .         ,ikp,(istop(j,ikp),j=1,5)
          write(message,'(a,2i4,3i3,a,a4)') 'Closing gap at '
     .         ,(istart(j,ikp),j=1,5),' for ',sitcod1  
          call report_stat('WARNING','MSTINF','htoglb/mstinf',' '
     .       ,message,0)
        endif
      enddo


CD DEBUG -- print the current list  
      if( debug ) then
        print *,'List after insertions (i site_indx(i) )'
        do i = 1, num_str
          ik = abs(site_indx(i)) 
          write(*,'(2i4,1x,1x,a4,1x,l1,1x,a4,2(1x,i4,1x,i3)
     .        ,1x,a20,1x,a15,1x,3f8.4,1x,f5.0)')
     .        i,site_indx(i),sitcod(ik),newflag(ik),sitcod1
     .      ,istart(1,ik),istart(2,ik),istop(1,ik),istop(2,ik)
     .      ,rctype(ik),anttyp(ik),(dUNE(j,ik),j=1,3),
     .       dAntAZ(ik)
        enddo 
      endif
CD END DEBUG       
       
*     For the new/old sorting logic in merge_stinf to work, we need to 
*     combine the new entries here before exiting.  The entries are time-ordered,
*     so we can assume that if the start times are the same within dup_toler 
*     (default 120s but setable in input), these are duplicates and set the site 
*     index negative to have the first value commented out.  
                  
      do i=num_save+1,num_str-1 
        sod=istart(3,i)*3600.d0 + istart(4,i)*60.d0 + istart(5,i)
        tk =  decyrs(istart(1,i),istart(2,i),sod)
        sod=istart(3,i+1)*3600.d0 + istart(4,i+1)*60.d0 + istart(5,i+1)
        tkp = decyrs(istart(1,i+1),istart(2,i+1),sod)
        if( dabs(tk-tkp).lt.dup_tol ) then  
          if(debug) print *,'duplicate set negative i tk tkp ',i,tk,tkp
          site_indx(i) = -abs(site_indx(i))
        endif
      enddo
                
*     Set the site names and GAMIT codes for all the new entries      
      
cd      print *,'Still in MERGE num_save+1, num_str ',num_save+1,num_str
      do i=num_save+1,num_str
        sitcod(i) = sitcod1
        sname(i) =  sname1      
        htcod(i) = 'DHARP' 
        comment(i)=
     .       'mstinf: '//fname(1:nblen(fname))   
        call read_rcvant( 2,2,char6,char20,radome(i)
     .                  ,rcvcod(i),rctype(i),pcncod )    
        anttyp20 = anttyp(i)//'     '
        call read_rcvant( 2,1,antcod(i),anttyp20,radome(i)
     .                  , char6,char20,pcncod )      
        if( radome(num_str).eq.'     ' ) radome(num_str)='UNKN'
c       use guess_rcvant.dat only for the firmware code 
        call guess_rcvant(2,fname,char20,char20,char5,rcvers(i)
     .                  , char6, swver(i),char6,char132 )
c       remove ^s added by guess_rcvant 
        call sub_char ( rcvers(i),'^',' ')   
*       --now get the GAMIT codes for the receiver and antenna reassign according to new IGS standards 
        if( debug ) then
          print *,'DEBUG  sitcod: ',sitcod(i) 
          print *,'       sname : ',sname(i)
          print *,'       rctype: ',rctype(i)
          print *,'       rcvcod: ',rcvcod(i) 
          print *,'       rcvers: ',rcvers(i)
          print *,'       swver : ',swver(i) 
          print *,'       anttyp: ',anttyp(i)
          print *,'       radome: ',radome(i)
        endif
        if( radome(num_str).eq.'     ' ) radome(num_str)='UNKN '
      enddo
                 
      if( debug ) then 
        write(*,*) 'Here is the full current list '
        do i=1,num_str
          write(*,'(1x,i3,1x,a4,1x,a16,2(1x,5i5),3(1x,a20),1x,a15
     .          ,1x,a20,1x,a4,4f10.4,1x,f5.0)')
     .        i,sitcod(i),sname(i),(istart(j,i),j=1,5)
     .       , (istop(j,i),j=1,5)
     .       , rctype(i),rcvrsn(i),rcvers(i),anttyp(i),antsn(i)
     .       , radome(i),(dUNE(j,i),j=1,3), dAntAZ(i)
        enddo
      endif        

CD DEBUG -- print the current list  
      if( debug ) then
        print *,'List after removing duplicates(i site_indx(i) )'
        do i = 1, num_str
          ik = abs(site_indx(i))
          write(*,'(2i4,1x,l1,1x,a4,2(1x,i4,1x,i3),1x,a20,1x,a15)') 
     .      i,site_indx(i),newflag(ik),sitcod(ik)               
     .      ,istart(1,ik),istart(2,ik),istop(1,ik),istop(2,ik)
     .      ,rctype(ik),anttyp(ik)
        enddo 
      endif
CD END DEBUG

      return
      end 
                
CTITLE REV_RCVLIST
          
*     Reverse the time order of the GIPSY receiver list

      Subroutine rev_rcvlist(max_log,nr,rcvr_start1,rcvr_stop1
     .                       ,rctype1,rcvrsn1,rcvers1)
                   
      implicit none     
                
*       Arguments
      integer*4 max_log,nr,rcvr_start1(5,max_log),rcvr_stop1(5,max_log)
      character*20 rctype1(max_log),rcvrsn1(max_log),rcvers1(max_log)


*       Local                                       
      integer*4  rcvr_start2(5,max_log),rcvr_stop2(5,max_log),i,j
      character*20  rctype2(max_log),rcvers2(max_log)
     .           , rcvrsn2(max_log)      
                   
      do i=1,nr  
        do j=1,5
          rcvr_start2(j,i) = rcvr_start1(j,nr+1-i)
          rcvr_stop2(j,i) = rcvr_stop1(j,nr+1-i)
        enddo
        rctype2(i) = rctype1(nr+1-i)
        rcvers2(i) = rcvers1(nr+1-i)
        rcvrsn2(i) = rcvrsn1(nr+1-i)
      enddo   
      do i=1,nr        
        do j=1,5
          rcvr_start1(j,i) = rcvr_start2(j,i)
          rcvr_stop1(j,i) = rcvr_stop2(j,i)
        enddo
        rctype1(i) = rctype2(i)
        rcvers1(i) = rcvers2(i)
        rcvrsn1(i) = rcvrsn2(i)
      enddo
      return
      end
                
CTITLE REV_ANTLIST
          
*     Reverse the time order of the GIPSY antenna list

      Subroutine rev_antlist( max_log,na,ant_start1,ant_stop1
     .                       , anttyp1,antsn1,radome1,dUNE1)

      integer*4 ant_start1(5,max_log),ant_stop1(5,max_log)
     .        , ant_start2(5,max_log),ant_stop2(5,max_log)
      character*15 anttyp1(max_log),anttyp2(max_log)
      character*20 antsn1(max_log),antsn2(max_log)
      character*5  radome1(max_log),radome2(max_log)
      real*8 dUNE1(3,max_log),dUNE2(3,max_log)
                   
      do i=1,na  
        do j=1,5
          ant_start2(j,i) = ant_start1(j,na+1-i)
          ant_stop2(j,i) = ant_stop1(j,na+1-i)
        enddo  
        anttyp2(i) = anttyp1(na+1-i)
        antsn2(i) = antsn1(na+1-i)
        radome2(i) = radome1(na+1-i)  
        do j=1,3
          dUNE2(j,i) = dUNE1(j,na+1-i)
        enddo
      enddo 
      do i=1,na        
        do j=1,5
          ant_start1(j,i) = ant_start2(j,i)
          ant_stop1(j,i) = ant_stop2(j,i)
        enddo
        anttyp1(i) = anttyp2(i)
        antsn1(i) = antsn2(i)
        radome1(i) = radome2(i)     
        do j=1,3
          dUNE1(j,i) = dUNE2(j,i)
        enddo
      enddo
      return
      end

CTITLE  INSERT_RCVR_INFO

      Subroutine insert_rcvr_info( rctype1,rcvrsn1,rcvers1 )    
                   
      implicit none

*       Primary variables and index num_str in common
      include 'mstinf.h' 
       
*       Local
* MOD TAH 200213: Changed from  C*20 and C*10 to character*(*) for consistency
      character*(*) rctype1,rcvers1
      character*(*) rcvrsn1      

      rctype(num_str) = rctype1
      rcvrsn(num_str) = rcvrsn1
      rcvers(num_str) = rcvers1              
      return
      end


CTITLE  INSERT_ANT_INFO                                      

      Subroutine insert_ant_info( anttyp1,antsn1,radome1,dUNE1,dAntAz1)
                     
      implicit none

*       Primary variables and index num_str in common
      include 'mstinf.h' 

*       Local                  
      character*20 antsn1
      character*15 anttyp1
      character*5  radome1
      real*8 dUNE1(3)
      real*8 dAntAz1  ! Antenna aligment from True N (deg).
      integer*4 i
      logical debug

      debug = dump_log
              
cd      print *,'In INSERT_ANT_INFO dUNE1 ',(dUNE1(i),i=1,3)
      anttyp(num_str) = anttyp1              
      antsn(num_str) = antsn1
      radome(num_str) = radome1
      do i=1,3
         dUNE(i,num_str) = dUNE1(i)
      enddo   
* MOD TAH 200203: Save antenna Azimuth.
      dAntAz(num_str) = dAntAz1       
      if( debug )
     .  write(*,'(a,i4,1x,a15,1x,a20,1x,a5,2(3f10.4),1x,F5.1)') 
     .     'INSERT_ANT_INFO num_str anttyp1 antsn1 radome1 dUNE1 dUNE '
     .    , num_str,anttyp1,antsn1,radome1,(dUNE1(i),i=1,3)
     .    , (dUNE(i,num_str),i=1,3), dAntAz(num_str)
      return
      end
     
CTITLE INSERT_TIMES
 
      Subroutine insert_times( epoch,time )

      implicit none

*       Primary variables and index num_str in common
      include 'mstinf.h' 
  
      integer*4 time(5),i
      character*5 epoch  
                        
      if( epoch.eq.'start' ) then
        do i=1,5
          istart(i,num_str) = time(i)
        enddo
      elseif( epoch.eq.'stop ') then
        do i=1,5
          istop(i,num_str) = time(i)
        enddo    
      endif
      return
      end
                                             
CCONVERT_IGST

      Subroutine convert_igst( A_or_R,S_or_E,atime,itime )    

      implicit none

*   Convert epochs from IGS log files to yy, doy, hr, min sec
*   If non-integer values encountered, set start to 1980 1 0 0 0 0 
*   and end to o 9999 999 0 0 0 
                                     
* Character indicating antenna (A) or receivers (R)  (used for error print
      character*3 A_or_R                

*  Character indicating start (S) or end (E) (affects defaults)
      character*1 S_or_E

*  Time as read from log CCYY-MM-DDThh:mmz
      character*17 atime                  

*  Times converted to integers  (yr,doy,hr,min,sec)
      integer*4 itime(5)  

*  External function
      integer*4 idoy
                                                      
*  Local                  
      character*256 message
      integer*4 imon,iday,ihr,min,ioerr  

      logical debug/.false./           

*  Seconds not specified on IGS log file
      itime(5) = 0
       
*  Check if date is template      
      if(debug) print *,'atime ',atime
      if( atime(2:5).eq.'CCYY'.or.atime(1:4).eq.'CCYY'.or.
     .    atime(1:4).eq.'    ' ) then
        if( S_or_E.eq.'S' ) then
          write(message,'(a,a3,a)') 'Start date missing for ',a_or_r
     .          ,' in IGS log file'
          call report_stat('WARNING','MSTINF','htoglb/mstinf',' '
     .       ,message,0)
          itime(1) = 1980
          itime(2) = 1
        elseif( S_or_E.eq.'E') then
          itime(1) = 9999
          itime(2) = 999
        endif
      else
        read(atime,'(i4,1x,i2,1x,i2)',iostat=ioerr) itime(1),imon,iday
        if(debug) print *,'itime(1),imon,iday ',itime(1),imon,iday
        if( ioerr.ne.0 ) then 
          if( S_or_E.eq.'S' ) then  
             write(message,'(a,a3,a,a17)') 'IGS log error decoding '
     .         ,a_or_r,' Date Installed: ',atime
             call report_stat('FATAL','MSTINF','htoglb/mstinf',' '
     .                        ,message,ioerr)
          else   
             write(message,'(a,a3,a,a17)') 'IGS log error decoding '
     .        ,a_or_r,' Date Removed: ',atime
            call report_stat('FATAL','MSTINF','htoglb/mstinf',' '
     .                      ,message,ioerr)
          endif
        else
          itime(2) = idoy(itime(1),imon,iday)
        endif
        if( atime(12:13).eq.'hh'.or.atime(12:13).eq.'  '.or.
     .      atime(13:14).eq.'hh'  ) then
          itime(3) = 0
          itime(4) = 0
          itime(5) = 0
        else
          read(atime(12:16),'(i2,1x,i2)',iostat=ioerr) 
     .         itime(3),itime(4)
          if( ioerr.ne.0 ) then
            if( S_or_E.eq.'S' ) then
              call report_stat('FATAL','MSTINF','htoglb/mstinf',' '
     .       ,'Error decoding hh:mm of Date Installed on IGS log',ioerr)
            else
              call report_stat('FATAL','MSTINF','htoglb/mstinf',' '
     .         ,'Error reading hh:mm of Date Removed on IGS log',ioerr)
            endif
          endif
        endif
      endif
 
      return
      end

CTITLE UPD_FROM_SINEX

      Subroutine upd_from_sinex ( usnx,isnxf )

      implicit none
     
*     Routine to get station.info entries from a SINEX file. 
*     rwk October 2008 from surboutine update_stinf (now upd_from_rinex). 
                                                              
*       iusnx is the unit number, set in the main program
*       ilogf  is the index in the list of SINEX files requested   
*       sitcod1 is the 4-character ID of the requested site; if blank, update all

*     Routine adds all log entries to the global array, relying on 
*     merge_stinf to replace old entries if requested
  
* COMMON VARIABLES

      include 'mstinf.h' 
      include '../../gamit/includes/makex.h'  
              
* LOCAL VARIABLES
              
      integer*4 usnx,isnxf
                                
      integer rcvr_start(5,max_snx_site),rcvr_stop(5,max_snx_site)
     .      , ant_start(5,max_snx_site),ant_stop(5,max_snx_site)
     .      , ecc_start(5,max_snx_site),ecc_stop(5,max_snx_site)
     .      , iy1,iday1,isod1,iy2,iday2,isod2,ihr,imin,isec
     .      , ioerr,num_save,nss,nsr,nsa,nse,i,j

      real*8 anth1(max_snx_site),antn1(max_snx_site),ante1(max_snx_site)
     .     , anthx,antnx,antex,sod,sec,coord(6),coordx,apr_epoch
      real*8 antdaz  ! Antenna aligment from True N (deg).
      integer*4 iantdaz  ! Integer version for reading from SINEX file
      real*8 antdaz1(max_snx_site) ! Save array values
                  
      character*1 pcncod                    
      character*3 ecctypx
      character*4 sitcod1(max_snx_site),sitcodx,radomex,crdtyp
      character*5 radome1(max_snx_site),rcvrsnx,antsnx,char5 
      character*6 char6      
      character*11 rcversx 
      character*15 anttypx
      character*16 sname1(max_snx_site),snamex                          
      character*20 rctype1(max_snx_site),rcvrsn1(max_snx_site)
     .            ,rcvers1(max_snx_site),anttyp1(max_snx_site)
     .            ,antsn1(max_snx_site),rctypex,anttyp20,char20
      character*24 char24
      character*80 line                                
      character*132 char132
      character*256 message
     
      logical found,block_end,finished,debug
      data debug/.false./
              
* FUNCTION               
      integer*4 nblen,itimdif
      real*8 decyrs

*      Report what we're doing
                 
      write(message,'(2a)') 
     . 'Updating station.info entries from SINEX file ',sx_file(isnxf)
      call report_stat('STATUS','MSTINF','htoglb/mstinf',' '
     .                ,message,0)  
      call report_stat('WARNING','MSTINF','htoglb/mstinf',' '
     .                ,message,0)

*     Initialize the times for facilate checking
      do i=1,max_snx_site
        rcvr_start(1,i) = 0
        rcvr_start(2,i) = 0
        rcvr_start(3,i) = 0
        rcvr_start(4,i) = 0
        rcvr_start(5,i) = 0
        ant_start(1,i) = 9999
        ant_start(2,i) = 999
        ant_start(3,i) = 0
        ant_start(4,i) = 0
        ant_start(5,i) = 0
      enddo                

      if( debug ) write(*,'(a,a4)') 'snx_site= ',snx_site
      
*     Read the file until all useful blocks searched
                                    
      finished = .false. 
      do while( .not.finished )
               
        read(usnx,'(a)', iostat=ioerr) line  
        if( debug )   print *,'read line ioerr ',line,ioerr
        if( ioerr.eq.-1 ) then
          finished = .true.
        elseif( ioerr.ne.0 ) then
          call report_stat('FATAL','MSTINF','htoglb/mstinf',' '
     .      ,'Error reading SINEX line',ioerr)
        elseif( line(1:8).eq. '+SITE/ID' ) then
          block_end = .false.  
          found = .false.    
          nss = 0
          do while (. not.block_end .and. .not.found ) 
            read(usnx,'(a)', iostat=ioerr) line
            if( ioerr.ne.0 ) then
              call report_stat('FATAL','MSTINF',
     .           'htoglb/mstinf',' ','Error reading SITE ID line',ioerr)    
            elseif( line(1:8).eq.'-SITE/ID' ) then 
              if( debug )  print *,'found -SITE/ID'
              block_end = .true.
            elseif( line(1:1).eq.' ' ) then
              read(line,'(1x,a4,16x,a16)',iostat=ioerr) sitcodx,snamex
              if( ioerr.ne.0 )  then
                call report_stat('FATAL','MSTINF','htoglb/mstinf',' '
     .                   ,'Error decoding SINEX site name entry',ioerr)
              else      
                if( snx_site.eq.'    ') then   
                  nss = nss + 1   
                  if( nss.gt.max_snx_site ) then
                     write(message,'(a,i4,a)') 
     .                 '# SINEX sites > max_snx_site (',max_snx_site,')'
                     call report_stat('FATAL','MSTINF','htoglb/mstinf'
     .                     ,' ',message,ioerr)
                  endif
                elseif( sitcodx.eq.snx_site ) then
                  nss = 1
                  found = .true.
                endif               
                sitcod1(nss) = sitcodx
                sname1(nss) = snamex 
                if( debug ) write(*,'(a,i4,1x,a4)') 
     .               'nss sitcod1 ',nss,sitcod1
              endif
            endif
          enddo

        elseif( line(1:14).eq. '+SITE/RECEIVER' ) then
          block_end = .false.  
          found = .false.  
          nsr = 0              
          if( debug )  print *,'found SITE/RECEIVER'
          do while (. not.block_end .and. .not.found ) 
            read(usnx,'(a)', iostat=ioerr) line
            if( ioerr.ne.0 ) then
              call report_stat('FATAL','MSTINF'
     .                        ,'htoglb/mstinf',' '
     .                        ,'Error reading SITE/RECEIVER line',ioerr)
            elseif( line(1:14).eq.'-SITE/RECEIVER' ) then 
              block_end = .true.
              if(debug) print *,'end of SITE/RECEIVER block'
            elseif( line(1:1).eq.' ' ) then
              read(line
     .         ,'(1x,a4,10x,2(1x,i2,1x,i3,1x,i5),1x,a20,1x,a5,1x,a11)'
     .         , iostat=ioerr) sitcodx,iy1,iday1,isod1,iy2,iday2,isod2
     .         , rctypex,rcvrsnx,rcversx
              if( debug ) write(*,'(a,a4)') 'read rcv sitcodx ',sitcodx
              if( ioerr.ne.0 )  then
                call report_stat('FATAL','MSTINF','htoglb/mstinf',' '
     .                     ,'Error decoding SINEX receiver entry',ioerr)
              else                   
                if( snx_site.eq.'    ') then   
                  nsr = nsr + 1    
                  if(debug) write(*,'(a,i4,2(1x,a4))') 
     .                 'nsr sitcodx sitcod1 ',nsr,sitcodx,sitcod1
                  if( sitcodx.ne.sitcod1(nsr) ) then
                    write(message,'(a,i4,a,a4,a,a4,a )') 'SINEX sitcod '
     .                 ,nsr,' for receiver (',sitcodx
     .               ,') does not match code for site id (',sitcod1(nsr)
     .                ,')'
                     call report_stat('FATAL','MSTINF','htoglb/mstinf'
     .                 ,' ',message,0) 
                  endif 
                else
                  nsr = 1 
                  found = .true.
                endif  
                call fix_y2k(iy1)   
                rcvr_start(1,nsr) = iy1
                rcvr_start(2,nsr) = iday1 
                sod = dfloat(isod1)
                call ds2hms(iy1,iday1,sod,ihr,imin,sec)
                rcvr_start(3,nsr) = ihr
                rcvr_start(4,nsr) = imin
                rcvr_start(5,nsr) = idint(sec)
                call fix_y2k(iy2)   
                rcvr_stop(1,nsr) = iy2
                rcvr_stop(2,nsr) = iday2    
                sod = dfloat(isod2)
                call ds2hms(iy2,iday2,sod,ihr,imin,sec) 
                rcvr_stop(3,nsr) = ihr
                rcvr_stop(4,nsr) = imin
                rcvr_stop(5,nsr) = idint(sec)  
                rctype1(nsr) = rctypex    
                rcvrsn1(nsr) = ' ' 
                rcvrsn1(nsr)(1:5) = rcvrsnx  
                rcvers1(nsr) = ' ' 
                rcvers1(nsr)(1:11) = rcversx
                if( debug )   write(*,*) 'nsr ',nsr
     .            , (rcvr_start(j,nsr),j=1,5),(rcvr_stop(j,nsr),j=1,5)
     .            , rcvrsn1(nsr),rcvers1(nsr)
              endif
            endif
          enddo          
          if( nsr.ne.nss ) then
            write(message,'(a,i3,a,i3,a)') '# SITE/RECEIVER entries ('
     .      ,nsr,') not equal # SITE/ID entries (',nss,') on SINEX file'
            call report_stat('FATAL','MSTINF','htoglb/mstinf',' '
     .                       ,message,0)
          endif
                                                      

        elseif( line(1:13).eq. '+SITE/ANTENNA' ) then
          block_end = .false.  
          found = .false.  
          nsa = 0   
         if( debug )   print *,'found SITE/ANTENNA '
          do while (. not.block_end .and. .not.found ) 
            read(usnx,'(a)', iostat=ioerr) line
            if( ioerr.ne.0 ) then
              call report_stat('FATAL','MSTINF'
     .                        ,'htoglb/mstinf',' '
     .                        ,'Error reading SITE/ANTENNA line',ioerr)
            elseif( line(1:13).eq.'-SITE/ANTENNA' ) then 
              block_end = .true.
            elseif( line(1:1).eq.' ' ) then

* MOD TAH 200203: Added AntDAZ to list of values from station.info
              read(line,310,iostat=ioerr) 
     .            sitcodx,iy1,iday1,isod1,iy2,iday2,isod2
     .         , anttypx,radomex,antsnx, iantdaz
 310          format(1x,a4,10x,2(1x,i2,1x,i3,1x,i5),1x,
     .                 a15,1x,a4,1x,a5,1x,I4)
              antdaz = iantdaz
              if( ioerr.ne.0 )  then
                call report_stat('FATAL','MSTINF','htoglb/mstinf',' '
     .                      ,'Error decoding SINEX antenna entry',ioerr)
              else                   
                if( snx_site.eq.'    ') then   
                  nsa = nsa + 1
                  if(sitcodx.ne.sitcod1(nsa)) then
                    write(message,'(a,i4,a,a4,a,a4,a )') 'SINEX sitcod '
     .                 ,nsa,' for antenna (',sitcodx
     .               ,') does not match code for site id (',sitcod1(nsa)
     .                ,')'
                     call report_stat('FATAL','MSTINF','htoglb/mstinf'
     .                 ,' ',message,0) 
                  endif 
                else
                  nsa = 1 
                  found = .true.
                endif  
                call fix_y2k(iy1)   
                ant_start(1,nsa) = iy1
                ant_start(2,nsa) = iday1   
                sod = dfloat(isod1)
                call ds2hms(iy1,iday1,sod,ihr,imin,sec)
                ant_start(3,nsa) = ihr
                ant_start(4,nsa) = imin
                ant_start(5,nsa) = idint(sec) 
                call fix_y2k(iy2)   
                ant_stop(1,nsa) = iy2
                ant_stop(2,nsa) = iday2      
                sod = dfloat(isod2)
                call ds2hms(iy2,iday2,sod,ihr,imin,sec)
                ant_stop(3,nsa) = ihr
                ant_stop(4,nsa) = imin
                ant_stop(5,nsa) = idint(sec)      
                anttyp1(nsa) = ' ' 
                anttyp1(nsa)(1:15) = anttypx   
                radome1(nsa) = ' '
                radome1(nsa)(1:4) = radomex 
                antsn1(nsa) = ' ' 
                antsn1(nsa)(1:5) = antsnx 
                antdaz1(nsa) = antdaz
                if( debug ) 
     .            print *,'nsa ',nsa
     .             ,(ant_start(j,nsa),j=1,5),(ant_stop(j,nsa),j=1,5)
     .             ,anttyp1(nsa),radome1(nsr),antsn1(nsr)
              endif
            endif
          enddo          
          if( nsa.ne.nss ) then
            write(message,'(a,i3,a,i3,a)') '# SITE/ANTENNA entries ('
     .      ,nsa,') not equal # SITE/ID entries (',nss,') on SINEX file'
            call report_stat('FATAL','MSTINF','htoglb/mstinf',' '
     .                       ,message,0)
          endif
      
        elseif( line(1:18).eq. '+SITE/ECCENTRICITY' ) then
          block_end = .false.  
          found = .false.  
          nse = 0         
          if( debug )  print *,'found SITE/ECCENTRICITY'
          do while (. not.block_end .and. .not.found ) 
            read(usnx,'(a)', iostat=ioerr) line
            if( ioerr.ne.0 ) then
              call report_stat('FATAL','MSTINF'
     .                        ,'htoglb/mstinf',' '
     .                    ,'Error reading SITE/ECCENTRICITY line',ioerr)
            elseif( line(1:18).eq.'-SITE/ECCENTRICITY' ) then 
              block_end = .true.
            elseif( line(1:1).eq.' ' ) then
              read(line
     .         ,'(1x,a4,10x,2(1x,i2,1x,i3,1x,i5),1x,a3,3f9.0)'
     .         , iostat=ioerr) sitcodx,iy1,iday1,isod1,iy2,iday2,isod2
     .         , ecctypx,anthx,antnx,antex
              if( ioerr.ne.0 )  then
                call report_stat('FATAL','MSTINF','htoglb/mstinf',' '
     .                       ,'Error decoding SINEX eccentricity',ioerr)
              else                   
                if( snx_site.eq.'    ') then   
                  nse = nse + 1
                  if(sitcodx.ne.sitcod1(nse)) then
                    write(message,'(a,i4,a,a4,a,a4,a )') 'SINEX sitcod '
     .                 ,nse,' for eccentricity (',sitcodx
     .               ,') does not match code for site id (',sitcod1(nse)
     .                ,')'
                     call report_stat('FATAL','MSTINF','htoglb/mstinf'
     .                 ,' ',message,0) 
                  endif 
                else
                  nse = 1 
                  found = .true.
                endif   
                if( ecctypx.ne.'UNE' ) call report_stat('FATAL','MSTINF'
     .            ,'htoglb/mstinf',' '
     .             ,'Cannot handle SINEX eccentricity type not UNE',0)
                call fix_y2k(iy1)   
                ecc_start(1,nse) = iy1
                ecc_start(2,nse) = iday1
                sod = dfloat(isod1)
                call ds2hms(iy1,iday1,sod,ihr,imin,sec)
                ecc_start(3,nse) = ihr
                ecc_start(4,nse) = imin
                ecc_start(5,nse) = idint(sec) 
                call fix_y2k(iy2)   
                ecc_stop(1,nse) = iy2
                ecc_stop(2,nse) = iday2      
                sod = dfloat(isod2)
                call ds2hms(iy2,iday2,sod,ihr,imin,sec)
                ecc_stop(3,nse) = ihr
                ecc_stop(4,nse) = imin
                ecc_stop(5,nse) = idint(sec)                                 
*               check if times match within 5 minutes
                if( abs(itimdif(ecc_start(1,nse),ant_start(1,nse)))
     .              .gt.300 .or. 
     .             abs(itimdif(ecc_start(1,nse),ant_start(1,nse)))
     .              .gt.300) then
                  write(message,'(3a,a4)') 'SINEX eccentrity times '
     .                 ,'do not match antenna times ','for site '
     .                 ,sitcod1(nse)
                  call report_stat('FATAL','MSTINF','htoglb/mstinf'
     .                             ,' ',message,0)
                endif
                anth1(nse) = anthx
                antn1(nse) = antnx
                ante1(nse) = antex   
                if( debug ) 
     .             print *,'nse ',nse
     .             , (ecc_start(j,nse),j=1,5),(ecc_stop(j,nse),j=1,5)
     .             , anth1(nse),antn1(nse),ante1(nse)
              endif
            endif
          enddo          
          if( nse.ne.nss ) then
            write(message,'(a,i3,a,i3,a)') 
     .           '# SITE/ECCENTRICITY entries (',nse
     .          ,') not equal # SITE/ID entries (',nss,') on SINEX file'
            call report_stat('FATAL','MSTINF','htoglb/mstinf',' '
     .                       ,message,0)
          endif  

        elseif(write_apr .and. line(1:17).eq. '+SOLUTION/APRIORI' ) then
          block_end = .false.  
          nse = 0         
          if( debug ) print *,'found SOLUTION/APRIORI'
          do while (. not.block_end )  
            read(usnx,'(a)', iostat=ioerr) line
            if( debug ) print *,'line ',line
            if( ioerr.ne.0 ) then
              call report_stat('FATAL','MSTINF'
     .                        ,'htoglb/mstinf',' '
     .                    ,'Error reading SOLUTION/APRIORI line',ioerr)
            elseif( line(1:17).eq.'-SOLUTION/APRIORI' ) then 
              block_end = .true.
            elseif( line(1:1).eq.' ' ) then  
              read(line
     .         ,'(7x,a4,3x,a4,9x,i2,1x,i3,1x,i5,8x,d21.15)'
     .         ,iostat=ioerr) crdtyp,sitcodx,iy1,iday1,isod1,coordx
              call fix_y2k(iy1)
              if( ioerr.ne.0 ) then  
                 call report_stat('FATAL','MSTINF','htoglb/mstinf',' '
     .               ,'Error decoding SINEX apriori coordinate',ioerr)   
              else                                   
                if( debug ) print *,'sitcodx crdtyp coordx '
     .              ,sitcodx,crdtyp,coordx
                if( crdtyp.eq.'STAX' ) then
                  coord(1) = coordx   
                  if( debug ) print *,'coord(1) ',coord(1)
                elseif( crdtyp.eq.'STAY' ) then
                  coord(2) = coordx        
                 if( debug ) print *,'coord(2) ',coord(2)
                elseif( crdtyp.eq.'STAZ' ) then
c                 see if next entry is velocity for this site or postion for another site
                  read(usnx,'(a)') line
                  if( line(8:11).eq.'STAX' ) then 
                    coord(3) = coordx  
                    if( debug ) print *,'coord(3) ',coord(3)
                    do i=4,6
                      coord(i) = 0.d0  
                    enddo     
c                   end of site, write apr entry
                    apr_epoch = decyrs(iy1,iday1,dfloat(isod1))
                    if( debug ) print *, 'coord ',coord
                   write(uapr,'(2x,2a4,3(2x,f13.4),3(2x,f8.4),2x,f8.3)')
     .              sitcodx,'_GPS', (coord(i),i=1,6),apr_epoch 
                    backspace(usnx)
                  elseif( line(8:11).eq.'VELX' ) then
                    backspace(usnx)
                  endif
                elseif( crdtyp.eq.'VELX' ) then
                  coord(4) = coordx
                elseif( crdtyp.eq.'VELY' ) then
                  coord(5) = coordx
                elseif( crdtyp.eq.'VELZ' ) then
                  coord(6) = coordx
                  apr_epoch = decyrs(iy1,iday1,dfloat(isod1))
                  write(uapr,'(2x,2a4,3(2x,f13.4),3(2x,f8.4),2x,f8.3)')
     .              sitcodx,'_GPS', (coord(i),i=1,6),apr_epoch 
*               end if on crdtyp
                endif    
*             end if on error reading values from line
              endif 
*           end if on checking for block end
            endif
*         end do on reading through block
          enddo
*       end if on checking for beginning of block
        endif
*     end do on reading through the file
      enddo
                 

*     Now insert the values into the global arrays
             
      num_save = num_str
      do i=1,nss 
         num_str = num_str + 1 
         site_indx(num_str) = num_str
         newflag(num_str) = .true.
         sitcod(num_str) = sitcod1(i)
         sname(num_str) = sname1(i)  
         do j=1,5
           istart(j,num_str) = rcvr_start(j,i)
           istop(j,num_str) =  rcvr_stop(j,i)
         enddo
         rctype(num_str) = rctype1(i)
         rcvrsn(num_str) = rcvrsn1(i)
         rcvers(num_str) = rcvers1(i) 
         anttyp(num_str) = anttyp1(i)
         antsn(num_str)  = antsn1(i)
         radome(num_str) = radome1(i)
         dUNE(1,num_str) = anth1(i)
         dUNE(2,num_str) = antn1(i)
         dUNE(3,num_str) = ante1(i)
         dAntAZ(num_str) = antdaz1(i)
      enddo   

*     Get the GAMIT codes

      do i=num_save+1,num_str
        htcod(i) = 'DHARP'
        comment(i)=
     .      ' mstinf: '//sx_file(isnxf)(1:nblen(sx_file(isnxf)))   
        call read_rcvant( 2,2,char6,char20,radome(i)
     .                  ,rcvcod(i),rctype(i),pcncod )    
        anttyp20 = anttyp(i)//'     '
        call read_rcvant( 2,1,antcod(i),anttyp20,radome(i)
     .                  , char6,char20,pcncod )      
        if( radome(num_str).eq.'     ' ) radome(num_str)='UNKN'
c       use guess_rcvant.dat only for the firmware code 
        call guess_rcvant(2,ig_file(isnxf),char20,char20,char5,rcvers(i)
     .                  , char6, swver(i),char6,char132 )
c       remove ^s added by guess_rcvant 
        call sub_char ( rcvers(i),'^',' ')   
        if( debug ) then
          print *,'DEBUG  sitcod: ',sitcod(i) 
          print *,'       sname : ',sname(i)
          print *,'       rctype: ',rctype(i)
          print *,'       rcvcod: ',rcvcod(i) 
          print *,'       rcvers: ',rcvers(i)
          print *,'       swver : ',swver(i) 
          print *,'       anttyp: ',anttyp(i)
          print *,'       radome: ',radome(i)
          print *,'      comment: ',comment(i)
        endif
      enddo

*     Write an entry to the output apr file if requested
      if( write_apr) 
     .   write(uapr,'(2x,2a4,3(2x,f13.4),3(2x,f8.4),2x,f8.3)')
     .      sitcod(i),'_GPS', (coord(i),i=1,6),2008.00

      if( debug ) then  
        do i=1,num_str
          write(*,*) i,sitcod(i),sname(i),(istart(j,i),j=1,5)
     .       , (istop(j,i),j=1,5)
     .       , rctype(i),rcvrsn(i),rcvers(i),anttyp(i),antsn(i)
     .       , radome(i),(dUNE(j,i),j=1,3)
        enddo  
      endif

      return
      end      


