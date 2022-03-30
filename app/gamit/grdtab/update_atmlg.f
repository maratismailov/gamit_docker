       Program UPDATE_ATMLG

c      Program to update the ATML grid file header and to add additional data
c      beyond the end of the current file, both to extend the file of the 
c      current year (with incrementals from Tonie Van Dam) and to add extra
c      dupliate records for one day to allow interpolation. 
c      It currently converts Version 1.0 (10)  to Version 2.0 (20):
c                                                             
c          1. Modify Record 1 to replace 2-byte spare by number of header records
c
c          2. Modify Record 2 from I*4 ngval nglon nglat 
c                             to   R*4 interval + I*2 ngval nglon nglat
c
c          3. Add Record 3 w/ start/end times (actual, decimal, e.g. 365.75 for orginal files)
c
c          4. Optionally add extra records at the end

c      R. King 4 Oct 2006, from P. Tregoning programs fixatlg and extend_atm.
c      Last modified by R. King 13 July 2007 to read and add a supplementary file
                                                  


c  Command line:  update_atmlg [input-file] [add-file] [output-file] [yr_start] [doy_start] [yr_end] [doy_end]

c            where [input-file] is the input grid file (Version 1.0 or 2.0) 
c                  [add-file]   is a file whose values are added to the input file (Version 2.0)
c                  [output-file] is the output grid file (Version 2.)
c                  [yr_start]   is 4-digit start-year of input and desired output file
c                  [doy-start]  is  start-day of input and desired output file
c                  [yr_end]     is 4-digit year of desired last value]
c                  [doy_end]    is  stop-day of desired last value  
c                    If the yr_end, doy_end is omitted or equal to 9999, 999, then
c                    set the end epoch from the last epoch on the add-file.
c
c     Examples:  update_atmlg ATMDISP_cm.03 ' ' atmdisp_cm.2003 2003 1 2004 1 
c                update_atmlg ATMDISP_cm.07 ATMDISP_cm_suppl.188_07 atmdisp_cm.2007 2007 1 
c    

c   Notes: 1. File will be byte-swapped if necessary
c          2. If input is Version 2.0, no conversion of headers except change of end time
c          3. End day is now interpreted to mean a complete day, ending with 24:00 UTC,
c             so an old file ending with 365.75 will have appended one epoch; the
c             usual extension will be one full day (epochs)
                                     
      implicit none                
                                        
      
c Input/output files
      character*24 file_in,file_add,file_out
        
c Start/stop epochs for output file
      integer*4 iyrout1,iyrout2,idoyout1,idoyout2
          
c Record counters
      integer*4 irec1,irec2,irec3,nepoch_in,nepoch_add,nrec_in,nrec_add
     .        , nepoch_pad,nrec_pad
      real*4 days_in,days_add,days_pad,days_out
                
c  Variables for grid-file headers
      integer*2 version_in,version_add,version_out,nhead_in,nhead_out
     .      , i2val,i2lon,i2lat,recl
     .      , start_yr_in,end_yr_in,start_yr_add,end_yr_add
     .      , start_yr_out,end_yr_out
     .      , nhead_add,i2vala,i2lona,i2lata
      integer*4 ngval,nglon,nglat,ngvala,nglona,nglata   
      real*4 interval,start_doy_in,end_doy_in
     .     , intervala,start_doy_add,end_doy_add
     .     , start_doy_out,end_doy_out,values(3)
     .     , debugstart,debugstop        
      character*2 spare
      character*5 press_src,press_srca
      character*3 ref_frame,ref_framea

c  Other variables
      integer*4 ngrid,ioerr,iy1,iy2,i4y,i,j
      real*4 fversion,tmpdoy
      character*20 arg   
      logical swapped,swappeda,addfile,debug

c  External functions
      integer*4 nydays
      real*4 diffdays

      debug = .true.
      debugstart = 1.0
      debugstop = 367.0
      swapped = .false.  
      swappeda = .false.                
         
c   define the record length. For now, all atm files have a record length of 12
      recl=12

c  Read the command line

      call getarg(1,file_in)
      if(file_in(1:1).eq." ") then
        write(*,'(//,2a,//,2a,//,2a,/,a)')      
     .  "Run:  update_atmlg [in-file] [add-file] [out-file] [start yr] "
     .        ,"[start doy] [end yr] [end doy]"
     .   ," e.g. : update_atmlg ATMDISP_cm.03 ' ' atmdisp_cm.2003 "
     .        ,"2003 1 2004 1"
     .  ,"        update_atmlg ATMDISP_cm.166_07 ATMDISP_cm.166_188_07 "
     .           ,"atmdisp_cm.2007_188  2007 1 2007 188" 
     .  ,"If the dates are omitted, the dates on the files will be used"
        stop  
      endif                  
      file_add = '  ' 
      call getarg(2,file_add)   
      addfile = .false.
      if( file_add(1:1).ne.' ') addfile = .true.
c     trap case of old input with out-file second argument
      call getarg(3,file_out)
      read(file_out,'(i4)',iostat=ioerr) iyrout1
      if( iyrout1.ne.0.and.ioerr.eq.0) then
c        third argument is year instead of out-file, stop the user
         write(*,'(a)') 'Old command-line, need blank for add-file'
         stop
      endif                   
      iyrout1 = 9999
      idoyout1 = 999
      call getarg(4,arg)
      read(arg,'(i4)') iyrout1
      call getarg(5,arg)
      read(arg,'(i3)') idoyout1
c     initialize out epochs (overridden by command-line or add-file
c     header, then padded to an even day
      iyrout2 = 9999
      idoyout2 = 999
      call getarg(6,arg)
      if( arg(1:1).ne." ") then
        read(arg,'(i4)') iyrout2
        call getarg(7,arg)
        read(arg,'(i3)') idoyout2
      endif           


c  Open the original input file

      open(10,file=file_in,status='old',access='direct',recl=recl
     .        ,form='unformatted',iostat=ioerr)
      if(ioerr.ne.0)then
          print*,'Error opening file ',file_in,'. Does it exist?'
          stop "Program ended with errors"
      else              
         write(*,'(/,2a)') 'Opened input file ',file_in
      endif 

c  Read Record 1 of the original file:   
c
c  Version 1.0 or 2.0       version  spare pressure-source  reference-frame
c                             I*2     C*2     C*5              C*3
c  
c  Version 2.1 same except           nhead 
c                                     I*2

      read(10,rec=1,iostat=ioerr) version_in,spare,press_src,ref_frame 
      if(version_in.gt.100.or.version_in.lt.0)then
*         If the grid size does not seem be correct, swap the bytes
*         and see if OK.   
        write(*,'(a,i4,a)') '**version =',version_in,' swapping**'
        call swap_bytes(2,version_in,1)
*         Now see if OK size.
        if(version_in.gt.100.or.version_in.lt.0) then  
           write(*,'(a,i4,a)') 'Invalid version (',version_in
     .            ,') of ATML grid file '
*             Still a problem so kill.    
              stop
        else
*             File seems to be swapped.  Set the swapped status and continue.
*             (We now need to swap every record read)
              swapped = .true.    
        end if
      end if 
      if( version_in.eq.10 ) then    
        read(10,rec=1,iostat=ioerr) version_in,spare,press_src,ref_frame 
        if( swapped ) call swap_bytes(2,version_in,1)
        nhead_in= 2    
        write(*,'(/,a,1x,i4,1x,a5,1x,a3)')       
     .    'Record 1 Version 1.0  version pressure-source frame '
     .            ,version_in,press_src,ref_frame  
c*** temporary to fix Tonie's error
c           ref_frame = ' CM'
c           print *,"**WARNING: changing CF to CM in header to fix error"
c****************
        write(*,'(a)') '--assume nhead = 2'   
      elseif( version_in.eq.20 ) then
        read(10,rec=1,iostat=ioerr) version_in,spare,press_src,ref_frame 
        if( swapped ) call swap_bytes(2,version_in,1)
        nhead_in = 3      
        write(*,'(/,a,1x,i4,1x,a5,1x,a3)') 
     .    'Read Record 1 Version 2.0  version pressure-source frame '
     .            ,version_in,press_src,ref_frame 
        write(*,'(a)') '--assume nhead = 3'
      elseif( version_in.eq.21 ) then    
        read(10,rec=1,iostat=ioerr) 
     .      version_in,nhead_in,press_src,ref_frame
        if( swapped ) then
          call swap_bytes(2,version_in,1)
          call swap_bytes(2,nhead_in,1)  
        endif          
        write(*,'(a,2i3,1x,a5,1x,a3)')  
     . 'Read Record 1 Version 2.1  version nhead pressure-source frame '
     .            ,version_in,nhead_in,press_src,ref_frame  
      else
         write(*,'(a,i2)') 'Unknown version ',version_in
         stop
      endif
       

c  Read Record 2:  

c     Version 1.0 :  ngval   nglon   nglat
c                     i4       i4     i4

c     Version 2.0ff :  interval  ngval   nglon   nglat
c                        r4        i2      i2      i2

           
      if( version_in.eq.10 ) then 
        read(10,rec=2,iostat=ioerr) ngval,nglon,nglat
        if(swapped)then
            call swap_bytes(4,ngval,1)
            call swap_bytes(4,nglon,1)
            call swap_bytes(4,nglat,1)
        endif                                             
        interval = 0.25  
        i2val =  ngval  
        i2lon =  nglon  
        i2lat =  nglat         
        write(*,'(/,a,3i5)')
     .      'Read Record 2 Version 1.0 ngval nglon nglat'
     .                     , ngval,nglon,nglat        
        write(*,'(a)') '--assume interval = 0.25 d'

      elseif( version_in.ge.20 ) then  
        read(10,rec=2,iostat=ioerr) interval,i2val,i2lon,i2lat 
        if( ioerr.ne.0 ) then
           write(*,'(a)') 'Error reading record 2 of grid file '
           stop              
        endif          
cd        print *,'DEBUG swapped i2val i2lon i2lat ',swapped,
cd     .    i2val,i2lon,i2lat
        if(swapped)then
            call swap_bytes(4,interval,1)
            call swap_bytes(2,i2val,1)
            call swap_bytes(2,i2lon,1)
            call swap_bytes(2,i2lat,1)   
        endif     
        ngval = i2val
        nglon = i2lon
        nglat = i2lat
        fversion = float(version_in)/10.d0                 
        write(*,'(/,a,f3.1,a,1x,f4.2,1x,3i5)')
     .   'Read Record 2 Version ',fversion,' interval ngval nglon nglat'
     .                     ,interval,ngval,nglon,nglat        
      endif   


c  Read Record 3 (absent in Version 1.0)
c
c     Version 2.0     yr_start  doy_start   yr_end   doy_end
c                       i2        r4          i2       r4 
                    
       if( version_in.eq.10 ) then     
c        assign time variables: assume day 1.0 thru day 325.75 
         start_yr_in = mod(iyrout1,100)
         start_doy_in = idoyout1
         end_yr_in = start_yr_in   
         days_in = nydays(iyrout1)
         end_doy_in = days_in + 0.75   
         interval = 0.25
         write(*,'(a,/,a,i2,1x,f6.2,/,a,1x,i2,1x,f6.2,/a,f4.2)')
     .     'Version 1.0 assume input:'
     .      , '  Start day = ',start_yr_in,start_doy_in
     .      , '  End day = ',end_yr_in,end_doy_in
     .      , '  Interval = ',interval   
 
      elseif( version_in.ge.20 ) then
        read(10,rec=3,iostat=ioerr) 
     .        start_yr_in,start_doy_in,end_yr_in,end_doy_in
c**Temporary to fix Tonie's error for the first file of 2011
c        start_yr_in = 11 
c ***end temp
        if( ioerr.ne.0 ) then
          write(*,'(a)') 'Error reading record 3 of grid file'
          stop            
        endif
        if( swapped ) then
           call swap_bytes(2,start_yr_in,1)
           call swap_bytes(4,start_doy_in,1)
           call swap_bytes(2,end_yr_in,1)
           call swap_bytes(4,end_doy_in,1) 
        endif  
        write(*,'(a,/,a,i2,1x,f6.2,/,a,1x,i2,1x,f6.2,/a,f4.2)')
     .     'Read Record 3 Version 2.0 :'
     .      , '  Start day = ',start_yr_in,start_doy_in
     .      , '  End day = ',end_yr_in,end_doy_in
     .      , '  Interval = ',interval   
       else
        write(*,'(a,i7)') 'Unknown version ',version_in
        stop  
      endif      


c  Open the supplementary file and read the header to check consistency  

      if( addfile ) then
        open(15,file=file_add,status='old',access='direct',recl=recl
     .        ,form='unformatted',iostat=ioerr)
        if(ioerr.ne.0)then
          print*,'Error opening add-file ',file_add,'. Does it exist?'
          stop "Program ended with errors"
        else              
          write(*,'(/,2a,/)') 'Opened add file ',file_add
        endif 
        read(15,rec=1,iostat=ioerr) 
     .     version_add,spare,press_srca,ref_framea 
        if(version_add.gt.100.or.version_in.lt.0) then
*         If the grid size does not seem be correct, swap the bytes
*         and see if OK.   
          write(*,'(a,i4,a)') 'Add-file **version =',version_add
     .              ,' swapping**'
          call swap_bytes(2,version_add,1)
*         Now see if OK size.
          if(version_add.gt.100.or.version_add.lt.0) then  
             write(*,'(a,i4,a)') 'Invalid add-file version ('
     .         ,version_add,') of ATML grid file '
*               Still a problem so kill.    
               stop
          else
*             File seems to be swapped.  Set the swapped status and continue.
*             (We now need to swap every record read)
              swappeda = .true.    
          endif
        endif             
        if( version_add.lt.20 ) then
          write(*,'(a)') 
     .    'Supplementary file must be version 2.0 or higher'
         stop
        elseif( version_add.eq.20.or.version_add.eq.21 ) then    
          read(15,rec=1,iostat=ioerr) 
     .      version_add,nhead_add,press_srca,ref_framea
          if( swappeda ) then
            call swap_bytes(2,version_add,1)
            if( version_add.ge.21 ) then
              call swap_bytes(2,nhead_add,1)   
            else
              nhead_add = 3
            endif
          endif  
          if(press_srca.ne.press_src .or. ref_framea.ne.ref_frame) then
            write(*,'(a,2(a5,1x,a3,1x,a5,1x,a3))') 
     .       'In-file and add-file models different: '
     .       ,press_src,ref_frame,press_srca,ref_framea
            stop
          endif        
          write(*,'(2a,2i3,1x,a5,1x,a3)')  
     .      'Read Add-file Record 1 Version 2.1 '
     .      , 'version nhead pressure-source frame '
     .      ,version_add,nhead_add,press_srca,ref_framea  
        else
          write(*,'(a,i2)') 'Unknown add-file version ',version_add
          stop
        endif
c       Record 2 (Ver 2.1):  interval  ngval   nglon   nglat
c                            r4        i2      i2      i2
        read(15,rec=2,iostat=ioerr) intervala,i2vala,i2lona,i2lata 
        if( ioerr.ne.0 ) then
          write(*,'(a)') 'Error reading record 2 of add-file '
          stop              
        endif
        if( swappeda )then
          call swap_bytes(4,intervala,1)
          call swap_bytes(2,i2vala,1)
          call swap_bytes(2,i2lona,1)         
          call swap_bytes(2,i2lata,1)  
        endif  
        ngvala = i2vala
        nglona = i2lona
        nglata = i2lata
        fversion = float(version_add)/10.d0                 
        write(*,'(/,a,f3.1,a,1x,f4.2,1x,3i5)')
     .  'Read Record 2 Version',fversion,' interval ngval nglon nglat'
     .                   ,intervala,ngvala,nglona,nglata        
c       Record 3:  Version 2.0     yr_start  doy_start   yr_end   doy_end
c                         i2        r4          i2       r4 
        read(15,rec=3,iostat=ioerr) 
     .        start_yr_add,start_doy_add,end_yr_add,end_doy_add
        if( ioerr.ne.0 ) then
          write(*,'(a)') 'Error reading record 3 of add-file'
          stop            
        endif
        if( swappeda ) then
          call swap_bytes(2,start_yr_add,1)
          call swap_bytes(4,start_doy_add,1)
          call swap_bytes(2,end_yr_add,1)
          call swap_bytes(4,end_doy_add,1) 
        endif   
        write(*,'(a,/,a,i2,1x,f6.2,/,a,1x,i2,1x,f6.2,/a,f4.2)')
     .     'Read Record 3 Version 2 :'
     .      , '  Start day = ',start_yr_add,start_doy_add
     .      , '  End day = ',end_yr_add,end_doy_add
     .      , '  Interval = ',intervala   
        if( nglata.ne.nglat.or.nglona.ne.nglon ) then
          write(*,'(a,2i4,a,2i4)') 'Add-file grid size (',nglata,nglona
     .        ,') different from in-file (',nglat,nglon
          stop   
        endif  
        if( end_doy_in.lt.start_doy_add ) then
          write(*,'(/,a,f6.2,a,f6.2,a)') '**End day of input ('
     .      ,end_doy_in,') < start day of add-file (',start_doy_add
     .      ,')  stop'
          stop
        endif
      endif


c Open the output file

      open(20,file=file_out,status='unknown',access='direct',recl=recl
     .        ,form='unformatted',iostat=ioerr)
      if(ioerr.ne.0)then
          print*,'Error opening output file ',file_out
          stop "Program ended with errors"
      else                
         write(*,'(/,2a)') 'Opened output file ',file_out
      endif

c  Define the grid size (# lat/lon pts per epoch)

      ngrid = nglat*nglon


c  Set the output epochs

      start_yr_out = start_yr_in
      start_doy_out = start_doy_in 
c     if there is an add-file, set the end-day from its header
      if( addfile ) then
        end_yr_out = end_yr_add
        end_doy_out = end_doy_add    
cd        print *,'DEBUG end_doy_out 2 from add-header ',end_doy_out
c     otherwise, use the input values
      else
        end_yr_out = end_yr_in
        end_doy_out = end_doy_in   
      endif
c     in either case, pad to the end of the day
cd      print *,'DEBUG end_doy_out amod nrec_pad '
cd     .       , end_doy_out,amod(end_doy_out,1.0),nrec_pad
      if( amod(end_doy_out,1.0).ne.0. ) then
        tmpdoy = end_doy_out
        end_doy_out = end_doy_out + (1. - amod(end_doy_out,1.0))
        days_pad = end_doy_out - tmpdoy 
cd        print *,'end_doy_out days_pad ',end_doy_out,days_pad                                                                       
        nepoch_pad = nint(days_pad/interval) 
        nrec_pad = nepoch_pad*ngrid
cd        print *,'days_pad nepoch_pad nrec_pad '
cd     .    ,days_pad,nepoch_pad,nrec_pad
      else
        nepoch_pad = 0
      endif     
c     override with input if specified and within the available grids
      if( iyrout1.ne.9999 ) then 
        i4y = start_yr_out 
        if( mod(iyrout1,100).eq.i4y .and. 
     .      idoyout1.gt.nint(start_doy_out) ) then
          start_doy_out = idoyout1
          write(*,'(a,i3,f6.1,a)') '** Output start epoch set to '
     .        ,start_yr_out,start_doy_out,' by command-line value'
        endif
      endif     
      if( iyrout2.ne.9999 ) then  
         i4y = end_yr_out             
         if( mod(iyrout2,100).eq.i4y .and. 
     .       idoyout2.lt.nint(end_doy_out) ) then
           end_doy_out = idoyout2 + 1.0
           write(*,'(a,i3,f6.1,a)') '** Output end epoch set to '
     .        ,end_yr_out,end_doy_out,' by command-line value'
         endif
       endif


c  Write the output records
                  
c   Record 1:
      version_out = 21     
      nhead_out = 3          
      write(*,'(/,a)') 'Headers of output file: '
      write(*,'(a,2i3,1x,a5,1x,a3)')
     .     '  Record 1: version nhead press_src frame '
     .    , version_out,nhead_out,press_src,ref_frame   
      write(20,rec=1) version_out,nhead_out,press_src,ref_frame
 

c   Record 2:
      write(*,'(a,1x,f4.2,1x,3i5)')
     .   '  Record 2: interval ngval nglon nglat'
     .                ,interval,ngval,nglon,nglat 
      write(20,rec=2) interval,i2val,i2lon,i2lat
                                                                
c   Record 3:
      write(*,'(a,2(1x,i2,1x,f6.2))') '  Record 3: '
     .       ,start_yr_out,start_doy_out,end_yr_out,end_doy_out
      write(20,rec=3) start_yr_out,start_doy_out,end_yr_out,end_doy_out

c  Calculate the records for the input and output file      
         
      iy1 = start_yr_in
      call fix_y2k(iy1)
      iy2 = end_yr_in
      call fix_y2k(iy2)
      days_in= diffdays(iy1,start_doy_in,iy2,end_doy_in)
cd      print *,' '
cd      print *,'DEBUG In: yr doy yr doy days '
cd     .    ,iy1,start_doy_in,iy2,end_doy_in,days_in      
      iy1 = start_yr_out
      call fix_y2k(iy1)
      iy2 = end_yr_out
      call fix_y2k(iy2)
      days_out= diffdays(iy1,start_doy_out,iy2,end_doy_out)
cd      print *,'DEBUG Out: yr doy yr doy days '
cd     .    ,iy1,start_doy_out,iy2,end_doy_out,days_out  
      nepoch_in = nint(days_in/interval) + 1 
      nrec_in  = nhead_in +  ngrid*nepoch_in  
cd      print *,'nepoch_in nrec_in ',nepoch_in,nrec_in
      if( addfile ) then  
cd        print *,'start_yr_add,start_doy_add,end_yr_add,end_doy_add '
cd    .         , start_yr_add,start_doy_add,end_yr_add,end_doy_add   
        if( start_doy_add.gt.(end_doy_in+interval) ) then
          write(*,'(a,f8.2,a,f8.2,a,f5.2,a)') 'First added epoch ('
     .       ,start_doy_add,') > last in-file epoch (',end_doy_in
     .       ,') + interval (',interval,')'
          stop
        endif
        iy1 = start_yr_add
        call fix_y2k(iy1)
        iy2 = end_yr_add
        call fix_y2k(iy2)
        days_add = 
     .    diffdays( iy1,start_doy_add,iy2,end_doy_add )
        nepoch_add = nint(days_add/interval) + 1
        nrec_add = ngrid*nepoch_add 
cd        print *,'days_add,nepoch_add,nrec_add '
cd     .      ,days_add,nepoch_add,nrec_add    
      endif


c  Read and write the values records from the in_file

      do i = nhead_in+1,nrec_in                 
        irec1 = i
        read(10,rec=irec1,iostat=ioerr) values   
        if( ioerr.ne.0) call grderr(10,irec1,ioerr)  
        if( swapped ) call swap_bytes(4,values,3)   
cd        if( irec1.lt.10 ) print *,'READ irec1 ',irec1,values
        irec3= irec1 + (nhead_out-nhead_in)
        write(20,rec=irec3) values           
cd        if( irec3.lt.10) print *,'WROTE irec3 ',irec3,values   
cd        if( i.gt.10 ) stop
      enddo    
cd      print *,'DEBUG copied in-file, irec1 irec3 ',irec1,irec3
      
c  Add records from the supplementary file and/or 
c  pad with last values for interpolation
             
      if( addfile ) then     
        irec2 = nhead_add + 1
        call grid_time('T',start_doy_in,interval,ngrid,nhead_in
     .    ,start_doy_add,irec3)
cd        print *,'start_doy_add irec2 irec3 ',start_doy_add,irec2,irec3
        do i = 1, nrec_add  
          read(15,rec=irec2,iostat=ioerr) values 
	       if( ioerr.ne.0 ) call grderr(15,irec2,ioerr)          
	       if( swappeda ) call swap_bytes(4,values,3)
cd          if( irec2.eq.4.or.
cd     .        irec2.eq.10516.or.
cd     .        irec2.eq.21028 ) then
cd              print *,'READ irec2',irec2,values
cd          endif
          write(20,rec=irec3) values     
cd          if( irec3.eq.8283460.or.
cd     .      irec3.eq.8293972.or.
cd     .      irec3.eq.8304484 ) then
cd            print *,'WROTE irec3 ',irec3,values
cd          endif 
          irec2 = irec2 + 1  
          irec3 = irec3 + 1
        enddo    
      endif
      print *,'Suppl records added, irec2 irec3 ',irec2,irec3
      do i = 1, nepoch_pad  
c       reset input counters to first record of epoch
        irec1 = irec1 - ngrid
        irec2 = irec2 - ngrid
        print *,'pad epoch irec1 irec2 ',i,irec1,irec2
        do j=1,ngrid 
          if( addfile ) then
            read(15,rec=irec2,iostat=ioerr) values
	         if( ioerr.ne.0 ) call grderr(15,irec2,ioerr)
	         if( swappeda ) call swap_bytes(4,values,3)  
    	     irec2 = irec2 + 1
          else    
            read(10,rec=irec1,iostat=ioerr) values  
            if( ioerr.ne.0 ) call grderr(10,irec1,ioerr)  
            if( swapped ) call swap_bytes(4,values,3)
          endif
          write(20,rec=irec3) values
          irec3 = irec3 + 1 
          if(j.eq.1) print *,'wrote pad i irec3 values ',i,irec3,values
        enddo
      enddo
      	

c  Close the files and write what we did

      close(10,iostat=ioerr)
      if( ioerr.ne. 0)  
     .  write(*,'(a,i4)') 'Error closing input file, iostat = ',ioerr  
      if( addfile) then
        close(15,iostat=ioerr)
         if( ioerr.ne. 0)  
     .    write(*,'(a,i4)') 'Error closing input file, iostat = ',ioerr  
      endif  
      close(20,iostat=ioerr)
      if( ioerr.ne. 0)  
     .  write(*,'(a,i4)') 'Error closing output file, iostat = ',ioerr 
      if( addfile ) write(*,'(a,i10,a)')
     .  'Wrote ',nrec_add,' data records from add_file to end '
      write(*,'(a,i10,a)') 
     .  'Padded with ',nrec_pad,' duplicate records'          
    
      stop
      end       
     
c--------------------------------------------------------

      function diffdays(iyr1,doy1,iyr2,doy2)

      integer*4 iyr1,iyr2,idoy1,idoy2
     .        , im1,id1,im2,id2,jd1,jd2,julday 
      real*4 doy1,doy2,fjd1,fjd2,diffdays     

c      print *,'DIFFDAYS iyr1 doy1 iyr2 doy2 '
c     .       , iyr1,doy1,iyr2,doy2
                   
      idoy1 = doy1
      call monday(idoy1,im1,id1,iyr1)
      jd1 = julday(im1,id1,iyr1) 
      fjd1 = float(jd1) + (doy1-float(idoy1)) 
c      print *,'DIFFDAYS im1 id1 jd1 fjd1 '
c     .              ,   im1,id1,jd1,fjd1
      idoy2 = doy2
      call monday(idoy2,im2,id2,iyr2)
      jd2 = julday(im2,id2,iyr2) 
      fjd2 = float(jd2) + (doy2-float(idoy2))   
c      print *,'DIFFDAYS im2 id2 jd2 fjd2 '
c     .              ,   im2,id2,jd2,fjd2
      diffdays = fjd2 -fjd1 
      return
      end
      
c----------------------------------------------------------------------

      subroutine grderr(unit,record,err)
                       
      integer*4 unit,record,err
      
      write(*,'(a,2i3,i10)') 'Read error (unit,err,record)'     
     .     , unit,err,record
      stop
      end
             
c---------------------------------------------------------------------

      Subroutine fix_y2k(iyear)

c     Check for 2 or 4-digit years.  If 2-digits, convert to 4 by guessing
c     the date.   Same as check_y2k except no warning printed.  R. King July 1999

      implicit none

      integer*4 iyear

      if( iyear.le.1900) then 
                  
c        earliest GPS launch date is 1978; earliest space-geodetic data about 1960 
         if( iyear.gt.60.and.iyear.le.99 ) then
            iyear = iyear + 1900
         else
            iyear = iyear + 2000
         endif  

      endif

      return
      end

c------------------------------------------------------------------------------------









