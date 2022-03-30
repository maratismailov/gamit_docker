       Program UPDATE_VMF3G

c      Program to create a VMF3 grid file or update it with  additional data
c      (and eventually to update the header for new versions.
c
c      M. Floyd  October 2020 from update_vmfg; last updated 28 October 2020


c  Command line:  update_vmf3g [input-file] [add-file] [output-file] [yr_start] [doy_start] [yr_end] [doy_end]

c            where [input-file] is the input binary, direct-access grid file
c                               for a new file, use 'new' (blank will give help-file)
c                  [add-file]   two options:
c                     if [add-file] begins with 'vmf', then add this binary, direct-access
c                        file to the input file
c                     if [add-file] is a number, interpret as a [yyyymmdd], and add all ascii
c                        files of the form all VMF3_[yyyymmdd].H[HH], nominally 4 (epochs)
c                  [output-file] is the output grid file
c                  [yr_start]   is 4-digit start-year of input and desired output file
c                  [doy-start]  is  start-day of input and desired output file
c                  [yr_end]     is 4-digit year of desired last value]
c                  [doy_end]    is  stop-day of desired last value
c
c                Note: For ascii files, the program currently works for only one day
c                      at a time, so that the date is taken from the [add-file] entry
c                      and and the last four arguments (start/stop) should be omitted.
c
c                      Currently untested for binary add-files.
c
c     Examples:  update_vmf3g vmf3grd.2007_155 vmf3grd.2007_154_160 vmf3grd.2007_160
c                update_vmf3g vmf3grd.2007_155 20070605 vmf3grd.2007_156
c
c

      implicit none

c Input/output files
      character*24 file_in,file_add,file_out

c Start/stop epochs for output file
      integer*4 iyrout1,iyrout2,idoyout1,idoyout2

c Record counters
      integer*4 irec1,irec2,irec3,nepoch_in,nepoch_add,nrec_in,nrec_add
     .        , nepoch_pad,nrec_pad,nrec_orostr
      real*4 days_in,days_add,days_pad,days_out

c  Maximum grid and file sizes

      integer*4 maxlat,maxlon,maxfil
      parameter(maxlat=180,maxlon=360,maxfil=4)

c  Variables from the VMF3 files
      integer*4 ah,ahval(maxfil,maxlon,maxlat)
     .         ,aw,awval(maxfil,maxlon,maxlat)
     .         ,iorog
      real*4    zh,zhval(maxfil,maxlon,maxlat)
     .         ,zw,zwval(maxfil,maxlon,maxlat)


c  Variables from the header records
      character*2 source,sourcea
      character*4 mapf,mapfa
      integer*4 nglon,nglat,nglona,nglata
      integer*2 version,intervalhr,intervalhra,ngval,ngvala
      real*4 start_time,stop_time,start_timea,stop_timea

c  Variables for grid-file headers
      integer*2 version_in,version_add,version_out,nhead_in,nhead_out
     .      , recl
     .      , start_yr_in,stop_yr_in,start_yr_add,stop_yr_add
     .      , start_yr_out,stop_yr_out
     .      , nhead_add
      real*4 intervald,start_doy_in,stop_doy_in
     .     , start_doy_add,stop_doy_add
     .     , start_doy_out,stop_doy_out
     .     , start_mmdd_add,stop_mmdd_add
     .     , debugstart,debugstop

c  Variables for the ascii orography file
      integer*4 maxrow
      parameter( maxrow=145)
      integer*4 ival3(maxrow),lrow

c  Unit numbers
      integer*4 iuin,iuorog,iuadd,iuout

c  Other variables
      integer*4 ngrid,ioerr,iy1,iy2,i4y,i,j,k
     .        , ifile,it,ilat,ilon,yr4
      real*4 alat1,alat2,alon1,alon2,adlat,adlon,tmp
CD    real*4 tmpdoy
      character*1 addtype
      character*5 tmpname
      character*20 arg
      character*80 line
      logical swapped,swappeda,addfile,debug,newfile


c  External functions
      integer*4 nydays
      real*4 diffdays

c  Initialization

      iuin = 10
      iuorog = 11
      iuadd = 15
      iuout = 20
      debug = .true.
      debugstart = 1.0
      debugstop = 367.0
      swapped = .false.
      swappeda = .false.
c     dummy to avoid compiler warning:
      it=1
      nhead_add = 3
      addfile = .false.
      newfile = .false.


c   define the record length. For now, all VMF3 files have a record length of 52
      recl=52

c  Read the command line

      call getarg(1,file_in)
      if(file_in(1:1).eq." ") then
        write(*,'(//,2a,//,2a,//,2a,/,a)')
     .  "Run:  update_vmf3g [in-file] [add-file] [out-file] [start yr] "
     .  ,       "[start doy] [end yr] [end doy]"
     .  ,     "where [add-file] can be a binary file name or the date "
     .  ,       "(YYDDD) of Vienna ascii files "
     .  ,       "(e.g., VMF3_20070605.H00)"
     .  , " e.g. : update_vmfg vmf1grd.2007_155 vmf1grd.2007_154_160 "
     .  ,     " vmf1grd.2007_160"
     .  , "        update_vmfg vmf1grd.2007_155 07156 vmf1grd.2007_156"
     .  ,"If the dates are omitted, the dates on the files will be used"
        stop
      endif
      if( file_in(1:3).eq.'new' ) then
        newfile = .true.
      endif
      file_add = '  '
      call getarg(2,file_add)
      addfile = .false.
      if( file_add(1:1).eq.' ') then
         addfile = .false.
      elseif ( file_add(1:3).eq.'vmf'.or.file_add(1:3).eq.'VMF') then
         addfile = .true.
         addtype = 'B'
         write(*,'(2a)') 'Adding values from grid file ',file_add
      else
         read(file_add(1:2),'(i2)',iostat=ioerr) start_yr_add
         if(ioerr.ne.0 ) then
           write(*,'(a)') 'Error reading stop_yr_out from file_add name'
           stop
         endif
         yr4 = start_yr_add
         call fix_y2k(yr4)
         start_yr_add = yr4
         read(file_add(3:6),'(f4.0)',iostat=ioerr) start_mmdd_add
         call ymd_to_doy(start_mmdd_add,start_doy_add)
cd         print *,'DEBUG file_add read start_yr_add start_doy_add '
cd     .      ,file_add,start_yr_add,start_doy_add
          if(ioerr.ne.0 ) then
           write(*,'(a)')'Error reading stop_doy_out from file_add name'
           stop
         endif
         if( ioerr.eq.0 ) then
           addfile = .true.
           addtype = 'A'
           stop_yr_add = start_yr_add
           stop_doy_add = start_doy_add + 0.75
           write(*,'(/,a,2i4)')
     .        'Adding ascii files for 0h 6h 12h 18h on day '
     .         ,stop_yr_add,int(start_doy_add)
c          last epoch will be 18h on day of add-file
cd           print *,'DEBUG  stop_doy_add ',stop_doy_add
         else
           write(*,'(3a)') 'Add-file name (',file_add
     .         ,') not vmf/VMF or integer yyddd'
           stop
         endif
      endif
      call getarg(3,file_out)
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

      if( .not.newfile ) then
        open(iuin,file=file_in,status='old',access='direct',recl=recl
     .        ,form='unformatted',iostat=ioerr)
        if(ioerr.ne.0)then
          write(*,'(2a)') 'Error opening file ',file_in
          stop "Program ended with errors"
        else
          write(*,'(/,2a)') 'Opened input file ',file_in
        endif

c  Read the headers of the original file:
c
c    Record 1:
c
c    ! Version:            1.0
c                      R*4     C*4     C*2       C*4
c
c             e.g  10VMF1JB (version 1.0, VMF1, source: JB = Johannes Boehm)

        read(iuin,rec=1,iostat=ioerr) version_in,mapf,source
        if(version_in.lt.10.or.version_in.gt.100) then
          call swap_bytes(2,version,1)
          if(version.lt.10.or.version.gt.100) then
            write(*,'(2a)') 'Byte-swapping input grid file ',file_in
            swapped = .true.
          else
            write(*,'(a)') 'Invalid values read from input file header'
            stop
          endif
         else
           swapped = .false.
         endif
         write(*,'(/,a,1x,i2,1x,a4,1x,a2)')
     .      'Record 1 Version / Source ',version_in,mapf,source
         nhead_in = 3

c    Record 2:
c
c    Version 1.0   interval  #values  #lon   #lat
c                    I*2       I*2     I*4    I*4
c            e.g      6         3      144    91
c                 (6h intervals, 3 values, long, lat at 2.5-deg spacing)
        read(iuin,rec=2,iostat=ioerr) intervalhr,ngval,nglon,nglat
        if(swapped)then
          call swap_bytes(2,intervalhr,1)
          call swap_bytes(2,ngval,1)
          call swap_bytes(4,nglon,1)
          call swap_bytes(4,nglat,1)
        endif
        write(*,'(a,4i4)')
     .      'Read Record 2 interval #val #lon #lat '
     .        ,intervalhr,ngval,nglon,nglat
         intervald = float(intervalhr)/24.

c     Record 3:

c     Version 1.0    start       stop
c                     R*4         R*4
c             e.g.   2006.0      2006.2157
c              ( 2007 day 1.0 thru day 78.75)
        read(iuin,rec=3) start_time,stop_time
        if(swapped)then
          call swap_bytes(4,start_time,1)
          call swap_bytes(4,stop_time,1)
        endif
        write(*,'(a,2f12.5)')
     .      'Read Record 3 start stop  ',start_time,stop_time
        start_yr_in = int(start_time)
        start_doy_in=amod(start_time,1.0)*float(nydays(start_yr_out))+1.0
        call round6h(start_doy_in)
        stop_yr_in = int(stop_time)
        stop_doy_in=amod(stop_time,1.0)*float(nydays(stop_yr_out)) + 1.0
        call round6h(stop_doy_in)
        write(*,'(22x,a,2(i6,f7.2),a)') '( = ',start_yr_in,start_doy_in,
     .      stop_yr_in,stop_doy_in,' )'
      else
        write(*,'(a)') 'No grid input file, create a new one'
      endif


c  If the supplementary file is a binary grid, open it and check the headers

      if( addfile.and.addtype.eq.'B' ) then
        open(iuadd,file=file_add,status='old',access='direct',recl=recl
     .        ,form='unformatted',iostat=ioerr)
        if(ioerr.ne.0)then
          write(*,'(2a)') 'Error opening binary add-file ',file_add
          stop "Program ended with errors"
        else
          write(*,'(/,2a,/)') 'Opened binary add file ',file_add
        endif
        read(iuadd,rec=1,iostat=ioerr) version_add,mapfa,sourcea
        if(version_add.lt.10.or.version_add.gt.100)then
        call swap_bytes(2,version_add,1)
        if(version_add.ge.10.or.version_add.lt.100) then
           write(*,'(a)') 'Byte-swapping supplementary grid file'
           swappeda = .true.
        else
          write(*,'(a)') 'Invalid values read from add-file'
          stop
        endif
       else
         swappeda = .false.
       endif
       write(*,'(a,1x,i4,1x,a4,1x,a2)')
     .    'Record 1 Version / Source ',version_add,mapfa,sourcea
       read(iuadd,rec=2,iostat=ioerr) intervalhra,ngvala,nglona,nglata
       if(swapped)then
         call swap_bytes(2,intervalhra,1)
         call swap_bytes(2,ngvala,1)
         call swap_bytes(4,nglona,1)
         call swap_bytes(4,nglata,1)
       endif
       if( version_add.eq.10 ) nhead_add = 3
       write(*,'(/,a,4i4)')
     .      'Read Record 2 interval #val #lon #lat '
     .        ,intervalhra,ngvala,nglona,nglata
       if( intervalhra.ne.intervalhr ) then
         write(*,'(a,2i3)') 'Interval of binary add-file NE input file '
     .        ,intervalhra,intervalhr
         stop
       endif
       read(iuin,rec=3) start_timea,stop_timea
       if( swapped ) then
         call swap_bytes(4,start_timea,1)
         call swap_bytes(4,stop_timea,1)
        endif
        write(*,'(/,a,4i4)')
     .      'Read Record 3 start stop  ',start_timea,stop_timea
        start_yr_add = int(start_timea)
        start_doy_add =
     .      amod(start_timea,1.0)*float(nydays(start_yr_out)) + 1.0
        call round6h(start_doy_add)
        stop_yr_add = int(stop_timea)
        stop_doy_add =
     .      amod(stop_timea,1.0)*float(nydays(stop_yr_out)) + 1.0
        call round6h(stop_doy_add)
        write(*,'(22x,a,2(i6,f7.2,a))') '( = '
     .     ,start_yr_add,start_doy_add,stop_yr_add,stop_doy_add,' )'
      endif
cd      print *,'DEBUG 2 start_yr_add ',start_yr_add


c If the supplementary files are ascii, open them sequentially and store the values

c     ahval(it,ilon,ilat)
c     awval(it,ilon,ilat)
c     zhval(it,ilon,ilat)
c       ifile =1-4  ah  (i*4)
c              5-8  aw  (i*4
c              9-12 zh  (r*4)
c       it=1-4  0h 6h 12h 18h
      if( addfile.and.addtype.eq.'A' ) then
        intervalhra = 6
        ngvala = 3
        tmpname = file_add(1:5)
        file_add(3:7)= tmpname
        file_add(12:24) = ' '
        read(tmpname(1:2),'(i2)') start_yr_add
        yr4 = start_yr_add
        call fix_y2k(yr4)
        start_yr_add = yr4
        read(tmpname(3:5),'(f3.0)') start_doy_add
        start_timea = float(start_yr_add)
     ,       + (start_doy_add -1.)/float(nydays(start_yr_add))
cd        print *,'tmpname start_yr_add start_doy_add start_timea'
cd     .         , tmpname,start_yr_add,start_doy_add,start_timea
        mapfa = "VMF3"
        sourcea = "DL"
        do ifile=1,12
          if( ifile.le.4 ) then
            file_add(1:2) = 'ah'
          elseif( ifile.le.8 ) then
            file_add(1:2) = 'aw'
          else
            file_add(1:2) = 'zh'
          endif
          it = mod(ifile,4)
          if( it.eq.1 ) file_add(8:11) = '.h00'
          if( it.eq.2 ) file_add(8:11) = '.h06'
          if( it.eq.3 ) file_add(8:11) = '.h12'
          if( it.eq.0 ) file_add(8:11) = '.h18'
          open(iuadd,file=file_add,status='old',iostat=ioerr)
          if(ioerr.ne.0)then
            write(*,'(2a)') 'Error opening add-file ',file_add
            stop
          else
           write(*,'(/,2a)') 'Opened ascii add-file ',file_add
          endif
          read(iuadd,'(1x,6f8.2)',iostat=ioerr)
     .      alat1,alat2,alon1,alon2,adlat,adlon
          write(*,'(a,6f8.2)') 'Header lat1 lat2 lon1 lon2 dlat dlon '
     .                , alat1,alat2,alon1,alon2,adlat,adlon
          nglata = nint((alat1-alat2)/adlat) + 1
          nglona = nint((alon2-alon1)/adlon) + 1
          do ilat = 1,nglata
            if( file_add(1:2).eq.'ah' ) then
              read(iuadd,'((12i7))',iostat=ioerr)
     .            (ahval(ifile,ilon,ilat),ilon=1,nglona)
cd              if( ilat.le.2.or.ilat.eq.91 )
cd     .          print *,'ifile ilat ahval(1-2) ',ifile,ilat
cd     .               ,ahval(ifile,1,ilat),ahval(ifile,2,ilat)
            elseif( file_add(1:2).eq.'aw') then
              read(iuadd,'((12i7))',iostat=ioerr)
     .           (awval(ifile,ilon,ilat),ilon=1,nglona)
cd              if( ilat.le.2.or.ilat.eq.91 )
cd     .          print *,'ifile ilat awval(1-2) ',ifile,ilat
cd     .               ,awval(ifile,1,ilat),awval(ifile,2,ilat)
            elseif( file_add(1:2).eq.'zh' ) then
              read(iuadd,'((12f7.4))',iostat=ioerr)
     .           (zhval(ifile,ilon,ilat),ilon=1,nglona)
cd              if( ilat.le.1.or.ilat.eq.91 )
cd     .           print *,'ifile ilat zhval(1-2) ',ifile,ilat
cd     .               ,zhval(ifile,1,ilat),zhval(ifile,2,ilat)
            else
              write(*,'(a,a2,a)') 'Value type ',file_add(1:2)
     .             ,' not recognized '
              stop
            endif
            if( ioerr.ne.0 ) then
              write(*,'(3a)') 'Error reading values from '
     .             ,file_add,'(ioerr ilat ilon )'
              write(*,*) ioerr,ilat,ilon
              stop
            endif
          enddo
          close(iuadd)
        enddo
      endif


c Open the output file

      open(iuout,file=file_out,status='unknown',access='direct'
     .     ,recl=recl,form='unformatted',iostat=ioerr)
      if(ioerr.ne.0)then
          write(*,'(2a)') 'Error opening output file ',file_out
          stop "Program ended with errors"
      else
         write(*,'(/,2a)') 'Opened output file ',file_out
      endif

c  If no input file, set the header values from the add-file

      if( newfile ) then
        if( .not.addfile ) then
          write(*,'(a)') 'No old file or add-file, stop'
          stop
        else
          mapf = mapfa
          source = sourcea
          intervalhr = intervalhra
          intervald = float(intervalhr)/24.
          nglat = nglata
          nglon = nglona
c         ascii files have both 0 and 360 lon value, but grid only 0
          if( addtype.eq.'A') nglon = nglon -1
          ngval = ngvala
          start_time = start_timea
        endif
      else
c       check compatibility of input and add-file
        if( nglata.ne.nglat ) then
          write(*,'(a,2i5)') 'Grid-node nlat mismatch ',nglata,nglat
          stop
        endif
        if( nglona.ne.(nglon+1)) then
c         TUV file has both 0E and 360E so one more value than
c         GAMIT grid files; ignore the value at 360E
          write(*,'(a,2i5)') 'Grid-node nlon mismatch',nglona,nglon
          stop
        elseif( nglat.gt.maxlat.or.nglon.gt.maxlon) then
          write(*,'(a,4i4)') '# lat/lon > max lat/lon '
     .         , nglat,nglon,maxlat,maxlon
          stop
        endif
      endif

c  Define the grid size (# lat/lon pts per epoch)

      ngrid = nglat*nglon


c  Set the output epochs

      start_yr_out = int(start_time)
      start_doy_out =
     .       amod(start_time,1.0)*float(nydays(start_yr_out)) + 1.0
c     if there is an add-file, set the stop-day from its header or file name
      if( addfile ) then
        stop_yr_out = stop_yr_add
        stop_doy_out = stop_doy_add
cd        print *,'DEBUG stop_doy_out 2 from add-header ',stop_doy_out
c     otherwise, use the input values
      else
        stop_yr_out = stop_yr_in
        stop_doy_out = stop_doy_in
      endif
c     in either case, pad to the end of the day
c**  rwk/pt 070814: No, remove the padding, end the day at 0.75
cd      print *,'DEBUG stop_doy_out amod  '
cd     .       , stop_doy_out,amod(stop_doy_out,1.0)
      if( amod(stop_doy_out,1.0).ne.0. ) then
c*        tmpdoy = stop_doy_out
c*        stop_doy_out = stop_doy_out + (1. - amod(stop_doy_out,1.0))
c*        days_pad = stop_doy_out - tmpdoy
          days_pad = 0
cd        print *,'stop_doy_out days_pad intervald '
cd     .      ,stop_doy_out,days_pad,intervald
        nepoch_pad = nint(days_pad/intervald)
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
         i4y = stop_yr_out
         if( mod(iyrout2,100).eq.i4y .and.
     .       idoyout2.lt.nint(stop_doy_out) ) then
           stop_doy_out = idoyout2 + 1.0
           write(*,'(a,i3,f6.1,a)') '** Output end epoch set to '
     .        ,stop_yr_out,stop_doy_out,' by command-line value'
         endif
       endif


c  Write the output header records

c   Record 1:
      version_out = 10
      nhead_out = 3
      write(*,'(a)') 'Headers of output file: '
      write(*,'(a,1x,i4,1x,a4,1x,a2)')
     .    'Record 1 Version / Source ',version_out,mapf,source
      write(iuout,rec=1,iostat=ioerr) version_out,mapf,source
c   Record 2:
      write(iuout,rec=2,iostat=ioerr) intervalhr,ngval,nglon,nglat
      write(*,'(a,4i4)')
     .      'Record 2 interval #val #lon #lat '
     .        ,intervalhr,ngval,nglon,nglat
c   Record 3:
      stop_time = float(stop_yr_out)
     .      + (stop_doy_out-1.)/float(nydays(stop_yr_out))
      write(iuout,rec=3) start_time,stop_time
      write(*,'(a,2f12.6)')
     .      'Record 3 start stop  ',start_time,stop_time
       write(*,'(22x,a,2(i6,f7.2),a)') '( = ',start_yr_out,start_doy_out
     .       ,stop_yr_out,stop_doy_out,' )'


c  Calculate the records for the input and output file

      iy1 = start_yr_in
cd      print *,'DEBUG start_yr_in iy1 ',start_yr_in,iy1
      call fix_y2k(iy1)
      iy2 = stop_yr_in
      call fix_y2k(iy2)
cd      print *,'DEBUG start_yr_out iy2 ',start_yr_out,iy2
      days_in= diffdays(iy1,start_doy_in,iy2,stop_doy_in)
cd      print *,' '
cd      print *,'DEBUG In: yr doy yr doy days '
cd     .    ,iy1,start_doy_in,iy2,stop_doy_in,days_in
      iy1 = start_yr_out
      call fix_y2k(iy1)
      iy2 = stop_yr_out
      call fix_y2k(iy2)
      days_out= diffdays(iy1,start_doy_out,iy2,stop_doy_out)
cd      print *,'DEBUG Out: yr doy yr doy days '
cd     .    ,iy1,start_doy_out,iy2,stop_doy_out,days_out
      nepoch_in = nint(days_in/intervald) + 1
      nrec_in  = nhead_in +  ngrid*nepoch_in
cd      print *,'nepoch_in nrec_in ',nepoch_in,nrec_in
      if( addfile ) then
cd        print *,'start_yr_add,start_doy_add,stop_yr_add,stop_doy_add '
cd     .         , start_yr_add,start_doy_add,stop_yr_add,stop_doy_add
       if(.not.newfile.and.start_doy_add.gt.(stop_doy_in+intervald))then
          write(*,'(a,f8.2,a,f8.2,a,f5.2,a)') 'First added epoch ('
     .       ,start_doy_add,') > last in-file epoch (',stop_doy_in
     .       ,') + interval (',intervald,')'
          stop
        endif
        iy1 = start_yr_add
        call fix_y2k(iy1)
        iy2 = stop_yr_add
        call fix_y2k(iy2)
        days_add =
     .    diffdays( iy1,start_doy_add,iy2,stop_doy_add )
        nepoch_add = nint(days_add/intervald) + 1
        nrec_add = ngrid*nepoch_add
cd        print *,'days_add,nepoch_add,nrec_add '
cd     .      ,days_add,nepoch_add,nrec_add
      endif


c  Write the orography at the end of the file

c     The intent was to write the values at the first record of day
c     368 ( 367*91*144*4 + 3[head] +1 = 19 236 676 ); in fact in
c     Version 1.0 they are written 4 records prior to the beginning
c     of day 369 (368*91*144*4 + 3 +1 - 4 = 19 289 088).
      if( version_out.ne.10 ) then
        write(*,'(a)')
     .     'Output file not version 1.0, orography not defined'
        stop
      endif
      if( start_doy_out.ne.1.0 ) then
         write(*,'(a)') 'Start doy not 1.0--orography may be wrong'
         stop
      endif
      if( newfile ) then
c       open the orography file
        open(iuorog,file='orography.zhd',status='old',iostat=ioerr)
        if(ioerr.ne.0)then
          write(*,'(a)') 'Topography file for ZHD not found'
          stop "Cannot continue"
        else
          write(*,'(a)') 'Opened orography file orography.zhd'
        endif
c       read and print the header line of the orography file.  File has 10 columns of ht values
        read(iuorog,'(a)',iostat=ioerr) line
        if(  ioerr.eq.0 ) then
          write(*,'(2a)') 'Orography header: ',line
        else
          write(*,'(a,i3)') 'Error reading header of orography file '
     .                    ,ioerr
          stop
        endif
c       now read 144 longitude values (10 per line) for each latitude
c       and write them beginning at record 19 289 088 of the output file
c       (VMF orography file has 145 values, but the last is redundant)
c       grid_time accounts for header
        call grid_time( 'T',start_doy_out,intervald,ngrid,nhead_out
     .                , 369.,nrec_orostr )
c       subtract 4 for mistake (see above),
        nrec_orostr = nrec_orostr  - 4
cd        print *,'nhead_out nrec_orostr ',nhead_out,nrec_orostr
c       subtract 1 since increment is before write in loop below
        irec3 = nrec_orostr - 1
        do lrow = 1,nglat
          do i=1,14
            read(iuorog,'(10i8)',iostat=ioerr)
     .                          (ival3((i-1)*10+j),j=1,10)
            if( ioerr.ne.0 ) then
             write(*,'(2a,i3,a,i2,a,i3)')'Error reading orography value'
     .           ,' latrow =',lrow,' i=',i,' iostat=',ioerr
              stop
            endif
          enddo
          read(iuorog,'(4i8)',iostat=ioerr) (ival3(140+j),j=1,4)
          if( ioerr.ne.0 ) then
            write(*,'(2a,i3,a,i3)') 'Error reading 15th orography value'
     .              ,' latrow=',lrow,' iostat=',ioerr
            stop
          endif
          tmp = 0.0
          do j = 1,nglon
            irec3 = irec3 + 1
            iorog = ival3(j)
            write(iuout,rec=irec3) tmp,tmp,iorog
cd            if(j.eq.1) print*,ival3(1),row
cd             if(lrow.eq.48.and.j.eq.74) then
cd               print *,'Canberra ht (rec # )',ival3(j),irec3
cd             endif
cd             if( irec3.ge.19289080.and.irec3.le.19289235 ) then
cd               print *,'irec3 iorog ',irec3,iorog
cd             endif
          enddo
        enddo

      else
c       copy from old file
        if( version_in.ne.10 ) then
          write(*,'(a)')
     .     'Input file not version 1.0, orography may be wrong'
          stop
        endif
        call grid_time( 'T',start_doy_in,intervald,ngrid,nhead_in
     .                , 369.,nrec_orostr )
c       version 1 has start accidentally misplaced (day 368)
        nrec_orostr = nrec_orostr -4
        do irec1 = nrec_orostr, nrec_orostr+ngrid-1
          read(iuin,rec=irec1,iostat=ioerr) tmp,tmp,iorog
          write(iuout,rec=irec1) tmp,tmp,iorog
        enddo
      endif


c  Read and write the values records from the in_file

      if( .not.newfile ) then
cd      print *,'Re-reading infile nhead,nrec_in ',nhead_in,nrec_in
        do i = nhead_in+1,nrec_in
          irec1 = i
          read(iuin,rec=irec1,iostat=ioerr) ah,aw,zh
          if( ioerr.ne.0) call grderr(iuin,irec1,ioerr)
          if( swapped ) call swap_bytes(4,ah,1)
          if( swapped ) call swap_bytes(4,aw,1)
          if( swapped ) call swap_bytes(4,zh,1)
          irec3= irec1 + (nhead_out-nhead_in)
          write(iuout,rec=irec3) ah,aw,zh
        enddo
cd        print *,'DEBUG copied in-file, irec1 irec3 ',irec1,irec3
      else
        irec3 = nhead_out + 1
cd        print *,'NEWFILE nhead_out irec3 ',nhead_out,irec3
      endif



c  Add records from the supplementary file and/or
c  pad with last values for interpolation

      if( addfile ) then
        irec2 = nhead_add + 1
        call grid_time('T',start_doy_out,intervald,ngrid,nhead_out
     .    ,start_doy_add,irec3)
cd       print *,'start_doy_out start_doy_add irec2 irec3 '
cd     .        , start_doy_out,start_doy_add,irec2,irec3
        ifile = 0
        ilat = 1
        ilon = 0
cd        print *,'Start addition nrec_add ',nrec_add
        do i = 1, nrec_add
           if( addtype.eq.'B' ) then
             read(iuadd,rec=irec2,iostat=ioerr) ah,aw,zh
             if( ioerr.ne.0 ) call grderr(iuadd,irec2,ioerr)
             if( swappeda ) call swap_bytes(4,ah,1)
             if( swappeda ) call swap_bytes(4,aw,1)
             if( swappeda ) call swap_bytes(4,zh,1)
             irec2 = irec2 + 1
           elseif( addtype.eq.'A') then
cd             if( i.le.5 ) print *,'nglon ilon nglat ilat ifile ah',
cd     .            nglon,ilon,nglat,ilat,ifile,ahval(ifile,ilon,ilat)
             if( ilon.eq.nglon ) then
               ilon = 1
               if( ilat.eq.nglat ) then
                 ilat = 1
               else
                 ilat = ilat + 1
               endif
             else
               ilon = ilon + 1
             endif
             if( ilat.eq.1.and.ilon.eq.1 ) then
                ifile = ifile + 1
             endif
             ah = ahval(ifile,ilon,ilat)
             aw = awval(ifile+4,ilon,ilat)
             zh = zhval(ifile+8,ilon,ilat)
           endif
           write(iuout,rec=irec3) ah,aw,zh
c* debug
cd           if( irec3.le.5.or.(irec3.ge.13106.and.irec3.le.13108))
cd     .         print *,'irec3 ah aw zh ',irec3,ah,aw,zh
cd           if( (ilat.le.2.and.ilon.le.2).or.
cd     .         (ilat.eq.90.and.ilon.eq.144).or.
cd     .         (ilat.eq.91.and.ilon.eq.143 ) )
cd     .      print *,'i irec3 ifile ilat ilon ah aw zh '
cd     .              ,i,irec3,ifile,ilat,ilon,ah,aw,zh
           irec3 = irec3 + 1
        enddo
      endif
      write(*,'(a,2i9)')'Suppl records added, irec2 irec3 ',irec2,irec3

      do i = 1, nepoch_pad
c       reset input counters to first record of epoch
        irec1 = irec1 - ngrid
        irec2 = irec2 - ngrid
        write(*,'(a,3i9)') 'pad epoch irec1 irec2 ',i,irec1,irec2
        do k=1,ngrid
          if( addfile ) then
            if( addtype.eq.'B' ) then
              read(iuadd,rec=irec2,iostat=ioerr) ah,aw,zh
	           if( ioerr.ne.0 ) call grderr(iuadd,irec2,ioerr)
	           if( swappeda ) call swap_bytes(4,ah,1)
	           if( swappeda ) call swap_bytes(4,aw,1)
 	           if( swappeda ) call swap_bytes(4,zh,1)
    	       irec2 = irec2 + 1
            elseif( addtype.eq.'A' ) then
              ah = ahval(it,ilon,ilat)
              aw = awval(it,ilon,ilat)
              zh = zhval(it,ilon,ilat)
            endif
            write(iuout,rec=irec3) ah,aw,zh
            irec3 = irec3 + 1
          else
            read(iuin,rec=irec1,iostat=ioerr) ah,aw,zh
            if( ioerr.ne.0 ) call grderr(iuin,irec1,ioerr)
            if( swapped ) call swap_bytes(4,ah,1)
            if( swapped ) call swap_bytes(4,aw,1)
            if( swapped ) call swap_bytes(4,zh,1)
          endif
          write(iuout,rec=irec3) ah,aw,zh
          irec3 = irec3 + 1
cd          if(k.eq.1) print *,'wrote pad i irec3 ah aw zh '
cd     .       ,i,irec3,ah,aw,zh
        enddo
      enddo
      	

c  Close the files and write what we did

      if( .not.newfile ) then
         close(iuin,iostat=ioerr)
         if( ioerr.ne. 0)
     .   write(*,'(a,i4)') 'Error closing input file, iostat = ',ioerr
      endif
      if( addfile) then
        close(iuadd,iostat=ioerr)
         if( ioerr.ne. 0)
     .    write(*,'(a,i4)') 'Error closing input file, iostat = ',ioerr
      endif
      close(iuout,iostat=ioerr)
      if( ioerr.ne. 0)
     .  write(*,'(a,i4)') 'Error closing output file, iostat = ',ioerr
      if( addfile ) write(*,'(a,i10,a)')
     .  'Wrote ',nrec_add,' data records from add file to end '
      write(*,'(a,i10,a,/)')
     .  'Padded with ',nrec_pad,' duplicate records'

      stop
      end

c--------------------------------------------------------

      function diffdays(iyr1,doy1,iyr2,doy2)

      integer*4 iyr1,iyr2,idoy1,idoy2
     .        , isec,julday
      real*4 doy1,doy2,diffdays
      real*8 fjd1,fjd2,sec

      idoy1 = int(doy1)
      isec = int((doy1-float(idoy1))*86400.)
      call yds_to_jd( iyr1,idoy1,isec,fjd1)
cd      print *,' iyr1 doy1 sec fjd1 ',iyr1,idoy1,sec,fjd1
      idoy2 = int(doy2)
      isec = int((doy2-float(idoy2))*86400.)
      call yds_to_jd( iyr2,idoy2,isec,fjd2)
cd      print *,' iyr2 doy2 sec fjd2 ',iyr2,idoy2,sec,fjd2
      diffdays = fjd2 -fjd1
cd      print *,'  fjd1 fjd2 diffdays ',fjd1,fjd2,diffdays
      return
      end

c----------------------------------------------------------------------

      subroutine grderr(unit,record,err)

      integer*4 unit,record,err

      write(*,'(a,2i5,i10)') 'Read error (unit,err,record)'
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

