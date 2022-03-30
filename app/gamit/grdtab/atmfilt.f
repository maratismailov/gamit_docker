      Program ATMFILT

c   Program to apply a 20th-order Butterworth filter to atm grid files
c   to remove any sub-daily signals with any power. In other words, to
c   turn the "partial-tidal" atm binary grid files into proper "non-tidal"
c   grid files.
c
c   The filter requires 50 epochs at the start/end of each year to avoid
c   edge effects, so we read the last 13 days of the previous year's file,
c   and either read the 1st 13 days of the next year's file (if a complete
c   year) or 13 days beyond the end-time of the output file (if a partial
c   year).
c
c   The process is
c
c    a) for year X, open the binary grid files for years X-1, X, and
c       if a complete year, X+1
c    b) for each grid node, do the following:
c       - read in a time series of (NEU) deformation values, starting at
c         day 353 (354 if a leap year), through to 13 days beyond the 
c         last epoch of the file to be filtered, if necessary through day
c         13 of year X+1
c       - run the Butterworth filter over the timeseries for N, E, U
c       - output the filtered values for year X for the node
c
c   The header of the new grid file will indicate "B" for Butterworth filtering
c   for the type of file (rather than the current "P" for partial tidal information).
c   Otherwise, there should be no need for any changes to the header or the
c   the format of the file. The Butterworth filter was constructed by Chris Watson (UTAS)

c   (The code is based on that of atmtoasc.f, P. Tregoning 2004, R. King 2006)

c   P. Tregoning 22 September 2008
c   R. King 30 August 2010, mods to handle partial years

      implicit none

c  Variables from the input headers and output header
c  (grid spacing, interval, and models must be the same for all files, 
c   but not #header-lines or start/stop times)

      integer*2 nhead1,nhead2,nhead3,nhead,i2val,i2lon,i2lat
     .        , yr_start1,yr_end1,yr_start2,yr_end2,yr_start3,yr_end3
     .        , yr_start,yr_end,version

      real*4 interval,interval1,doy_start1,doy_end1,doy_start2,doy_end2
     .     , doy_start3,doy_end3,doy_start,doy_end,fversion,values(3)

      character*5 press_src,press_src1
      character*3 ref_frame,ref_frame1          
      character*2 spare

c  Other variables
      logical swapped1,swapped2,swapped3

      integer*4 maxepoch,nepoch,ngrid,ngrid1,ioerr,ts,tend    
     .         , ngval,ngval1,nglon,nglon1,nglat,nglat1
     .         ,irecbox(2,2),nrec,ilat1,ilon1,i,iepoch,start_epoch
     .         ,end_epoch,recl,irec
      parameter (maxepoch = 2000)

      real*4 dlat,dlon,n(4),e(4),u(4),dx,dx1,dy,dy1,t,zero
            
      real*8 displ(maxepoch,3),displfilt(maxepoch,3),end_time

      character atmfile1*26,atmfile2*26,atmfile3*26,arg1*30,arg2*3
     .          ,filtfile*50 

      logical partial_year,leap_year2                       
             
     
c  Initialization

      atmfile1 = ' '
      atmfile2 = ' ' 
      atmfile3 = ' ' 
      swapped1 = .false.
      swapped2 = .false.
      swapped3 = .false.
      leap_year2 = .false.
      partial_year = .false.

c  Define the record length. For now, all atm files have a record length of 12
      recl=12

c   the input atmospheric loading file contains displacements
c   (not corrections) of the Earth as a result of the pressure
c   loading. The first line of the file contains three values:
c
c   number of values per line, number of grid nodes in long, num in lat
c
c   The file then contains the displacements in N, E, U (in mm) for
c   each grid node in rows of longitude, starting at 90N, 00E then
c   running across to 357.5E before going to 87.5N, 00E. The first
c   set of grid entries correspond to decimal day 001.00, followed
c   by 001.25 etc through until the last day of the year. Thus, the
c   last row of the file will contain NEU values for grid point
c   90S, 357.5E at decimal day 365.75 (or 366.75 if it is a leap year).
c
c   Clear as mud? 

c  each line has 3 x real*4 values. Therefore, the record length is 12
              
c------------------------------------------------------------------------
c Decode the command lines and open the files
c------------------------------------------------------------------------
      call getarg(1,atmfile1)
      if(atmfile1(1:1).eq.' ')then
         print*,'Runstring: atmfilt atmfile1 atmfile2 atmfile3 '
     .         ,'filtfile'                
         print *,'atmfilt1 is the earliest (extra beginning) file'
         print *,'atmfilt2 is the middle (primary) file'
         print *,'atmfilt3 is the latest (extra-end) file'
         print *,"For a partial year, replace atmfile3 by ''"
         stop 'No command line arguments found. Program stopped'
      else
         open(unit=11,file=atmfile1,form='unformatted',access='direct'
     .        ,status='old',recl=recl,iostat=ioerr)
         if(ioerr.ne.0) then 
           write(*,'(3a,i4)') 'Error opening ',atmfile1,' iostat=',ioerr  
           stop
         else
            write(*,'(/,2a)') 'Opened previous year ',atmfile1
         endif             
        
      endif

      call getarg(2,atmfile2)
      if(atmfile2(1:1).eq.' ')then
         print*,'Runstring: atmfilt atmfile1 atmfile2 atmfile3 '
     .         ,'filtfile'
         stop 'Missing file name for atmfile2'
      else
         open(unit=12,file=atmfile2,form='unformatted',access='direct'
     .        ,status='old',recl=recl,iostat=ioerr)
         if(ioerr.ne.0) then 
           write(*,'(3a,i4)') 'Error opening ',atmfile2,' iostat=',ioerr  
           stop                                             
         else
           write(*,'(2a)') 'Opened current year ',atmfile2
         endif
      endif

      call getarg(3,atmfile3)
      if(atmfile3(1:1).eq.' ')then   
         partial_year = .true.
      else
         open(unit=13,file=atmfile3,form='unformatted',access='direct'
     .        ,status='old',recl=recl,iostat=ioerr)
         if(ioerr.ne.0) then 
           write(*,'(3a,i4)') 'Error opening ',atmfile3,' iostat=',ioerr  
           stop                           
         else
           write(*,'(2a)') 'Opened successive year ',atmfile3
         endif
      endif

      call getarg(4,filtfile)
      if(filtfile(1:1).eq.' ')then
         print*,'Runstring: atmfilt atmfile1 atmfile2 atmfile3 '
     .         ,'filtfile'
         stop 'Missing file name for output filtfile'
      else
         open(unit=14,file=filtfile,form='unformatted',access='direct'
     .        ,status='unknown',recl=recl,iostat=ioerr)
         if(ioerr.ne.0) then 
           write(*,'(3a,i4)') 'Error opening ',filtfile,' iostat=',ioerr  
           stop                                              
         else      
           write(*,'(2a)') 'Opened file for filtered output ',filtfile
         endif

      endif
 
c-----------------------------------------------------------------
c Read the headers of the three files
c------------------------------------------------------------------

      call read_atml_headers( 12,nhead2,swapped2,press_src,ref_frame
     .                      , interval,ngval,nglon,nglat
     .                      , yr_start2,doy_start2,yr_end2,doy_end2 )
      ngrid = nglon*nglat
 
      call read_atml_headers( 11,nhead1,swapped1,press_src1,ref_frame1
     .                      , interval1,ngval1,nglon1,nglat1
     .                      , yr_start1,doy_start1,yr_end1,doy_end1 )
      ngrid1 = nglon1*nglat1
      if(press_src1.ne.press_src .or. 
     .   ref_frame1.ne.ref_frame .or.
     .   ngval1.ne.ngval .or.
     .   ngrid1.ne.ngrid .or.
     .   interval1.ne.interval ) then
         print *,'File for preceeding year incompatible'
         print *,'            File 1    File2 '
         print *,'press_src ',press_src1,press_src
         print *,'ref_frame ',ref_frame1,ref_frame
         print *,'interval  ',interval1,interval
         print *,'ngval     ',ngval1,ngval
         print *,'ngrid     ',ngrid1,ngrid
         stop
      endif
                              
      if( .not.partial_year ) then       
        call read_atml_headers( 13,nhead3,swapped3,press_src1,ref_frame1
     .                      , interval1,ngval1,nglon1,nglat1
     .                      , yr_start3,doy_start3,yr_end3,doy_end3 )
        if(press_src1.ne.press_src .or. 
     .    ref_frame1.ne.ref_frame .or.     
     .    ngval1.ne.ngval .or.
     .    ngrid1.ne.ngrid .or.
     .    interval1.ne.interval ) then
          print *,'File for succeeding year incompatible'
          print *,'            File 3    File2 '
          print *,'press_src ',press_src1,press_src
          print *,'ref_frame ',ref_frame1,ref_frame
          print *,'interval  ',interval1,interval           
          print *,'ngval     ',ngval1,ngval
          print *,'ngrid     ',ngrid,ngrid1
          stop
        endif     
      endif
      
c PT040106: NOTE that the ATM grids go from 90N to 90S 
c PT040106: divide by nglat-1 because the files include 90S and 90N
       dlat = 180./float(nglat-1)
       dlon = 360./float(nglon)

c -------------------------------------------------------------------
c         Setup   output   grid   file   header
c -------------------------------------------------------------------
             
c  Write the file in the current output version 2.1.
      version = 21
      nhead = 3
c  Change the pressure source to reflect that it is c  now non-tidal. 
c  Currently all GAMIT files have press_src = "NCEP ".  We will set 
c  the 1st character to "F" to show that it is filtered, that is, "FNCEP"
      press_src(1:1) = "F"
c  Reference frame must be the same for all the input and output files (checked above)
      write(14,rec=1,iostat=ioerr) version,nhead,press_src,ref_frame 

c Record 2                                                               
c   Grid spacing must be the same for all input and output files (checked above)  
      spare = "  "                                
      i2val = ngval
      i2lon = nglon
      i2lat = nglat 
      write(14,rec=2,iostat=ioerr) interval,i2val,i2lon,i2lat,spare 

c Record 3. The new file goes either to 13 days before the end of the current
c           file or through the first epoch of the next year (for interpolation)
      yr_start = yr_start2                  
      doy_start = doy_start2
      if( partial_year ) then 
         yr_end = yr_start2
         doy_end = doy_end2 - 13.0
      else      
        yr_end = yr_start2 + 1.0    
        doy_end = 1.75 
      endif
      write(14,rec=3,iostat=ioerr) yr_start,doy_start,yr_end,doy_end  
      write(*,'(/,a,/,a,i2,1x,i2,1x,a5,1x,a3,1x,/,a,f4.2,1x,3i3,/
     ,          ,a,2(1x,i2,1x,f6.2))') 
     .   'Wrote header of output file '
     .  , ' Record 1: ',version,nhead,press_src,ref_frame 
     .  , ' Record 2: ',interval,i2val,i2lon,i2lat  
     .  , ' Record 3: ',yr_start,doy_start,yr_end,doy_end

c -------------------------------------------------------------------
c         Start   of   loop   over   each   grid   node
c -------------------------------------------------------------------

      do ilat1 = 1,nglat
        Print*,"Latitude ",90 - (ilat1-1) * dlat
        do ilon1 = 1,nglon

c         reset the epoch counter for each node
          nepoch = 0
cd        print*,'DEBUG: LL grid lat long:',ilat1,ilon1
  
c         Calculate the record-number offsets from the 1st point of the epoch
c            (i.e., upper left = 0)
c RWK100902: Fix count so 90N is record 0 
c          irecbox(1,1) = (ilat1-1)*nglon+ilon1
          irecbox(1,1) = (ilat1-1)*nglon+ilon1 - 1

c PT081008: DEBUG  hardwire it for coords -25, 132.5E to test the output
c          irecbox(1,1) = 6677

c---------Read the first file from day 353 (or 354) to through the end (13 days, 52 epochs)

c RWK not needed?   doy_start = 1.0
          if(mod(float(yr_start-1),4.).eq.0)then  ! The first file is a leap year file
            end_time = 366.75  
            t = 354.0                                        
c            print*,'leap year'
          else
            end_time = 365.75                                          
            t = 353.0
c            print*,'no leap year'
          endif     
c         set the time counter for the grid calls    
c         print*,'Reading period ',t,' to',end_time,' from file ',atmfile1
          do while(t.lt.(end_time+interval))  
            nepoch = nepoch + 1  
            call grid_time('T',doy_start1,interval,ngrid,nhead1,t,nrec)
cd            print*,'DEBUG file 1: nepoch,t and nrec:',nepoch,t,nrec
c            read values for the single node at this epoch
             read(11,rec=irecbox(1,1)+nrec,iostat=ioerr) n(1),e(1),u(1)
cd            print*,nepoch,-13.0 + float(nepoch)/4.0,t,u(1)
             if( ioerr.ne.0 ) then  
               print *,'ioerr rec ',ioerr,irecbox(1,1)+nrec
               write(*,'(a,i8,a,f6.2)') 
     .          'Attempt to read beyond EOF of file 1, record = '
     .            ,nrec,'   t=',t
               stop
             else
               if( swapped1 ) then
                 call swap_bytes(4,n,4)
                 call swap_bytes(4,e,4)
                 call swap_bytes(4,u,4)
               endif  
             endif
             displ(nepoch,1) = n(1)
             displ(nepoch,2) = e(1)
             displ(nepoch,3) = u(1)
cDEBUG      write out the time and displacements
cd            write(*,'(i6,f7.2,3f11.4)') nepoch,t,(displ(nepoch,i),i=1,3)
              t = t + interval
          enddo  ! end of read of first file days 353-365 (or 366)
          

c---------Read the second file from day 001 through the end 

          t = 1.0         
          if( partial_year ) then   
            end_time = doy_end2
          else 
c           check for leap year
            if(mod(float(yr_start),4.).eq.0) then  
              end_time = 366.5  ! it therefore reads to 366.75 
              leap_year2 = .true.                                    
            else
              end_time = 365.5   ! it therefore reads to 365.75  
            endif
c             print*,'Reading period ',t,' to',end_time,' from file '
c     .         ,atmfile2
          endif
c          do while(t.le.(end_time+interval))  
          do while(t.le.end_time)
            nepoch = nepoch + 1  
            call grid_time('T',doy_start2,interval,ngrid,nhead2,t,nrec)
cd            print*,'DEBUG file 2: nepoch,t and nrec:',nepoch,t,nrec
c           read values for the single node at this epoch
            read(12,rec=irecbox(1,1)+nrec,iostat=ioerr) n(1),e(1),u(1)
            if( ioerr.ne.0 ) then  
              print *,'ioerr rec ',ioerr,irecbox(1,1)+nrec
              write(*,'(a,i8,a,f6.2)') 
     .         'Attempt to read beyond EOF, record = '
     .           ,nrec,'   t=',t
              stop
            else
              if( swapped2 ) then
                call swap_bytes(4,n,4)
                call swap_bytes(4,e,4)
                call swap_bytes(4,u,4)
              endif  
            endif
            displ(nepoch,1) = n(1)
            displ(nepoch,2) = e(1)
            displ(nepoch,3) = u(1)      
cd DEBUG 
cd            write(*,'(a,4i8)') 'ilat1 ilon1 irecbox11 nrec '
cd     .        ,ilat1,ilon1,irecbox(1,1),nrec 
cd           write out the time and displacements
cd            write(*,'(i6,f7.2,3f11.4)') nepoch,t,(displ(nepoch,i),i=1,3)
cd end DEBUG
            t = t + interval
          enddo  ! end of read of second file 

c-------- Read the third file through to the end of day 13
c                        
          if( .not.partial_year ) then
            end_time = 13.75                                          
            t = 1.0
c            print*,'Reading period ',t,' to',end_time,' from file '
c     .         ,atmfile3
            do while(t.le.(end_time+interval))  
              nepoch = nepoch + 1  
              call grid_time('T',doy_start,interval,ngrid,nhead,t,nrec)      
cd              print*,'DEBUG file 3: nepoch,t and nrec:',nepoch,t,nrec
c             read values for the single node at this epoch
              read(13,rec=irecbox(1,1)+nrec,iostat=ioerr) n(1),e(1),u(1)
              if( ioerr.ne.0 ) then  
                print *,'ioerr rec ',ioerr,irecbox(1,1)+nrec
                write(*,'(a,i8,a,f6.2)') 
     .            'Attempt to read beyond EOF, record = ',nrec,'   t=',t
                stop
              else
                if( swapped3 ) then
                  call swap_bytes(4,n,4)
                  call swap_bytes(4,e,4)
                  call swap_bytes(4,u,4)
                endif  
              endif
              displ(nepoch,1) = n(1)
              displ(nepoch,2) = e(1)
              displ(nepoch,3) = u(1)
c             write out the time and displacements
cd DEBUG
cd            write(*,'(i6,f7.2,3f11.4)') nepoch,t,(displ(nepoch,i),i=1,3)
              t = t + interval
            enddo  ! end of read of third file days 001-050
         endif  ! endif for need to read third file
                          
c--------------------------------------------------------------------------
c DEBUG: print out the time series for each file. The subtraction of 13.25
c        means that the first epoch of the second file comes out with a
c        time stamp of 1.00, ie equivalent to the epoch that comes out of
c        atmtoasc. This represents YYYY.0 or the first epoch of the first day of year
c
c PT08008: tested on files atmdisp_cm.200[5-7] and it was correct.
c          print*,'Raw atml from binary grid files'
c          do iepoch = 1,nepoch
c            write(*,'(4f10.4)'),-13.25 + float(iepoch)/4.0
c     .                   , (displ(iepoch,i),i=1,3)
c          enddo
c-------------------------------------------------------------------------------------                
                                                          
c-------------------------------------------------------------------------------------
c   Filter the complete time series
c-------------------------------------------------------------------------------------

           call atm_butter(maxepoch,nepoch,displ,interval,displfilt)
           

c--------------------------------------------------------------------------------------
c   Write the filtered values into a new file for the year of interest
c-------------------------------------------------------------------------------------

c  facile! Now we write out the records for the year of interest into the 
c  output  binary grid file. The first record for the year starts at 53 
c  (i.e 13 days * 4/day + 1). The last epoch will be either the end of the 
c  year + 4 (for interpolation) or (for partial years) 13 days before the
c  end of the input file.  
                                                  
          start_epoch = 53
          if( partial_year ) then
            end_epoch = start_epoch + doy_end*4 
          else
            if (leap_year2) then
              end_epoch = start_epoch + 366*4 + 4
            else
              end_epoch = start_epoch + 365*4 + 4  
            endif 
          endif
cd          print *,'For output file, end_epoch ',end_epoch                                                  

c---DEBUG----------------------------------------------------------
c  PT081008: just output to the screen the raw and the filtered time series
c            so that I can plot them and see whether the filter is working
c      print*,'writefiltnode: filtered time series'
c      do i=start_epoch,end_epoch
c        write(*,'(7f10.4)')float(i)/4.0 - 50.25,(displfilt(i,j),j=1,3)
c      enddo
c----------------------------------------------------------------------

c         now we simply loop through the whole year and write out the values into
c         the appropriate records for this node. The first record is going to
c         be irecbox + 3 header lines     
c PT/RWK 100902: With irecbox(1,1) now zero for upper left, need to add 1 more in write
c         irec = irecbox(1,1) + 3                                                      
          irec = irecbox(1,1) + 4
          do i = start_epoch,end_epoch
            values(1) = displfilt(i,1)
            values(2) = displfilt(i,2)
            values(3) = displfilt(i,3)
cDEBUG use the input values to check epoch alignment
c           values(1) = displ(i,1)
c           values(2) = displ(i,2)
c           values(3) = displ(i,3)
c           print*,recl,values
c           recl = 4
            write(14,rec=irec) values  
cd            write(*,'(a,2i8,f11.4)') 'writing epoch rec Up '
cd     .           ,i,irec,values(3)
            irec = irec + ngrid
          enddo
 
cd          print *,'DEBUG stop'
cd          stop 

        enddo   ! end of longitude loop
      enddo   ! end of latitude loop
c -------------------------------------------------------------------------------------

c eh voila, c'est finit. Sante!
 
      close(11)
      close(12)
      close(13)
      close(14)

      print*,"C''est tout finit!"

      end
