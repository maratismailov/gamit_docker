      program vmfasc2bin

c  program to convert the input data for the VMF1 mapping
c  function from two separate, ascii files into one
c  single, direct access binary file.
c
c  The two coefficients per node are Ah and Aw. 
c
c  The input files are coming from Johannes Boehm (Tuviens) and
c  contain values on a global grid whose dimensions are given in
c  the first line of the file. As at 26 May, 2004, this grid is
c  2 x 2.5 deg in lat, long, respectively. 
c
c  The file starts at 90N, running from 0 E to 360 E (including 360 E)
c  then steps down to 88 N, 0-360E etc
c
c  To be compatible with the atmospheric loading grid file
c  used in gamit (atmdisp.YYYY), I will write out the IMF grid
c  from 0 E to (360-dlon) E. That is, one less grid point per
c  row than the input data. Also, both parameters will be
c  written into the same file.
c
c  I will also take the liberty of writing the sea-level
c  atmopsheric pressure value as the third column of this
c  file. We don't yet have an input source of data for this
c  value; however, once we do it will be easy to regenerate
c  these files. The code will then exist within gamit to read
c  it directly at the same time as reading in the VMF1 coefficients.
c
c  P. Tregoning
c  30 May 2006.
c
c  PT060829: remove the atm pressure from this file. Output only
c            the Ah and Aw values. The pressure and/or ZHD will
c            be output into different file(s).
c 
c            Also, add header lines so that it is compatible
c            with other binary grid files.
c
c Header line 1: 10VMF1JB  (I*2 C*4 C*2 then 4 blanks)
c Header line 2: 062144091 (06 hour epochs, 2 values/record, nglon=144, nglat=091)
c Header line 3: start_time,stop_time (R*4 R*4)
c
c PT060830: Bob King wants the ZHD written to this grid file as well as the
c           mapping function coefficients.
c
c PT070329: modified to output only the start/stop time if argument "head"
c           is requested.

      implicit none

c  maximum dimensions for a 0.5 x 0.5 deg grid
      integer maxrow,maxcol
      parameter (maxrow=721,maxcol=361)

      character outfile*13,infile*12,cyear*4,line*100,head*4
      integer*4 i,j,lu1,lu2,lu3,lu4,lu5
     .         ,start_yr,stop_yr,ioerr,nrec,doy,doy_last
     .         ,row,hr,hr_last

      real*4 minlat,minlon,maxlat,maxlon,dlat,dlon
     .      ,header(6),tmp(6)
      
      integer*2 version,interval,ngval
      character mapf*4,source*2,spare*4
      integer*4 Ahval(maxrow),Awval(maxrow)
     .         ,ival1(maxrow),ival2(maxrow),ival3(maxrow)
     .         ,nglon,nglat
      real*4 start_time,stop_time,Zhval(maxrow),start_doy,stop_doy
      logical foundfile 

c  function
      integer*4 nydays
      

c unit numbers
c lu1: Ah files
c lu2: Aw files
c lu3: Zh files
c lu4: orography.zhd
c lu5: output grid file

      lu1 = 1
      lu2 = 2
      lu3 = 3
      lu4 = 4
      lu5 = 10

c Get the year from the command line argument
      call getarg(1,cyear)
      if(cyear(1:1).eq.' ')then
        print*,"Runstring: vmfasc2bin YYYY"
        stop
      elseif(cyear.eq."head")then
        print*,'Will read only the start/stop time'
c  open the output file. This will be a binary, direct access file
        head = cyear
        call getarg(2,cyear)
        read(cyear,'(i4.4)') start_yr
        write(outfile,'(a8,i4.4)')'vmf1grd.',start_yr
        open(lu5,file=outfile,status='old',access='direct'
     .     ,recl=12,iostat=ioerr)
        read(lu5,rec=3) start_time,stop_time   
        print*,"start:",start_time,' Stop:',stop_time  
        start_doy=amod(start_time,1.0)*float(nydays(start_yr))+1.0
        call round6h(start_doy)  
        stop_yr = nint(stop_time)
        stop_doy=amod(stop_time,1.0)*float(nydays(stop_yr))+1.0   
        call round6h(stop_doy)             
        write(*,'(22x,a,2(i6,f7.2),a)') '( = ',start_yr,start_doy
     .     , stop_yr,stop_doy,' )'
        stop
      else
        read(cyear,'(i4.4)') start_yr
        write(*,'(a,i4)') 'Year = ',start_yr
      endif                

c  RWK 070816: For now require the file to stop within the same year
c              (though it may go through day 367)
                     
      start_doy = 1.0
      stop_yr = start_yr     
      stop_doy = float(nydays(stop_yr)) + 0.75
                        

c  PT060530: input files have the format:
c            aw05193.h18  (ie wet coefficients on 12 July 2005 at 18 UT)
                        
c  open the output file. This will be a binary, direct access file
c  with one grid node (three values per node) written per line. Therefore,
c  the record length will be 3 x R*4 = 12
      write(outfile,'(a8,i4)')'vmf1grd.', start_yr

c  PT060530: the record length required will be:
c            I*4  I*4  = 8  
c  PT060830: adding the ZHD makes it I*4 I*4 R*4  = 12

      open(lu5,file=outfile,status='unknown',access='direct',recl=12
     .     ,iostat=ioerr)
      if(ioerr.ne.0)then
        print*,'Error opening output file ',outfile
        stop 'stopped with errors'
      else
        print*,'Opened output file :',outfile
      endif



c PT060530: for the VMF1, we are going to read each 6-hourly
c           file one at a time (rather than read a concatenated
c           file as per the imfasc2bin program). Therefore,
c           from here down is a big loop from day of year 001
c           until the file doesn't exist to be opened

      doy = 001    


c  open the orography file
      open(lu4,file='orography.zhd',status='old',iostat=ioerr)
      if(ioerr.ne.0)then
        print*,'Topography file for ZHD not found'
        stop "Cannot continue"
      else
        print*,'Opened orography file orography.zhd'
      endif

      ioerr = 0
      foundfile = .true.
      do while ( foundfile )
        do hr = 1,4

c  open the input files
          write(infile,'(a2,a2,i3.3,a2,i2.2)')'ah',cyear(3:4),doy
     .           ,'.h',(hr-1)*6
          open(lu1,file=infile,status='old',iostat=ioerr)
          if( ioerr.ne.0.and.foundfile ) then 
            write(*,'(a)') 'Error opening input file ',infile
            write(*,'(a,2i3,a)') 'Assume last available epoch is '
     .          ,doy_last,hr_last,'h'   
            stop_doy = float(doy_last) + float(hr_last)/24.
            stop_time = float(start_yr) 
     .         + (stop_doy-1.)/float(nydays(start_yr))
c* RWK 070816: Don't write the third header record until we've read the input files
c*            write(lu5,rec=3) start_time,stop_time   
            foundfile = .false.
          endif
          if( foundfile ) then
            write(*,'(2a)') 'Opened file ',infile  
c           save the doy and hr to determine the stop epoch
            doy_last = doy
            hr_last = (hr-1)*6
            write(infile,'(a2,a2,i3.3,a2,i2.2)')'aw',cyear(3:4),doy
     .            ,'.h',(hr-1)*6
            open(lu2,file=infile,status='old',iostat=ioerr)
            if(ioerr.ne.0) then   
            write(*,'(3a)') 'Tried to open ',infile
     .           ,'--something wrong, stop'  
              stop
            endif
            write(infile,'(a2,a2,i3.3,a2,i2.2)')'zh',cyear(3:4),doy
     .            ,'.h',(hr-1)*6
            open(lu3,file=infile,status='old',iostat=ioerr)
            if(ioerr.ne.0) then 
            write(*,'(3a)') 'Tried to open  ',infile
     .           ,'--something wrong, stop'  
              stop
            endif

c  read the header lines of the Ah-input file                
            read(lu1,'(1x,6f8.2)')(header(i),i=1,6)                     

c check that the Aw file is the same:
            read(lu2,'(1x,6f8.2)')(tmp(i),i=1,6)
            do i=1,6
              if(tmp(i).ne.header(i)) then
                print*,'header of Aw file does not match Ah file'
                print*,'Ah-file: ',header
                print*,'Aw-file: ',tmp
                stop 'cannot proceed'
              endif
            enddo

c check that the Zh file is the same:
            read(lu3,'(1x,6f8.2)')(tmp(i),i=1,6)
            do i=1,6
              if(tmp(i).ne.header(i)) then
                print*,'header of Zh file does not match Ah file'
                print*,'Ah-file: ',header
                print*,'Zh-file: ',tmp
                stop 'cannot proceed'
              endif
            enddo

c  assign the header values to meaningful names
c PT060530: Johannes has changed the order of the header
c           values to be:   maxlat  minlat  minlong  maxlong  dlat  dlon
            minlat = header(2)
            maxlat = header(1)
            minlon = header(3)
            maxlon = header(4)
            dlat   = header(5)
            dlon   = header(6)
c  calculate the number of grid nodes in latitude
            nglat = (maxlat - minlat)/dlat + 1
            nglon = (maxlon - minlon)/dlon
            if(nglon.gt.maxrow)then
              stop 'grid too big. Redimension maxrow!'
            endif

c  set the number of parameters per node for the output file
            ngval = 3

c  write out the header for the binary grid file
            if(doy.eq.001.and.hr.eq.1)then

c  PT060829: write out a header line indicating information about the VMF1 version etc
c  col 1-2: version
c  col 3-6: mapping function
c  col 7-8: source (JB=Johannes Boehm)
              version = 10
              mapf = "VMF1"
              source = "JB"
              spare = "    "
              write(lu5,rec=1)version,mapf,source,spare

c  PT060829: the time interval between epochs is 6 hours
              interval = 6
              print*,interval,ngval,nglon,nglat
              write(lu5,rec=2)interval,ngval,nglon,nglat

c  PT060829: third header line contains start and stop times (decimal year)
c            This will be written out once the end time is known from
c            reading the input files.
              nrec = 3
            endif

c  ok. Now read the global grid one time-step at a time and
c  write out the three parameters into the binary file. The
c  input files are written one band of latitude at a time,
c  with 12 values/line. That is, (12 rows of 12 values + one
c  single value on the 13th row) per band of latitude.

            ioerr = 0
            do row = 1,nglat

              do i=1,12
                read(lu1,'(12i7)',iostat=ioerr)
     .                       (ival1((i-1)*12+j),j=1,12)
                if(ioerr.eq.0) read(lu2,'(12i7)',iostat=ioerr)
     .                           (ival2((i-1)*12+j),j=1,12)
                if(ioerr.eq.0) read(lu3,'(12f7.5)',iostat=ioerr)
     .                          (zhval((i-1)*12+j),j=1,12)
c       print*,'ival1(j),',i,row,(ival1(j),j=1,6)
              enddo

c  read the 145th value that sits on its own line
              if(ioerr.eq.0)
     .           read(lu1,'(i7)',iostat=ioerr)ival1(nglon+1)
              if(ioerr.eq.0)
     .           read(lu2,'(i7)',iostat=ioerr)ival2(nglon+1)
              if(ioerr.eq.0)
     .           read(lu3,'(f7.5)',iostat=ioerr)zhval(nglon+1)

              do i = 1,nglon
                Ahval(i) = ival1(i) 
                Awval(i) = ival2(i) 
              enddo

c  debug 
c  work out what latitude band we are up to
c        print*,'Latitude: ', 90.0 - (float(nrec-3)/(nglon-1))*dlat
c        print*,'nrec = ',nrec,ioerr
c        print*,(Ahval(i),i=1,nglon)
c        print*,' '
c        print*,(Awval(i),i=1,nglon)

              if(ioerr.eq.0)then
c  input data is 0 to 360 deg (ie last column repeats first). Write it out
c  as only 0 to (360 - dlon). This makes it a) smaller ; b) compatible with
c  the atm loading grid files.
                do j = 1,nglon
                  nrec = nrec + 1
                  write(lu5,rec=nrec)Ahval(j),Awval(j),Zhval(j)
                enddo
c  print the day of year to show what is going on
c  there are 91 x 144 x 4 records /day (=52416) in the binary output files
                if(mod(float(nrec)/52416.0,1.0).lt.0.001)
     .            print*,'Day of year ',1+nrec/52416,' rec #',nrec
              endif

c  end of the loop over latitudes for each particular epoch file
            enddo
            close(lu1)
            close(lu2)
            close(lu3)

c  end of if statement for found files
          endif       

c  end of loop for epoch (ie the 4 files per day)
        enddo

c  increment the day of year
        doy = doy + 1

c  end of day-of-year loop
      enddo

c'est tout!
      write(*,'(a,i10,a)') 'Wrote ',nrec,' records to output file ' 
            
c  Now write out the start, stop epochs
             
      start_time = float(start_yr) 
     .      + (start_doy-1.)/float(nydays(start_yr))
      stop_time = float(stop_yr) 
     .      + (stop_doy-1.)/float(nydays(stop_yr))
 
      write(*,'(a,11x,2f10.4)') 
     .  'Writing start and stop time to header:',start_time,stop_time
      write(lu5,rec=3) start_time,stop_time    
      write(*,'(22x,a,2(i6,f7.2),a)') '( = '
     .      ,start_yr,start_doy,stop_yr,stop_doy,' )'

c PT060607: now, after the 367th day, read and write a grid of
c           topography so that the ZHD can be corrected in GRDTAB
c           to the height of the requested station.  
c RWK070816: This code actually writes the orthography beginning 
c           at the 369th day - 4 (unintentional miscalculation),
c           but don't change it for Version 1.0.

c  the orography is read in for each epoch. We then
c  write it out. Leave one day blank so that I can add day 001 of
c  the next year ....
      nrec = 368 * nglon*nglat*24/interval - 1
      write(*,'(a,i10)') 'Writing orography starting at record ',nrec+1

c read a header line of the orography file. File has 10 columns of height values
      read(lu4,'(a)')line
c      print*,'line = ',line
      ioerr = 0
      do row = 1,nglat
        do i=1,14
          if(ioerr.eq.0) read(lu4,'(10i8)',iostat=ioerr)
     .                         (ival3((i-1)*10+j),j=1,10)
c          print*,(ival3((i-1)*10+j),j=1,10)
        enddo
        if(ioerr.eq.0)read(lu4,'(4i8)',iostat=ioerr)
     .         (ival3(140+j),j=1,4)
        tmp(1) = 0.0
        do j = 1,nglon
          nrec = nrec + 1
          write(lu5,rec=nrec)tmp(1),tmp(1),ival3(j)
c          if(j.eq.1)print*,ival3(1),row 
        if(row.eq.48.and.j.eq.74)then
          write(*,'(a,i4)') 
     .       '**Test: topographic heights for Canberra:',ival3(j)
        endif
        enddo
      enddo


c'est tout!


      close(lu4)
      close(lu5)
      end



