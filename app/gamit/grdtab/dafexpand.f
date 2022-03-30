c     Program to convert a set of gridded ascii ocean tide files produced by
c     Hans-Georg Scherneck (Onsala) from 1-degree spacing to half-degree spacing.
c     Used to allow easier conversion of multiple ascii files to a single GAMIT 
c     direct-access file in program /gamit/utils/daf2da.f.

c     Written by R. King  14 July 2000
                                      
c     The Scherneck global models come as 11 files, one for each tidal
c     constituent, with the values (complex pairs) for the diurnal and 
c     semi-dirunal tides tabulated at 0.5-degree intervals.  There are 180 
c     latitude rings, ordered from 90 S to 89.5 N; within each ring are 720 
c     longitude pairs, ordered from 0.5 to 360.0 E.  The long-period tides
c     are tabulated at 1-degree intervals from 89.5 S to 89.5 N and 0.5 to
c     359.5 E.  This program interpolates between 1-degree values to obtain
c     0.5-degree values.  At the south pole, the values for 89.5 S are simply
c     copied.  

      implicit none
              
c   Runstring:  dafexpand [in-file]  [out-file]

c      where in-file is a 1-degree grid and out-file is a half-degree grid.


      logical fcheck
                   
      character*1 ftyp
      character*2 atide
      character*4 vartyp           
      character*9 module
      character*20 infile,outfile

      integer*4 luin,luout,ilat,iarg,iclarg,nlon,nlat,nblen
     .        , irec, nval,ioerr,i

      data luin/1/,luout/2/

      real*4 inval(720),outval1(1440),outval2(1440),outval3(1440)
        
***** Remove old versions of the status, warning, and error files
                         
      module = 'DAFEXPAND'
      call report_stat('CLEAR',module,' ',' ', ' ',0)
      call report_stat('CLEAR','LIB',' ',' ', ' ',0)


c  Get the name of the input file from the command-line argument 
                      
      iarg = iclarg(1,infile)
      if  (iarg .le. 0) then    
        call report_stat('FATAL',module,'utils/dafexpand',' '
     .                   ,'Missing input file in command-line',0)
      else
c       stop if the input file given on the command line does not exist
        if( .not.fcheck(infile) ) 
     .   call report_stat('FATAL',module,'utils/dafexpand',infile
     .                   ,'Input file not found:',0)
      endif                      
      open(luin,file=infile,iostat=ioerr,status='old')   
      if( ioerr.eq.0 )  then
        call report_stat('STATUS',module,'utils/dafexpand',infile
     .                     ,'Opened input grid file:',ioerr)
      else
        call report_stat('FATAL',module,'utils/dafexpand',infile
     .                     ,'Error opening input file:',ioerr)
      endif
   
c  Get the name of the output file from the command-line argument 
                      
      iarg = iclarg(2,outfile)
      if  (iarg .le. 0) then 
        outfile(1:nblen(infile)+4) = infile(1:nblen(infile))//'.exp'
c        call report_stat('WARNING',module,'utils/dafexpand',outfile
c     .                   ,'Assigning outfile name: ',0)
      endif
      open(luout,file=outfile,iostat=ioerr,status='unknown')   
      if( ioerr.eq.0 )  then
        call report_stat('STATUS',module,'utils/dafexpand',outfile
     .                  ,'Opened output grid file:',ioerr)
      else
        call report_stat('FATAL',module,'utils/dafexpand',outfile
     .                  ,'Error opening output file:',ioerr)
      endif

c  Read and write the header line 
                             
c     get grid spacing and tide type from first line of file
      read(luin,'(2i5,1x,a1,1x,a4,1x,a2,10x,i1)',iostat=ioerr)  
     .         nlon,nlat,ftyp,vartyp,atide,nval
      if( ioerr.ne.0 ) 
     .  call report_stat('FATAL',module,'utils/dafexpand',infile
     .               ,'Error reading header of input file',ioerr)
      if( nlon.ne.360 .or. nlat.ne.180 ) 
     .      call report_stat('FATAL',module,'utils/dafexpand',' '
     .                      ,'Input file not 1 deg x 1 deg ',0)
      nlon = 720
      nlat = 360  

      write(luout,'(2i5,1x,a1,1x,a4,1x,a2,10x,i1)',iostat=ioerr)
     .         nlon,nlat,ftyp,vartyp,atide,nval

                  
c    each input latitude band contains 720 1-deg long values (360 complex pairs),
c      spread over 72 lines (10 values per line)
c    each output laitutde band contains 1440 half-deg long values (720 complex pairs),
c      srpead over 144 lines (10 values per line)
                
c     read the first latitude band (89.5 S) into storage

        read(luin,'((10e12.4))') (inval(i),i=1,720)
c       expand to 1440 values by interpolation
        call expand_long(inval,outval1)   
c       copy this band into the 90 S  and 89.5 S bands in the output
        write(luout,'((1p10e12.4))') (outval1(i),i=1,1440)
        write(luout,'((1p10e12.4))') (outval1(i),i=1,1440)
        irec = 2

c    loop over all the other latitude bands, copying if an even index
c    corresponding to a non-even latitude (i.e. 2=89.5, 3=89.0, etc)
c    and interpolating if an odd index.   
                                                         
      do ilat = 2,180

c       read a new 1-degree band (e.g. 88.5) and expand it by interpolation to fill all .5-deg longitudes
        read(luin,'((10e12.4))') (inval(i),i=1,720)
        call expand_long(inval,outval3)   
c       interpolate the intervening latitude band (e.g. 89.0)
        call expand_lat(outval1,outval3,outval2)
c       write the two latitude bands
        write(luout,'((1p10e12.4))') (outval2(i),i=1,1440)
        write(luout,'((1p10e12.4))') (outval3(i),i=1,1440)
        irec = irec + 2
c       shift storage to read the next 1-degree band
        do i=1,1440
          outval1(i) = outval3(i)
        enddo
        
      enddo 
          
      write(6,'(/,a,i7,a)') 'Output file written with ',irec
     .                       ,' latitude bands'
      call report_stat('STATUS',module,'utils/dafexpand',outfile
     .       ,'Wrote half-deg grid file',0)
                         
      stop
      end
                                  

      subroutine expand_long(inval,outval)

c     expand a latitude band of 1-deg longtidue complex pairs to half-degree complex pairs
c     Scherneck input values are from 0.5 to 359.5 deg E, output from 0.5 to 360.
                  
      implicit none

      integer*4 i
      
      real*4 inval(720),outval(1440)

c     loop through the input from 0.5E to 359.5E, interpolating and copying
      do i=1,719
        if( mod(i,2).ne.0 ) then
c         odd index in out-longitude: copy from input values
          outval(2*i-1) = inval(i)
          outval(2*i)   = inval(i+1) 
        else
c         even index in out-longitude: average the adjacent input values
          outval(2*i-1) = 0.5 * ( inval(i-1) + inval(i+1) )
          outval(2*i)  =  0.5 * ( inval(i)   + inval(i+2) )
        endif
      enddo

c     interpolate the last pair (0 or 360)
      outval(1439) = 0.5*( inval(719) + inval(1) )
      outval(1440) = 0.5*( inval(720) + inval(2) )

      return
      end


      subroutine expand_lat(outval1,outval3,outval2)

c     compute a latitude band half-way between two others, averaging
c     the values at each longitude

      implicit none

      integer*4 i

      real*4 outval1(1440),outval2(1440),outval3(1440)

      do i=1,1440
        outval2(i) = 0.5*(outval1(i)+outval3(i))
      enddo

      return                            
      end

