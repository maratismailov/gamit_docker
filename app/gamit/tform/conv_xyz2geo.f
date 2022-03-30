      Program CONV_XYZ2GEO

c     Program to convert Cartesian coordinates for a single station to geodetic using command-line 
c     arguments.  Based on putgeo.f, R. King  15 October 2002
c
c       Format:  conv_xyz2geo site x y z datum  format
c
c          where site is up to 16 charactes
c                 x,y,z are in meters   
c                 epoch is in decimal years
c                 datum is 5-characters
c                 format is 'dec' (decimal) or 'dms' (deg/min/sec)   
c
c       Example: conv_xyz2geo algo_gps 5378678.345  2434123.987  5267897.541  2002.42 nad83 dec    
               
c     Output written to (1-line) file  conv_xyz2geo.out

      implicit none

      include '../includes/tform.h'

      integer*4 ioerr,idms,numsit,ofile,iarg,iclarg

      real*8 x(3,nsdim),epoch
                        
      character*3 format
      character*5 datum
      character*16 sitnam,arg 

       
c     Get the input values
             
      iarg = iclarg(1,arg)
      if( iarg.le.0 ) then 
         write(*,'(6(/,a))')
     .'Program CONV_XYZ2GEO',
     .'Runstring: conv_xyz2geo <site> <x> <y> <z> <epoch> <datum> <forma
     .t>',
     .'  where <site> is the site-name up to 16 characters',
     .'        <x> <y> <z>  are Cartesian coordintes in meters', 
     .'        <epoch> is in decimal years',
     .'        <datum> is a 5-character datum name (e.g NAD83)',
     .'        <format> is dec for decimal degrees, dms for deg/min/sec'
     .,' ',
     .'Output is one-line file named conv_xyz2geo.out'
        stop
      endif
      
      read(arg,'(a16)',iostat=ioerr) sitnam 
        if( ioerr.ne.0 .or. sitnam(1:1).eq.' ' )
     .      call report_stat('FATAL','TFORM','conv_xyz2geo'
     .          ,' ','Error reading site name from command line', ioerr)
      iarg = iclarg(2,arg) 
      if( iarg.eq.0 ) call report_stat('FATAL','TFORM','conv_xyz2geo'
     .                  ,' ','Missing x coordinate in command line',0)
      read(arg,*,iostat=ioerr) x(1,1)  
      if( ioerr.ne.0 ) call report_stat('FATAL','TFORM','conv_xyz2geo'
     .               ,' ','Error reading x(1) from command line',ioerr )
      iarg = iclarg(3,arg)    
      if( iarg.eq.0 ) call report_stat('FATAL','TFORM','conv_xyz2geo'
     .                 ,' ','Missing y coordinate in command line',0)
      read(arg,*,iostat=ioerr) x(2,1) 
      if( ioerr.ne.0 ) call report_stat('FATAL','TFORM','conv_xyz2geo'
     .               ,' ','Error reading x(2) from command line',ioerr )
      iarg = iclarg(4,arg)  
      if( iarg.eq.0 ) call report_stat('FATAL','TFORM','conv_xyz2geo'
     .                  ,' ','Missing z coordinate in command line',0)
      read(arg,*,iostat=ioerr) x(3,1) 
      if( ioerr.ne.0 ) call report_stat('FATAL','TFORM','conv_xyz2geo'
     .               ,' ','Error reading x(3) from command line',ioerr )
      iarg = iclarg(5,arg)
      if( iarg.eq.0 ) call report_stat('FATAL','TFORM','conv_xyz2geo'
     .                           ,' ','Missing epoch in command line',0)
       read(arg,*,iostat=ioerr) epoch 
       if( ioerr.ne.0 ) call report_stat('FATAL','TFORM','conv_xyz2geo'
     .              ,' ','Error reading epoch from command line',ioerr )   
      iarg = iclarg(6,arg)   
      if( iarg.eq.0 ) call report_stat('FATAL','TFORM','conv_xyz2geo'
     .                     ,' ','Missing datum command line',0)
      read(arg(1:5),'(a5)',iostat=ioerr) datum  
      if( ioerr.ne.0 ) call report_stat('FATAL','TFORM','conv_xyz2geo'
     .              ,' ','Error reading datum from command line',ioerr )
      call uppers(datum) 
      iarg = iclarg(7,arg) 
      if( iarg.eq.0 ) call report_stat('FATAL','TFORM','conv_xyz2geo'
     .                        ,' ','Missing format type command line',0)
      read(arg(1:3),'(a3)',iostat=ioerr) format   
      if( ioerr.ne.0 ) call report_stat('FATAL','TFORM','conv_xyz2geo'
     .             ,' ','Error reading format from command line',ioerr )
   
      
c     Set the parameters for the putgeo call
                   
c     must have L-file format to use deg/min/sec  
      if( format.eq.'dec' ) then
        idms = 1
      else
        idms = 2  
      endif   
      numsit = 1
      ofile = 1                           
      iprnt = 0
      iterm = 0
c     optional
c*      iscrn = 6  
      iscrn = 0 

c     Open the output file
      
      open(unit=ofile,file='conv_xyz2geo.out',status='unknown')


c     Output the coordinates

      call putgeo( idms,numsit,x,sitnam,datum,epoch,ofile )

      stop
      end


