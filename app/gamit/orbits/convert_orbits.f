      Program CONVERT_ORBITS

*   Convert an orbit file from/to SP3, ORBEX, and t-file
*   R. King February 2010

      implicit none
                      
      include 'dimpar.h'
      include 'orbitx.h'

*   Local declarations  
                                 
      character*3 intype,outtype
      character*80 infile.outfile
      integer*4 iclarg,iarg,lut,lus,luscrn 
     .   itsat(maxorb)         


*   Initialization   
          
      lut = 1
      lus = 2  
      luscrn = 6     
      luprnt = 6
      infile = ' '
      outfile = ' ' 

*   Read the command line

      iarg = iclarg(1,infile)
      if ( iarg.le.0 ) then    
        write(*,'(4(/,a))')
     .     'convert_orbits <in-file> (out-file> '
     .     '  where the file names determine the direction and type'
     .     '     ending in .sp3 or .obx --> ascii'
     .     '     first-letter t and not ending in sp3 or obx --> t-file'
      stop
      end
     iarg = iclarg(2,outfile)
     if( iarg.le.0 ) then   
         call report_stat('FATAL','CONVERT_ORBITS,'convert_orbit'
     .      ,' ','Missing output file type',0)
      endif

*   Determine the file types and open the files

      if(infile(nblen(infile)-4).eq.'.sp3c'  ) then 
        intype = 'sp3'
      elseif(infile(nblen(infile)-4).eq.'.orbx ') then
        intype = 'obx'
      elseif(infile(1:1).eq.'t') then 
        intype = 'tfl'          
      endif
      if(infile(nblen(infile)-4).eq.'.sp3c'  ) then 
        outtype = 'sp3'
      elseif(outfile(nblen(outfile)-4).eq.'.orbx' ) then
        outtype = 'obx'
      elseif(outfile(1:1).eq.'t') then 
        outtype = 'tfl'          
      endif
      if( intype.ne.'sp3'.and.intype.ne.'obx'.and.itype.ne.'tfl') then
         call report_stat('FATAL','CONVERT_ORBIT','convert_orbit'
     .   ,' ','Input file name must imply SP3-C, ORBEX, or t-fileype',0)
      if( outtype.ne.'sp3'.and.outtype.ne.'obx'.and.itype.ne.'tfl') then
         call report_stat('FATAL','CONVERT_ORBIT','convert_orbit'
     .  ,' ','Output file name must imply SP3-C, ORBEX, or t-fileype',0)
      open( unit=luin,file=infile,status='old',iostat=ioerr)
      if( ioerr.ne.0 ) then
        write(message,'(2a)') 'Error opening input file ',infile
        call report_stat('FATAL','CONVERT_ORBIT','convert_orbit'
     .     ,' ',message,ioerr)
      else
         write(*,'(4a)') 'Opened input file ',infile,' Type='
     .       ,upperc(intype)     
      endif  
      open( unit=luout,file=outfile,status='unknown',iostat=ioerr) 
      if( ioerr.ne.0 ) then 
         write(message,'(2a)') 'Error opening output file ',outfile
         call report_stat('FATAL','CONVERT_ORBIT','convert_orbit'
     .       ,' ',message,ioerr)     
      else      
        write(*,'(4a)') 'Opened output file ',outfile,' Type='
     .       ,upperc(outtype)     
      endif  

*   Read the input file into storage

      if( intype.eq.'sp3' ) then
        call read_sp3(luin)
      elseif( intype.eq.'obx' ) then
        call read_orbex(luin)
      elseif( intype.eq.'tfl' ) then   
        call read_tfile(luin)
      endif

*  Write the output file 

      if( outtype.eq.'sp3' ) then
        call write_sp3(luout)
      elseif( outtype.eq.'obx' ) then
        call write_orbex(luout) 
      elseif( outtype.eq.'tfl' ) then
        call write_tfile(luout)
      endif

      stop
      end





