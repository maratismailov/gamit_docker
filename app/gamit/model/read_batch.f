      Subroutine read_batch( scratch_dir,trash ) 

c     Read the batch file controls   
c     R. King 26 July 2012 

      implicit none

      include '../includes/dimpar.h'   
      include '../includes/units.h' 
      include '../includes/global.h'
      include '../includes/model.h' 

      integer*4 ioerr,idatum,indx
                               
      real*8 value
       
      character*1 trash   
      character*6 etid6
      character*80 line,cvalue,which_line,scratch_dir 
      character*256 message  

c     function
      integer*4 trimlen
      logical fcheck
        
c  Read the GNSS code (replaces obsolete SKD code)
      which_line = 'gnss line'
      read(iterm,'(a)',iostat=ioerr) gnss
      if( gnss.eq.'S' ) gnss = 'G' 
 
c  Read the name of the archive (p-) file
      which_line = 'printfile'
      read(iterm,'(a16)',err=1000,iostat=ioerr) pfiln

c  Read the name of the station clock (i-) file
      which_line ='ifile'
      read(iterm,'(a16)',err=1000,iostat=ioerr) ifiln

c  Read the name of the coordinate (l-) file
      which_line ='lfile'
      read(iterm,'(a16)',err=1000,iostat=ioerr) lfiln

c Read the name of the observation (x- or c-) or simulation file
      which_line ='input x, c, or s-file'
      read(iterm,'(a16)',err=1000,iostat=ioerr) obfiln
   
c  Read the name of the output observation (c-) file and scratch directory
      which_line ='Output C-file and scratch directory'
      read(iterm,'(a)',iostat=ioerr,err=1000 ) line
      call multiread(line,indx,'CH',ioerr,value,cvalue,1) 
      if( ioerr.ne.0 ) goto 1000
      cfiln = cvalue(1:trimlen(cvalue)) 
      if( cfiln.eq.obfiln ) call report_stat('FATAL','MODEL'
     .     ,'read_batch',cfiln,'Output C-file name same as input: ',0)
      call multiread(line,indx,'CH',ioerr,value,cvalue,1) 
      if( ioerr.ne.0 ) goto 1000
      if ( cvalue(1:trimlen(cvalue)) .ne. 'output' .and.
     .     cvalue(1:trimlen(cvalue)) .ne. 'Output' ) then
          scratch_dir = cvalue(1:trimlen(cvalue))
      else
        scratch_dir = 'NONE'
      endif   

c  Read the control for deleting c-files
      which_line ='delete cfiles'
      read(iterm,'(a1)',err=1000,iostat=ioerr) trash
 
c  Read the name of the orbit (t-) file
      which_line ='tfile'
      read(iterm,'(a16)',err=1000,iostat=ioerr) tfiln
                                                 
c  Read the source of the higher-order ionospheric corrections 
c     Only a global IONEX (f-file) map is currently coded
      which_line = 'ion source'
      read(iterm,'(a10,1x,a6)',err=1000,iostat=ioerr) ionfiln,magfield

c  Read the name of the RINEX met file
      which_line = 'RINEX met file'
      read(iterm,'(a16)',err=1000,iostat=ioerr) metfiln

c  Read the name of the output met values (z-) file
      which_line = 'RINEX met file'
      read(iterm,'(a16)',err=1000,iostat=ioerr) zfiln

c  Read the name of the satellite clock (j-) file
      which_line = 'jfile'
      read(iterm,'(a16)',err=1000,iostat=ioerr) jfiln
                                                 
c  Read the geodetic datum, tide, EOP, and loading models requested
      which_line = 'tide controls'
* MOD TAH 200219: Changed to free format read to accomodate 
*     extra digit in ietide optios (IERS20 mean pole).
C     read(iterm,'(i1,1x,i2,1x,i2,1x,a6,1x,a1,1x,a1)',err=1000
C    .          ,iostat=ioerr)
C    .     idatum,ietide,isptide,etid6,atmlflg,hydrlflg     
      read(iterm,*,err=1000,iostat=ioerr)
     .     idatum,ietide,isptide,etid6,atmlflg,hydrlflg  
   
      etidemod = etid6

c  Read the antenna models
      which_line = 'antenna models'
      read(iterm,'(a1,5x,a4,1x,a4)',err=1000,iostat=ioerr) 
     .     epcvflg,antmod_in,svantmod_in                 
c     old batch files may have a numerical value of the elevation
c     cutoff angle (no longer used) in the first 4 spaces
      if( epcvflg.ne.'Y' .and. epcvflg.ne.'N' ) epcvflg = 'N'
                
c  Read the clock model and yaw (y-) file name
      which_line = 'clock and yaw'                    
      read(iterm,'(i1,1x,a10)',iostat=ioerr) klock,yfiln

c  Read the met model options
      which_line = 'met options'
      read(iterm,'(a18)',err=1000,iostat=ioerr) metopts

c  Read the met models
      which_line = 'met models'
      read(iterm,'(a4,3(1x,a4))',err=1000,iostat=ioerr) 
     .         dryzen,wetzen,drymap,wetmap  
c     if GMF requested by GPT2 available, use GPT2
      if( fcheck('gpt.grid') ) then
         call report_stat('STATUS','MODEL','read_batch',' '
     .     , 'Replacing GMF with GPT for mapping functions' ,ioerr)
         if( drymap.eq.'GMFH' ) drymap = 'GPT '
         if( wetmap.eq.'GMFW' ) wetmap = 'GPT '
      endif
    
 1000 if (ioerr .ne. 0) then
         if( ioerr.eq.-1 ) then   
           write(message,1005) which_line
 1005      format('Check FIXDRV run, EOF found at input B-file line: '
     .     ,a80)
           write(iprnt,1010)message
 1010      format(//,A80)
           call report_stat('FATAL','MODEL','read_batch',' '
     .                     ,message,ioerr)
         endif 
         write(message,1015)which_line
 1015    format('Check FIXDRV run: Batch file error at B-file line: '
     .   ,a80)
         call report_stat('FATAL','MODEL','read_batch',' '
     .                   ,message,ioerr)
      endif 

      end


