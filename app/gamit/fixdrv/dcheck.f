Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego. 1997.   All rights reserved.

      Subroutine dcheck(dfile,ldfile)

c     Check the D-file list input to FIXDRV and remove any stations
c     for which the X-file is missing (not created by MAKEX after
c     the D-file was created by MAKEXP).  This feature facilitates
c     more graceful batch processing. 

c     R. King   18 December 1997
c     Last modified by M. Floyd (2021/01/13): Changed variable "nsoln"
c       (integer*4) to "gnss" (character*1) to ensure accurate rewriting
c       of d-file header when "edit_needed" is true, following makexp
c       scheme, and replaced skipping of GNSS code and # sessions
c       information with definition of "gnss" and "nsess", respectively.

      Implicit none

      include '../includes/dimpar.h'

c     input and output D-file names
      character*16 dfile,dfile1
                  
c     unit numbers for input and back-up D-files
      integer*4 ldfile,ldfile1  
c     use scratch file number
      parameter(ldfile1=8)

c     number of sessions and stations
      integer*4 nsess,nstat,nstat1
             
c     array of xfile names
      character*16 xfile(maxsit),xcheck1,xcheck2

c     flag and logical array for presence/absence of X-file
      logical edit_needed,lxcheck(maxsit)
               
c     external functions
      logical fcheck    
      integer nblen

c     local variables
c     integer*4 ioerr,nsoln,istat
      integer*4 ioerr,istat
      character*1 gnss
      character*16 lfile,tfile,ifile,jfile 
      character*20 buff20
      character*256 message  

c     open the D-file  
      open( ldfile, file=dfile, iostat=ioerr,status='old')
      if( ioerr.ne.0 ) 
     .   call report_stat('FATAL','FIXDRV','dcheck',dfile
     .                   ,'Error opening D-file',ioerr)

c     Skip the GNSS code and # sessions 
c     read(ldfile,'(1x)')
c     read(ldfile,'(1x)')
      read(ldfile,'(a)') gnss
      read(ldfile,'(i2)') nsess

c     Read the L, T, and I file names
      read(ldfile,'(a16)',iostat=ioerr) lfile  
      if( ioerr.ne.0 ) call report_stat('FATAL','FIXDRV','dcheck',dfile
     .  ,'Error reading L-file name from D-file',ioerr)
      read(ldfile,'(a16)',iostat=ioerr) tfile  
      if( ioerr.ne.0 ) call report_stat('FATAL','FIXDRV','dcheck',dfile
     .  ,'Error reading T-file name from D-file',ioerr)
      read(ldfile,'(a16)',iostat=ioerr) ifile 
      if( ioerr.ne.0 ) call report_stat('FATAL','FIXDRV','dcheck',dfile
     .  ,'Error reading I-file name from D-file',ioerr)

c     initialize the 'edit-needed' flag
      edit_needed = .false.
  
c     read the J file name
       read(ldfile,'(a16)',iostat=ioerr) jfile 
       if( ioerr.ne.0 ) call report_stat('FATAL','FIXDRV','dcheck',dfile
     .     ,'Error reading J-file name from D-file',ioerr) 

c     read the number of stations -- 99 current limit from parameter pointers
        read(ldfile,*,iostat=ioerr) nstat 
        if( ioerr.ne.0 ) then
          call report_stat('FATAL','FIXDRV','dcheck'
     .     ,dfile,'Error reading # stations from D-file',ioerr)  
        elseif( nstat.gt.maxsit ) then   
          write(message,'(a,i2,a)') 
     .       'Number of stations on d-file > maxsit (',maxsit,')'
          call report_stat('FATAL','FIXDRV','dcheck',' ',message,0)
        elseif( nstat.gt.99 ) then
          call report_stat('FATAL','FIXDRV','dcheck',' '
     .     ,'Number of stations on d-file > 99 ',0) 
        endif  

c       initialize the output station count
        nstat1 = nstat
            
c       station loop to check X-files
        do istat = 1,nstat 

          read(ldfile,'(a16)',iostat=ioerr) 
     .                  xfile(istat)
          if( ioerr.ne.0 ) then
            write(message,'(a,i2,a)') 
     .        'Error reading X-file name from D-file (istats='
     .        ,istat,')'
            call report_stat('FATAL','FIXDRV','dcheck',dfile
     .                      ,message,ioerr)   
          endif 
          call lowers(xfile(istat))
          xcheck1 = xfile(istat)
          xcheck2 = xcheck1(1:nblen(xcheck1))//'.Z'
          if( .not.fcheck(xcheck1) .and. fcheck(xcheck2)) then
             call report_stat('WARNING','FIXDRV','dcheck',xcheck2
     .             ,'X-file compressed: ',0)
          endif  
          if( .not.fcheck(xcheck1) .and. .not.fcheck(xcheck2)) then
             call report_stat('WARNING','FIXDRV','dcheck'
     .                ,xfile(istat)
     .                ,'X-file does not exist, remove from D-file: ',0)
             lxcheck(istat) = .false.  
             nstat1 = nstat1  - 1
             edit_needed = .true. 
          else
             lxcheck(istat) = .true.
          endif
        enddo
        if( nstat1.le.1 ) call report_stat('FATAL','FIXDRV'
     .      ,'dcheck',' ','Only one or no existing X-files',0)

      rewind (ldfile) 

c     now write out a new d-file if necessary
      if( edit_needed ) then
c       copy the D-file to a new file-name so that we can overwrite the original 
        dfile1 = " " 
        dfile1 = dfile(1:10) // ".orig"   
        open( ldfile1, file=dfile1, iostat=ioerr,status='unknown') 
        if( ioerr.ne.0 ) 
     .   call report_stat('FATAL','FIXDRV','fixdrv',dfile
     .                   ,'Error opening backup D-file: ',ioerr)  
        do while ( ioerr.eq.0 ) 
          read(ldfile,'(a20)',iostat=ioerr) buff20 
          if( ioerr.ne.0.and.ioerr.ne.-1) call report_stat('FATAL',
     .           'FIXDRV','dcheck','Error copying D-file',' ',ioerr)
          if( ioerr.eq.0) write(ldfile1,'(a20)') buff20  
        enddo  
c       close and reopen the original as a new file--write the header info 
        close(ldfile1)
        close(ldfile)                    
        open( ldfile, file=dfile, iostat=ioerr,status='unknown')  
        if( ioerr.ne.0 ) 
     .   call report_stat('FATAL','FIXDRV','fixdrv',dfile
     .                   ,'Error re-opening D-file',ioerr) 
c       write(ldfile,'(i1)') nsoln 
        write(ldfile,'(a1)') gnss
c       write(ldfile,'(i1)') nsess    
        write(ldfile,'(i2)') nsess    
        write(ldfile,'(a16)') lfile
        write(ldfile,'(a16)') tfile
        write(ldfile,'(a16)') ifile
        write(ldfile,'(a16)') jfile
c       write out only the valid x-files  
        write(ldfile,'(i2)') nstat1
        do istat = 1,nstat
          if( lxcheck(istat) ) then  
             write(ldfile,'(a16)') xfile(istat)
          endif
        enddo
        rewind (ldfile)

      endif         

      return
      end
