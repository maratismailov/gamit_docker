      Program YTOASC

c  Program to dump a binary yaw table to an ascii file
c
c  P Tregoning      8 Dec 1997
c  Mods:  R. King  13 Jan 1998
c         R. King  12 Nov 2015

      implicit none

      include '../includes/dimpar.h'

                           
      character*16 yfile,afile,tfile
      character*20 antbody(maxsat)
      character*256 message

      integer*4 iscrn,iterm,ibin,iascii,ii,isat,nsat
     .         ,itsat(maxsat),isvn(maxsat),iblk(maxsat)
     .         ,iclarg,ioerr,iys,imonth,iday,idoys,ihs,imins
     .         ,iye,idoye,ihe,imine,iepoch,inter,nepoch,idoy,nversn
     .         ,ievent(maxsat),i
      real*8 t_start,t_end,t_send,secs,sece,attit(maxsat)

      ibin = 10
      iascii = 11
      iscrn = 6
      iterm = 5

c     Print the version and machine

      write(message,5)
    5 format('Started YTOASC ',a80)
      call report_stat('STATUS','YTOASC','orbits/ytoasc',' ',message,0)
C

c     get input file name from command line
c     ask if it is not there.
      ii = iclarg(1,yfile)
      if (ii .le. 0) then
         WRITE(ISCRN,11)
11       FORMAT(/,' Enter binary yaw table file name : ')
         READ(ITERM,'(A)') yfile
      endif

c     open the file
      open(ibin,file=yfile,form='unformatted'
     .     ,access='sequential',status='unknown',iostat=ioerr)

c     output file name
      afile = yfile
      afile(1:1)='A'

      OPEN (UNIT=Iascii,FILE=afile,STATUS='UNKNOWN',iostat=ioerr)
      if (ioerr .ne. 0) then
        call report_stat('FATAL','YTOASC','orbits/ytoasc',afile,
     .  'Error opening ASCII Y-table dump: ',ioerr)
      endif
           
c write a title for the ASCII dump
      write(iascii,'(a,a16)') 'ASCII dump of tabular yaw file ',yfile

c read the header record of the yaw file
       read(ibin,iostat=ioerr) nversn  
       backspace(ibin) 
       if( nversn.lt.1051.or.nversn.gt.1061) then
         write(message,'(a,i5,a)') 'Unknown version of y-file',nversn
         call report_stat('FATAL','YTOASC','orbits/ytoasc',afile
     .                   ,message,ioerr)
       endif    
       if( nversn.eq.1051 ) then
         read(ibin,iostat=ioerr) nversn,tfile,t_start,t_end,inter,nepoch
     .          , nsat,(itsat(isat),isvn(isat),iblk(isat),isat=1,nsat)  
       elseif( nversn.eq.1061) then    
         read(ibin,iostat=ioerr) nversn,tfile,t_start,t_end,inter,nepoch
     .         , nsat,(itsat(isat),isvn(isat),antbody(isat),isat=1,nsat)  
       endif
       if(ioerr.ne.0 ) call report_stat('FATAL','YTOASC','orbits/ytoasc'
     .    ,yfile,'Error reading header of binary y-file',ioerr)

c write the first record to the ascii file
       write(iascii,'(a,i4,a,a16,a,f16.7,a,f16.7,a,i4,a,i6)')    
     . ' Version: ',nversn,'  Input T-file: ',tfile,'  Start PJD:'
     .  ,t_start,'  Stop PJD:', t_end,'  Interval:',inter
     .  ,'  No. epochs:',nepoch,'  No. SVs:',nsat 
c      convert to y m d h m s for comparison with X-file    
c      add .001s to JD to avoid roundoff leading to N min 60. sec
       t_start = t_start + 1.157d-8
       call jul2hms(t_start,iys,imonth,iday,ihs,imins,secs) 
       idoys = idoy( iys,imonth,iday) 
c      add .001s to JD to avoid roundoff leading to N min 60. sec
       t_end = t_end + 1.157d-8
       call jul2hms(t_end,iye,imonth,iday,ihe,imine,sece)
       idoye = idoy( iye,imonth,iday) 
       if( sece.eq.60.d0 ) then
         sece = 0.d0              
         imine = imine + 1 
       endif
       
       write(iascii,'(41x,a,2i4,1x,2i3,f5.1,7x,2i4,1x,2i3,f5.1,a)')
     . '( Converted: ',iys,idoys,ihs,imins,secs,iye,idoye,ihe,imine,sece
     . ,' )' 
       write(iascii,'(/,a,9x,a,40i12)')
     . ' Epoch','Yaw angle (deg) for PRN : ',(itsat(i),i=1,nsat)
       write(iascii,'(34x,a,40i12)') ' SVN : ',(isvn(i),i=1,nsat) 
       if( nversn.eq.1051) then
         write(iascii,'(34x,a,40i12)') ' BLK : ',(iblk(i),i=1,nsat)
       elseif( nversn.eq.1061) then
         write(iascii,'(31x,a,40(1x,a11))') 'ANTBODY :     '
     .        ,(antbody(i)(1:11),i=1,nsat)
       endif
              
c now read all the epochs and write them out until the end of file is reached
       ioerr = 0 
       iepoch = 0
       do while (ioerr.eq.0)
         if( nversn.le.970 ) then 
           read(ibin,iostat=ioerr) t_send,(attit(i),i=1,nsat)   
         elseif( nversn.eq.1051.or.nversn.eq.1061 ) then
           read(ibin,iostat=ioerr) t_send,(attit(i),ievent(i),i=1,nsat)
         else
             write(message,'(a,i4)') 'Invalid y-file version ',nversn
             call report_stat('FATAL','YTOASC','orbits/ytoasc',yfile
     .                       ,message,ioerr)
         endif
         if(ioerr.eq.0) then 
           iepoch = iepoch + 1   
c          add .001s to JD to avoid roundoff leading to N min 60. sec
           t_send = t_send + 1.157d-8
           call jul2hms(t_send,iys,imonth,iday,ihs,imins,secs) 
           idoys = idoy( iys,imonth,iday)                     
           if( nversn.le.970 ) then 
              write(iascii,'(i4,f20.6,i5,i4,1x,2i3,f4.0,40f9.2)')
     .       iepoch,t_send,iys,idoys,ihs,imins,secs,(attit(i),i=1,nsat)  
           elseif( nversn.eq.1051.or.nversn.eq.1061 ) then        
              write(iascii,'(i4,f20.6,i5,i4,1x,2i3,f4.0,40(f9.2,i3))')
     .       iepoch,t_send,iys,idoys,ihs,imins,secs
     .        ,(attit(i),ievent(i),i=1,nsat)  
           endif
         elseif (ioerr.eq.-1 ) then
           if( iepoch.eq.nepoch ) then  
             write(message,'(a,a)')' Created file: ',afile
             call report_stat('STATUS','YTOASC','orbits/ytoasc',' '
     .                       ,message,0)
             write(message,'(a)')' Normal stop in YTOASC'
             call report_stat('STATUS','YTOASC','orbits/ytoasc',' '
     .                       ,message,0)
           else
             write(message,'(a,i5)') 'Premature end of y-file at epoch'
     .                              , iepoch+1
             call report_stat('FATAL','YTOASC','orbits/ytoasc',yfile
     .                       ,message,ioerr)
           endif
         else
           write(message,'(a,i5,a)') 'Error reading epoch ',iepoch+1
     .                              ,' of y-file'
           call report_stat('FATAL','YTOASC','orbits/ytoasc',yfile
     .                     ,message,ioerr)
         endif
       enddo
       end

