       SUBROUTINE GETCAR(NUMSIT,X,SITNAM,IFILE)

C       Read in a set or file of Cartesian coordinates.

      implicit none

      include '../includes/tform.h'

      integer*4 numsit,ifile,ioerr,i
      real*8 x(3,nsdim)
                       
      character*8 blank8,siteid
      character*12 sname
      character*16 sitnam(nsdim)
      character*256 line    

      data blank8/'       '/
      
      logical eof


      IF( IFILE.EQ.0 ) THEN
         NUMSIT=1
         WRITE(ISCRN,10)
   10    FORMAT(/,1X,'Enter X,Y,Z in meters')
         READ(ITERM,*,iostat=ioerr) X(1,1),X(2,1),X(3,1)  
         if(ioerr.ne.0 ) call report_stat('FATAL','TFORM','getcar'
     .            ,' ','Error reading input coordinates',ioerr) 
         sitnam(1) = ' '
      ELSE
         numsit = 0  
         eof = .false.
         do while (.not.eof)
           read(ifile,'(a)',iostat=ioerr) line    
           if( ioerr.eq.0 ) then   
             if( line(1:1).eq.' ') then 
              numsit = numsit + 1 
              sitnam(numsit) = ' ' 
              if( numsit.gt.nsdim ) call report_stat('FATAL','TFORM'
     .               ,'getcar',' '
     .                ,'Number of sites in file exceeds dimensions',0)
              read(line(1:9),'(1x,a8)',iostat=ioerr) siteid
              if(ioerr.ne.0 ) call report_stat('FATAL','TFORM','getcar'
     .            ,' ','Error reading site name from input file',ioerr) 
              call getsname(siteid,sname,sitnam(numsit),-2)  
              read(line(10:256),*,iostat=ioerr) (x(i,numsit),i=1,3) 
              if(ioerr.gt.0) then
                 call report_stat('FATAL','TFORM','getcar'
     .          ,' ','Error reading coordinates from input file',ioerr)  
              elseif(ioerr.eq.-1.or.sitnam(numsit)(1:8).eq.blank8 ) then
                 numsit= numsit - 1
              else
                if(iprnt.gt.0) write(iprnt,'(1x,a8,1x,3f15.4)') 
     .              siteid,(x(i,numsit),i=1,3)
              endif
            endif
           elseif( ioerr.eq.-1 ) then
              eof = .true. 
           else
              call report_stat('FATAL','TFORM','getcar'
     .            ,' ','Error reading line from input file',ioerr)  
           endif
         enddo   
      endif

      return
      end
