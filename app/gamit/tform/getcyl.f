       SUBROUTINE GETCYL(NUMSIT,X,SITNAM,IFILE)
                   
C       Read in a set or file of cylindrical coordinates.

      implicit none

      include '../includes/tform.h'

      integer*4 numsit,ifile,ioerr
      real*8 x(3,nsdim),long,pi,rho,z
                           
      character*8 siteid
      character*12 sname
      character*16 sitnam(nsdim)
      character*256 line    

      logical eof

      IF( IFILE.EQ.0 ) THEN
        NUMSIT=1
        WRITE(ISCRN,10)     
   10   FORMAT(/,1X,'Enter RHO, LONG(E), Z in deg or dms, m')    
        READ(ITERM,*) RHO,LONG,Z
        X(1,1)= RHO*DCOS(LONG*PI/180.D0)
        X(2,1)= RHO*DSIN(LONG*PI/180.D0)
        X(3,1)= Z
        sitnam(1) = '                '
      ELSE
         numsit = 0  
         eof = .false.
         do while (.not.eof)
           read(ifile,'(a)',iostat=ioerr) line
           if( ioerr.eq.0 ) then   
             if( line(1:1).eq.' ') then 
               numsit = numsit + 1  
               if( numsit.gt.nsdim ) call report_stat('FATAL','TFORM'
     .               ,'getcyl',' '
     .                ,'Number of sites in file exceeds dimensions',0)
               read(line(1:9),'(1x,a8)',iostat=ioerr) siteid  
               call getsname(siteid,sname,sitnam(numsit),-2)
               if(ioerr.ne.0 ) call report_stat('FATAL','TFORM','getcyl'
     .            ,' ','Error reading site name from input file',ioerr) 
               read(line(10:256),*,iostat=ioerr) rho,long,z
               if(ioerr.ne.0) call report_stat('FATAL','TFORM','getcyl'
     .           ,' ','Error reading coordinates from input file',ioerr)   
               if(iprnt.gt.0) write(iprnt,'(1x,a8,1x,3f15.4)') 
     .               sitnam(numsit)(1:8),rho,long,z
               X(1,numsit)= RHO*DCOS(LONG*PI/180.D0)
               X(2,numsit)= RHO*DSIN(LONG*PI/180.D0)
               X(3,numsit)= Z
             endif
           elseif( ioerr.eq.-1 ) then
              eof = .true. 
           else
              call report_stat('FATAL','TFORM','getcyl'
     .            ,' ','Error reading line from input file',ioerr)  
           endif
         enddo   
      endif

      return
      end
