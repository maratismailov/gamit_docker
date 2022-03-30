      Subroutine simred( iy,im,id,ihr,min,sec )
c
c      Read the simulation-control (S-file) and session.info to setup parameters
c      for dummy observations.  
c       R. King November 1998   

c  Format example (non-blank first column indicates comments)
c   spanid: 2009 125     1     (year  day-of-year  session, free format)
c   elev_cutoff: 15.           (degrees)
c   data_types: 4 L1 L2 P1 P2  (# observables, observables)
c   l1_fact: 1                 (wavelength factor for L1)
c   l2_fact: 1                 (wavelength factor for L2)
c   random noise: .01 .02 .2 .3  (std deviation of random noise for each observable, cycles) 

      implicit none
           
      include '../includes/dimpar.h'
      include '../includes/model.h'

      integer*4 check_flg,iy,im,id,ihr,min,idoy
     .        , count_arg,lift_arg,lcmd,ic,ib,iudisp,trimlen,ioerr,i,j
      real*8 sec,sod,fact
                   
      
      character*2 code,buf2
      character*4 upperc      
      character*60 dispfiln
      character*120 wcmd     
      character*256 line,message,buf60
                               
      logical debug,found
      
c Initialization

      debug = .false.  
      iudisp = 38  
      do i=1,3
       dispneu(i) = 0.d0
      enddo

c  Read the simulation-control file 

      call getcmd(iobs,'spanid',wcmd,lcmd,1)
      if ( lcmd.le.0 )  call report_stat('FATAL','MODEL','simred',' '
     .    ,'Missing spanid in simulation-control file',0)
      ic =count_arg(wcmd)
      if( ic.lt.3) call report_stat('FATAL','MODEL','simred',' '
     .     ,'Spanid error in simulation-control file',0 )
      read(wcmd,*) iy,idoy,isessn
      call check_y2k(iy)


c  Read the span variables from session.info
           
      check_flg = 2.       
      call rsesfo(isesfo,debug,check_flg,iy,idoy,isessn
     .           ,ihr,min,inter,nepoch,nchan,ischan )  
c*      print *,'SIMRED iwkn sow inter nepoch nchan ischan '
c*     .         ,iwkn,sow,inter,nepoch,nchan,ischan
      sec = 0.d0
      sod = ihr*3600.d0 + min*60.d0 + sec
c     need y/m/d h/m/s to match expected epoch from X- or C-files
      call monday(idoy,im,id,iy) 
      call ds2hms(iy,idoy,sod,ihr,min,sec) 
            
c       Get the elevation cutoff   (read this from the MODEL batch file?)

      call getcmd(iobs,'elev_cut',wcmd,lcmd,1)
      if( lcmd.le.0 ) then
c      set the default
        elvcut = 15.d0
      else
        ic = count_arg(wcmd)
        if( ic.eq. 0 ) call report_stat('FATAL','MODEL','simred',' '
     .   ,'Too few arguments for elev_cut in simulation-control file',0)
        read(wcmd,*) elvcut
      endif               
c*      print *,'SIMRED elvcut ',elvcut
               
c Get the observables

      call getcmd(iobs,'data_typ',wcmd,lcmd,1)
      if( lcmd.le.0 ) then
         call report_stat('WARNING','MODEL','simred',' ' 
     .         ,'Data types not input, assume L1 L2 P1 P2',0)
         ndat = 4
         dattyp(1) = 1
         dattyp(2) = 2
         dattyp(3) = 3
         dattyp(4) = 4
      else
         ic = count_arg(wcmd) 
cd         print *,'wcmd ic ',wcmd,ic
         if( ic.eq.0 ) call report_stat('FATAL','MODEL','simred',' '
     .     ,'No arguments for data_type in simulation-control file',0)
         ib = lift_arg(wcmd,code,1) 
cd         print *,'code ib ',code,ib 
         read(code,'(i2)') ndat
c         read(code,*) ndat 
         do i=1,ndat
           ib = lift_arg(wcmd,code,i+1)      
           if( ib.le.0 ) call report_stat('FATAL','MODEL','simred',' '
     . ,'Too few arguments for data_type in simulation-control file,',0)  
           read(code,'(a2)') buf2
           if( buf2.eq.'L1') dattyp(i) = 1
           if( buf2.eq.'L2') dattyp(i) = 2
           if( buf2.eq.'P1'.or.buf2.eq.'C1') dattyp(i) = 3
           if( buf2.eq.'P2') dattyp(i) = 4
         enddo                         
      endif
      call getcmd(iobs,'l1_fact',wcmd,lcmd,1)
      if( lcmd.le.0 ) then      
         call report_stat('WARNING','MODEL','simred',' ' 
     .                   ,'L1 factor not input, assume = 1',0)
         fact = 1
      else
         read(wcmd,*) fact
      endif               
      do i=1,nchan
        lambda(i,1) = -fact
      enddo
      call getcmd(iobs,'l2_fact',wcmd,lcmd,1)
      if( lcmd.le.0 ) then      
         call report_stat('WARNING','MODEL','simred',' ' 
     .                   ,'L2 factor not input, assume = 1',0)
         fact = -1
      else
         read(wcmd,*) fact
      endif               
      do i=1,nchan
        lambda(i,2) = -fact
      enddo   
      if( ndat.ge.3.and.dattyp(3).eq.3 ) then
         do i=1,nchan
           lambda(i,3) = 1
         enddo 
      endif
      if( ndat.ge.4.and.dattyp(4).eq.4 ) then
         do i=1,nchan
           lambda(i,4) = 1
         enddo
       endif
                
c Get the standard deviation of Gaussan random noise for L1, L2, P1, P2   
             
      noise(1) = .03
      noise(2) = .05
      noise(3) = .3
      noise(4) = .5
      call getcmd(iobs,'random',wcmd,lcmd,1)
      if( lcmd.gt.0 ) then
        ic = count_arg(wcmd) 
        if( ic.eq.0 ) call report_stat('FATAL','MODEL','simred',' '
     .    ,'No arguments for random noise in simulation-control file',0)
        read(wcmd,*) (noise(i),i=1,ic)
      endif
     
                 
c Get the displacement file and read it
      
      call getcmd(iobs,'disp',wcmd,lcmd,1)
      if( lcmd.gt.0 ) then
         ic = count_arg(wcmd) 
         if( ic.eq.0 ) call report_stat('FATAL','MODEL','simred',' '
     .  ,'No arguments for disp_file name in simulation-control file',0) 
         dispfiln = wcmd(1:lcmd)
         open( unit=iudisp,file=dispfiln,form='formatted',status='old'
     .       , iostat=ioerr ) 
         if( ioerr.ne.0 ) then
           call report_stat('FATAL','MODEL','simred',dispfiln
     .         ,'Error opening displacement file for simulation',ioerr)  
         else        
           write(message,'(a,a20)') 'Simulation displacement file : '
     .                  ,dispfiln
           call report_stat('STATUS','MODEL','simred',' ',message,0)
         endif
         found = .false.
         do while (.not.found .and. ioerr.eq.0 )
           read(iudisp,'(a)',iostat=ioerr) line 
           if( ioerr.eq.0 .and. line(1:1).eq.' '.and.
     .           trimlen(line).gt.0) then 
              if( upperc(line(2:5)).eq.upperc(sitecd) ) then 
                 buf60 = line(10:69)
                 read(buf60,*,iostat=ioerr) (dispneu(j),j=1,3) 
                 if( ioerr.gt.0) call report_stat('FATAL','MODEL'
     .              ,'simred',' ','Error reading displacements',ioerr)
                 found = .true. 
              endif 
           endif
         enddo  
         if( ioerr.gt.0 ) then
           write(message,'(2a)') 'Bad line in displacement file: '
     .                        ,  line(1:trimlen(line))
           call report_stat('FATAL','MODEL','simred',' ',message,ioerr)
         endif               
      endif

c For now, cancel epoch-by-epoch clocks to avoid complication of no observed pseudorange
     
      klock = 2    
                
c Set mtime to GPST by default
     
      mtime = 2

      return
      end
         

          



                                                      

      
