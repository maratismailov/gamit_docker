Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995.  All rights reserved.

C     Get apriori constraints for Earth orientation parameters

      Subroutine GET_EOP_APR

      implicit none
      include '../includes/dimpar.h'
      include 'solve.h'

      real*8 temp(2),temp1(2)
      character*256 message
      character*120 wcmd
      integer i,j,ic
      integer lcmd,count_arg

c
c      Output (in common /constr/ in solve.h)
c          eop_apr(maxsat,maxorb)   real*8      array of tight constaints
c          eop_apr2(maxsat,maxorb)  real*8      array of loose constraints


c--- Initialize eop_apr, eop_apr2, temp and temp1 arrays

      do i = 1,6
        eop_apr(i) = 0.0d0
        eop_apr2(i) = 0.0d0
      enddo
      do i = 1,2
        temp(i)  = 0.0d0
        temp1(i) = 0.0d0
      enddo

c----Pole parameters

      call getcmd(5,'aprior',wcmd,lcmd,2)
      call getcmd(5,'tight_apr_wo',wcmd,lcmd,3)  
c     decompose the command line
      ic = count_arg(wcmd)
      if (ic.eq.2) then
c       read tight pole constraints from command line
        read (wcmd(1:lcmd),*) (temp(i),i=1,ic) 
        if( temp(1).eq.0.d0 .or. temp(2).eq.0.d0 ) 
     .    call report_stat('FATAL','SOLVE','get_eop_apr',' ',
     .    'Pole tight constraints are zero in batch file',0)
      else  
c       not enough arguments--set defaults
        temp(1) = 3.0d0
        temp(2) = 0.3d0
        write(message,'(2a)') 
     .  'Missing or incomplete a priori constraints for pole parameters'
     .     ,'--seting defaults: 3.0 arc, 0.3 arcs/day'
         call report_stat('WARNING','SOLVE','get_eop_apr',' ',message,0)
      endif
      eop_apr(1) = temp(1)
      eop_apr(2) = temp(2)
      eop_apr(3) = temp(1)
      eop_apr(4) = temp(2)   
      if( logprt ) write(6,160)
      write(10,160)
160   format(/,
     1    2x,'A priori pole position errors in arcs and arcs/day',/,
     2    '    Xp      Xp_rate      Yp       Yp_rate ',/)
      if( logprt ) write( 6,161) (eop_apr(j),j=1,4)
      write(10,161) (eop_apr(j),j=1,4)
161   format(1x,4(f8.6,3x))
      call getcmd(5,'aprior',wcmd,lcmd,2)
      call getcmd(5,'loose_apr_wo',wcmd,lcmd,3)  
      if (ic.eq.2) then
c       read loose pole constraints from command line
        read (wcmd(1:lcmd),*) (temp(i),i=1,ic) 
        if( temp(1).eq.0.d0 .or. temp(2).eq.0.d0 ) 
     .    call report_stat('FATAL','SOLVE','get_eop_apr',' ',
     .    'Pole loose constraints are zero in batch file',0)
      else  
c       not enough arguments--set defaults
        temp(1) = 3.0d0
        temp(2) = 0.3d0
      endif
      eop_apr2(1) = temp(1)
      eop_apr2(2) = temp(2)
      eop_apr2(3) = temp(1)
      eop_apr2(4) = temp(2)

c-----UT1 parameters

      call getcmd(5,'aprior',wcmd,lcmd,2)
      call getcmd(5,'tight_apr_ut',wcmd,lcmd,3)    
c     decompose the command line
      ic = count_arg(wcmd)  
      if (ic.eq.2) then
c       read tight ut1 constraints from command line
        read (wcmd(1:lcmd),*) (temp(i),i=1,ic)   
        if( temp(1).eq.0.d0 .or. temp(2).eq.0.d0 ) 
     .     call report_stat('FATAL','SOLVE','get_eop_apr',' ',
     .    'UT1 tight constraints are zero in batch file',0)
      else
c       not enough arguments--set defaults
        temp(1) = 0.00002d0
        temp(2) = 0.02d0
        write(message,'(2a)') 
     .   'Missing or incomplete a priori constraints for UT1 parameters'
     .   ,'--seting defaults: 3.0 arc, 0.3 arcs/day'
         call report_stat('WARNING','SOLVE','get_eop_apr',' ',message,0)
      endif
      eop_apr(5) = temp(1)
      eop_apr(6) = temp(2)
      if( logprt ) write(6,180)
      write(10,180)
180   format(/,
     1    2x,'A priori earth rotation errors in sec and sec/day',/,
     2    '    UT1      UT1_rate',/)
      if( logprt ) write( 6,182) (eop_apr(j),j=5,6)
      write(10,182) (eop_apr(j),j=5,6)
182   format(1x,2(f8.6,3x))
      continue    
      call getcmd(5,'aprior',wcmd,lcmd,2)
      call getcmd(5,'loose_apr_wo',wcmd,lcmd,3) 
      if (ic.eq.2) then
c       read loose UT1 constraints from command line
        read (wcmd(1:lcmd),*) (temp(i),i=1,ic)  
        if( temp(1).eq.0.d0 .or.temp(2).eq.0.d0 ) 
     .    call report_stat('FATAL','SOLVE','get_eop_apr',' ',
     .    'UT1 loose constraints are zero in batch file',0)
      else  
c       not enough arguments--set defaults
        temp(1) = 0.2d0
        temp(2) = 0.02d0
      endif
      eop_apr2(5) = temp(1)
      eop_apr2(6) = temp(2)

      return
      end
