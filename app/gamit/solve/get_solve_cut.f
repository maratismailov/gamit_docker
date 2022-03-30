c Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
c
      Subroutine GET_SOLVE_CUT

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'  
      include 'parameters.h'

      character*4 snam,upperc,buf4
      character*16 code
      character*120 wcmd
      integer i,i1,i2,i0,ib,ic,type,ioer
      integer lcmd,count_arg,lift_arg
      real*8 tmp_cutoff
c
c     Initialise data arrays
      call zero1d(1,maxsit,elvcut_solve)

c     Get a priori data error for sites
      i2 = 0
      call getcmd(5,'cutoff_elev',wcmd,lcmd,2)
      do 150 i0 = 1,1000
         call getcmd(5,'cutoff',wcmd,lcmd,3)
         if (lcmd.le.0) goto 155
         i2 = i2+1
c        decompose the command line
         ic = count_arg(wcmd)
c        not enough arguments
         if (ic.le.1) goto 150
c        pointer to first argument
         ib = lift_arg(wcmd,code,1)
         if (ib.le.0) goto 150
         type = 1
         if (upperc(code(1:4)).eq.'ALL_') type = 2
         if (type.eq.1) snam = upperc(code(1:4))
         read (wcmd(ib+1:lcmd),*,iostat=ioer)tmp_cutoff
         if( ioer .gt. 0 )  call report_stat('FATAL','SOLVE'
     .       ,'solve_cut',' ','Error reading cutoff values',0)
         do 140 i=1,nsite
            i1=3*(i-1)+1
            if(free(i1).eq.0) go to 140
            if (type.eq.1) then
               buf4=upperc(rlabel(i1)(1:4))
               if (buf4.ne.snam) goto 140
            endif
               elvcut_solve(i)=tmp_cutoff
  140    continue
  150 continue
  155 continue
        if( logprt ) write(6,160)
        write(10,160)
  160   format(/,
     1    4x,'Cutoff elevation angle in SOLVE batch file (degrees): ',/,
     2    'Station                   Cutoff angle   ',/)
        do 165 i=1,nsite
          if( logprt ) write( 6,162) i,rlabel((i-1)*3+1)(1:4),sitnam(i)
     .                  ,elvcut_solve(i) 
          write(10,162) i,rlabel((i-1)*3+1)(1:4),sitnam(i)
     .                  ,elvcut_solve(i)
  162     format(i3,2x,a4,1x,a12,5x,f7.2)
  165   continue
c
      return
      end
