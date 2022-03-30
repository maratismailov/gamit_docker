      subroutine open_eop(inut,iut1,ipole)

c Purpose: to open the spatial files (pole, nuttab etc)
c
c NOTE: code borrowed from FILOPN (from ARC) but it was simpler just to
c       use the bits required rather than make GTOG more complicated in
c       order to use the existing routine
c
c 25th July, 1995

      implicit none

       integer*4 inut,iut1,ipole,ioerr,len,rcpar
       character*80 prog_name
       logical fcheck 

c.....get program name calling this routine
      len = rcpar(0,prog_name)

c.....open tai-ut1 table file
      iut1 = 30
      open (unit=iut1,file='ut1.',status='old',iostat=ioerr)
      if (ioerr .ne. 0 ) then
          call report_stat('FATAL',prog_name,'orbits/open_eop','ut1.',
     .    'Error opening UT1 table: ',ioerr)
      endif
c
c.....open nutation table file
c     open only if needed (no nbody file)
      if( .not.fcheck('nbody') ) then
        inut = 31
        open (unit=inut,file='nutabl.',status='old',iostat=ioerr)
        if (ioerr .ne. 0 ) then
          call report_stat('FATAL',prog_name,'orbits/open_eop','nutabl.'
     .                    ,  'Error opening nutation table: ',ioerr)
        endif
      endif 
c
c.....open pole table file    
      ipole = 32
      open(unit=ipole,file='pole.',status='old',iostat=ioerr)
      if (ioerr .ne. 0 ) then
          call report_stat('FATAL',prog_name,'orbits/open_eop','pole.',
     .    'Error opening pole table: ',ioerr)
      endif

      return
      end

