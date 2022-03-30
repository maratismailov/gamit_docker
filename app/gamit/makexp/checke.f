      Subroutine checke ( iue,nprns,gnss_sel,iprns )

c     Check the E-file and return the available satellites

c     R. King   23 February 1992

      implicit none

      include '../includes/dimpar.h'

      logical          debug

      integer*4        iwkn,jprn,iprns(maxsat),nprns
     .               , iflag,iue,icall,i

c     satellite to debug

      real*8 trans_sow, bephem(16), bclock(6), subfr1(8)      

c     move this eventually to includes/orbits
      integer maxhed
      parameter(maxhed=20)   
                         
      character*1 gnss_sel,gnss
      character*256 message

      debug = .false.
                     
c     Initialization

      nprns=0
      do 10 i=1,maxsat
        iprns(i)=0
   10 continue

c-----Begin loop over all records of the ephemeris file

      icall = 0                             
      gnss = gnss_sel
  100 continue
c     Read and decode a record of the E-file 
      call reade(iue,icall,gnss
     .          , iflag,trans_sow,jprn,iwkn,bephem,bclock,subfr1)
      if( iflag.eq..0 ) then
        continue  
      elseif( iflag.gt.0 ) then
        call report_stat('WARNING','MAKEXP','checke',' '
     .          ,'Bad record on nav-file, try reading more ',0)
      endif
      icall = 1

c      if(debug .and. jprn .eq. iprndb) then
c          write(uscren,101)
c     .    ff(20),ff(6),bclock(1),bclock(2),bclock(3),bclock(4)
c  101     format(' Read from E-file: ',6d12.5)
c      endif

      do 120 i=1,maxsat
cd          write(*,*) 'i,jprn,iprns(i):',i,jprn,iprns(i)
          if( gnss.ne.gnss_sel ) goto 130 
          if(  iprns(i).eq.jprn ) goto 130
          if( iprns(i).gt.0 .and. iprns(i).ne.jprn ) goto 120
          if( iprns(i).eq.0 ) then
             nprns = nprns + 1
             if( nprns.gt.maxsat ) then
                write(message,'(  )') 'Number of SVs on nav file ('
     .            ,nprns,') > maxsat (',maxsat,')'
                call report_stat('FATAL','MAKEXP','checke',' '
     .            ,message,0)
             endif
             iprns(nprns) = jprn
             goto 130
          endif
  120 continue
  130 continue

c     loop for another record
      if (iflag .ne. -1) goto 100

c-----end of loop over ephemeris file records

      return
      end
