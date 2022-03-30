      Subroutine j_from_nav( batch,afmt,gnss,nprns,iprns,debug,iprndb )

c     write a J-file based on a broadcast navigation file

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'

      logical          debug,first,batch
                  
      character*1      gnss,sys 
      character*60     afmt
      character*256    message

      integer*4        iwkn,iyr,idoy,ihr,min,jprn,iflag
     .               , iyrstart,idoystart,ihrstart,minstart
     .               , iyrend,idoyend,ihrend,minend
     .               , iprns(maxsat),nprns
     .               , icall,ioerr,len,i

c     satellite to debug
      integer iprndb

c     return non-blank length of string

      real*8           sow,sec,utcoff,delta,bephem(16), bclock(6)
     .               , xeaf0,xeaf1,xeaf2,subfr1(8)
     .               , time_end,time_end_save,trans_sow   


c     move this to includes/orbits.h eventually
      integer maxhed
      parameter(maxhed=20)

c
c    RUN INITIALIZATION
c    ******************

c     write a record every half hour
      delta = 1800
      first=.true.
      nprns=0
      time_end_save = 0.d0
      do i=1,maxsat
        iprns(i)=0
      enddo
      iyrstart= 0
      idoystart= 0  
      ihrstart= 0  
      minstart= 0   
c***  This premature since need to print the first epoch in the summary
c      first=.false.               
c** can't do this here:  iwkn, sow, iyr undefined; why done like this? rwk 010314
c      time_end_save = iwkn*605800.d0 + sow
c      time_end = iwkn*605800.d0 + sow 
c      time_end_save = time_end
c      iyrend = iyr
      idoyend = 0
      ihrend = 0
      minend = 0

c-----Begin loop over all records of the ephemeris file

      icall = 0
      sys = gnss
  100 continue

c     Read and decode a record of the navigation file  
      call reade( unav,icall,sys
     .          , iflag,trans_sow,jprn,iwkn,bephem,bclock,subfr1 )
      icall = 1   
      if( sys.eq.gnss ) then 
cd      if(debug .and. jprn .eq. iprndb)
cd     .   print *,jprn,iwkn,bclock(1),bclock(2),bclock(3),bclock(4)
        sow  = bclock(1)
        xeaf0= bclock(2)
        xeaf1= bclock(3) 
        if( gnss.eq.'R' .or. gnss.eq.'S') then
          xeaf2 = 0.d0
        else
          xeaf2= bclock(4)  
        endif   
c       For all systems, week number and sow have been put in GPST; 
c       need UTC for hr/min/sec on j-file      
        call timcon(1,iwkn,sow,iyr,idoy,ihr,min,sec,utcoff)
c       For Glonass, subtract the leap-second for the clock offset from GPST: No, do not!
c        if( gnss.eq.'R' ) then   
c           xeaf0 = xeaf0 - utcoff
c         endif         
c       Write a record of the J-file
c         skip if bad nav-file record
        if( iflag.eq.0 ) then
          write(usvclk,afmt) iyr,idoy,ihr,min,sec,iwkn,sow,gnss,jprn
     .                    , xeaf0,xeaf1,xeaf2                       
        endif

c       Save the first and last epochs and all of the PRN numbers

        if(first) then
           iyrstart=iyr
           idoystart=idoy
           ihrstart=ihr
           minstart=min
           first=.false.
           time_end_save = iwkn*605800.d0 + sow   
           first = .false.
        endif
        time_end = iwkn*605800.d0 + sow
        if( time_end.gt.time_end_save ) then
           iyrend = iyr
           idoyend = idoy
           ihrend = ihr
           minend = min
           time_end_save = time_end
        endif
        do 120 i=1,maxsat
CD           write(6,*) 'i,iprns(i),nprn:',i,iprns(i),nprn
            if( iprns(i).eq.jprn ) goto 130
            if( iprns(i).gt.0 .and. iprns(i).ne.jprn ) goto 120
            if( iprns(i).eq.0 ) then
               nprns = nprns + 1
               iprns(nprns) = jprn
               goto 130
            endif
  120   continue
  130   continue   
      endif

c     loop for another record
      if (iflag .ne. -1) goto 100

c-----end of loop over ephemeris file records

c     End-of-file reach, print summary
      write(message,510) nprns,iyrstart,idoystart,ihrstart,minstart
     .                , iyrend,idoyend,ihrend,minend
  510 format('J-File written for ',i2,' satellites',
     .' Start: ',i4,2x,i4,1x,2i3,' Stop : ',i4,2x,i4,1x,2i3)
      call report_stat('STATUS','MAKEJ','j_from_e',' ',message,0)

      return
      end

