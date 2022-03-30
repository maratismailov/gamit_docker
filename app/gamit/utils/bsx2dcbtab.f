      Program BSX2DCBTAB
                     
c     Create or update a Version 2 GAMIT dcb.dat file from a SINEX DCB file (.bsx)

c     Command-line arguments control the mode:
c       First argument is the input file, either a multi-month BSX file to be 
c        translated to GAMIT version 2 dcb.dat, or an existing version 2 dcb.dat 
c        to be updated
c       Second argument, if given, is the BSX file to be added to get dcb.dat.new
c          if omitted, the program will translate the GSX file, with the output 
c          file named dcb.dat.new.
c     

c     R. King 30 December 2015

      implicit none
                   
* nl    - number of total values in the arrays/dcb.dat
* axnl - dimension of the arrays -- change also in get_times, insert_arrays, sortsvn

      integer*4 maxnl,nl,ng,nr,nc,ne
      parameter(maxnl= 30000) 

* Array values associated with each entry
      
      character*1 gnss(maxnl)
      integer*4 svn(maxnl),prn(maxnl),start(5,maxnl),stop(5,maxnl)
     .        , svng(maxnl),prng(maxnl),startg(5,maxnl),stopg(5,maxnl)
     .        , svnr(maxnl),prnr(maxnl),startr(5,maxnl),stopr(5,maxnl)
     .        , svnc(maxnl),prnc(maxnl),startc(5,maxnl),stopc(5,maxnl)
     .        , svne(maxnl),prne(maxnl),starte(5,maxnl),stope(5,maxnl)
*     Although we write only yr doy hr min, make the time arrays 5 for 
*     for consistency with subroutine itimdif (and other GAMIT routines)

      real*4 dcb(maxnl),rms(maxnl) 
     .     , dcbg(maxnl),rmsg(maxnl)
     .     , dcbr(maxnl),rmsr(maxnl)
     .     , dcbc(maxnl),rmsc(maxnl)
     .     , dcbe(maxnl),rmse(maxnl)

   
* Pointers to SVN within the arrays

      integer*4 indx(maxnl)
                           
* Input and output files
        
* infile - current dcb.dat file, either version 1 of version 2
* incfile - incremental BSX file with values to be added to the input
* outfile - new version 2 dcb.dat file, always named dcb.dat.new
      character*20 infile,incfile,outfile   
      integer*4 luin,luinc,luout
      parameter(luin=1,luinc=2,luout=3)

* Other variables

      character*1 ag
      character*128 line 
      character*256 message
      integer*4 iprn,isvn,year,month,day,doy
     .        , startyr,startdoy,dcbstart(5)
     .        , stopyr,stopdoy,dcbstop(5),svnlast
     .        , irunt(3),iarg,iclarg,ioerr,i,j
      logical eof,endsvs

* Function 
      integer*4 idoy,nblen
      integer*8 itimdif
      logical leapyr 

                    
c  Read the command-line to get the input file name

      iarg = iclarg(1,infile)
      if( iarg.le.0 ) then  
         write(*,'(a)') 'Missing arguments for bsx2dcbtab '
         write(*,'(a)') '  bsx2dcbtab dcb.dat [BSX-file] '
         write(*,'(a)') '     ---output file is dcb.dat.new'
         stop
      endif 
      incfile = ' '   
      iarg = iclarg(2,incfile)
      if( iarg.le.0 ) then
         write(*,'(a )') 'Translating a BSX file to v2 dcb.dat'
      endif       
            
*  Open the input and output files

      open(unit=luin,file=infile,status='old',iostat=ioerr)
      if( ioerr.ne.0 )  then
         call report_stat('FATAL','BSX2DCBTAB','Main',infile
     .                                  ,'Cannot open input file',ioerr) 
      else
          write(*,'(2a)') 'Opened input file ',infile
      endif
      if( incfile(1:1).ne.' ') then 
        open(unit=luinc,file=incfile,status='old',iostat=ioerr)
        if( ioerr.ne.0 ) then
            call report_stat('FATAL','BSX2DCBTAB','Main',' '
     .              ,'Cannot open BSX dcb file',ioerr) 
        else
           write(*,'(2a)') 'Opened monthly file ',incfile
        endif
      endif
      open(unit=luout,file='dcb.dat.newgnss',status='unknown'
     .      ,iostat=ioerr)
      if( ioerr.ne.0 )  call report_stat('FATAL','BSX2DCBTAB','Main',' '
     .      ,'Failure in opening the new file dcb.dat.newgnss',ioerr)

*  Read all the values from the input dcb file into storage

      if( incfile(1:1).eq.' ' ) then 
c-----translating a BSX file 
        eof = .false.
        nl = 0                                       
        do while( .not.eof)
          read(luin,'(a)',iostat=ioerr) line   
cd        print *,'line ',line
          if( ioerr.eq.-1.or.line(1:3).eq.'   '  ) then
            eof = .true.
          elseif( ioerr.ne.0 ) then
            call report_stat('FATAL','BSX2DCBTAB','Main','dcb.dat',
     .            'Error reading BSX DCB file ',ioerr)
          elseif( line(1:4).eq.' DCB'.and.
     .            (line(31:38).eq.'C1C  C1W' .or.
     .             line(31:38).eq.'C1C  C5Q' .or.
     .             line(31:38).eq.'C1C  C1P' .or.   
     ,             line(31:38).eq.'C2I  C6I') ) then
cd          print *,'DCB line ',line
            nl = nl + 1
            call checkmax(nl,maxnl)
            read(line,
     .        '(6x,a1,i3,2x,i2,26x,i2,1x,i3,7x,i2,1x,i3,21x,2f12.4)'
     .        ,iostat=ioerr) gnss(nl),svn(nl),prn(nl)
     .        ,startyr,startdoy,stopyr,stopdoy,dcb(nl),rms(nl)
cd           print *
cd   .       ,'nl gnss svn prn startyr startdoy stopyr stopdoy dcb rms '
cd   .       , nl,gnss(nl),svn(nl),prn(nl)
cd   .       , startyr,startdoy,stopyr,stopdoy,dcb(nl),rms(nl) 
            if(ioerr.ne.0 ) then
               write(message,'(a,i6)') 
     .           'Error reading input file dcb line ',nl
               call report_stat('FATAL','BSX2DCBTAB','Main',infile
     .                            , message,ioerr)
            endif  
            call fix_y2k(startyr)      
            start(1,nl) = startyr
            start(2,nl) = startdoy   
            do i=3,5
              start(i,nl) = 0 
            enddo 
c           SINEX is C1C-C1W, AIUB and GAMIT are P1-C1
            dcb(nl) = -dcb(nl)  
            if( stopyr.eq.0 ) then  
c             0 entry in SINEX denotes no stop date
              stop(1,nl) = 2100
              stop(2,nl) = 1 
            else 
              call fix_y2k(stopyr)
              stop(1,nl) = stopyr
              stop(2,nl) = stopdoy -1  
              stop(3,nl) = 23
              stop(4,nl) = 59
              stop(5,nl) = 0 
            endif
          endif
        enddo

      else
c-------incrementing an existing V2 file
c       read the input dcb.dat values into storage
        eof = .false.
        do while( .not.eof ) 
          read(luin,'(a)',iostat=ioerr) line
          if( ioerr.eq.-1 ) then
             eof = .true.
          elseif( ioerr.ne.0 ) then
            call report_stat('FATAL','BSX2DCBTAB','Main','dcb.dat'
     .             , 'Error reading original DCB file ',ioerr)
          elseif( line(1:1).ne.' ' ) then
c           comment
            continue
          else                
            nl = nl + 1   
            call checkmax(nl,maxnl)
            read(line,'(1x,a1,i6,i4,2(i6,i4,2i3),2f10.3)')  
     .         gnss(nl),svn(nl),prn(nl) 
     .        ,(start(i,nl),i=1,4),(stop(i,nl),i=1,4),dcb(nl),rms(nl)
              start(5,nl) = 0
              stop(5,nl) = 0 
          endif
        enddo                                     
cd      print *,'after reading input file, nl = ',nl 
c       read the new values and add them to the arrays     
        do while( .not.eof)
          read(luinc,'(a)',iostat=ioerr) line   
cd          print *,'line ',line
          if( ioerr.eq.-1.or.line(1:3).eq.'   '  ) then
            eof = .true.
          elseif( ioerr.ne.0 ) then
            call report_stat('FATAL','BSX2DCBTAB','Main','dcb.dat',
     .            'Error reading BSX DCB file ',ioerr)  
          elseif(line(1:4).eq.' DCB'.and.
     .           (line(31:38).eq.'C1C  C1W' .or.
     .            line(31:38).eq.'C1C  C5Q' .or.
     .            line(31:38).eq.'C1C  C1P' .or.   
     ,            line(31:38).eq.'C2I  C6I') ) then
cd            print *,'DCB line ',line
            nl = nl + 1
            call checkmax(nl,maxnl)
            read(line,
     .        '(6x,a1,i3,2x,i2,26x,i2,1x,i3,7x,i2,1x,i3,21x,2f12.0)'
     .        ,iostat=ioerr) gnss(nl),svn(nl),prn(nl)
     .         ,startyr,startdoy,stopyr,stopdoy,dcb(nl),rms(nl)
            if(ioerr.ne.0 ) then
               write(message,'(a,i6)') 
     .           'Error reading input file dcb line ',nl
               call report_stat('FATAL','BSX2DCBTAB','Main',infile
     .                            , message,ioerr)
            endif    
c           change the sign to match the AIUB convention
            dcb(nl) = -dcb(nl)
            call fix_y2k(startyr)
            call fix_y2k(stopyr)
c           Set the stop date of the last entry to be the start date of the new one,
c           and make the stop date open-ended
            stop(1,nl-1) = startyr
            stop(2,nl-1) = startdoy -1
            if(stop(2,nl-1).le.0 ) then
              stop(1,nl-1) = startyr - 1
              if( leapyr(stop(1,nl-1)) ) then
                stop(2,nl-1) = 366
              else
                stop(2,nl-1) = 365  
              endif
            endif                              
            start(1,nl) = startyr
            start(2,nl) = startdoy
            do i=3,5
              start(i,nl) = 0 
            enddo   
            if( stopyr.eq.0 ) then
c             0 stop date denotes open-ended entry
              stop(1,nl) = 2100
              stop(2,nl) = 1 
            else
              stop(1,nl) = stopyr
              stop(2,nl) = stopdoy
              stop(3,nl) = 23
              stop(4,nl) = 59
              stop(5,nl) = 0  
            endif
          endif
        enddo

c     endif on whether translation or incrementing
      endif
                 
* Populate sub-arrays for each system 

      ng = 0 
      nr = 0 
      nc = 0 
      ne = 0         
      do i=1,nl
        if( gnss(i).eq.'G') then
           ng = ng+1 
           call insert_arrays( 
     .          gnss(i),prn(i),svn(i),dcb(i),rms(i),start(1,i),stop(1,i)
     .        , ng,prng,svng,dcbg,rmsg,startg,stopg )
        elseif( gnss(i).eq.'R') then
          nr = nr + 1
          call insert_arrays( 
     .          gnss(i),prn(i),svn(i),dcb(i),rms(i),start(1,i),stop(1,i)
     .        , nr,prnr,svnr,dcbr,rmsr,startr,stopr )
        elseif( gnss(i).eq.'C') then
          nc = nc + 1
          call insert_arrays( 
     .          gnss(i),prn(i),svn(i),dcb(i),rms(i),start(1,i),stop(1,i)
     .        ,  nc,prnc,svnc,dcbc,rmsc,startc,stopc )
        elseif( gnss(i).eq.'E') then
          ne = ne + 1
          call insert_arrays( 
     .          gnss(i),prn(i),svn(i),dcb(i),rms(i),start(1,i),stop(1,i)
     .        , ne,prne,svne,dcbe,rmse,starte,stope ) 
        else
          call report_stat('FATAL','BSX2DCBTAB','Main','dcb.dat'
     .                  , 'Only gnss G R C E coded',0)
        endif 
      enddo
cd    print *,'E values ne 'ne
cd    do i=1,ne
cd      print *,svne(i),prne(i),dcbe(i)
cd   .      ,(starte(j,i),j=1,2),(stope(j,i),j=1,2)
cd    enddo

*  Within each system, sort by SVN
                      
      call sort_svn( ng,svng,prng,dcbg,rmsg,startg,stopg )
      call sort_svn( nr,svnr,prnr,dcbr,rmsr,startr,stopr )
      call sort_svn( nc,svnc,prnc,dcbc,rmsc,startc,stopc )
      call sort_svn( ne,svne,prne,dcbe,rmse,starte,stope )
  
cd    print *,'Sorted E values'
cd    do i=1,ne
cd      print *,svne(i),prne(i),dcbe(i)
cd   .  ,(starte(j,i),j=1,2),(stope(j,i),j=1,2)
cd    enddo
        
* Write out the merged file
                                              
c     get the date of the update and write the headers
      call getdat( irunt(1),irunt(2),irunt(3) )
      write(luout,'(a)')  '* dcb.dat Version 2.0 - units are ns '
      write(luout,'(a,i4,a,i2,a,i2 )')  
     .   '* Last updated ',irunt(1),'-',irunt(2),'-',irunt(3)   
      write(luout,'(a)') 
     .  '*SYS SVN PRN  Start           Stop               P1-C1    rms'
c     write the values
      if( ng.ne.0 ) then
        do i=1,ng
          write(luout,'(1x,a1,i6,i4,2(i6,i4,2i3),2f10.3)')
     .           'G',svng(i),prng(i)
     .           ,(startg(j,i),j=1,4),(stopg(j,i),j=1,4),dcbg(i),rmsg(i)
        enddo
      endif
      if( nr.ne.0 ) then
        do i=1,nr
          write(luout,'(1x,a1,i6,i4,2(i6,i4,2i3),2f10.3)')
     .           'R',svnr(i),prnr(i)
     .           ,(startr(j,i),j=1,4),(stopr(j,i),j=1,4),dcbr(i),rmsr(i)
        enddo
      endif                                                        
      if( nc.ne.0 ) then
        do i=1,nc
          write(luout,'(1x,a1,i6,i4,2(i6,i4,2i3),2f10.3)')
     .           'C',svnc(i),prnc(i)
     .           ,(startc(j,i),j=1,4),(stopc(j,i),j=1,4),dcbc(i),rmsc(i)
        enddo
      endif                                                        
      if( ne.ne.0 ) then
        do i=1,ne
          write(luout,'(1x,a1,i6,i4,2(i6,i4,2i3),2f10.3)')
     .           'E',svne(i),prne(i)
     .           ,(starte(j,i),j=1,4),(stope(j,i),j=1,4),dcbe(i),rmse(i)
        enddo
      endif                                                        

      stop
      end

c------------------------------------------------------------------------------

      Subroutine get_times( ag,iprn,dcbstart,dcbstop
     .                    , isvn,newstart,newstop )

*     Check the iprn for a change midmonth and reset the start and
*     stop times to match the correct SVN
                 
      integer*4 maxnl
      parameter(maxnl= 30000) 

      integer*4 iprn,dcbstart(5),dcbstop(5),isvn,svnstart(5),svnstop(5)
     .        , newstart(5),newstop(5),startdoy,idum
      real*8 rdum
      character*1 ag  
      character*256 message
                  
cd    print *,'ag iprn dcbstart dcbstop ',ag,iprn,dcbstart,dcbstop
      if( iprn.le.32 ) then
* MOD TAH 190702: Added place holder for antpwr to snav_read call
        call svnav_read( -1,dcbstart(1),dcbstart(2),0,0,ag,iprn,isvn
     .                 , idum,rdum,rdum,rdum,rdum, rdum
     .                 , svnstart,svnstop )  
c       if svn is returned as zero, PRN not active at start, search through
c       the month until it's valid:
        if( isvn.eq.0 ) then
          startdoy = dcbstart(2)
          do while( isvn.eq.0 ) 
            startdoy = startdoy + 1
* MOD TAH 190702: Added place holder for antpwr to snav_read call
            call svnav_read( -1,dcbstart(1),startdoy,0,0,ag,iprn,isvn
     .                 , idum,rdum,rdum,rdum,rdum, rdum
     .                 , svnstart,svnstop )  
cd            print *,'DEBUG iprn isvn,yr doy '
cd     .            ,iprn,isvn,dcbstart(1),startdoy 
            if( (startdoy-dcbstart(2)).gt.32 ) then    
cd              print *,'DEBUG iprn isvn dcbstart startdoy '
cd     .           , iprn,isvn,dcbstart(2),startdoy
              write(message,'(a,i3,a,4i5)') 'No valid SVN entry for PRN'
     .             , iprn,' in span '
     .             , dcbstart(1),dcbstart(2),dcbstop(1),dcbstop(2)
              call report_stat('FATAL','BSX2DCBTAB','get_times',' '
     .                  ,message,0)
            endif
          enddo 
          do i=1,5
            newstart(i) = svnstart(i)
            newstop(i) =  dcbstop(i)
          enddo
        else
          do i=1,5   
            newstart(i) = dcbstart(i)
            if( svnstop(1).ne.0 ) then
              if( itimdif(svnstop,dcbstop).lt.0 ) then
                 newstop(i) = svnstop(i)
              else
                newstop(i) = dcbstop(i)    
              endif 
            else
              newstop(i) = dcbstop(i) 
            endif
          enddo             
        endif 
cd  DEBUG
cd       write(*,'(a,2i4,4i5)') 'iprn isvn svnstart svnstop '
cd     .    ,iprn,isvn,(svnstart(i),i=1,2),(svnstop(i),i=1,2)
      else
        print *,'** NEW iprn ',iprn
        iprn = iprn - 50    
* MOD TAH 190702: Added place holder for antpwr to snav_read call
        call svnav_read( -1,dcbstop(1),dcbstop(2),0,0,ag,iprn,isvn
     .                 , idum,rdum,rdum,rdum,rdum, rdum 
     .                 , svnstart,svnstop )  
        do i=1,5
          newstart(i) = svnstart(i)
          newstop(i) = dcbstop(i)
        enddo    
cd        write(*,'(a,2i4,4i5)') 'iprn isvn svnstart svnstop '
cd     .    ,iprn,isvn,(svnstart(i),i=1,2),(svnstop(i),i=1,2)
      endif 

      return
      end

             
c---------------------------------------------------------------------

      Subroutine checkmax(nl,maxnl)
           
      integer*4 nl,maxnl
      character*80 message

      if( nl.gt.maxnl ) then
         write(message,'(a,i7)') 'Number of dcb entries > maxnl ',maxnl
         call report_stat('FATAL','BSX2DCBTAB','checkmax',' ',message,0) 
      endif   
      end                                                                

c-----------------------------------------------------------------------
                                                                        
      Subroutine insert_arrays( cgnss,iprn,isvn,rdcb,rrms,istart,istop 
     .                        , ng,prng,svng,dcbg,rmsg,startg,stopg )

*     Insert one set of values from the master arrays into the GNSS-
*     specific arrays. 
                 
      integer*4 maxnl
      parameter(maxnl= 30000) 
          
*     Counter for the GNSS-specific array
      integer*4 ng

*     Values to be inserted 
      integer*4 iprn,isvn,istart(5),istop(5)
      real*4 rdcb,rrms
      character*1 cgnss

*     GNSS-specific arrays
      integer*4 prng(maxnl),svng(maxnl),startg(5,maxnl),stopg(5,maxnl)
      real*4 dcbg(maxnl),rmsg(maxnl)

*     Local
      integer*4 i
                 
cd      print *,'INSERT ng iprn isvn rdcb rrms istart istop '
cd     .        ,ng,iprn,isvn,rdcb,rrms,istart,istop  
      prng(ng) = iprn
      svng(ng) = isvn
      dcbg(ng) = rdcb
      rmsg(ng) = rrms
      do i=1,5
        startg(i,ng) = istart(i)                                   
        stopg(i,ng) = istop(i)
      enddo      

      return
      end

c-----------------------------------------------------------------------
      Subroutine sort_svn( n,svn,prn,dcb,rms,start,stop )

*     Sort the DCB arrays by SVN 
      
      implicit none
                 
      integer*4 maxnl
      parameter(maxnl= 30000) 


      integer*4 n,i,j,k
     .        , svn(maxnl),prn(maxnl),start(5,maxnl),stop(5,maxnl)
     .        , svnbuf1,svnbuf2,prnbuf1,prnbuf2
     .        , startbuf1(5),stopbuf1(5),startbuf2(5),stopbuf2(5)              

      real*4 dcb(maxnl),rms(maxnl),dcbbuf1,dcbbuf2,rmsbuf1,rmsbuf2
                 
c     Sort by SVN
               
      do i = 1,n-1 
        do j = 1,n-i
           svnbuf1 = svn(j)
           svnbuf2 = svn(j+1)  
           prnbuf1 = prn(j)
           prnbuf2 = prn(j+1)
           dcbbuf1 = dcb(j)
           dcbbuf2 = dcb(j+1)
           rmsbuf1 = rms(j)
           rmsbuf2=  rms(j+1) 
           do k=1,5
             startbuf1(k) = start(k,j)      
             startbuf2(k) = start(k,j+1)
             stopbuf1(k) = stop(k,j)      
             stopbuf2(k) = stop(k,j+1)   
           enddo                                   
           if( svnbuf1.le.svnbuf2 ) then  
             svn(j) = svnbuf1
             svn(j+1) = svnbuf2  
             prn(j) = prnbuf1
             prn(j+1) = prnbuf2
             dcb(j) = dcbbuf1
             dcb(j+1) = dcbbuf2 
             rms(j) = rmsbuf1
             rms(j+1) = rmsbuf2     
             do k=1,5
               start(k,j)   = startbuf1(k)
               start(k,j+1) = startbuf2(k)
               stop(k,j)   = stopbuf1(k)
               stop(k,j+1) = stopbuf2(k)
             enddo
           else    
             svn(j) = svnbuf2
             svn(j+1) = svnbuf1  
             prn(j) = prnbuf2
             prn(j+1) = prnbuf1
             dcb(j) = dcbbuf2
             dcb(j+1) = dcbbuf1 
             rms(j) = rmsbuf2
             rms(j+1) = rmsbuf1     
             do k=1,5
               start(k,j)   = startbuf2(k)
               start(k,j+1) = startbuf1(k)
               stop(k,j)   = stopbuf2(k)
               stop(k,j+1) = stopbuf1(k)
             enddo                   
           endif 
         enddo
      enddo   

      return
      end
                 
c-----------------------------------------------------------------------

      Subroutine fix_y2k(iyear)

c     Check for 2 or 4-digit years.  If 2-digits, convert to 4 by guessing
c     the date.   R. King December 2015 from gamit/lib/check_y2k

      implicit none

      integer*4 iyear,len,rcpar

       character*80 prog_name

      if( iyear.le.1900) then 
                  
c        earliest GPS launch date is 1978; earliest space-geodetic data about 1960 
         if( iyear.gt.60 ) then
            iyear = iyear + 1900
         else
            iyear = iyear + 2000
         endif  

      endif
      return
      end
