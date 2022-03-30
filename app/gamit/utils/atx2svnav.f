      Program ATX2SVNAV
                       
c     Read an ANTEX file to generate an svnav.dat file

c     Two command-line arguments: input file and output file

c     R. King 28 August 2015

      implicit none
                   
* nl    - number of total values in the arrays/svnav.dat
* maxnl - dimension of the arrays -- change also in adj_times, insert_arrays, sortsvn, sort_times

      integer*4 maxnl,nl,ng,nr,nc,ne
      parameter(maxnl= 30000) 

* Array values associated with each entry
      
      character*1 gnss(maxnl)
      character*20 antbody(maxnl),antbodyg(maxnl),antbodyr(maxnl)
     .           , antbodyc(maxnl),antbodye(maxnl)
      integer*4 svn(maxnl),prn(maxnl),start(5,maxnl),stop(5,maxnl)
     .        , svng(maxnl),prng(maxnl),startg(5,maxnl),stopg(5,maxnl)
     .        , svnr(maxnl),prnr(maxnl),startr(5,maxnl),stopr(5,maxnl)
     .        , svnc(maxnl),prnc(maxnl),startc(5,maxnl),stopc(5,maxnl)
     .        , svne(maxnl),prne(maxnl),starte(5,maxnl),stope(5,maxnl)  
*     Although we write only yr doy hr min, make the time arrays 5 for 
*     for consistency with subroutine itimdif (and other GAMIT routines)

   
* Pointers to SVN within the arrays

      integer*4 indx(maxnl)
                           
* Input and output files
        
* infile - ANTEX file name
* outfile - svnav.dat file name
      character*20 infile,outfile   
      integer*4 luin,luout
      parameter(luin=1,luout=2)

* Other variables

      character*1 yawbias
      character*80 line 
      character*256 message
      integer*4 iprn,isvn,year,month,day,doy,hr,min
     .        , chan,irunt(3),iarg,iclarg,ioerr,i,j
      real*4 mass,yawrate 
      real*8 sec 
      logical eof,foundant,founddate               
     

* Function 
      integer*4 idoy  
      integer*8 itimdif

                    
c  Read the command-line to get the input file name

      iarg = iclarg(1,infile)
      if( iarg.le.0 ) then  
         write(*,'(a)') 'Missing arguments for atx2svnav '
         write(*,'(a)') '  atx2svnav [ANTEX-file] [svnav.dat file]'
         stop
      endif 
      iarg = iclarg(2,outfile)
      if( iarg.le.0 ) then
         write(*,'(a )') 'Missing output svnav.dat file name'
      endif       
            
*  Open the input and output files

      open(unit=luin,file=infile,status='old',iostat=ioerr)
      if( ioerr.ne.0 )  then
         call report_stat('FATAL','ATX2SVNAV','Main',infile
     .                                  ,'Cannot open input file',ioerr) 
      else
          write(*,'(2a)') 'Opened input file ',infile
      endif
      open(unit=luout,file=outfile,status='unknown'
     .      ,iostat=ioerr)
      if( ioerr.ne.0 ) then
          call report_stat('FATAL','ATX2SVNAV','Main'
     .      ,outfile,'Failure in opening the output file',ioerr)
      else
         write(*,'(2a)') 'Opened output file ',outfile
      endif
                                   
*  Initialize the seconds field of the times
    
      do i=1,maxnl
        start(5,i) = 0
        stop(5,i) = 0               
      enddo                                                               

*  Read all the values from the input ANTEX file into storage   
                                   
      eof = .false.
      nl = 0 
      do while( .not.eof )   
        foundant = .false.
c       look for 'START OF ANTENNA'
        do while( .not.foundant )
           read(luin,'(a)',iostat=ioerr) line       
           if( ioerr.eq.-1 ) goto 100
           if( line(61:76).eq.'START OF ANTENNA' ) then
             read(luin,'(a)') line
             if( line(61:64).eq.'TYPE') then
c               columns 21-44 are empty for ground antennas
               if(line(21:21).ne.' ' ) then
                 foundant = .true.
                 nl = nl+1 
                 read(line,'(a20,a1,i2,18x,i3)') 
     .                 antbody(nl),gnss(nl),prn(nl),svn(nl)
cd                 print *,'nl ',nl,antbody(nl),gnss(nl),prn(nl),svn(nl)
               endif
             endif
           endif
        enddo
        if( foundant ) then
          founddate = .false.
          do while( .not.founddate ) 
            read(luin,'(a)',iostat=ioerr) line
            if( line(61:70).eq.'VALID FROM' ) then
cd              print *,'VALID FROM',line
              founddate = .true.
              read(line,'(5i6,f13.0)',iostat=ioerr) 
     .                        year,month,day,hr,min,sec
              doy = idoy(year,month,day)
              start(1,nl) = year
              start(2,nl) = doy
              start(3,nl) = hr
              start(4,nl) = min   
cd              print *,'nl start ',nl,(start(i,nl),i=1,4)
c             assume the next line is stop date, if present
              read(luin,'(a)',iostat=ioerr ) line
              if( line(61:71).eq.'VALID UNTIL') then  
cd                print *,'VALID UNTIL',line
                read(line,'(5i6,f13.0)',iostat=ioerr) 
     .                        year,month,day,hr,min,sec
                doy = idoy(year,month,day)
                stop(1,nl) = year
                stop(2,nl) = doy
                stop(3,nl) = hr
                stop(4,nl) = min        
              else
                stop(1,nl) = 2100
                stop(2,nl) = 1
                stop(3,nl) = 0 
                stop(4,nl) = 0 
              endif
            endif
          enddo                                 
        endif
      enddo      

* Populate sub-arrays for each system 

 100  ng = 0 
      nr = 0 
      nc = 0 
      ne = 0 
      do i=1,nl   
        if( gnss(i).eq.'G') then
           ng = ng+1 
           call insert_arrays( 
     .          gnss(i),prn(i),svn(i),antbody(i),start(1,i),stop(1,i)
     .        , ng,prng,svng,antbodyg,startg,stopg )
        elseif( gnss(i).eq.'R') then
          nr = nr + 1
          call insert_arrays( 
     .          gnss(i),prn(i),svn(i),antbody(i),start(1,i),stop(1,i)
     .        , nr,prnr,svnr,antbodyr,startr,stopr )
        elseif( gnss(i).eq.'C') then
          nc = nc + 1
          call insert_arrays( 
     .          gnss(i),prn(i),svn(i),antbody(i),start(1,i),stop(1,i)
     .        ,  nc,prnc,svnc,antbodyc,startc,stopc )
        elseif( gnss(i).eq.'E') then
          ne = ne + 1
          call insert_arrays( 
     .          gnss(i),prn(i),svn(i),antbody(i),start(1,i),stop(1,i)
     .        , ne,prne,svne,antbodye,starte,stope ) 
        else
          call report_stat('WARNING','ATX2SVNAV','Main',' '
     .                  , 'Only gnss G R C E coded',0)
        endif 
      enddo


*  Within each system, sort by SVN and time
                      
      if( ng.gt.0 ) then
        write(*,*) 'GPS unsorted (PRN-priotity)'
        print *,'ng ',ng 
        do i=1,ng
          write(*,*) prng(i),svng(i)
     .       ,(startg(j,i),j=1,5),(stopg(j,i),j=1,5)
        enddo
        call sort_svn( ng,svng,prng,antbodyg,startg,stopg )
      endif
      if( nr.gt.0 ) then
        write(*,*) 'Glonass unsorted (PRN-priotity)'
        do i=1,ng
         write(*,*) prng(i),svng(i)
     .      ,(startg(j,i),j=1,5),(stopg(j,i),j=1,5)
        enddo  
        call sort_svn( nr,svnr,prnr,antbodyr,startr,stopr )
      endif
      if( nc.gt.0 ) then
        call sort_svn( nc,svnc,prnc,antbodyc,startc,stopc )
      endif
      if( ne.gt.0 ) then
        call sort_svn( ne,svne,prne,antbodye,starte,stope )
      endif
  
cd      print *,'Sorted G values'
cd      do i=1,ng
cd        print *,svng(i),prng(i),dcbg(i)
cd     .  ,(startg(j,i),j=1,2),(stopg(j,i),j=1,2)
cd      enddo
                                                                                   
*  In order to allow ARC integrations onto adjacent days, where there may
*  be no valid transmission, decrement the start times by one day and
*  increment the stop times by one day, checking for overlaps
* RWK/TAH 150831: We've decided not to do this

      if( ng.gt.0 ) call adj_times( ng,'G',svng,startg,stopg ) 
      if( nr.gt.0 ) call adj_times( nr,'R',svnr,startr,stopr ) 
      if( nc.gt.0 ) call adj_times( nc,'C',svnc,startc,stopc ) 
      if( ne.gt.0 ) call adj_times( ne,'E',svne,starte,stope ) 
                                     

* Write out the merged file
                                              
      write(luout,'(a)') 'svnav.dat  Version  2.0'  
      write(luout,'(3a)') 'SYS SVN  PRN CHAN ANT/BODY               '
     .   , 'MASS(G) YAW BIAS  YAW RATE  START           STOP       '
     .   , '      COMMENTS'       
      if( ng.ne.0 ) then  
         do i=1,ng
c**        read mass and yaw values from a previous svnav.dat
           mass = 1000.
           yawbias = 'Y'
           yawrate = 0.1
           chan = 0 
           write(luout,'(1x,a1,2x,i4,2x,i2,2x,i2,2x,a20,f11.0,5x
     .                 , a1,f12.4,2(i6,i4,2i3),3x,a40)')  
     .         'G',svng(i),prng(i),chan,antbodyg(i),mass,yawbias,yawrate
     .       , (startg(j,i),j=1,4),(stopg(j,i),j=1,4)  
         enddo
       endif
            
      if( nr.ne.0 ) then  
         do i=1,nr           
c**        compute channel and add mass and yaw values  
           mass = 2000.
           yawbias = 'Y'
           yawrate = 0.2
           chan = 1
           write(luout,'(1x,a1,2x,i4,2x,i2,2x,i2,2x,a20,f11.0,5x
     .                 , a1,f12.4,2(i6,i4,2i3),3x,a40)')  
     .         'R',svnr(i),prnr(i),chan,antbodyr(i),mass,yawbias,yawrate
     .       , (startr(j,i),j=1,4),(stopr(j,i),j=1,4)  
         enddo
       endif
            
      if( nc.ne.0 ) then  
        do i=1,nc           
c**       mass and yaw values  
          mass = 3000.
          yawbias = 'Y'
          yawrate = 0.3
          chan = 0 
          write(luout,'(1x,a1,2x,i4,2x,i2,2x,i2,2x,a20,f11.0,5x
     .                , a1,f12.4,2(i6,i4,2i3),3x,a40)')  
     .        'C',svnc(i),prnc(i),chan,antbodyc(i),mass,yawbias,yawrate
     .      , (startc(j,i),j=1,4),(stopc(j,i),j=1,4)  
        enddo
      endif

      if( ne.ne.0 ) then
        do i=1,ne           
c**       add mass and yaw values   
          mass = 4000.
          yawbias = 'Y'
          yawrate = 0.4
          chan = 0 
          write(luout,'(1x,a1,2x,i4,2x,i2,2x,i2,2x,a20,f11.0,5x
     .                , a1,f12.4,2(i6,i4,2i3),3x,a40)')  
     .        'E',svne(i),prne(i),chan,antbodye(i),mass,yawbias,yawrate
     .      , (starte(j,i),j=1,4),(stope(j,i),j=1,4)  
        enddo
      endif

      stop
      end

c---------------------------------------------------------------------

      Subroutine checkmax(nl,maxnl)
           
      integer*4 nl,maxnl
      character*80 message

      if( nl.gt.maxnl ) then
         write(message,'(a,i7)') 'Number of dcb entries > maxnl ',maxnl
         call report_stat('FATAL','ATX2SVNAV','checkmax',' ',message,0) 
      endif   
      end                                                                

c-----------------------------------------------------------------------
                                                                        
      Subroutine insert_arrays( cgnss,iprn,isvn,antbody,istart,istop 
     .                        , ng,prng,svng,antbodyg,startg,stopg )

*     Insert one set of values from the master arrays into the GNSS-
*     specific arrays. 
                 
      integer*4 maxnl
      parameter(maxnl= 30000) 
          
*     Counter for the GNSS-specific array
      integer*4 ng

*     Values to be inserted 
      integer*4 iprn,isvn,istart(5),istop(5)
      character*1 cgnss

*     GNSS-specific arrays
      integer*4 prng(maxnl),svng(maxnl),startg(5,maxnl),stopg(5,maxnl)
      character*20 antbody,antbodyg(maxnl)

*     Local
      integer*4 i
                 
cd      print *,'INSERT ng iprn isvn rdcb rrms istart istop '
cd     .        ,ng,iprn,isvn,rdcb,rrms,istart,istop  
      prng(ng) = iprn
      svng(ng) = isvn             
      antbodyg(ng) = antbody
      do i=1,5
        startg(i,ng) = istart(i)                                   
        stopg(i,ng) = istop(i)
      enddo      

      return
      end

c-----------------------------------------------------------------------
      Subroutine sort_svn( n,svn,prn,antbody,start,stop )

*     Sort the arrays by SVN and time 
      
      implicit none
                 
      integer*4 maxnl
      parameter(maxnl= 30000) 


      integer*4 n,i,j,k,ns
     .        , svn(maxnl),prn(maxnl),start(5,maxnl),stop(5,maxnl)
     .        , svnbuf1,svnbuf2,prnbuf1,prnbuf2
     .        , startbuf1(5),stopbuf1(5),startbuf2(5),stopbuf2(5)              
     .        , istrt,iend,svnlast
      character*20 antbody(maxnl),antbuf1,antbuf2

c     Sort by SVN
               
      do i = 1,n-1 
        do j = 1,n-i
           svnbuf1 = svn(j)
           svnbuf2 = svn(j+1)  
           prnbuf1 = prn(j)
           prnbuf2 = prn(j+1)                    
           antbuf1 = antbody(j)
           antbuf2 = antbody(j+1)
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
             antbody(j) = antbuf1
             antbody(j+1) = antbuf2
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
             antbody(j) = antbuf2
             antbody(j+1) = antbuf1
             do k=1,5
               start(k,j)   = startbuf2(k)
               start(k,j+1) = startbuf1(k)
               stop(k,j)   = stopbuf2(k)
               stop(k,j+1) = stopbuf1(k)
             enddo                   
           endif 
         enddo
      enddo      
                 
cd      print *,'SVN-sorted values '
cd      do i=1,n
cd        write(*,*) 'i svn prn start ',i,svn(i),prn(i),(start(j,i),j=1,5)
cd      enddo
    
* Now sort by time within each SVN

      istrt = 1
      iend = n         
      ns = 1
      svnlast = svn(1)
      do i=1,n-1 
cd        print *,'MAIN i svn svnlast ',i,svn(i),svnlast 
        if( svn(i+1).eq.svnlast ) then
          ns = ns+1 
          iend = i+1  
cd          print *,'ns iend ',ns,iend
        else        
          if( ns.gt.1 ) then     
            istrt = iend - ns + 1
cd            print *,'ns istrt iend ',ns,istrt,iend
            call sort_times(istrt,iend,svn,prn,antbody,start,stop)  
            ns = 1 
          endif
        endif
        svnlast = svn(i+1)
      enddo

      return
      end
c------------------------------------------------------------------

      Subroutine sort_times(istart,iend,svn,prn,antbody,start,stop)

c     Sort an SVN subgroup of the array by start time of the entry
                                                                        
      implicit none
                                 
      integer*4 maxnl
      parameter(maxnl=30000)

      integer*4 istart,iend,svn(maxnl),prn(maxnl)
     .        , start(5,maxnl),stop(5,maxnl)
     .        , startbuf1(5),startbuf2(5),stopbuf1(5),stopbuf2(5)
     .        , svnbuf1,svnbuf2,prnbuf1,prnbuf2,i,j,k
      character*20 antbody(maxnl),antbuf1,antbuf2 

c     Function
      integer*8 itimdif
                 
cd      print *,'svn istart iend ',svn(istart),istart,iend
cd      do i=istart,iend
cd        write(*,*) i,svn(i),(start(j,i),j=1,5)
cd      enddo
      do i = istart,iend-1        
        do j = istart,istart+iend-i-1 
cd           print *,'DEBUG j ',j
           svnbuf1 = svn(j)
           svnbuf2 = svn(j+1)  
           prnbuf1 = prn(j)
           prnbuf2 = prn(j+1)                    
           antbuf1 = antbody(j)
           antbuf2 = antbody(j+1)
           do k=1,5
             startbuf1(k) = start(k,j)      
             startbuf2(k) = start(k,j+1)
             stopbuf1(k) = stop(k,j)      
             stopbuf2(k) = stop(k,j+1)   
           enddo  
cd           print *,'startbuf1 startbuf2 '
cd     .          ,(startbuf1(k),k=1,5),(startbuf2(k),k=1,5)
           if( itimdif(startbuf1,startbuf2).lt.0 ) then
cd             print *,'2 > 1 no swap '
             svn(j) = svnbuf1
             svn(j+1) = svnbuf2  
             prn(j) = prnbuf1
             prn(j+1) = prnbuf2  
             antbody(j) = antbuf1
             antbody(j+1) = antbuf2
             do k=1,5
               start(k,j)   = startbuf1(k)
               start(k,j+1) = startbuf2(k)
               stop(k,j)   = stopbuf1(k)
               stop(k,j+1) = stopbuf2(k)
             enddo
           else                      
cd             print *,'2 < 1 swap '
             svn(j) = svnbuf2
             svn(j+1) = svnbuf1  
             prn(j) = prnbuf2
             prn(j+1) = prnbuf1    
             antbody(j) = antbuf2
             antbody(j+1) = antbuf1
             do k=1,5
               start(k,j)   = startbuf2(k)
               start(k,j+1) = startbuf1(k)
               stop(k,j)   = stopbuf2(k)
               stop(k,j+1) = stopbuf1(k)
             enddo                   
           endif 
         enddo
      enddo    

cd      print *,'sorted'              
cd      do i=istart,iend
cd        write(*,*) i,svn(i),(start(j,i),j=1,5)
cd      enddo

      return
      end  

c-----------------------------------------------------------------

      Subroutine adj_times( n,gnss,svn,start,stop )

c     Decrement the start time by 6 hours and increment the stop time
c     by 6 hours to allow for interpolation in ARC; check that no overlap

      implicit none
                                              
      integer maxnl
      parameter(maxnl=30000) 
                          
      character*1 gnss
      integer*4 n,svn(maxnl),start(5,maxnl),stop(5,maxnl),jd,mon,day,i,j
      real*8 sod,sec      

c     function                  
      integer*4 julday
      integer*8 itimdif

c       First decrement the start times increment the stop times
c*** No, we want do this step, but rather pad the search in svnav_read
c     do i= 1,n                                      
c       call monday(start(2,i),mon,day,start(1,i))
c       jd = julday(mon,day,start(1,i))
c       sod = start(3,i)*3600.d0 + start(4,i)*60.d0
c       call timinc(jd,sod,-21600.d0) 
c       call dayjul(jd,start(1,i),start(2,i))
c       call ds2hms(start(1,i),start(2,i),sod,start(3,i),start(4,i),sec)
c       call monday(stop(2,i),mon,day,stop(1,i))
c       jd = julday(mon,day,stop(1,i)) 
c       sod = stop(3,i)*3600.d0 + stop(4,i)*60.d0
c       call timinc(jd,sod,21600.d0)  
c       call dayjul(jd,stop(1,i),stop(2,i))
c       call ds2hms(stop(1,i),stop(2,i),sod,stop(3,i),stop(4,i),sec)
c     enddo
                                                                      
c       However, we will round the stop times to the day boundary if 23 59 
      print *,'in adjust times n ',n
      do i=1,n           
        if(i.eq.1 ) then
          print *,'stop ',(stop(j,i),j=1,5)
        endif
        if( stop(3,i).eq.23.and.stop(4,i).eq.59 ) then
          stop(2,i) = stop(2,i) + 1
          stop(3,i) = 0
          stop(4,i) = 0 
        endif
        if(i.eq.1 ) then
          print *,'adj stop ',(stop(j,i),j=1,5)
        endif
      enddo  
      

c       Now check for overlaps 
      do i=2,n
         if( itimdif(stop(1,i-1),start(1,i)).gt.0.d0 ) then
c          write a warning to the screen for a manual fix-up
           write(*,'(2a,a1,i3,5i6)') '**WARNING: '
     .       ,'Time overlaps for SVN ',gnss,svn(i),(start(j,i),j=1,4)
         endif
      enddo
      return
      end
                           
 

