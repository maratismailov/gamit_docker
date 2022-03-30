Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1996. All rights reserved.

      Program BCCHECK

C    Written by R King November 1996 from BCTOT

C    Read a RINEX navigation file and evaluate
c    the consistency of the elements at the midpoint.

      implicit none

      include '../includes/dimpar.h'   
      include '../includes/orbits.h'

      logical found
                 
      integer maxhed
c     put this in ../includes/orbits.h eventually
      parameter(maxhed=20)                       


      integer*4 iye,idne,ihe,imine,iy,idn,ih,imin
     2        , iwkne,iwkn,iewkn
     3        , iesat,nsat,nprn,iflag
     4        , ibrd,nbrd,ibrd0,isat,icall
     5        , iprnt,iubc,iubcout,namelen,nblen
     6        , igood
     7        , ibcerr,i,j
     8        , ioerr,iarg,iarg2,iclarg
     9        , iyear,imonth,iday,ihr,imn,isec,ihnsec

      real*8 tsow,trans_sow,xsow,clock,ephem,bclock,bephem
     1     , utcoff,sece,sec,sowe,bsubfr1,subfr1
     2     , acheck,offset_hrs,x,x0,rtol,dr
                             
      character*1  asterisk,gnss
      character*2  hr,buf2
      character*3  doy  
      character*4  yr,buf4
      character*16 bcfile
      character*26 bcout
      character*120 version,head(maxhed)
      character*256 message

      dimension trans_sow(maxbrd,maxsat)
     .        , bephem(16,maxbrd,maxsat), ephem(16)
     .        , bclock(6,maxbrd,maxsat ), clock(6)
     .        , iesat(maxsat)
     .        , iewkn(maxbrd,maxsat)
     .        , nbrd(maxsat)
     .        , bsubfr1(8,maxbrd,maxsat),subfr1(8)
     .        , igood(maxsat,maxbrd)    
     .        , x(6),x0(6)

      data rtol/1000.d0/

c        Get version number
      call oversn(version)

      write(6,'(a)')' '
      write(message,5) version
    5 format('Started BCCHECK ',a80)
      call report_stat('STATUS','BCCHECK','orbits/bccheck',' '
     .                 ,message,0)
c        
c       Initialize the quality flag array
      do i = 1,maxsat
       do j = 1,maxbrd
         igood(i,j) = 0
       enddo
      enddo
      offset_hrs = 0.d0
                                

c** TEMPORARY: Force GPS
      gnss = 'G'

c--------Input ------------------------------------------------------------
       
c        Get the navigation file name and the reference epoch from the
c        command line

      iarg = iclarg(1,bcfile)   
      if ( iarg.le.0 ) then
        call report_stat('FATAL','BCCHECK','orbits/bccheck',' '
     .     ,'Navigation file name required on command line',0)
      else  
         iarg = iclarg(2,yr)  
         iarg2 = iclarg(3,doy)  
         if( iarg.le..0 .or. iarg2.le.0 ) then
            call report_stat('FATAL','BCCHECK','orbits/bccheck',' ' 
     .       ,'Need 3 arguments:  bccheck <file> <yr> <doy>',0) 
         else
             read(yr,'(i4)') iye 
             call check_y2k(iye)
             read(doy,'(i3)') idne 
         endif
         iarg = iclarg(4,hr)
         if( iarg.le.0 ) then
            ihe = 12
         else
            read(hr,'(i2)') ihe
         endif
         imine = 0
         sece = 0.d0 
      endif         

C        Open the files
           
      iprnt = 30 
      open( unit=iprnt,file='bccheck.out',form='formatted'
     .      ,status='unknown',iostat=ioerr)
      if (ioerr .ne. 0 ) 
     .    call report_stat('FATAL','BCCHECK','orbits/bccheck'
     .            ,'bccheck.out','Error opening print file: ',ioerr)
      write(iprnt,31) bcfile
31    FORMAT(/,' Input Broadcast (E-) File :   ',A16)  
      iubc = 31
      open(unit=iubc,file=bcfile,status='old',form='formatted'
     .    ,iostat=ioerr)
      if ( ioerr .ne. 0 ) then
          call report_stat('FATAL','BCCHECK','orbits/bccheck',bcfile,
     .    'Error opening input E-file: ',ioerr)
      endif  
      iubcout = 32 
      namelen = nblen(bcfile)          
      bcout = bcfile(1:namelen)//'.bcchecked'
      open(unit=iubcout,file=bcout,status='unknown',iostat=ioerr)
      if (ioerr .ne. 0 ) 
     .    call report_stat('FATAL','BCCHECK','orbits/bccheck',bcout,
     .    'Error opening print file: ',ioerr)



C        Convert the times to GPS week, seconds-of-week for BC evaluation

c     convert Julian day of ephemeris epoch to year and day number
c     convert ICS GPST yr/day/hms to GPS week number and seconds
      call timcon( -4,iwkne,sowe,iye,idne,ihe,imine,sece,utcoff )
      write(iprnt,'(/,a,i4,1x,i3,2x,2i3,f4.0,2x,i4,1x,f10.2)')
     .        ' Reference epoch for testing broacast elements: '
     .       , iye,idne,ihe,imine,sece,iwkne,sowe


c-------- Read all ephemeris values into storage-------------------------
            
      nsat = 0
      do i=1,maxsat
        iesat(i) = 0
        nbrd(i)  = 0
      enddo

c      read one record of the navigation file
       icall = 0
   20  call reade( iubc,icall,gnss
     .           ( iflag,xsow,gnss,nprn,iwkn,ephem,clock,subfr1 )
c      quit on eof
       if( iflag.eq.-1 ) goto 50 
c      set icall to indicate that the header has been read
       icall = 1
               
c      crude check for bad semi-major axis
       acheck = ephem(4)*ephem(4)/1.d3
       if ( nprn.eq.0 .or. nprn.gt.37 .or.
     .     dabs(acheck) .lt. 1.d4 .or. dabs(acheck) .gt. 1.d6) then
           write (message,'(a,i3,a,f7.0,a,d12.4)')
     .       'Bad BC ephemeris record   prn=',nprn,' sow='
     .       ,ephem(1),' a=',acheck
           call report_stat('WARNING','BCCHECK','orbits/bccheck'
     .                     ,' ',message,0)
       else

c        if the PRN is not in the existing list, add it
         found = .false.
         do  isat= 1, maxsat
           if( nprn.eq.iesat(isat) ) found = .true.
         enddo
         if ( .not.found ) then
           nsat = nsat + 1
           if( nsat.gt.maxsat) then
              write(message,30) nsat,maxsat
30            format('Number of satellites NSAT: ',i4,
     .        ' Exceeds program dimensions, MAXSAT:',i4)
              call report_stat('FATAL','BCCHECK','orbits/bccheck'
     .                         ,' ',message,0)
           endif
           iesat(nsat) = nprn
         endif

c        store the ephemeris message by satellite;
         do  isat = 1, nsat
           if( nprn.eq.iesat(isat) ) then
             nbrd(isat) = nbrd(isat) + 1
             if( nbrd(isat).gt.maxbrd ) then
               write(message,40)nbrd(isat),iesat(isat),maxbrd
40             format('Number of broadcast values (',i3,
     .         ') for satellite PRN ',i3,' exceeds MAXBRD (',i4,')')
              call report_stat('FATAL','BCCHECK','orbits/bccheck',' '
     .                        ,message,0)
              endif                          
              trans_sow(nbrd(isat),isat) = xsow
              iewkn(nbrd(isat),isat) = iwkn
              do  i=1,6
                 bclock(i,nbrd(isat),isat) = clock(i)
              enddo
              do  i=1,16
                 bephem(i,nbrd(isat),isat) = ephem(i)
              enddo                      
              do  i= 1,8
                 bsubfr1(i,nbrd(isat),isat) = subfr1(i)
              enddo
           endif
         enddo

c     end if on valid navigation record; go read another
      endif
      goto 20

   50 continue  
      if( nsat.le.0 ) then
        call report_stat('FATAL','BCCHECK','orbits/bccheck',' '
     .      ,'No satellites found ',0)
      else
         call report_stat('STATUS','BCCHECK','orbits/bccheck',bcfile
     .,'Successfully read broadcast navigation file: ',0)
      endif

c------ Display the satellites selected and their reference times at the IC epoch-
                  
       write(iprnt,'(/,3a)')
     .    ' PRN   Message reference epoch at test epoch '
     .   ,'   Offset from test epoch (hrs) '  
       write(iprnt,'(/,a,f5.0,a)') 
     .     ' Reject a message if its position is greater than '
     .     ,rtol,' m off at the reference epoch'

      do isat=1,nsat
           
         if( nbrd(isat).gt.0 )  then
           do ibrd = 1,nbrd(isat)

           call timcon( 4,iewkn(ibrd,isat),bephem(1,ibrd,isat)
     .                , iy,idn,ih,imin,sec,utcoff ) 
           offset_hrs= (iy-iye)*8760 + (idn-idne)*24 
     .                + (ih-ihe) + (imin-imine)/60. + (sec-sece)/3600.d0
c            print *,'debug isat,ibrd,iewkn,bephem offset_hrs '
c     .            ,isat,ibrd,iewkn(ibrd,isat),bephem(1,ibrd,isat)
c     .            , offset_hrs
c            write(iscrn
c     .           ,'(1x,i3,4x,i4,2x,f8.0,3x,2i4,2x,2i3,f6.2,6x,f8.2 )')
c     .             iesat(isat),iewkn(ibrd,isat),bephem(1,ibrd,isat)
c     .           , iy,idn,ih,imin,sec
c     .           , offset_hrs
c            write(iprnt
c     .           ,'(1x,i3,4x,i4,2x,f8.0,3x,2i4,2x,2i3,f6.2,6x,f8.2)')
c     .             iesat(isat),iewkn(ibrd,isat),bephem(1,ibrd,isat)
c     .           , iy,idn,ih,imin,sec
c     .           , offset_hrs
            enddo
          endif
      enddo
               

C------- For each set of elements, evaluate x y z at the reference time----------
          
c     Loop over satellites   
      
      do isat = 1,nsat   
        
        
c       Choose the middle epoch as reference for the quality check
        if( nbrd(isat).gt.1 ) then
           ibrd0 = nbrd(isat)/2
   60      call brdxyz ( iwkne,sowe,iewkn(ibrd0,isat)
     .                 ,bephem(1,ibrd0,isat),x0,ibcerr )  
c          crude check for bad semi-major axis 
           acheck = bephem(4,ibrd0,isat)*bephem(4,ibrd0,isat)/1.d3
           if( dabs(acheck) .lt. 1.d4 .or. dabs(acheck) .gt. 1.d6) then 
             write(message,'(a,i3,a)') 
     .        ' **Bad elements at midpoint (ibrd=',ibrd0
     .        ,' choosing another reference' 
             call report_stat('WARNING','BCCHECK','orbits/bccheck',' '       
     .                       ,message,0)
             igood(isat,ibrd0) = -1
             ibrd0 = ibrd0 + 1
             if(ibrd0.gt.nbrd(isat)) then
               call report_stat('FATAL','BCCHECK','orbits/bccheck',' '
     .                      ,'Cannot find a good reference epoch',0)
              endif
             goto 60  
           else          
                continue
           endif 
        else if ( nbrd(isat).gt.1 ) then
           ibrd0 = 1
        else
c          ibrd0 = 0
        endif
            
          
                        
c       Loop over the messages and set the quality flag

c          igood = -1    failed check
c                =  0    initial value (missing or not checked)
c                =  1    passed check
                 
        if( nbrd(isat).gt.0 ) then
          do ibrd = 1, nbrd(isat)  

             call brdxyz ( iwkne,sowe,iewkn(ibrd,isat)
     .                   ,bephem(1,ibrd,isat),x,ibcerr )  
                
c             print *,'isat prn ibrd t xyz ',isat,iesat(isat),ibrd
c     .          , iewkn(ibrd,isat),bephem(1,ibrd,isat),(x(i),i=1,6) 

             call timcon( 4,iewkn(ibrd,isat),bephem(1,ibrd,isat)
     .                  , iy,idn,ih,imin,sec,utcoff )  
             offset_hrs= (iy-iye)*8760 + (idn-idne)*24 
     .                + (ih-ihe) + (imin-imine)/60. + (sec-sece)/3600.d0


c           crude check for bad semi-major axis
            acheck = bephem(4,ibrd,isat)*bephem(4,ibrd,isat)/1.d3
            if ( iesat(isat).eq.0 .or. iesat(isat).gt.37 .or.
     .         dabs(acheck) .lt. 1.d4 .or. dabs(acheck) .gt. 1.d6) then
               write (message,'(a,i3,a,f7.0,a,d12.4)')
     .           'Bad prn or semi-major axis:  prn=',iesat(isat)
     .        ,' sow=',bephem(1,ibrd,isat),' a=',acheck
            call report_stat('WARNING','BCCHECK','orbits/bccheck',' '
     .                      ,message,0)
              igood(isat,ibrd) = -1
           endif  
           if( igood(isat,ibrd).ne.-1 ) then
               dr = dsqrt((x(1)-x0(1))**2 + (x(2)-x0(2))**2 
     .                 + (x(3)-x0(3))**2)*1.d3
              if( dr.le.rtol ) then
                  igood(isat,ibrd) = 1 
                  asterisk = ' '   
               else 
                  write(message,'(a,i3,a,i5,i4,2i3,a,f6.1,f8.0,a)')   
     .           'Elements for PRN ',iesat(isat),' exceed tolerance at '
     .             ,iy,idn,ih,imin,' offset hrs= ',offset_hrs,dr,' m'  
                  call report_stat('WARNING','BCCHECK','orbits/bccheck'
     .                             ,' ',message,0)
                  igood(isat,ibrd) = -1
                  asterisk = '*'
               endif

               write(iprnt,
     . '(1x,i3,4x,i4,f10.0,3x,2i4,2x,2i3,f6.2,4x,f8.2,1x,f10.1,1x,a1)')
     .             iesat(isat),iewkn(ibrd,isat),bephem(1,ibrd,isat)
     .           , iy,idn,ih,imin,sec
     .           , offset_hrs,dr,asterisk     
             endif
          enddo
        endif

      enddo   


c-------- Write out a navigation file with only good elements 

c     write the header
              
      nhead = nhead + 1
      if( nhead.gt.maxhed )
     .   call report_stat('FATAL','BCCHECK','orbits/bccheck',' '
     .        ,'Too many header records on output',0)  
c     shift 'END OF HEADER' record down one to make room for new comment
      head(nhead) = head(nhead-1)
      head(nhead-1)(1:31) = 'File filtered by GAMIT BCCHECK ' 
      call getdat(iyear,imonth,iday)
      call gettim(ihr,imn,isec,ihnsec)   
      write (buf4,'(i4)') iyear  
      read  (buf4,'(a4)') head(nhead-1)(34:37)
      head(nhead-1)(38:38)='-'
      write (buf2,'(i2)') imonth   
      read  (buf2,'(a2)') head(nhead-1)(39:40)
      head(nhead-1)(41:41)='-'
      write (buf2,'(i2)') iday
      read  (buf2,'(a2)') head(nhead-1)(42:43)
      write (buf2,'(i2)') ihr 
      read  (buf2,'(a2)') head(nhead-1)(45:46)
      head(nhead-1)(47:47)=':'
      write (buf2,'(i2)') imn 
      read  (buf2,'(a2)') head(nhead-1)(48:49)
      head(nhead-1)(50:50)=':'
      write (buf2,'(i2)') isec 
      read  (buf2,'(a2)') head(nhead-1)(51:52)
      head(nhead-1)(53:60) = '        '
      head(nhead-1)(61:80) = 'COMMENT             '       
      icall = 0  
c     all arguments dummy except for first four
      call writee( iubcout,icall,nhead,head,tsow
     .           , nprn,iewkn,ephem,clock,subfr1 )

c     write the valid records
    
      do isat = 1,nsat
         
        if( nbrd(isat).gt.0 ) then
          do ibrd = 1, nbrd(isat)

            if( igood(isat,ibrd).eq.1 ) then

              tsow = trans_sow(ibrd,isat)
              nprn = iesat(isat)
              iwkn = iewkn(ibrd,isat)
              do i = 1,16
                ephem(i) = bephem(i,ibrd,isat)
              enddo
              do i = 1,6
                clock(i) = bclock(i,ibrd,isat)
              enddo
              do i = 1,8
                subfr1(i) = bsubfr1(i,ibrd,isat)
              enddo  
              icall = 1
              call writee( iubcout,icall,nhead,head,tsow
     .                   , nprn,iwkn,ephem,clock,subfr1 )

            endif

          enddo
        
        endif

      enddo 
                                                     
      call report_stat('STATUS','BCCHECK','orbits/bccheck',bcout
     .                ,'Successfully wrote filtered navigation file',0)
   
      stop
      end
