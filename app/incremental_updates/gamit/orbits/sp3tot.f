Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1994. All rights reserved.
      Program SP3TOT
         
c  Read an IGS Standard Product (SP-3) earth-fixed orbit ascii orbit tile and create 
c  a binary inertial T-file for GAMIT

c  R. King October 2009 rewrite of program ngstot (R King & Y Bock  March 1988) 
c  This version stores coordinates in memory, does not create an intermediate 
c  earth-fixed t-file, pre-scans the values to detect missing satellites or 
c  epochs, cannot be run interactively, and cannot read SP-1.

      implicit none

      include '../includes/dimpar.h'  

      integer*4 maxepc
      parameter (maxepc=300)   
                
      integer*4 jds,jdf,jde,mjds,iyear,imon,iday,ihr,imin,nblen 
     .        , jd(maxepc),sigx(4)
     .        , iusp3,iut,inut,iut1,ipole,iscrn
     .        , issat,itsat,nsvs,nsvt,nepchs,iepchs,nepcht
     .        , nintrs,idir,icall,nics,nepctic
     .        , iy,iy1,idoy,idoy1,iweek,idow,julday,iarg,iclarg
     .        , iut1pol,notl,id,ioerr,i,j,k  
C     iut1_pol - Bit mapping setting of short period EOP terms.  

c     function
      integer*4 mchkey

      real*8 sec, delts, frcts, dt, te, tf, ts, a, t(maxepc)
      real*8 satics(maxorb,maxsat),accsat(maxsat),pvsigb,clksigb
     .     , docmc(3),x(3,maxsat,maxepc),clock
      real*8  stveci(6),stveco(6),rot(3,3),rotdot(3,3),sidtm,tdtgpst

      character*1 lowerc,upperc,spver,eflg,cflg,mflg,pflg,abatch
      character*2 linesym
      character*4 icstfl(15),icsnam(maxorb)
      character*5 precmod,nutmod,gravmod,frame,srpmod,eradmod,antradmod
      character*6 prog
      character*8 otlmod
      character*16 spfile,tfile,gfile,xfile,tfilef,tmpnam
      character*80 nut_header,line
      character*120 version
      character*256 message

      logical batch,sat_ok(maxsat),eof,truncate

      dimension itsat(maxsat),issat(maxsat)

      data icstfl/'X   ','Y   ','Z   ','XDOT','YDOT','ZDOT'
     .           ,'RAD1','RAD2','RAD3','    ','    ','    '
     .           ,'    ','    ','    '/

      data nutmod/'     '/,gravmod/'EGM96'/,srpmod/'ECOM1'/

c Initialization
c     Screen
      iscrn = 6
      iusp3 = 1
c     T-file
      iut = 2  
      spfile = " "
      tfile  = " "  
      inut = 3 
      iut1 = 4
      ipole = 5 
      do j=1,maxsat   
        do i=1,maxorb
          satics(i,j)=0.d0 
        enddo
      enddo
      do k=1,maxepc   
        jd(k) = 0
        t(k) = 0.    
        do j=1,maxorb
          do i=1,3
            x(i,j,k) = 0.d0 
          enddo
        enddo 
      enddo      
      batch = .true.
      abatch = 'b' 
                     

c Get the version number

      call oversn(version)
c*      write(iscrn,'(a)')' '
      write(message,5) version
    5 format('Started SP3TOT ',a120)
      call report_stat('STATUS','SP3TOT','orbits/sp3tot',' ',message,0)
                        

c Read the command-line input

c  if there are command-line arguments, use them and skip the interactive questions
c
c     sp3tot  [sp3file]  [t-file]  [-trun]
c            (required) (required  (optional)

      iarg = iclarg(1,spfile)    
      if( iarg.eq.0 ) then
         write(*,'(a)')
         write(*,'(a)') 'sp3tot [sp3file] [t-file] [-trun]'
         write(*,'(a)') 
         write(*,'(a)') '  First two arguments required'
         write(*,'(a)') ' -trun truncates file at first bad entry'
         write(*,'(a)') '  default is to remove the SV entirely'
         stop
      else
c       file name must be standard format, without full path for logic in SP3TOT
        if( nblen(spfile).gt.16 ) call report_stat('FATAL','SP3TOT'
     .       , 'orbits/sp3tot',' '
     .       , 'SP3 file name too long --cannot use full path',0)
        iarg = iclarg(2,tfile)                                  
        if( iarg.le.0 )  
     .            call report_stat('FATAL','SP3TOT','orbits/sp3tot',' '
     .               ,'Missing command-line argument for T-file name',0)
        truncate = .false.
        iarg = iclarg(3,tmpnam)
        if( iarg.gt.0 ) then    
          if( tmpnam(1:2).eq.'-t' ) then
            truncate = .true.
          else            
            call report_stat('FATAL','SP3TOT','orbits/sp3tot',' '
     .         ,'Third argument must be -trun',0)
          endif
        endif
      endif   
                  
          
c Open the files and set the frame parameters for the t-file

      open(unit=iusp3,file=spfile,status='old',form='formatted',
     .    iostat=ioerr)
      if (ioerr .ne. 0 ) then
          call report_stat('FATAL','SP3TOT','orbits/sp3tot',spfile
     .                    ,'Error opening SP3 file: ',ioerr)
      endif
      open(unit=iut,file=tfile,status='unknown',form='unformatted'
     .     ,iostat=ioerr)
      if (ioerr .ne. 0) then
        call report_stat('FATAL','SP3TOT','orbits/sp3tot',tfile
     .                  ,'Error opening t-file',0)
          stop
      endif 
c     open the nutation table if needed (no 'nbody' for Sun and Moon)
      if( .not.fcheck('nbody') then                
        open(unit=inut,file='nutabl.',status='old',form='formatted'
     .       ,iostat=ioerr )
        if (ioerr .ne. 0) then
          call report_stat('FATAL','SP3TOT','orbits/sp3tot','nutabl.'
     .                  ,'Error opening nutation file',0)
        endif  
        read(inut,'(a)',iostat=ioerr) nut_header
        if( ioerr.ne.0) call report_stat('FATAL','SP3TOT','orbits/sp3tot'
     .   ,' ','Error reading first line of nutation file to get model'
     .   ,ioerr)
        rewind(inut)
        if(mchkey(nut_header,'IAU20',80,5).gt.0) then
          nutmod='IAU00'
        else
          nutmod='IAU80'
        endif  
      else
c SCM changed default to new sofa IAU00A model 190801
c        nutmod = 'IAU00'
        nutmod = 'IAU0A'
      endif
      frame = 'J2000'
c      precmod = 'IAU76'      
      precmod = 'IAU0A'      
      open(unit=iut1,file='ut1.',status='old',form='formatted'
     .    ,iostat=ioerr )
      if (ioerr .ne. 0) then
        call report_stat('FATAL','SP3TOT','orbits/sp3tot','ut1.'
     .                  , 'Error opening UT1 file',0)
          stop
      endif    
      open(unit=ipole,file='pole.',status='old',form='formatted'
     .    ,iostat=ioerr )
      if (ioerr .ne. 0) then
        call report_stat('FATAL','SP3TOT','orbits/sp3tot','pole.'
     .                  ,'Error opening pole-postion file',0)
          stop
      endif                      


c Determine the file format using the first two characters of the first line

c   SP1    ' # '  --- no longer supported
c   SP3-a  '# '  or '# a'
c   SP3-b  '#b'
c   SP3-c  '#c'
      read(iusp3,'(a2)') linesym
      if ( linesym.eq." #") then  
        call report_stat('FATAL','SP3TOT','orbits/sp3tot',tfile
     .                  , 'SP-1 files no longer supported',0)
        spver = '1'                
      elseif ( linesym.eq.'# '.or.linesym.eq.'#a' ) then
        spver = 'a'
      elseif ( linesym.eq.'#b' ) then
        spver = 'b'
      elseif ( linesym.eq.'#c' ) then
        spver = 'c'
      else
        call report_stat('FATAL','ORBITS','sp3tot',' '
     .                  ,'Unrecognized SP orbit format',0 )
      endif
      rewind(iusp3)
       

c Read the SP file header to get times and satellites

      call rsp3hd( batch,spver,iusp3,iyear,imon,iday,ihr,imin,sec
     .           , delts,mjds,frcts,nepchs,nsvs,issat,accsat
     .           , pvsigb, clksigb, otlmod )
c     convert SP times to PEP JD
      jds = julday(imon,iday,iyear)
      ts = dble(ihr)*3600.d0 + dble(imin)*60.d0 + sec  
      ts = anint(ts)               

c Read the values into storage

      iepchs = 0    
      eof = .false.
      do while ( .not.eof )  
        read(iusp3,'(a)',iostat=ioerr) line
        if( ioerr.eq.-1.or.line(1:3).eq.'EOF' ) then
          eof = .true.
        elseif( ioerr.ne.0 ) then
          call report_stat('FATAL','SP3TOT','orbits/sp3tot'
     .         ,' ','Error reading epoch line of SP3 file',ioerr)
        else
          read(line,'(3x,i4,4(1x,i2),1x,f10.7)',iostat=ioerr) 
     .        iyear,imon,iday,ihr,imin,sec    
          if( ioerr.ne.0 ) then
            call report_stat('FATAL','SP3TOT','orbits/sp3tot'
     .         ,' ','Error decoding  epoch line of SP3 file',ioerr)
          else       
            iepchs = iepchs + 1
            jd(iepchs)= julday( imon,iday,iyear )
            t(iepchs)= dble(ihr)*3600.d0 + dble(imin)*60.d0 + sec
            do j=1,nsvs
              read(iusp3,'(a)',iostat=ioerr) line   
              if( ioerr.ne.0 ) then
                 call report_stat('FATAL','SP3TOT','orbits/sp3tot',' '
     .               ,'Error reading coordinate line of sp3 file',ioerr)
              endif
c             assume that sp3-c has the full 80 columns, though some can be blank
              if( spver.eq.'c') then      
                read(line,'(2x,i2,4f14.6,3i3,i4,1x,2a1,2x,2a1)'
     .                       ,iostat=ioerr) 
     .              id,(x(i,j,iepchs),i=1,3),clock,(sigx(i),i=1,4)
     .            , eflg,cflg,mflg,pflg
              else                     
                read(line,'(2x,i2,3(1x,f13.6))',iostat=ioerr) 
     .                  id,(x(i,j,iepchs),i=1,3)
              endif  
              if( ioerr.ne.0 ) then                  
                 write(message,'(a)') line     
                 call report_stat('WARNING','SP3TOT','orbits/tdtrit',' '
     .               ,message,ioerr)
                 call report_stat('FATAL','SP3TOT','orbits/tdtrit',' '
     .             ,'Error decoding coordinates line in SP3 file',ioerr)
              endif
              if( id.ne.issat(j) ) then  
                 write(message,'(a,i3,a,i3,a,i4,a)') 
     .               'PRN = ',id,' but expected '
     .              ,issat(j),' at epoch ',iepchs,' of SP3 file'
                  call report_stat('FATAL','SP3TOT','orbits/sp3tot',' '
     .                          ,message,0 )
               endif
            enddo
          endif
        endif
      enddo                              
      if( iepchs.ne.nepchs ) then
        write(message,'(a,i4,a,i4,a)') '# epochs read (',iepchs
     .      ,' .ne. # expected (',nepchs,'), using # read '    
         call report_stat('WARNING','FATAL','SP3TOT','orbits/sp3tot'
     .                   ,' ',message,0)
         nepchs = iepchs
      endif
                                                       

c Scan the values to get eliminate satellites with bogus epochs (or shorten the
c span) and get the ICs at the center epoch

      call scansp3( nsvs,issat,nepchs,jd,t,x,truncate,nsvt,itsat,nepcht)
    
c Get the start, stop, and IC epoch for the t-file 

c     start epoch                               
      jds = jd(1)
      ts  = t(1)    
c     IC epoch -- middle of span
      dt = delts*(nepcht/2)
      jde= jds
      te=  ts
      call timinc( jde,te,dt ) 
c     if initial condition time is within 1 hour of midday make it exactly 12:00
      if ( dabs(te-43200.d0) .le. 3600.d0 ) then
        te = 43200.d0
      endif   
c     end epoch 
      jdf = jd(nepcht)  
      tf = t(nepcht)
c     avoid roundoff by assuming that the epochs are even seconds
      te = anint(te)
      tf = anint(tf)      
        
c Get the values at the IC epoch
            
       nepctic = 0 
       do i = (nepcht/2) - 5 , (nepcht/2) + 5
          if(  jd(i).eq.jde.and.dabs(t(i)-te).lt.1.d-6 ) nepctic = i
       enddo         
      if( nepctic.ne.0 ) then
         do j=1,nsvt
           do i=1,3
             satics(i,j) = x(i,j,nepctic)
           enddo
         enddo  
       else
         call report_stat('FATAL','SP3TOT','orbits/sp3tot',' ',
     .    'Cannot match epoch number with IC time',0)
       endif
             

c Positions only - no velocities or partials           
      nintrs = 3
      nics = 9            

c Set the non-gravitational force parameters to nominal values
      do j=1,nsvt
         satics(7,j)= 1.D0
         satics(8,j)= 0.D0
         satics(9,j)= 0.D0
      enddo
      do i=1,nics
        icsnam(i) = icstfl(i)
      enddo
                     

c If ocean-tidal loading correction used, read in the coefficients and 
c correct the ICs and values

      if( otlmod(1:1).ne.' ') then                         
c       hard-wire # components
        notl = 11    
c       1st call reads in the coefficients, most arguments dummy                            
        call otlcmc( jds,ts,otlmod,notl,1,docmc )  
c       2nd call computes the offsets at the specified time
        call otlcmc( jde,te,otlmod,notl,2,docmc )  
        do j=1,nsvt
          do i=1,3
            satics(i,j) = satics(i,j) + docmc(i)/1.d3
          enddo
        enddo   
        do k=1,nepcht
          do j=1,nsvt
            do i=1,3
              x(i,j,k) = x(i,j,k) + docmc(i)/1.d3
            enddo
          enddo
        enddo
      endif 
 

c Rotate the ICs and  orbit into the inertial frame

      idir = 1  
      tdtgpst = 32.184d0 + 19.d0 
* MOD TAH 200505: Read the sestbl. to get the correct value   
C     iut1pol = 7 
      call get_iut1pol( iut1pol )

      call rotsnp( idir,jde,te,tdtgpst,iut1pol
     .           , frame,precmod,iut1,ipole,inut
     .           , rot,rotdot,sidtm,xpole,ypole ) 
      do j=1,nsvt
        do i=1,6
          stveci(i) = satics(i,j)
        enddo                         
        call rotcrd( rot,rotdot,stveci,stveco,0,1 )
        do i=1,6
           satics(i,j) = stveco(i)
        enddo  
      enddo  
* MOD TAH 200505" Removed line below.  Not needed.           
C     iut1pol = 7    
      do k=1,nepcht 
        call rotsnp( idir,jd(k),t(k),tdtgpst,iut1pol,rot,rotdot,sidtm
     .           , iut1,ipole,inut,frame,precmod )
        do j=1,nsvt
          do i=1,3
            stveci(i) = x(i,j,k)
          enddo                         
          call rotcrd( rot,rotdot,stveci,stveco,0,1 )
          do i=1,3
             x(i,j,k) = stveco(i)
          enddo    
        enddo  
      enddo

       
c** Note that the velocity at the initial epoch has not yet been
c   dynamically calculated (it's 0 in the E-fixed frame and just 
c   has the precession, nutation, and sidereal rotation in the 
c   inertial frame).  A proper velocity is calculated for the g-file
c   gmake, below.

c Write the t-file header
                                     
      call thdrit( iut,jde,te,jds,ts,jdf,tf,delts,nepcht,nintrs
     .           , nsvt,itsat,satics,nics,spfile,icsnam
     .           , precmod,nutmod,gravmod,frame,srpmod
     .           , eradmod,antradmod )
       
         
c Write the t-file data records
   
      do k=1,nepcht
         write(iut) ((x(i,j,k),i=1,3),j=1,nsvt)  
      enddo
                                    

c Now compute dynamic ICs and write out a G-file

      abatch = 'b'
      icall = 1
      prog = 'SP3TOT'
      close(iut) 
      call gmake( prog,tfile,gfile,itsat )

      call report_stat('STATUS','SP3TOT','orbits/sp3tot',' ',
     .               'Normal end to SP3TOT ',0)

      stop
      end
