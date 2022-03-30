Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 2015. All rights reserved.

      Subroutine GETSP3( debug,icall,iwkn,sow,gnss_sel,iprn
     .                 , xsat,svclock,satok )
C
C Written by R.W. King December 2015 based on orbits/ngstot
C
C Read an NGS Standard Product emphemeris file at the requested epoch
c and return the Earth-fixed coordinates.  This routine requires an
c exact alignment of the observation epoch and the SP3 entry, which
c will be the case if we invoke it for even 15-minute intervals of GPST
       
c  Input
c      icall = 0   : Read all of the records for the selected GNSS into storage
c      icall > 0   : Return the coordinates at nearest requested time
c      iwkn,sow    : GPS week, sec-of-week requested 
c      gnss_sel    : GNSS system requested (G R C E J)
c      iprn        : PRN requested  
c  Output
c      xsat(3)     : Earth-fixed Cartesian coordinates of the SV
c      svclock     : SV clock offset 
c      satok       : T/F if SV found on sp3 file

      implicit none

      include '../includes/dimpar.h'  
      include '../includes/makex.h'

      integer*4 icall,iprn,id,iyr,imon,iday,idn,ihr,imin,jds,mjds
     .        , iwkn,iwkne(maxepc),numsp3sv,isp3sv(maxsp3sv)
     .        , numsat,itsat(maxsat),msat,nepochs,iepoch,iy
     .        , ioerr,jsat,isatcnt,notl,i,j,k

      real*8 pos(3),accsat(maxsat),pvsigb,clksigb,svclock
     .     , docmc(3),sec,sow,sowe(maxepc),ts,utcoff
     .     , xsat(3),xsatdot(3)
     .     , satics(maxorb,maxsat),frcts,delts,tdiff
     .     , clocks(maxsp3sv, maxepc)
     .     , y(3,maxsp3sv, maxepc)
* MOD TAH 180309: Fixed order of indices of clock and y to 
*                 be consistent with assignment
C    .     , clocks(maxepc,maxsp3sv)
C    .     , y(3,maxepc,maxsp3sv),

      character*1 spver,gnss_sel,asvid
      character*2 linesym
      character*8 otlmod
      character*80 line 
      character*256 message
                     
      logical debug,satok

c     Functions
      integer*4 idoy,julday
      real*8 secdif

      save y,clocks,iwkne,sowe,nepochs

c* First call to read all values into storage
  
      if( icall.eq.0 ) then
    
c Determine the file format using the first two characters of the first line
* MOD TAH 190703: Added support for SP3-d
c   SP1    ' # '
c   SP3-a  '# '  or '# a'
c   SP3-b  '#b'
c   SP3-c  '#c'
c   SP3-d  '#d'
            
      read(usp3,'(a2)',iostat=ioerr) linesym 
      if( ioerr.eq.0 ) then 
        if ( linesym.eq." #") then  
          spver = '1'
          call report_stat('FATAL','MAKEX','getsp3',' '
     .                  , 'SP1 format not supported',0)
        elseif ( linesym.eq.'# '.or.linesym.eq.'#a' ) then
          spver = 'a'
        elseif ( linesym.eq.'#b' ) then
          spver = 'b'
        elseif ( linesym.eq.'#c' ) then
          spver = 'c'
* MOD TAH 190703: Check for sp3 d version
        elseif ( linesym.eq.'#d' ) then
          spver = 'd'
        else
          call report_stat('FATAL','MAKEX','getsp3',' '
     .                  ,'Unrecognized SP orbit format',0 )
        endif
      else
        call report_stat('FATAL','MAKEX','getsp3',' '
     .                  , 'Error reading SP3 file type',ioerr)
      endif 
      rewind(usp3)

c Read the SP file header to get times and satellites
c MOD TAH: rsp3hd seems to handle sp3d even though version is not 
c     passed to routine?
      call rsp3hd( usp3,gnss_sel
     .             , iyr,imon,iday,ihr,imin,sec
     .             , delts,mjds,frcts,nepochs
     .             , numsp3sv,numsat,itsat,accsat
     .             , pvsigb,clksigb,otlmod )
      if( debug ) then 
        print *,'GETSP3 numsat gnss_sel,clksigb otlmod '
     .              , numsat,gnss_sel,clksigb,otlmod
        print *,' itsat',(itsat(i),i=1,numsat)
        print *,'       accsat  ',(accsat(i),i=1,numsat)
        write(*,'(a,i3,100i3)') 
     .   'In GETSP3 numsat itsat ', numsat,(itsat(i),i=1,numsat)
      endif

c If ocean-tidal loading correction used, read in the coefficients
c  (saved in /otlcmc; docmc not calculated during this call)

      if( otlmod(1:1).ne.' ') then                         
c      hard-wire # components
       notl = 11                             
       call otlcmc( jds,ts,otlmod,notl,1,docmc )
       write(message,'(a,a8)') 
     .   'Converting SP3 CE to CM using otlcmc.dat offsets for ',otlmod 
       call report_stat('STATUS','MAKEX','getsp3',' ',message,0)
      endif
                    
c Initialize the array to provide a check for a valid entry later

      do k=1,maxsat
        do j=1,maxepc
          do i=1,3
C            y(i,j,k) = 0.d0
* MOD TAH 180309: Fixed ordering to be consistent with deminensioning.
            y(i,k,j) = 0.d0
          enddo
         enddo
       enddo

c Read the values into storage

      do iepoch =1,nepochs

        read(usp3,'(3x,i4,4(1x,i2),1x,f10.7)',iostat=ioerr) 
     .           iyr,imon,iday,ihr,imin,sec
        if( ioerr.ne.0 ) then 
          write(message,'(a,f5.2,a)') 'Error reading epoch ',iepoch
     .        ,' for SP ',spver,'file'  
          call report_stat('FATAL','MAKEX','getsp3',' '
     .           ,message,ioerr)
        endif       
c       convert to GPS week,sec-of-week and store                 
        idn = idoy(iyr,imon,iday)
        call timcon( -4,iwkne(iepoch),sowe(iepoch),iyr,idn,ihr,imin,sec
     .             , utcoff )
        isatcnt = 0          
        do j=1,numsp3sv
          read(usp3,'(a)',iostat=ioerr) line   
c         assume that sp3-c has the full 80 columns, though some can be blank
          read(line,'(1x,a1,i2,4f14.6,3i3,i4)'
     .          ,iostat=ioerr) asvid,id,(pos(i),i=1,3),svclock
          if( asvid.eq.' ' ) asvid = 'G'
cd        if(debug) print *,'read asvid id x ',asvid,id,(pos(i),i=1,3)
          if( ioerr.ne.0 ) then
            write(message,'(a,i5,a,a1,a)') 
     .         'Error reading coordinates at epoch ',iepoch
     .            ,' for SP ',spver,' file'
            call report_stat('FATAL','MAKEX','getsp3',' ',message,0)
          endif
          if( asvid.eq.gnss_sel ) then
            isatcnt = isatcnt + 1 
            if(id.ne.itsat(isatcnt)) then
              write(message,'(a,i2,a,i2,a,i5)') 'Unexpected PRN (',id
     .          ,') at slot ',isatcnt,' for epoch ',iepoch
              call report_stat('FATAL','MAKEX','getsp3',' ',message,0)
            endif
            do i=1,3
              y(i,isatcnt,iepoch) = pos(i) 
            enddo          
            clocks(isatcnt,iepoch) = svclock*1.d-6 
          endif 
        enddo  
        if( isatcnt.ne.numsat ) then
          write(message,'(a,i4,a,i3,a,i3,a)') 
     .       'isatcn at epoch ',iepoch,'= ',isatcnt
     .      ,' not equal numsat from header (',numsat,')' 
          call report_stat('FATAL','MAKEX','getsp3',' ',message,0)
        endif
c       apply the correction from CE to CM for ocean loading 
        if( otlmod(1:1).ne.' ' ) then   
          jds= julday( imon,iday,iyr )
          ts= dble(ihr)*3600.d0 + dble(imin)*60.d0 + sec
c         hard-wire the # of components
          notl = 11
          call otlcmc( jds,ts,otlmod,notl,2,docmc ) 
cd           write(6,'(a,i8,f8.1,3f7.4)') 'CMC db ',jds,t
cd            ,(docmc(i),i=1,3)
          do j=1,numsat
            do i=1,3
              y(i,j,iepoch) = y(i,j,iepoch)  + docmc(i)/1.d3
            enddo                                          
          enddo
        endif
            
c     end loop on epochs and initial read (icall=0)
      enddo
      return  

c     Come here with icall > 1 to find the epoch closest to the one requested
      else 

      satok = .false.
      do j=1,numsat
        if(iprn.eq.itsat(j)) jsat=j 
      enddo          
      do iepoch=1,nepochs
        if( debug.and.iepoch.eq.1.and.iprn.eq.itsat(jsat)) then
          write(*,'(a,a1,i2)')'GETSP3: first SP3 value for PRN '
     .             ,gnss_sel,iprn
          call wtime (6,iwkne(iepoch),sowe(iepoch),'GPST'
     .                    ,'  Epoch read:   ')
        endif   
        tdiff = secdif( iwkn,sow,iwkne(iepoch),sowe(iepoch) ) 
cd        print *,'iepoch iwkn sow iwkne sowe tdiff '
cd     .         , iepoch,iwkn,sow,iwkne(iepoch),sowe(iepoch),tdiff
c       allow 1s slack in testing for match
        if( dabs(tdiff).le.1.d0 ) then
           do i=1,3
              xsat(i) = y(i,jsat,iepoch)*1.d3
c             no velocities available currently
              xsatdot(i) = 0.d0
           enddo                                
           svclock = clocks(jsat,iepoch)
           if (debug) then                
             write(*,'(a,a1,i2)')'GETSP3: epoch match for PRN '
     .                ,gnss_sel,iprn
             call wtime (6,iwkne(iepoch),sowe(iepoch),'GPST'
     .                    ,'  Epoch read:   ')
           endif    
c          make sure there was a valid position  
           if( xsat(1).ne.0.d0 ) then
             satok = .true.
           else    
             call timcon( 4,iwkne(iepoch),sowe(iepoch)
     .                  , iyr,idn,ihr,imin,sec,utcoff )
             write(message,'(a,a1,i2,a,3i3,f5.1)') 
     .          'Bad SP3 position for PRN ',gnss_sel,iprn,' at GPST '
     .           ,idn,ihr,imin,sec 
             call report_stat('WARNING','MAKEX','getsp3',' ',message,0)
           endif
c          exit when found
           return
        endif
      enddo
 
c     endif for ICALL
      endif
      return
      end

