      Program eqfupd

c     Create or update an eq_file from SINEX or station.info
c     R. King 070112

c*** This version quick and dirty to simply translate an IGS SINEX.
c     Next version will so comprehensive merges with checks and take
c     SINEX and station.info.


      implicit none
                     

      logical found,eod

      integer*4 isnx,ieq,yrb,yre,doyb,doye,sodb,sode
     .        , monb,mone,dayb,daye,hrb,hre,minb,mine,secb,sece
     .        , ioerr
               
      character*1 achr,pvflg      
      character*16 snxf,eqf  
      character*4 site4,ext
      character*2 pt   ! SINEX PT code (normally A)
      integer*4 occ    ! Sinex OCC code (increments by 1 for each offset)
      character*8 site1,site2,sitelast
      character*80 line  

      integer*4 rcpar, lenrun

      real*8 dsod,dsec

      data isnx/1/,ieq/2/

***   Decode runstring
      lenrun = rcpar(1,snxf)
      if( lenrun.le.0 ) then
          write(*,100)
 100      format('EQFUPD: Earthquake file update based on SINEX ',
     .           'discontinuity file',/,
     .           'Runstring: ',/,
     .           'eqfupd <sinex file> [output eqfile]',/,
     .           'Default [output eqfile] is eq_file.out',/)
          stop 'eqfupd: Incomplete runstring'
      endif

C      write(*,*) 'Enter SINEX file name'
C      read(*,*) snxf  
      open(unit=isnx,file=snxf,status='old',iostat=ioerr) 
      if( ioerr.ne.0 ) then
        write(*,'(a,a,i5)') 'Error opening input file ',snxf,ioerr   
        stop
      else  
        write(*,'(a,a)') ' Opened input file ',snxf
      endif
                    
      lenrun = rcpar(2,eqf)
      if( lenrun.le.0 ) then
         eqf = 'eq_file.out'
      end if

      open(unit=ieq,file=eqf,status='unknown',iostat=ioerr)
      if( ioerr.ne.0 ) then
        write(*,'(a,a,i5)') 'Error opening ouput rename  file '
     .     ,eqf,ioerr
        stop
      else  
        write(*,'(a,a)') ' Opened output rename file ',eqf  
      endif  

      
c  Position the SINEX file at the discontinuities

      found = .false.
      do while( .not.found )
         read(isnx,'(a)',iostat=ioerr) line
         if( ioerr.eq.-1 ) then  
           print *,'EOF on SINEX file '
           stop
         elseif ( ioerr.ne.0 ) then
            write(*,*) 'Error reading SINEX file ', ioerr 
            stop
         elseif (line(11:17).eq.'DISCONT' ) then
            found = .true.
         endif
      enddo

c  Skip a line

      read(isnx,'(a)') line

c  Read and translate the discontinuities

      eod = .false.  
      site4 = '    '
      do while(.not.eod ) 
        read(isnx,'(a)',iostat=ioerr) line
        if( ioerr.eq.-1 ) then  
          print *,'EOF on SINEX file '
          stop
        elseif ( ioerr.ne.0 ) then
          write(*,*) 'Error reading SINEX file ', ioerr 
          stop   
        elseif( line(11:17).eq.'DISCONT' ) then
          eod = .true.  
        else
c          print *,'site line ',line
CABER  A    1 P 00:000:00000 00:000:00000 P - 
          read(line,120,iostat=ioerr)
     .            site4,pt, occ, yrb,doyb,sodb,yre,doye,sode,pvflg
 120      format(1x,a4,1x,a2,1x,i4,3x,i2,1x,i3,1x,i5,1x,i2,1x,i3,
     .           1x,i5,1x,a1)         
c          print *,'site4 sitelast ',site4,sitelast 
c          print *,'doyb doye ',doyb,doye
          if( ioerr.ne.0 ) then
             print *,'Error decoding line ',ioerr
          elseif( doyb.eq.0 .and.doye.eq.0 ) then
c            print *,'doys=0 skip '
c           last or dummy entries have 0 for start and end, skip
            continue                       
          elseif( pvflg.eq.'V' ) then
c            ignore velocity steps, which should correspond to position steps
             continue
          elseif( site4.ne.sitelast(1:4) ) then
c           new site, reinitialize and skip to next entry
            ext = "_GPS"
            sitelast = site4//ext
            achr = '1'   
c            print *,'skip write'
          else
            site1 = sitelast 
            if( pvflg.eq.'X' ) then
c             X is exclude data flag. Write this time span with an XCL extension
              ext = "_XCL"
            else
              ext = "_GPS"
              site1(5:8) = "_GPS" 
              if( pt(2:2).ne.'A' ) then
                ext(3:3) = pt(2:2)
              end if
              if( occ.le.9 ) then
                 write(achr,'(I1)') occ
              else
                 achr = char(65 + occ - 10)
                 if( occ.gt.90 ) then
                   call report_stat('FATAL','GLOBK',
     .                  'eqfupd',occ,
     .                  'Too many BREAKS at site',0)
                 endif
              end if
              ext(2:2) = achr
            endif
            site2 = site4 // ext    
c           if end doy zero, must be open ended: set to 2100 1 1 0 
            if( doye.eq.0 ) then
              yre = 2100
              doye = 1
            endif  
            call fix_y2k(yrb)
            call fix_y2k(yre)
            call monday(doyb,monb,dayb,yrb)
            call monday(doye,mone,daye,yre)    
            dsod = sodb
            call ds2hms(yrb,doyb,dsod,hrb,minb,dsec) 
            secb= dsec
            if(secb.ge.30) minb = minb + 1    
            dsod = sode
            call ds2hms(yre,doye,dsod,hre,mine,dsec)
            sece = dsec
            if(sece.ge.30) mine = mine + 1
c            print *,'dates read ',yrb,doyb,sodb,yre,doye,sode  
c            print *,'dates converted ',yrb,monb,dayb,hrb,minb,secb
c     .                                ,yre,mone,daye,hre,mine,sece
             write(ieq,'(a,1x,a4,5x,a8,2(2x,i4,4i3))') 
     .          ' rename ',site4,site2
     .          ,yrb,monb,dayb,hrb,minb
     .          ,yre,mone,daye,hre,mine   
            sitelast = site2        
            achr = char(ichar(achr)+1)
            if( achr.eq.":" ) achr = "A"
          endif                                   
        endif
      enddo       

      print *,'Normal end of eqfupd '

      stop
      end


