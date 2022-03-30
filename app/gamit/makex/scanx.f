      Subroutine scanx ( istart,istop
     .                 , iy1,im1,id1,ih1,imn1,sec1
     .                 , iy2,im2,id2,ih2,imn2,sec2 
     .                 , nchan,numobsx )

c     Scan the X-file to get the start and stop epochs of actual data
c     Written by R. King 19 June 1990 for XTORX; separated 21 Mar 97

      implicit none
                  
c     --maximum dimensions for satellites in session
       include '../includes/dimpar.h'  

c     --file names and unit numbers
      include '../includes/makex.h'

      logical start,lgood

      character*3 x_time_flag
      integer iyr,jdoy,ihr,min,iwkn
     .      , iy1,im1,id1,ih1,imn1,istart
     .      , iy2,im2,id2,ih2,imn2,istop
     .      , iwkn0,iy0,jdoy0
     .      , msat,inter,iepoch,nepoch,ier,jsat,itflag,i,j
     .      , mtime 
     .      , nchan,numobsx(maxsat)

      real*8 sec,sow,sow0,sec1,sec2,delta,utcoff

      character*4  buf4
      mtime = 0

c  Read the header records
       
c     read until END
10    read(uxfile,'(a4)',end=990) buf4
      if( buf4.ne.'END ' ) goto 10

c     read 8 lines for site and antenna offset
      do  i=1,8
       read(uxfile,'(a4)',end=990) buf4  
      enddo

c     read number of channels and sats in channels
      read(uxfile,'(i3)',end=990) nchan
      do  i=1,nchan
       read(uxfile,'(a4)') buf4   
      enddo

C     Skip 1 line
      read(uxfile,'(a4)',end=990) buf4 

c     Read time type (UTC or GPST)
      read(uxfile,'(14x,a3)') x_time_flag 
      if( x_time_flag.eq.'   '.or.x_time_flag.eq.'UTC') then
        mtime = 1
      elseif (x_time_flag.eq.'GPS') then
        mtime = 2
      else
        call report_stat('FATAL','XTORX','scanx',' '
     .                  ,'Unknown time flag',0 )
      endif

      read(uxfile,'(a4)',end=990) buf4

C     Read the start time and interval  
      read(uxfile,50,end=990) iy0,jdoy0,ihr,min,sec,inter
50    format (1X,I2,1X, I3,1X, I2,1X,I2,1X,F6.3,1X, I6)
      call fix_y2k(iy0)

c     convert UTC yr,doy to GPS time wk,sow
      if(mtime.eq.1) itflag = -2
c     convert GPS yr,doy to GPS time wk,sow
      if(mtime.eq.2) itflag = -4
c     UTC yr,doy to GPST wk,sow
      call timcon ( itflag,iwkn0,sow0,iy0,jdoy0
     .            , ihr,min,sec,utcoff)

C     Skip one line and read the number of epochs 
      read(uxfile,'(a4)',end=990) buf4
      read (uxfile,'(i4)',end=990) nepoch  
c     skip 4 more lines
      do  i=1,4
        read(uxfile,'(a)',end=990) buf4  
      enddo

c  Read through the data records, saving the first and last times
c  counting the number of obs for each SV

      start= .false.
      do i = 1,nchan
          numobsx(i) = 0
      enddo

      do 200 i=1,nepoch

      read ( uxfile
     .,      fmt = '(/,2I4,I5,I4,2I3,F11.7)'
     .,      end = 992 )
     .       iepoch,msat,iyr,jdoy,ihr,min,sec
          if( jdoy.gt.0 ) call fix_y2k(iyr)
c          convert UTC to GPS time
           if( iyr.gt.0 ) then
c             UTC yr,doy,hms to GPST wk,sow
              if(mtime.eq.1) itflag = -2
c             GPS yr,doy,hms to GPST wk,sow
              if(mtime.eq.2) itflag = -4
              call timcon( itflag,iwkn,sow,iyr,jdoy,ihr,min,sec,utcoff )
           else
c             compute the observation time if it's not there
              delta = (iepoch-1)*dble(inter)
              call secsum( iwkn0,sow0,delta,iwkn,sow )
              endif
           itflag = +4
c          GPST wk,sow to GPST yr,doy,hms
           call timcon( itflag,iwkn,sow,iyr,jdoy
     .               , ihr,min,sec,utcoff)

      if( .not.start .and. msat.gt.0 ) then
           start = .true.
           iy1 = iyr                      
           call monday( jdoy,im1,id1,iy1 )
           ih1= ihr
           imn1 = min
           sec1 = sec
           istart = iepoch
           endif
      if( start .and. msat.gt.0 ) then
           iy2 = iyr                   
           call monday( jdoy,im2,id2,iy2 )
           ih2= ihr
           imn2 = min
           sec2 = sec
           istop = iepoch
           endif

      do  j=1,msat  
         read( uxfile,fmt = '(10X,2I2)',end=992) ier,jsat
        if( lgood(ier) ) numobsx(jsat) = numobsx(jsat) + 1
      enddo

200   continue 
      
      return
                    
990   write(uinfor,'(a)') 'Unexpected EOF in reading X-file header'
      write(uinfor,'(a)') 'Exit from SCANX but continue' 
      call report_stat('WARNING','XTORX','scanx',' '
     .   ,'Unexpected EOF in reading X-file header--continue',0)
      return

992   write(uinfor,'(a)') 'Unexpected EOF in reading X-file data'
      write(uinfor,'(a)') 'Exit from SCANX but continue'     
      call report_stat('WARNING','XTORX','scanx',' '
     .   ,'Unexpected EOF in X-file data--continue',0)

      return

      end

