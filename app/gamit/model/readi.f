      subroutine READI ( iui,iscrn,iprnt,iyr,idoy,isessn,sitecd,jd0,t0
     .                 , clkepc,clkrat,clkacc,clkcub )
C
C     Read receiver clock corrections from external data set (I-file)
C     R.W. King   30 December 1991 - from  READJ and SREAD

c     I-file specifications:
c              
c     time tags receiver clock in UTC at the beginning of the session

c   Input:
c     iui          Unit number for I-file
c     iscrn        Unit number for screen output
c     iprnt        Unit numbe for print output (P-file in MODEL)
c     iyr          year of session
c     idoy         day-of-year of session
c     isessn       session number
c     sited        4-character site code
c     jd0          PEP JD of session
c     t0           seconds-of-day (UTC) of session

c   Output:
c     clkepc       clock epoch correction at initial time (seconds)
c     clkrat       coefficient of t (rate) term (s/s/)
c     clkacc       coefficient of t**2 term (1/s)clock acceleration term
c     clkcub       coefficient of t**3 term (1/s**2)

      implicit none

      include '../includes/dimpar.h'

      character*3   i_time_flag
      character*4   sitecd,site4 
      character*8   atime
      character*60  afmt
      character*256 message

      integer*4 iui,iscrn,iprnt
     .        , isessn,iyr ,idoy
     .        , isessi,iyri,idoyi,imi,idi,ihri,mini
     .        , jd0,jdi,julday,nblen,ioerr,mchkey,i

      real*8 t0,ti,clkepc,clkrat,clkacc,clkcub,seci,utcoff,taiutc


C     Read the header records

      read(iui,'(1x)',iostat=ioerr)
      if (ioerr .ne. 0) then
        call report_stat('FATAL','MODEL','readi',' '
     .  ,'Error reading I-file: line 1',0)
      endif
c      replace the following to allow for shifted header tokens with y2k changes rwk 990805
c      read(iui,'(24x,a3)',iostat=ioerr) i_time_flag
      read(iui,'(24x,a8)') atime 
      i = mchkey(atime,'GPST',8,4) 
      if( i.gt.0 ) then 
        i_time_flag = 'GPS'  
      else
        i_time_flag = 'UTC'
      endif
c       end changes
      if (ioerr .ne. 0) then
        call report_stat('FATAL','MODEL','readi',' '
     .  ,'Error reading I-file: line 2',0)
      endif
      read(iui,'(a60)',iostat=ioerr) afmt
      if (ioerr .ne. 0) then
        call report_stat('FATAL','MODEL','readi',' '
     .  ,'Error reading I-file: line 3',0)
      endif
c     find the end of the format statement
      i = nblen(afmt)
c     choke if you can't find a format statement
      if (i .lt. 1) goto 993
      afmt = afmt(1:i)


c     Find the site code and day number that match the session

      ioerr = 0
      do 200 while (ioerr.eq.0)

        read(iui,fmt=afmt,err=993,END=990,iostat=ioerr) site4
     .      , iyri,idoyi,isessi,ihri,mini,seci
     .      , clkepc,clkrat,clkacc,clkcub
        call check_y2k(iyri)
        if (isessi.eq.0 ) isessi = 1

        if( site4.eq.sitecd ) then

c         first convert the times
          call monday( idoyi,imi,idi,iyri )
          jdi = julday( imi,idi,iyri )
          ti  = ihri*3600.d0 + mini*60.d0 + seci
c         convert from UTC to GPST if necessary
          if( i_time_flag.ne.'GPS' ) then
            utcoff = taiutc(jdi) - 19.d0
            print *,'READI iyri ti jdi utcoff ',iyri,ti,jdi,utcoff
            call timinc(jdi,ti,utcoff)
            call dayjul( jdi,iyri,idoyi ) 
            print *,'READI 2 iyri ti jdi utcoff ',iyri,ti,jdi,utcoff

          endif
c         now check for a match
          if( iyri.eq.iyr.and.idoyi.eq.idoy .and.isessi.eq.isessn ) then
c             check to make sure the times agree with X-file
              if( jdi.ne.jd0 .or. dabs(ti-t0).gt.1.d-10 ) goto 995
c             write the values to the P-file
              if( iprnt.gt.0 ) then
                write(iprnt,75)
75              format(/,'Input I-file clock ploynomial coefficients ')
                write(iprnt,100) clkepc,clkrat,clkacc,clkcub
100             format(  ' Epoch  (seconds)',     12X,D13.6
     .                ,/,' Rate  (seconds/second)',6X,D13.6
     .                ,/,' Acceleration  (1/SEC)', 7X,D13.6
     .                ,/,' Cubic term    (1/S/S)', 7X,D13.6)
              endif
              return
          endif
        endif

  200 continue


990   write(iscrn,992) sitecd,idoy
      write(iprnt,992) sitecd,idoy
      write(message,991) sitecd,idoy
991   format('Site ',a4,' or Day of year ', i3,' not found on I-File',
     .       ' Clock terms set to zero')
992   format(/,1x,'**Warning from READI:  Site ',a4,' or Day of year '
     .       , i3,' not found on I-File',/
     .       , 1x,'  Clock terms set to zero')
      call report_stat('WARNING','MODEL','readi',' ',message,0)
      clkepc = 0.d0
      clkrat = 0.d0
      clkacc = 0.d0  
      clkcub = 0.d0
      return

993   call report_stat('FATAL','MODEL','readi',' ',
     .'File error in I-file, incorrect format descriptor at line 3',0)
      return

995   write(iscrn,996) jdi,ti,jd0,t0
      write(message,997) jdi,ti,jd0,t0
996   format(/,1X,'JD, T =',I7,F14.6,' read from I-File do not agree wit
     1h',/   ,1X,'JD, T =',I7,F14.6,' read from observation file',/
     2        ,1X,'Stop in READI')
997   format('Error wrong DOY: JD, T =',I7,F14.6,' read from I-File do',
     1'not agree with JD, T =',I7,F14.6,' read from observation file')
      call report_stat('FATAL','MODEL','readi',' ',message,0)

      return
      end
