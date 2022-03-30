      Subroutine AVCLCK( iprnt,clksav,ischan,nchan,ier
     .                  , rclock0,itimel,bad)

c     Average receiver clock offsets to form RCLOCK0

      include '../includes/dimpar.h'
      include '../includes/errflg.h'

      character*256 message
      real*8 clksav(maxsat),rclock0
      integer*4 ischan(maxsat),ier(maxsat),nchan,itimel
      logical bad

      real*8 rsum,rms,clkerr
      integer ngood,mgood,iprt,i
      logical lgood
         
cd      print *,'AVCLCK clksav ',clksav
cd      print *,'AVCLCK rclock0 itime1 bad ',rclock0,itime1,bad 

c     perform a mean
      rsum = 0.0d0
      ngood = 0
      DO 830 I=1,NCHAN
c        Skip if no data, deleted, hardware error, or unweighted
c         print *,'epoch ichan ier clksav ',itimel,i,ier(i),clksav(i)
         if (lgood(ier(i))) then
            rsum = rsum + clksav(i)
            ngood = ngood + 1
         endif
 830  continue

      if (ngood .gt. 0) then
         rclock0 = rsum/dble(ngood)
      else
         rclock0 = 0.d0
      endif

c     get an RMS for outlier detection
      rms = 0.0d0
      ngood = 0
      do 832 i = 1,nchan
         if (lgood(ier(i))) then
c           deviation from mean in microseconds
            clkerr = (clksav(i) - rclock0) * 1.d6
            rms = rms + clkerr*clkerr
            ngood = ngood + 1
         endif
 832  continue

      if (ngood .gt. 0) then
         rms = dsqrt(rms)/dble(ngood)
         if (rms .lt. 0.5d0) then
c          0.5 microsecond is good enough
           rms = 0.5d0
         endif
      else
c        1 microsecond
         rms = 1.d0
      endif

c     check for outliers and make new mean
      rsum = 0.0d0
      ngood = 0
      mgood = 0
      DO 840 I=1,NCHAN
c        Skip if no data, deleted, hardware error, or unweighted
         if (lgood(ier(i))) then
c           deviation from mean in microseconds
            clkerr = (clksav(i) - rclock0) * 1.d6
            if (dabs(clkerr) .lt. 2.0d0*rms) then
c              only use if good to a microsecond
                  rsum = rsum + clksav(i)
                  ngood = ngood + 1
               mgood = mgood + 1
            else
c              warn if off by more than 5 microseconds
               write(iprnt,810) ischan(i),itimel,clkerr
 810           format(1x,'MODEL warning in AVCLCK: RCLOCK from ',
     .         'PRN',i2.2,' at epoch ',i4,1x,f10.1,
     .         ' micros from mean')
               write(message,820) ischan(i),itimel,clkerr
 820           format('RCLOCK from PRN',i2.2,' at epoch '
     .         ,i4,1x,f10.1,' micros from mean')
               call report_stat('WARNING','MODEL','avclck',' '
     .         ,message,0)
            endif
         endif
 840  continue

      if (ngood .gt. 0) then
         rclock0 = rsum/dble(ngood)
      else
         rclock0 = 0.d0
      endif

c     check final RMS
      rsum = 0.d0
      do 850 i = 1,nchan
         if (lgood(ier(i))) then
            clkerr = (clksav(i) - rclock0) * 1.d6
c       print*,' clock offset for ',ischan(i),' is ',clksav(i)

            if (dabs(clkerr) .lt. 2.0d0*rms) then

c  debug to check size of clock residuals
c       print*,' clock residual for ',ischan(i),' is ',clkerr
               rsum = clkerr*clkerr
            endif
         endif
 850  continue
c       print*,itimel,' rclock from pseudoranges is ',rclock0
c       stop

      if (mgood .gt. 0) then
         rms = dsqrt(rsum)/dble(mgood)
         if (rms .gt. 1.d0) then
            write (iprnt,860) itimel,rms
  860       format(1x,'MODEL warning in AVCLCK: '
     .       ,' RMS for RCLOCK at epoch ',i4,' is ',f8.1
     .       ,' microseconds')
            write (message,870) itimel,rms
  870       format('RMS for RCLOCK at epoch ',i4,' is ',f8.1
     .      ,' microseconds')
            call report_stat('WARNING','MODEL','avclck',' ',message,0)
         endif
         bad = .false.
      else
         bad = .true.
      endif

      return
      end
