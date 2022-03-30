      Program CONVERT_JPLYBIAS

c     Program to convert the times in the JPL yaw_bias_table from seconds-from J2000
c     to YMD HMS for use in svnav.dat

c     Command-line inputs are the old svnav.dat file and yaw_bias_table.

c     R. King 140217
                 
      implicit none

      integer*4 date(5),iyb(50),nval,svn,nblen,ioerr,i
      real*8 t(50),xjd,sec   
      character*1 bias   
      character*256 line
      logical eof                    


      open(unit=1,status='old',file='yaw_bias_table',iostat=ioerr)
      if( ioerr.ne.0 ) then
        write(*,'(a,a30)') 'Error opening yaw_bias_table'
        stop
       endif

       eof = .false.
       do while(.not.eof) 
                
          read(1,'(a)',iostat=ioerr ) line
          if( ioerr.eq.-1 ) then
            eof = .true.
            goto 999 
          endif
          if( line(1:3).eq.'GPS' ) then
            read(line(4:5),'(i2)') svn
            read(line(7:7),'(i1)') nval
            read(line(8:nblen(line)),*) (t(i),iyb(i),i=1,nval)
            write(*,'(a2,1x,i2)') 'SV',svn
            do i=1,nval
c             units of the JPL tables are seconds after J2000 (true JD 2451545.0),
c             convert to a calender date
              xjd = 2451545.0d0 + t(i)/86400.d0 
              call jd_to_ymdhms(xjd,date,sec)
              if(iyb(i).eq.0) then
                bias='U'
              elseif(iyb(i).eq.1) then 
                bias='Y'
              else if(iyb(i).eq.-1) then
                bias = 'A' 
              elseif(iyb(i).eq.2) then
                 bias = 'P'
              elseif(iyb(i).eq.-2) then
                 bias = 'N'
              else
                write(*,*) 'Unknown bias integer '
                bias = 'X'
              endif
              write(*,'(f14.1,i3,1x,5i6,1x,a1)') t(i),iyb(i),date,bias
            enddo
          endif
                                   
       enddo
999    write(*,'(a)') 'Normal stop'
       end

               

              
    




