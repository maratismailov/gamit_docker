      subroutine runtim( dattim )

c     Program to return the date and time of the run in a character string
c     R. King 940719

      character*2 buf2
      character*4 buf4
      character*19   dattim
      integer*4 iyear,imonth,iday,ihr,imn,isec,ihnsec

      dattim=' '
      call getdat(iyear,imonth,iday)
      call gettim(ihr,imn,isec,ihnsec)
      write (buf2,'(i2)') iday
      read  (buf2,'(a2)') dattim(1:2)
      dattim(3:3) = '-'
      write (buf2,'(i2)') imonth
      read  (buf2,'(a2)') dattim(4:5)
      dattim(6:6)= '-'
      write (buf4,'(i4)') iyear
      read  (buf4,'(a4)') dattim(7:10)
      write (buf2,'(i2)') ihr
      read  (buf2,'(a2)') dattim(12:13)
      dattim(14:14)=':'
      write (buf2,'(i2)') imn
      read  (buf2,'(a2)') dattim(15:16)
      dattim(17:17)=':'
      write (buf2,'(i2)') isec
      read  (buf2,'(a2)') dattim(18:19)
      return
      end

