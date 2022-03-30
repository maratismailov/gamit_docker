      subroutine wtime (lu,iwkn,sow,otype,label)
c
c     translate a time tag and write it to logical unit lu 
c
      implicit none
c
      integer*4      lu
      integer*4      iyear,iday,hh,mm,iwkn
      integer*4      iflag,len,rcpar
      real*8         sow
      real*8         dummy,seconds
      character*(*)  label,otype
      character*80  prog_name
      character*256 message

      iflag = 0

c     get calling program for report_stat
      len = rcpar(0,prog_name)

c     check if output yr,doy,hr,mm,sec should be UTC or GPS time  
c      call uppers(otype)
      if ( otype(1:3).eq. 'GPS' ) then
        iflag = +4
      elseif ( otype(1:3) .eq. 'UTC' ) then
        iflag = +1
      else
        write(message,'(3a)')'Only conversion from (GPS week,sow) to'
     .,'(yr,doy,hr,min,sec) in GPST or UTC allowed. You tried: ',otype 
        call report_stat('FATAL',prog_name,'lib/wtime',' ',message,0)
      endif

c     write out an epoch to the desired unit, but only if the
c     time is reasonable

      if (iwkn .gt. 0 .    and. iwkn .lt. 5000       .and.
     .    sow  .ge. 0.0d0 .and. sow  .le. 806400.0d0) then
c        convert from gps time back to utc
         if ( iflag .eq. 1 ) then
           call timcon (iflag,iwkn,sow,iyear,
     .                  iday,hh,mm,seconds,dummy)   
           write(lu,20)label,iyear,iday,hh,mm,seconds,iwkn,sow
  20       format(a,' UTC: ',
     .     i4,1x,i3.3,1x,i2.2,':',i2.2,1x,f6.3,'  GPS: ',i4,1x,f12.4)
         else  
           call timcon (iflag,iwkn,sow,iyear,
     .                  iday,hh,mm,seconds,dummy) 
           write(lu,25)label,iyear,iday,hh,mm,seconds,iwkn,sow
  25       format(a,' GPS: ',
     .     i4,1x,i3.3,1x,i2.2,':',i2.2,1x,f6.3,'  GPS: ',i4,1x,f12.4) 
         endif
c
      else
         write(lu,30)label,iwkn,sow
  30     format(a20,' UTC: ',
     .   4x,1x,3x,1x,2x,1x,2x,1x,6x,  '  GPS: ',i4,1x,g12.4)
      endif
c
      return

      end

