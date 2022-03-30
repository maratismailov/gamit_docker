      subroutine wclock (lu,afmt,gnss,nprn,iwkn,sow,
     .    prgl1,svclock,rcvclock,iwkne,sowe )
c
c     write the appropiate data to the clock file
c         
c       lu  unit number of k-file
c       afmt   format for k-file record
c       gnss    c*1  GNSS id
c       nprn    i*4  PRN of SV
c       iwkn    i*4  GPST week number of entry
c       sow     r*8  GPST seconds-of-week of entry  
c       prngl1  r*8  Upper-frequency pseudo-range 
c       svlock  r*8  SV clock offsets from j-file
c       iwkne   i*4  GPST week number of nav-file values for orbit (=0 if sp3)
c       sowe    r*8  GPST seconds-of-week of nav-file values



      implicit none

C     logical unit (assumed open)
      character*1 gnss 
      integer*4   lu,nprn,iwkn,iwkne,iy,hour,min,itflag,id

      real*8   sow,sowe,ltvel,dummy,
     .         prgl1,svclock,rcvclock,seconds
                                 
      character*(*) afmt 

      data ltvel   /2.99792458d+08/


c     case conversion functions

c     convert from gps (wknum,secof wk) to utc time (hms format)
      itflag = +1
      call timcon (itflag,iwkn,sow,iy,id,hour,min,
     1   seconds,dummy)

c** old format (no headeror version number, site-name included)
c      site_upper = site
c      call uppers(site_upper)
c      write(lu,516) site_upper,nprn, iy, id, hour, min, seconds,
c     .              prgl1/ltvel,svclock,rcvclock
c 516  format(1x,a4,2x,i2,2x,i4,1x,i3,1x,i2,1x,i2,f9.4,3(2x,f12.8))

      if( iwkne.gt.0 ) then
c       write an entry computed from a navigation file 
        write(lu,afmt) iy,id,hour,min,seconds,iwkn,sow,gnss,nprn
     .     ,prgl1/ltvel,svclock,rcvclock,iwkne,sowe
      else
c       write an etnry computed from an sp3 file
        write(lu,afmt) iy,id,hour,min,seconds,iwkn,sow,gnss,nprn
     .     ,prgl1/ltvel,svclock,rcvclock
      endif
 
      return
      end

