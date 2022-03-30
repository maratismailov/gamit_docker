      Subroutine grid_time(call,start,interval,ngrid,nhead,time,nrec )

c     Given decimal day-of-year or record number, return the other from
c     a gridded direct-access file

c     R. King 070813  

      implicit none

      character*1 call
c       'T' :  given time, return the number of the first record of that epoch
c       'R' :  given a record number, return the time

      real*4 start,interval,time              
c        start    : decimal day-of-year of first epoch
c        interval : interval in days between epochs (nominally 0.25)
c        time     : decimal day-of-year of requested or returned epoch

            
      integer*2 nhead
      integer*4 ngrid,nrec
c        ngrid: number of lat/lon values per epoch                    
c        nhead: number of header records 
c        nrec : requested record number or first record of requested epoch

      if( call.eq.'T' ) then
        nrec = nint((time-start)/interval)*ngrid + nhead + 1

      elseif( call.eq.'R') then
        time =  interval*float(nrec-nhead-1)/float(ngrid) + start 

      endif

      return
      end 

