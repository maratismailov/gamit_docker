      subroutine obsrep( nobs, GPSWeek, rtobs)

      implicit none 

*     Suboutine to report results

      include 'obsInternal.h'

      type ( obsInternal ) rtobs

      integer*4 nobs, GPSWeek 

!      integer*4 i   ! Loop counter

      print *,'Num Data, GPS Week ', nobs, GPS Week

      print *,'SAT ',rtobs%satNum,' ',rtobs%flags, ' ', rtobs%satSys
      print *,'Data L1 ',rtobs%P1, rtobs%P2, rtobs%L1, rtobs%L2

      return
      end
