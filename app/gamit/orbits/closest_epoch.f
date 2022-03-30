Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995. All rights reserved.

      Function closest_epoch( nepoch,iwkbc,sowbc,iwk,sow,offset_hrs )

C    R. King  7 October 1995

      implicit none

      include '../includes/dimpar.h'
      include '../includes/orbits.h'

      integer*4 closest_epoch,nepoch,iwkbc(maxbrd),iwk,i

      real*8 sowbc(maxbrd),sow,offset_hrs,tdiff,tdiff_save,secdif

      character*256 message

c     check array lengths
      if( nepoch.gt.maxbrd) then
        write(message,'(a,i2,a,i2)') 'Number of epochs ',nepoch
     .    ,' exceeds array dimension MAXBRD ',maxbrd  
         call report_stat('FATAL','BCTOT','orbits/closet_epoch',' '
     .                    ,message,0) 
      endif


c     store the time difference from the first epoch

      tdiff_save = secdif( iwk,sow,iwkbc(1),sowbc(1) )
      closest_epoch = 1

c     loop over the other epochs to find one closer

      do i = 2, nepoch

         tdiff = secdif( iwk,sow,iwkbc(i),sowbc(i) )

         if( dabs(tdiff).lt.dabs(tdiff_save) ) then
            closest_epoch = i
            tdiff_save = tdiff
         endif

      enddo

c     save the time difference in hours for warnings and
c     possible use in weighting the data.

      offset_hrs = tdiff_save/3600.d0

      return
      end
