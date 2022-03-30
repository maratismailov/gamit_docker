      subroutine openf( lunit,nam,how,frm,acs,forreal )
c     open a file

      implicit none

      include '../includes/makex.h'

      character*(*)   nam,how,frm,acs
      character*256   message
      integer*4      lunit,ioerr
      logical        forreal

      if (forreal) then
         open (unit   =  lunit,
     .         file   =  nam,
     .         status =  how,
     .         form   =  frm,
     .         access =  acs,
     .         iostat =  ioerr)

         if (ioerr .ne. 0) then
            if (lunit .eq. ubatch .or.
     .          lunit .eq. usceno .or. 
     .          lunit .eq. ucoord .or. 
     .          lunit .eq. usited ) then
                call report_stat('FATAL','MAKEX','openf',nam,
     .                       'Error opening file: ',ioerr)
            else if (lunit .eq. urinex .or.
     .          lunit .eq. uxfile .or.
     .          lunit .eq. usp3 .or.
     .          lunit .eq. usvclk .or.
     .          lunit .eq. uclock) then
C               go to next session    
                call report_stat('FATAL','MAKEX','openf',nam,
     .                           'Error opening file: ',ioerr)
               call closem
            else
               call closem
               call report_stat('FATAL','MAKEX','openf',nam,
     .                         'Unknown error: ',ioerr)
            endif
         else   
            write(message,'(a,a80)') 'Opened: ',nam
            write (uinfor,'(1x,a)') message
            call report_stat('STATUS','MAKEX','openf',' ',message,0)
         endif
      endif
      return
      end



