      program seek
      logical seekob
      integer ichan1,ichan2,isite1,isite2,nsat,ncfls

      ncfls = 4
      nsat  = 3

      read *,ichan1,ichan2,isite1,isite2

      seekob = .true.

      do while (seekob)
c        look for next series with data
         if (seekob) then
            if (ichan2 .ne. 0 .and. isite2 .ne. 0) then
               ichan2 = ichan2 + 1
               if (ichan2 .gt. nsat) then
                  ichan1 = ichan1 + 1
                  if (ichan1 .gt. nsat-1) then
                     isite2 = isite2 + 1
                     if (isite2 .gt. ncfls) then
                        isite1 = isite1 + 1
                        if (isite1 .gt. ncfls-1) then
                           goto 100
                        endif
                        isite2 = isite1 + 1
                     endif
                     ichan1 = 1
                  endif
                  ichan2 = ichan1 + 1
               endif
            else if (ichan2 .eq. 0 .and. isite2 .ne. 0) then
               ichan1 = ichan1 + 1
               if (ichan1 .gt. nsat-1) then
                  isite2 = isite2 + 1
                  if (isite2 .gt. ncfls) then
                     isite1 = isite1 + 1
                     if (isite1 .gt. ncfls-1) then
                        goto 100
                     endif
                     isite2 = isite1 + 1
                  endif
                  ichan1 = 1
               endif
            else if (isite2 .eq. 0 .and. ichan2 .ne. 0) then
               ichan2 = ichan2 + 1
               if (ichan2 .gt. nsat) then
                  ichan1 = ichan1 + 1
                  if (ichan1 .gt. nsat-1) then
                     isite1 = isite1 + 1
                     if (isite1 .gt. ncfls) then
                        goto 100
                     endif
                     ichan1 = 1
                  endif
                  ichan2 = ichan1 + 1
               endif
            else if (isite2 .eq. 0 .and. ichan2 .eq. 0) then
               ichan1 = ichan1 + 1
               if (ichan1 .gt. nsat) then
                  isite1 = isite1 + 1
                  if (isite1 .gt. ncfls) then
                      goto 100
                  endif
                  ichan1 = 1
               endif
            endif
         endif
         print '(4(i2,1x))',ichan1,ichan2,isite1,isite2
      enddo

 100  continue
      seekob = .false.
c     newobs = .false.
c     redraw = .true.
c     newbox = .false.
c     imenu = -1
      if (isite2 .ne. 0) then
         isite1 = ncfls-1
         isite2 = ncfls
      else
         isite1 = ncfls
      endif
      if (ichan2 .ne. 0) then
         ichan1 = nsat -1
         ichan2 = nsat
      else
         ichan1 = nsat
      endif
      call message (alarm,'Last series!')
      print '(4(i2,1x))',ichan1,ichan2,isite1,isite2

      stop
      end


      subroutine message (a,msg)
      character*(*) msg

      print *,msg

      return
      end
