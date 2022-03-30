c**********************************************************************
      subroutine  nextobserv (
     .           ichan1,ichan2,isite1,isite2,nsat,ncfls,
     .           iseek,llist,ilist,mwarn,seekob)

c     find the next observation 4-tuple

      implicit none

      logical         seekob
      integer         ichan1,ichan2,isite1,isite2,iseek
      integer         nsat,ncfls,llist,ilist
      integer*2       mwarn(4)


c     find the next observation 4-tuple
c     there are 3 modes of search: forward seek, backward seek, and list.

c     list search
      if (iseek .eq. 3) then
          call  seeklist (
     .           ichan1,ichan2,isite1,isite2,nsat,ncfls,
     .           llist,ilist,mwarn,seekob)

c     backward search
      else if (iseek .eq. 2) then
          call  seekback
     .       (ichan1,ichan2,nsat,isite1,isite2,ncfls,seekob)

c     forward search - default
      else
          call  nextob
     .       (ichan1,ichan2,nsat,isite1,isite2,ncfls,seekob)
      endif

      return
      end
c
c**********************************************************************
      subroutine nextob
     .       (ichan1,ichan2,nsat,isite1,isite2,ncfls,seekob)

c     find the next observation 4-tuple

c     K. Feigl May 89

      implicit none

      integer ichan1,ichan2,nsat,isite1,isite2,ncfls
      logical seekob

      integer jj0,jj1

c     come here if there is no data
 5    continue

      if (isite2 .eq. 0 .and. ichan2 .eq. 0) then
         ichan1 = ichan1 + 1
         if (ichan1 .gt. nsat) then
            isite1 = isite1 + 1
            if (isite1 .gt. ncfls) then
               seekob = .false.
               goto 10
            endif
            ichan1 = 1
         endif
      else if (ichan2 .eq. 0 .and. isite1 .ne. 0) then
         ichan1 = ichan1 + 1
         if (ichan1 .gt. nsat-1) then
            isite1 = isite1 + 1
            if (isite1 .gt. ncfls) then
               isite2 = isite2 + 1
               if (isite2 .gt. ncfls-1) then
                  seekob = .false.
                  goto 10
               endif
               isite1 = isite2 + 1
            endif
            ichan1 = 1
         endif
      else if (isite2 .eq. 0 .and. ichan1 .ne. 0) then
         ichan1 = ichan1 + 1
         if (ichan1 .gt. nsat) then
            ichan2 = ichan2 + 1
            if (ichan2 .gt. nsat-1) then
               isite1 = isite1 + 1
               if (isite1 .gt. ncfls) then
                  seekob = .false.
                  goto 10
               endif
               ichan2 = 1
            endif
            ichan1 = ichan2 + 1
         endif
      else if (ichan1 .ne. 0 .and. isite1 .ne. 0) then
         if (ichan1 .gt. ichan2 .or. isite1 .gt. isite2) then
            ichan1 = ichan1 + 1
            if (ichan1 .gt. nsat) then
               ichan2 = ichan2 + 1
               if (ichan2 .gt. nsat-1) then
                  isite1 = isite1 + 1
                  if (isite1 .gt. ncfls) then
                     isite2 = isite2 + 1
                     if (isite2 .gt. ncfls-1) then
                        seekob = .false.
                        goto 10
                     endif
                     isite1 = isite2 + 1
                  endif
                  ichan2 = 1
               endif
               ichan1 = ichan2 + 1
            endif
         else
            ichan2 = ichan2 + 1
            if (ichan2 .gt. nsat) then
               ichan1 = ichan1 + 1
               if (ichan1 .gt. nsat-1) then
                  isite2 = isite2 + 1
                  if (isite2 .gt. ncfls) then
                     isite1 = isite1 + 1
                     if (isite1 .gt. ncfls-1) then
                        seekob = .false.
                        goto 10
                     endif
                     isite2 = isite1 + 1
                  endif
                  ichan1 = 1
               endif
               ichan2 = ichan1 + 1
            endif
         endif
      endif

      call bound (ichan1,ichan2,isite1,isite2,jj0,jj1)

c     if the index of the last point is smaller than the index
c     of the index of the first point, then there is no data in
c     this n-tuple, so keep looking
      if (jj0 .gt. jj1) goto 5

c     come here after last series
 10   continue
      if (.not. seekob) then
         if (isite2 .eq. 0) then
            isite1 = 1
         else
            if (ichan1 .gt. ichan2 .or. isite1 .gt. isite2) then
               isite1 = 2
               isite2 = 1
            else
               isite1 = 1
               isite2 = 2
            endif
         endif
         if (ichan2 .eq. 0) then
            ichan1 = 1
         else
            if (ichan1 .gt. ichan2 .or. isite1 .gt. isite2) then
               ichan1 = 2
               ichan2 = 1
            else
               ichan1 = 1
               ichan2 = 2
            endif
         endif
      endif
      return
      end
c
c**********************************************************************
      subroutine seekback
     .       (ichan1,ichan2,nsat,isite1,isite2,ncfls,seekob)

c     find the next observation 4-tuple

c     S. Wdowinski 1991

      implicit none

      integer ichan1,ichan2,nsat,isite1,isite2,ncfls
      logical seekob

      integer jj0,jj1

c     come here if there is no data
 5    continue

      if (isite2 .eq. 0 .and. ichan2 .eq. 0) then
         ichan1 = ichan1 - 1
         if (ichan1 .lt. 1) then
            isite1 = isite1 - 1
            if (isite1 .lt. 1) then
               seekob = .false.
               goto 10
            endif
            ichan1 = nsat
         endif
      else if (ichan2 .eq. 0 .and. isite1 .ne. 0) then
         ichan1 = ichan1 - 1
         if (ichan1 .lt. 1) then
            isite1 = isite1 - 1
            if (isite1 .le. isite2) then
               isite2 = isite2 - 1
               if (isite2 .lt. 1) then
                  seekob = .false.
                  goto 10
               endif
               isite1 = ncfls
            endif
            ichan1 = nsat
         endif
      else if (isite2 .eq. 0 .and. ichan1 .ne. 0) then
         ichan1 = ichan1 - 1
         if (ichan1 .le. ichan2) then
            ichan2 = ichan2 - 1
            if (ichan2 .lt. 1) then
               isite1 = isite1 - 1
               if (isite1 .lt. 1) then
                  seekob = .false.
                  goto 10
               endif
               ichan2 = nsat-1
            endif
            ichan1 = nsat
         endif
      else if (ichan1 .ne. 0 .and. isite1 .ne. 0) then
         if (ichan1 .gt. ichan2 .or. isite1 .gt. isite2) then
            ichan1 = ichan1 - 1
            if (ichan1 .le. ichan2) then
               ichan2 = ichan2 - 1
               if (ichan2 .lt. 1) then
                  isite1 = isite1 - 1
                  if (isite1 .lt. isite2) then
                     isite2 = isite2 - 1
                     if (isite2 .lt. 1) then
                        seekob = .false.
                        goto 10
                     endif
                     isite1 = ncfls
                  endif
                  ichan2 = nsat - 1
               endif
               ichan1 = nsat
            endif
         else
            ichan2 = ichan2 - 1
            if (ichan2 .le. ichan1) then
               ichan1 = ichan1 - 1
               if (ichan1 .lt. 2) then
                  isite2 = isite2 - 1
                  if (isite2 .le. isite1) then
                     isite1 = isite1 - 1
                     if (isite1 .lt. 2) then
                        seekob = .false.
                        goto 10
                     endif
                     isite2 = nsat-1
                  endif
                  ichan1 = nsat
               endif
               ichan2 = ncfls - 1
            endif
         endif
      endif
      call bound (ichan1,ichan2,isite1,isite2,jj0,jj1)

c     if the index of the last point is smaller than the index
c     of the index of the first point, then there is no data in
c     this n-tuple, so keep looking
      if (jj0 .gt. jj1) goto 5

c     come here after last series
 10   continue
      if (.not. seekob) then
         if (isite2 .eq. 0) then
            isite1 = 1
         else
            if (ichan1 .gt. ichan2 .or. isite1 .gt. isite2) then
               isite1 = 2
               isite2 = 1
            else
               isite1 = 1
               isite2 = 2
            endif
         endif
         if (ichan2 .eq. 0) then
            ichan1 = 1
         else
            if (ichan1 .gt. ichan2 .or. isite1 .gt. isite2) then
               ichan1 = 2
               ichan2 = 1
            else
               ichan1 = 1
               ichan2 = 2
            endif
         endif
      endif
      return
      end
c
c**********************************************************************
      subroutine seeklist (
     .           ichan1,ichan2,isite1,isite2,nsat,ncfls,
     .           llist,ilist,mwarn,seekob)

c     find the next observation 4-tuple

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

      logical         seekob
      integer         ichan1,ichan2,isite1,isite2
      integer         kchan1,kchan2,ksite1,ksite2
      integer         nsat,ncfls,llist,ilist,nlist
      integer*2       mwarn(4)
      character*96    amsg      

      integer*4 len,rcpar

      character*80 prog_name
 
c     get calling program name to know which version of gmsg to call
      len = rcpar(0,prog_name)



      if(llist .ge. 1) then

          ilist = ilist + 1

          if (ilist .le. llist) then

              kchan2 = klist(ilist,1)
              kchan1 = klist(ilist,2)
              ksite2 = klist(ilist,3)
              ksite1 = klist(ilist,4)
              nlist  = klist(ilist,5)

              if ((kchan1 .ge. 1) .and. (kchan2 .le. nsat) .and.
     .           (kchan2 .ge. 0) .and. (kchan1 .le. nsat) .and.
     .           (ksite1 .ge. 1) .and. (ksite2 .le. ncfls).and.
     .           (ksite2 .ge. 0) .and. (ksite1 .le. ncfls)) then

                   ichan2 = kchan2
                   ichan1 = kchan1
                   isite2 = ksite2
                   isite1 = ksite1

                   write (amsg,'(''LIST #'',i4,'': '',4(i2,1x))')
     .                     nlist,kchan1,kchan2,ksite1,ksite2

              else
                   write (amsg,
     .                '(''LIST #'',i4,'': OUT OF RANGE '',4(i2,1x))')
     .                            nlist,kchan1,kchan2,ksite1,ksite2
              endif

          else
              write (amsg,'(''end of LIST '')')
              seekob = .false.
          endif
      else
          write (amsg,'(''LIST is empty '')')
          seekob = .false.
      endif
           
      if( prog_name(1:5).eq.'cview' ) then
          call gmsg(mwarn, amsg)
      else
          call gmsgd(mwarn, amsg) 
      endif

      return
      end
c
c**********************************************************************
