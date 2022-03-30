      subroutine read_list (llist)

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

      logical         readlist,wish,lask
      integer         kchan1,kchan2,ksite1,ksite2,nlist
      integer         il,llist,ipick

c     initialize parameters
      llist = 0
      readlist = .true.
      il = 47

      do 300 while (readlist)

c         ask about reading a list
          write (6,1000) 'Do you wish to read a LIST?'
          wish = lask()

          if (wish) then

              write(6,2000)
              write(6,3000)
              call imenu (ipick,5)

              if (ipick .eq. 5) then
                   open (unit=il,file='cview.list',status='unknown')
              else if (ipick .eq. 4) then
                   open (unit=il,file='dd.srt'    ,status='unknown')
              else if (ipick .eq. 3) then
                   open (unit=il,file='rms.tot',status='unknown')
              else if (ipick .eq. 2) then
                   open (unit=il,file='rms.ful',status='unknown')
              else
                   open (unit=il,file='rms.qui'   ,status='unknown')
              endif

              do 100 while (readlist .and. (llist .lt. maxlist))
                  llist = llist+1

                  if (ipick .eq. 5) then
                      read(il,*,err=100,end=200)
     .                       kchan1,kchan2,ksite1,ksite2
                      nlist = llist
                  else if (ipick .eq. 4)then
                      read(il,6000,err=100,end=200)
     .                       nlist,kchan1,kchan2,ksite1,ksite2
                  else
                      read(il,7000,err=100,end=200)
     .                       nlist,kchan1,kchan2,ksite1,ksite2
                  endif

                  if((llist .gt. 0) .and. (llist .le. maxlist))then

                      klist(llist,1) = kchan2
                      klist(llist,2) = kchan1
                      klist(llist,3) = ksite2
                      klist(llist,4) = ksite1
                      klist(llist,5) = nlist
                   else
                      klist(llist,1) = 0
                      klist(llist,2) = 0
                      klist(llist,3) = 0
                      klist(llist,4) = 0
                      klist(llist,5) = 0
                   endif


  100         continue

          endif

  200     continue
          readlist = .false.
          close(il)

          if((wish) .and. (llist .eq. 0)) then
               write(6,2000)
               write(6,4000)
               readlist = .true.
          endif
  300 continue


 1000 format (1x,a,1x,$)
 2000 format(//)
 3000    format (1x,'You may read one of the following files:'/,
     .   '    (these files are the output of SCANRMS, and SCANDD ',/,
     .   '                            that are sorted by SORTER) ',/,
     .   '  1. rms.qui    - LC-RMS of quick solutions.',/,
     .   '  2. rms.ful    - LC-RMS of full  solutions',/,
     .   '  3. rms.tot    - LC-RMS as in CVIEW',/,
     .   '  4. dd.srt     - a sorted version of "scan.dd".',/,
     .   '  5. cview.list  - a file with format of "C2 C1 S2 S1".',///)
 4000    format (1x,'LIST file is empty (or bad format) - ',
     .                                         'try again.',//)
 6000 format (i4,5x,4x,7x,3x,4i5,4x,6x,4x,6x)
 7000 format (i4,5x,6x,6x,8x,5x,8x,4x,5x,2x,4i6)

      return
      end
c***********************************************************************
      subroutine ed_list (
     .           ichan1,ichan2,isite1,isite2,nsat,ncfls,
     .           llist,ilist,key,newser,redraw,mwarn,mhelp)

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      logical         redraw,newser
      integer         ichan1,ichan2,isite1,isite2
      integer         kchan1,kchan2,ksite1,ksite2
      integer         nsat,ncfls,llist,ilist,nlist
      character*1     key
      integer*2       mwarn(4),mhelp(4)
      character*96    amsg


      newser = .false.
      if(llist .ge. 1) then

          if(key .eq. MOUSE1) ilist = ilist - 1
          if(key .eq. MOUSE2) ilist = ilist + 1
          if(key .eq. MOUSE3) ilist = ilist + 1

          if(ilist .lt. 1) ilist = llist
          if(ilist .gt. llist) ilist = 1

          kchan2 = klist(ilist,1)
          kchan1 = klist(ilist,2)
          ksite2 = klist(ilist,3)
          ksite1 = klist(ilist,4)
          nlist  = klist(ilist,5)

          if ((kchan1 .ge. 1) .and. (kchan2 .le. nsat) .and.
     .        (kchan2 .ge. 0) .and. (kchan1 .le. nsat) .and.
     .        (ksite1 .ge. 1) .and. (ksite2 .le. ncfls).and.
     .        (ksite2 .ge. 0) .and. (ksite1 .le. ncfls)) then

               ichan2 = kchan2
               ichan1 = kchan1
               isite2 = ksite2
               isite1 = ksite1

               if(key .eq. MOUSE2) newser = .true.

               write (amsg,'(''LIST #'',i4,'': '',4(i2,1x))')
     .                     nlist,kchan1,kchan2,ksite1,ksite2

          else
               write (amsg,
     .                '(''LIST #'',i4,'': OUT OF RANGE '',4(i2,1x))')
     .                  nlist,kchan1,kchan2,ksite1,ksite2
          endif
      else
          write (amsg,'(''LIST is empty '')')
      endif

      call gmsg (mwarn,amsg)
      call gmsg (mhelp,amsg)
      redraw = .true.

      return
      end
c***********************************************************************

