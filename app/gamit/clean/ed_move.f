c************************************************************************
c           **** MOVE ***
      subroutine ed_move (
     .        nbrack,nobs,iibr,iibl,l12e,imenu,
     .        redraw,newbox,replot,yadd,mwarn)


      integer         nbrack,nobs,iibl,iibr,l12e,imenu
      logical         replot,newbox,redraw
      real*8          yadd
      integer*2       mwarn(4)

c              assume user wants to move from iibl to iibr inclusive
c              and that this has been determined from the 2 brackets
               if (nbrack .eq. 1 .or. nbrack .eq. 2) then
c                 will occur on all points to the right of the bracket
                  if (nbrack .eq. 1) iibr = nobs

c                 This operation is only allowed on L1 or L2
                  if (l12e .eq. 1 .or. l12e .eq. 4) then
c                    set up manual MOVE
                     imenu = 12
                     redraw = .true.
                     newbox = .true.
                  else
                     call gmsg(mwarn,'May not move this observable.')
                  endif
               else if (nbrack .le. 0) then
                  call gmsg (mwarn,'Need at least 1 bracket for MOVE.')
                  replot = .true.
               else if (nbrack .gt. 2) then
                  replot = .true.
               endif
               yadd = 0.0d0
               nbrack = 0

      return
      end
c************************************************************************
c           **** MOVE MENU (menu12) ***
      subroutine ed_menu12 (
     .        nobs,iibl,iibr,l12e,icomm,itemp,isite1,ichan1,
     .        imenu,yadd,ymult,key,replot,redraw,newbox,mwarn,amsg)

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         nobs,iibl,iibr,l12e
      character*1     key
      integer         isite1,ichan1,imenu,icomm
      integer         i,itemp
      logical         replot,redraw,newbox,ltemp
      real*8          yadd,ymult
      character*96    amsg
      integer*2       mwarn(4),iabs2



c           MOVE menus are set up, choose amount to move

            if (icomm .ge. 11 .and. icomm .le. 20) then
c              subtract, reset, or add?
               if (key .eq. MOUSE1) then
                  ymult = -1.0d0
               else if  (key .eq. MOUSE2) then
                  yadd = 0.0d0
                  ymult = 1.0d0
               else
                  ymult = +1.0d0
               endif
            endif


c           Determine amount to add.
c           Smallest allowable cycle slip is 1/abs(2*lambda)
c           This is a factor of 2 smaller than what we understand, or
c           what is currently allowed by PATCH or REPAIR.
c           Allow it here to deal with those weirdo cases.
            ltemp = .true.
            if (icomm .eq. 11) then
c              Is a quarter-cycle slip allowed?
               if (l12e .eq. 1) then
c                  itemp = 1
                  itemp = 2
               else if (l12e .eq. 4) then
                  itemp = 2
               endif

               if (iabs2(lambds(isite1,ichan1,itemp)) .eq. 2) then
                  yadd = yadd + ymult * 0.25d0
               else
                  call gmsg
     .            (mwarn,'No quarter-cycle slips in this observable.')
                  ltemp = .false.
               endif
            else if  (icomm .eq. 12) then
               yadd = yadd + ymult * 0.50d0
            else if  (icomm .eq. 13) then
               yadd = yadd + ymult * 1.0d0
            else if  (icomm .eq. 14) then
               yadd = yadd + ymult * 2.0d0
            else if  (icomm .eq. 15) then
               yadd = yadd + ymult * 5.0d0
            else if  (icomm .eq. 16) then
               yadd = yadd + ymult * 10.0d0
            else if  (icomm .eq. 17) then
               yadd = yadd + ymult * 50.0d0
            else if  (icomm .eq. 18) then
               yadd = yadd + ymult * 100.0d0
            else if  (icomm .eq. 19) then
               yadd = yadd + ymult * 500.0d0
            else if  (icomm .eq. 20) then
               yadd = yadd + ymult * 1000.0d0

c              **** CANCEL **** don't really want to move
            else if (icomm .eq. 1 .or. key .eq. MOUSE2) then
               yadd = 0.
               redraw = .true.
               replot = .true.
               newbox = .true.
               imenu = 0

c              **** PERFORM **** the move previously determined
            else if  (icomm .eq. 2 .or. key .eq. MOUSE3) then
c              first, buffer the old values into b vectors
               if (l12e .eq. 1) then
                  call vcopy (kww,wl1,kbb,bl1,nobs)
               else if (l12e .eq. 4) then
                  call vcopy (kww,wl2,kbb,bl2,nobs)
               endif

c              add yadd to all points
               if (l12e .eq. 1) then
                  do 4250 i=iibl,iibr
                     wl1(i) = wl1(i) + yadd
 4250             continue
               else if (l12e .eq. 4) then
                  do 4260 i=iibl,iibr
                     wl2(i) = wl2(i) + yadd
 4260             continue
               endif
               redraw = .true.
               replot = .true.
               newbox = .true.
               imenu = 0
            end if


c           show the user the move amount
            if (ltemp) then
               write (amsg,4265) yadd
 4265           format (' Add ',f20.2)
               call gmsg (mwarn, amsg)
            endif

      return
      end
c************************************************************************

