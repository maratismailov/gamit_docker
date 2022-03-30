c************************************************************************
c           **** UNWEIGHT ***
      subroutine ed_unwe (
     .               replot,mwarn,nbrack,nobs,
     .               iibl,iibr,l12e)


      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'


      logical       replot
      logical       lgood,lmarg
      integer       nbrack,nobs
      integer       i,iibl,iibr,l12e
      integer*2     mwarn(4)


               if (nbrack .eq. 1 .or .nbrack .eq. 2) then
c                 one bracket is for one point
                  if (nbrack .eq. 1) iibr = iibl

c                 first, buffer the old values into b vectors
                  if (l12e .eq. 1) then
                     call vcopy (kww,wl1,kbb,bl1,nobs)
                  else if (l12e .eq. 4) then
                     call vcopy (kww,wl2,kbb,bl2,nobs)
                  else
                     call vcopy (kww,wl1,kbb,bl1,nobs)
                     call vcopy (kww,wl2,kbb,bl2,nobs)
                  endif

c                 flag the data as unweighted
                  do 3310 i = iibl, iibr
                     if(lgood(kww(i)).or.lmarg(kww(i)))kww(i)=igunwt
 3310             continue
                  replot = .true.
               else if (nbrack .le. 0) then
                  call gmsg (mwarn,'Need brackets for UNWEIGHT')
               else if (nbrack .gt. 2) then
                  call gmsg (mwarn,'Too many brackets for UNWEIGHT')
                  replot = .true.
               endif
               nbrack = 0


      return
      end
c************************************************************************
c           **** REWEIGHT ****
      subroutine ed_rewe (
     .               replot,mwarn,nbrack,nobs,
     .               iibl,iibr,l12e)


      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'


      logical       replot,lmarg
      integer       nbrack,nobs,iibl,iibr,l12e,i
      integer*2     mwarn(4)



               if (nbrack .eq. 1 .or. nbrack .eq. 2) then
c                 first, buffer the old values into b vectors
                  if (l12e .eq. 1) then
                     call vcopy (kww,wl1,kbb,bl1,nobs)
                  else if (l12e .eq. 4) then
                     call vcopy (kww,wl2,kbb,bl2,nobs)
                  else
                     call vcopy (kww,wl1,kbb,bl1,nobs)
                     call vcopy (kww,wl2,kbb,bl2,nobs)
                  endif

c                 one bracket is for one point
                  if (nbrack .eq. 1) iibr = iibl

c                 if the data were marginal before, flag them as good now
c                 do this on L1, L2 or both
                  do 3340 i = iibl, iibr
                     if (lmarg(kww(i))) kww(i) = iggood
 3340             continue
                  replot = .true.
               else if (nbrack .le. 0) then
                  call gmsg (mwarn,'REWEIGHT needs at least 1 bracket.')
               else if (nbrack .gt. 2) then
                  replot = .true.
                  call gmsg (mwarn,'Too many brackets for REWEIGHT.')
               endif
               nbrack = 0

      return
      end
c************************************************************************
c           **** UNDO ****
      subroutine ed_undo (replot,nobs,nelim,npatch,l12e)

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'


      logical       replot
      integer       l12e,nobs,nelim,npatch


c              recover the values from the buffered b-vectors
               if (l12e .eq. 1) then
                  call vcopy (kbb,bl1,kww,wl1,nobs)
               else if (l12e .eq. 4) then
                  call vcopy (kbb,bl2,kww,wl2,nobs)
               else
                  call vcopy (kbb,bl1,kww,wl1,nobs)
                  call vcopy (kbb,bl2,kww,wl2,nobs)
               endif
               nelim = 0
               npatch = 0
               replot = .true.


      return
      end
c************************************************************************
c           **** CANCEL ****
      subroutine ed_cncl (
     .          redraw,replot,nbrack,nelim,npatch,nsave)

         logical    redraw,replot
         integer    nbrack,nsave,nelim,npatch

               redraw = .true.
               replot = .true.
               nbrack = 0
               nelim = 0
               nsave = 0

      return
      end
c************************************************************************
c           **** NOTE ***
      subroutine ed_note (
     .               replot,mwarn,nbrack,iibl,luout)


      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'


      logical       replot
      integer       nbrack,iibl,luout,j
      integer*2     mwarn(4)
      integer*4     igerr
      character*12  alabel(6)
      character*32  vfile



               if (nbrack .eq. 1) then
c                 cycle slip to file name
                  vfile = 'vfile.out'
c%%               append is not ANSI standard, unfortunately
                  open (file    = vfile,
     .                  unit    = luout,
c**  .                  status  = 'append',
c**               status = 'append' is not ANSI standard, unfortunately
     .                  status  = 'unknown',
     .                  position  = 'append',
     .                  err     = 3075,
     .                  iostat  = igerr,
     .                  form    = 'formatted') 
                  print *,'DEBUG Opening vfile.out igerr ',igerr  
 3075             continue
                  if (igerr .eq. 0) then  
c**   this code to trap the undefined 'alabel' and give me an example to check it,
c**   while at the same time avoiding a distracting compiler warning for other users
c**    ---rwk 010731
                     do j=1,6
                        alabel(j) = ' ' 
                     enddo
                     print *,'DEBUG iibl alabel ',iibl,alabel 
                     print *,'DEBUG variable alabel undefined '
                     print *,
     .                  ' --Please report this to rwk@chandler.mit.edu'   
c**    end debug code
                     write (luout,3077) (alabel(j),j=1,4),iibl
 3077                format (4(1x,a12),1x,i5)
                     close (luout)
                  else
                     call gmsg (mwarn,'Could not open V file.')
                     print *,'Could not open V file ',vfile
                     call ferror (6,igerr)
                  endif
               else  if (nbrack .ne. 1) then
                  call gmsg(mwarn,'Need 1 bracket for NOTE')
                  replot = .true.
               endif
               nbrack = 0


      return
      end
c************************************************************************

