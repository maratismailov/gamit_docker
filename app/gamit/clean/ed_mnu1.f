c************************************************************************
      subroutine ed_seek (newser,seekob,nobs)
         logical    newser,seekob
         integer    nobs

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

c     first, buffer the old values into b vectors
      call vcopy (kww,wl1,kbb,bl1,nobs)
      call vcopy (kww,wl2,kbb,bl2,nobs)

      newser = .true.
      seekob = .true.

      return
      end
c************************************************************************
      subroutine ed_plot (newser,nobs)
      logical    newser
      integer    nobs

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

c     first, buffer the old values into b vectors
      call vcopy (kww,wl1,kbb,bl1,nobs)
      call vcopy (kww,wl2,kbb,bl2,nobs)

      newser = .true.

      return
      end
c************************************************************************
      subroutine ed_file (newser,penplt)
      logical    newser,penplt

      newser = .true.
      penplt = .true.

      return
      end
c************************************************************************
      subroutine ed_1way (ichan1,ichan2,isite1,isite2,
     .                         mwarn,amsg,newser,redraw)
      logical    redraw,newser
      integer    ichan1,ichan2,isite1,isite2
      integer*2       mwarn(4)
      character*96    amsg


      write (amsg,100) ichan1,ichan2,isite1,isite2
 100  format('1-WAY combination of: ',4i3)
      call gmsg(mwarn, amsg)

      ichan2 = 0
      isite2 = 0

      redraw = .true.
      newser = .true.

      return
      end
c************************************************************************
      subroutine ed_sat1 (MOUSE1,MOUSE2,MOUSE3,key,
     .                            ichan1,nsat,redraw)
      logical    redraw
      integer    ichan1,nsat
      character*1 MOUSE1,MOUSE2,MOUSE3,key


      if (key .eq. MOUSE1) then
         ichan1 = ichan1 - 1
      else if (key .eq. MOUSE2) then
         ichan1 = 1
      else if (key .eq. MOUSE3) then
         ichan1 = ichan1 + 1
      endif
      if (ichan1 .gt. nsat) then
         ichan1 = 1
      else if (ichan1 .lt. 1) then
         ichan1 = nsat
      endif
      redraw = .true.

      return
      end
c************************************************************************
      subroutine ed_sat2 (MOUSE1,MOUSE2,MOUSE3,key,
     .                            ichan2,nsat,redraw)
      logical    redraw
      integer    ichan2,nsat
      character*1 MOUSE1,MOUSE2,MOUSE3,key

      if (key .eq. MOUSE1) then
         ichan2 = ichan2 - 1
      else if (key .eq. MOUSE2) then
         ichan2 = 0
      else if (key .eq. MOUSE3) then
         ichan2 = ichan2 + 1
      endif
      if (ichan2 .gt. nsat) then
         ichan2 = 0
      else if (ichan2 .lt. 0) then
         ichan2 = nsat
      endif
      redraw = .true.

      return
      end
c************************************************************************
      subroutine ed_site1 (MOUSE1,MOUSE2,MOUSE3,key,
     .                            isite1,ncfls,redraw)

c     change choice of site 1

      logical    redraw
      integer    isite1,ncfls
      character*1 MOUSE1,MOUSE2,MOUSE3,key

      if (key .eq. MOUSE1) then
         isite1 = isite1 - 1
      else if (key .eq. MOUSE2) then
         isite1 = 1
      else if (key .eq. MOUSE3) then
         isite1 = isite1 + 1
      endif
      if (isite1 .gt. ncfls) then
         isite1 = 1
      else if (isite1 .lt. 1) then
         isite1 = ncfls
      endif
      redraw = .true.

      return
      end
c************************************************************************
      subroutine ed_site2 (MOUSE1,MOUSE2,MOUSE3,key,
     .                            isite2,ncfls,redraw)
c     change choice of site 2

      logical    redraw
      integer    isite2,ncfls
      character*1 MOUSE1,MOUSE2,MOUSE3,key

      if (key .eq. MOUSE1) then
         isite2 = isite2 - 1
      else if (key .eq. MOUSE2) then
         isite2 = 0
      else if (key .eq. MOUSE3) then
         isite2 = isite2 + 1
      endif
      if (isite2 .gt. ncfls) then
         isite2 = 0
      else if (isite2 .lt. 0) then
         isite2 = ncfls
      endif
      redraw = .true.

      return
      end
c************************************************************************
      subroutine ed_seekch (MOUSE1,MOUSE2,MOUSE3,key,
     .                     iseek,newbox,redraw)
      logical    redraw,newbox
      integer    iseek
      character*1 MOUSE1,MOUSE2,MOUSE3,key

      if (key .eq. MOUSE1) then
         iseek = iseek - 1
      else if (key .eq. MOUSE2) then
         iseek = 1
      else if (key .eq. MOUSE3) then
         iseek = iseek + 1
      else
         iseek = 1
      endif
      if (iseek .gt. 3) then
         iseek = 1
      else if (iseek .lt. 1) then
         iseek = 3
      endif
      redraw = .true.
      newbox = .true.

      return
      end
c************************************************************************
      subroutine ed_print (mwarn,ibmx0,ibmy0,nx,ny,
     .                                  lcolor,idump)
      integer*2     mwarn(4),ibmx0,ibmy0,nx,ny
      integer       idump
      logical       lcolor

c     make a screen dump

c     Does not work: June 1991

      call gmsg (mwarn,'PRINT does not work.')
      return

c     Comment the following until made to work, to avoid compiler warning
c     about unaccessible statements.
c      call gdump (ibmx0,ibmy0,nx,ny,lcolor,idump)
c      if (idump .eq. 0) then
c         call gmsg (mwarn,'PRINT of screen queued.')
c      else
c         call gmsg (mwarn,'PRINT failed.')
c      endif
c      return
      end
c************************************************************************
      subroutine ed_l12v (key,l12v,icomm,redraw)
      logical    redraw
      integer    l12v(*)
      integer    icomm,itemp
      character*1 key
      integer    i

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      itemp = icomm - nmenus/2
c     not too many choices
      if (itemp .lt. 1 .or. itemp .gt. maxwin) then
         return
      endif

      if (key .eq. MOUSE1) then
         l12v(itemp) = l12v(itemp) - 1
c        Multipath (# 13-16 not implemented yet!)
         if (l12v(itemp) .eq. 16) l12v(itemp) = 12
      else if (key .eq. MOUSE2) then
c        default should be L1,L2,LC,LG,CL
         if (icomm .eq. nmenus/2+1) then
            l12v(itemp) = 1
         else if(icomm .eq. nmenus/2+2) then
            l12v(itemp) = 4
         else if(icomm .eq. nmenus/2+3) then
            l12v(itemp) = 3
         else if(icomm .eq. nmenus/2+4) then
            l12v(itemp) = 2
         else if(icomm .eq. nmenus/2+5) then
            l12v(itemp) = ntypes
         endif
      else if (key .eq. MOUSE3) then
         l12v(itemp) = l12v(itemp) + 1
c        Multipath (# 13-16 not implemented yet!)
         if (l12v(itemp) .eq. 13) l12v(itemp) = ntypes
      endif

      if (l12v(itemp) .gt. ntypes) then
         l12v(itemp) = -2
      else if (l12v(itemp) .lt. -2) then
         l12v(itemp) = ntypes
      endif

c     make sure that slots are filled from left to right
      do i = 1,maxwin-1
         if (l12v(i) .eq. 0 .and.  l12v(i+1) .ne. 0) then
            l12v(i) = l12v(i+1)
         endif
      enddo
      if (l12v(maxwin-1) .eq. l12v(maxwin))  then
            l12v(maxwin) = 0
      endif


      redraw = .true.

      return
      end
c************************************************************************
      subroutine ed_ipoly (MOUSE1,MOUSE2,MOUSE3,key,
     .                     ipoly,maxply,newbox,redraw)
      logical    redraw,newbox
      integer    ipoly,maxply
      character*1 MOUSE1,MOUSE2,MOUSE3,key

c     choose the degree of the polynomial
c     -1 is first difference

      if (key .eq. MOUSE1) then
         ipoly = ipoly - 1
      else if (key .eq. MOUSE2) then
         ipoly = 0
      else if (key .eq. MOUSE3) then
         ipoly = ipoly + 1
      else
         ipoly = 0
      endif
      if (ipoly .gt. maxply) then
         ipoly = -1
      else if (ipoly .lt. -2) then
         ipoly = maxply
      endif
      redraw = .true.
      newbox = .false.

      return
      end
c************************************************************************
      subroutine ed_skymap (MOUSE1,MOUSE2,MOUSE3,key,
     .                     igraph,maxgph,newbox,redraw)
      logical    redraw,newbox
      integer    igraph,maxgph
      character*1 MOUSE1,MOUSE2,MOUSE3,key

      if (key .eq. MOUSE1) then
         igraph = igraph - 1
      else if (key .eq. MOUSE2) then
         igraph = 1
      else if (key .eq. MOUSE3) then
         igraph = igraph + 1
      else
         igraph = 1
      endif
      if (igraph .gt. maxgph) then
         igraph = 1
      else if (igraph .lt. 1) then
         igraph = maxgph
      endif
      redraw = .true.
      newbox = .true.

      return
      end
c************************************************************************
      subroutine ed_stack (lstack,redraw,newbox,
     .                            allpts,lspan,nstack)

c     toggle stacking of plots

      logical    lstack,redraw,newbox,allpts,lspan
      integer    nstack

      lstack = (.not. lstack)
      redraw = .true.
      newbox = .false.
      if (lstack) then
         allpts = .false.
         lspan  = .false.
      endif
      nstack = 0

      return
      end
c************************************************************************
      subroutine ed_movie (
     .            key,MOUSE3,lmovie,stopmv,redraw,seekob,newser)

      logical       lmovie,stopmv,redraw,seekob,newser
      character*1   key,MOUSE3

      lmovie = (.not. lmovie)
      redraw = .true.
      if (lmovie) then
         seekob = .true.
         newser = .true.
      endif
      return
      end
c************************************************************************
      subroutine ed_stop (imenu,quitit)

c     tired of editing, leave the menus

      logical    quitit
      integer    imenu

      imenu = 0
      quitit = .true.
      return
      end
c************************************************************************

