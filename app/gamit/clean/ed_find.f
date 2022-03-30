c************************************************************************
c           **** FIND ***
      subroutine ed_find (
     .       nbrack,ntemp,iibl,itemp,jsite,jchan,
     .       redraw,replot,newser,lfound,mwarn,
     .       amsg,isite1,isite2,ichan1,ichan2)


      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         nbrack,iibl,ntemp,itemp,jsite,jchan
      integer         isite1,ichan1,isite2,ichan2
      logical         replot,redraw,newser,lfound
      integer*2       mwarn(4)
      character*96    amsg

c     find the one-way in which a missing point occurs
      if (nbrack .eq. 1) then

          lfound = .false.

c         check number of bias at iibl
          ntemp = 0
          call check_bias(isite1,isite2,ichan1,ichan2,
     .                     jsite,jchan,ntemp,iibl)
c         needs only one bias
          if (ntemp .eq. 1) then
               write (amsg,3061) jchan,jsite
 3061          format ('Bias in chan ',i2,' site ',i2)
               call gmsg (mwarn,amsg)
               lfound = .true.

           else if (ntemp .eq. 0) then
c              if no bias check number of gap at iibl
               ntemp = 0
               call check_gap(isite1,isite2,ichan1,ichan2,
     .                     jsite,jchan,ntemp,iibl)
c              needs only one gap
               if (ntemp .eq. 1) then
                    write (amsg,3062) jchan,jsite
 3062               format ('Gap in chan ',i2,' site ',i2)
                    call gmsg (mwarn,amsg)
                    lfound = .true.
                else if (ntemp .gt. 1) then
                    call gmsg (mwarn,'More than 1 gap.')
                else if (ntemp .eq. 0) then
                    call gmsg (mwarn,'No bias or gap here.')
                endif

           else if (ntemp .gt. 1) then
c              if more then 1 bias check number of (bias .and. gap) at iibl
               ntemp = 0
               call check_bgap(isite1,isite2,ichan1,ichan2,
     .                     jsite,jchan,ntemp,iibl)
c              needs only one bias-gap
               if (ntemp .eq. 1) then
                    write (amsg,3063) jchan,jsite
 3063               format ('Bias-gap in chan ',i2,' site ',i2)
                    call gmsg (mwarn,amsg)
                    lfound = .true.
                else if (ntemp .gt. 1) then
                    call gmsg (mwarn,'More than 1 bias-gap.')
                else if (ntemp .eq. 0) then
                    call gmsg (mwarn,'more than 1 bias.')
                endif
           endif

      else if (nbrack .le. 0) then
           call gmsg (mwarn,'FIND needs 1 bracket')
           lfound = .false.
      else if (nbrack .gt. 1) then
           call gmsg (mwarn,'Too many brackets for FIND')
           lfound = .false.
      endif
      if(lfound) then
           if (ichan1 .ne. jchan) then
                itemp  = ichan1
                ichan1 = ichan2
                ichan2 = itemp
           endif
           if (isite1 .ne. jsite) then
                itemp  = isite1
                isite1 = isite2
                isite2 = itemp
           endif
           newser = .true.
           redraw = .true.
      endif
      replot = .true.
      nbrack = 0


      return
      end
c************************************************************************
c           **** FIND ***
      subroutine ed_find1 (
     .       nbrack,ntemp,iibl,itemp,jsite,jchan,
     .       redraw,replot,newser,lfound,mwarn,
     .       npatch,amsg,isite1,isite2,ichan1,ichan2,
     .       nprob,pbias,mbias)


      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         nbrack,iibl,ntemp,itemp,jsite,jchan
      integer         isite1,ichan1,isite2,ichan2,npatch
      integer         nprob,pbias(maxepc)
      logical         replot,redraw,newser,lfound
      logical         mbias
      integer*2       mwarn(4)
      character*96    amsg

          mbias  = .false.
          lfound = .false.

c         check number of bias at iibl
          ntemp = 0
          call check_bias(isite1,isite2,ichan1,ichan2,
     .                     jsite,jchan,ntemp,iibl)
c         needs only one bias
          if (ntemp .eq. 1) then
               write (amsg,3061) jchan,jsite
 3061          format ('Bias in chan ',i2,' site ',i2)
               call gmsg (mwarn,amsg)
               lfound = .true.

           else if (ntemp .eq. 0) then
c              if no bias check number of gap at iibl
               ntemp = 0
               call check_gap(isite1,isite2,ichan1,ichan2,
     .                     jsite,jchan,ntemp,iibl)
c              needs only one gap
               if (ntemp .eq. 1) then
                    write (amsg,3062) jchan,jsite
 3062               format ('Gap in chan ',i2,' site ',i2)
                    call gmsg (mwarn,amsg)
                    lfound = .true.
                else if (ntemp .gt. 1) then
                    call gmsg (mwarn,'More than 1 gap.')
                else if (ntemp .eq. 0) then
                    call gmsg (mwarn,'No bias or gap here.')
                endif

           else if (ntemp .gt. 1) then
c              if more then 1 bias check number of (bias .and. gap) at iibl
               ntemp = 0
               call check_bgap(isite1,isite2,ichan1,ichan2,
     .                     jsite,jchan,ntemp,iibl)
c              needs only one bias-gap
               if (ntemp .eq. 1) then
                    write (amsg,3063) jchan,jsite
 3063               format ('Bias-gap in chan ',i2,' site ',i2)
                    call gmsg (mwarn,amsg)
                    lfound = .true.
                else if (ntemp .gt. 1) then
                    call gmsg (mwarn,'More than 1 bias-gap.')
                    mbias = .true.
                else if (ntemp .eq. 0) then
                    call gmsg (mwarn,'more than 1 bias.')
                    mbias = .true.
                endif
           endif


      return
      end
c************************************************************************
c           **** FIND ***
      subroutine ed_find2 (
     .       iibl,ntemp,s1,s2,c1,c2,mwarn,
     .       isite1,isite2,ichan1,ichan2)


      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         iibl,ntemp
      integer         isite1,ichan1,isite2,ichan2
      logical         lbias,s1,s2,c1,c2
      integer*2       mwarn(4)
      character*96    amsg


      ntemp = 0
      s1 = .false.
      s2 = .false.
      c1 = .false.
      c2 = .false.

      if (isite1 .gt. 0 .and. ichan1 .gt. 0) then
          if (lbias(ierr(iibl,ichan1,isite1))) then
                 s1 = .true.
                 c1 = .true.
                 ntemp = ntemp + 1
          endif
      endif
      if (isite1 .gt. 0 .and. ichan2 .gt. 0) then
          if (lbias(ierr(iibl,ichan2,isite1))) then
                 s1 = .true.
                 c2 = .true.
                 ntemp = ntemp + 1
          endif
      endif
      if (isite2 .gt. 0 .and. ichan1 .gt. 0) then
          if (lbias(ierr(iibl,ichan1,isite2))) then
                 s2 = .true.
                 c1 = .true.
                 ntemp = ntemp + 1
          endif
      endif
      if (isite2 .gt. 0 .and. ichan2 .gt. 0) then
          if (lbias(ierr(iibl,ichan2,isite2))) then
                 s2 = .true.
                 c2 = .true.
                 ntemp = ntemp + 1
          endif
      endif

       if (ntemp .eq. 2) then
             if (s1 .and. .not. s2) then
                  write (amsg,'(''at '',2i4," bias at s1   c1+c2")')
     .                  iibl,ntemp
             else if (s2 .and. .not. s1) then
                  write (amsg,'(''at '',2i4," bias at s2   c1+c2")')
     .                  iibl,ntemp
             else if (c1 .and. .not. c2) then
                  write (amsg,'(''at '',2i4," bias at s1+s2 -  c1")')
     .                  iibl,ntemp
             else if (c2 .and. .not. c1) then
                  write (amsg,'(''at '',2i4," bias at s1+s2 -  c2")')
     .                  iibl,ntemp
             else if (c2 .and. .not. c1) then
                  write (amsg,'(''at '',2i4," bias at s1+s2 -  c2")')
     .                  iibl,ntemp
             else
                  write (amsg,'(''at '',2i4," bias at s1,s2,c1,c2")')
     .                  iibl,ntemp
             endif

       else if (ntemp .gt. 2) then
             write (amsg,'(''at '',2i4," bias")')iibl,ntemp

       else
             write (amsg,'(''at '',2i4," bias ?????????????")')
     .              iibl,ntemp

       endif

       call gmsg (mwarn,amsg)


      return
      end
c************************************************************************
      subroutine check_bias (isite1,isite2,ichan1,ichan2,
     .                     jsite,jchan,ntemp,iibl)


      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         iibl,ntemp,jsite,jchan
      integer         isite1,ichan1,isite2,ichan2
      logical         lbias

          if (isite1 .gt. 0 .and. ichan1 .gt. 0) then
              if (lbias(ierr(iibl,ichan1,isite1))) then
                     jchan = ichan1
                     jsite = isite1
                     ntemp = ntemp + 1
               endif
          endif
          if (isite1 .gt. 0 .and. ichan2 .gt. 0) then
              if (lbias(ierr(iibl,ichan2,isite1))) then
                     jchan = ichan2
                     jsite = isite1
                     ntemp = ntemp + 1
               endif
          endif
          if (isite2 .gt. 0 .and. ichan1 .gt. 0) then
              if (lbias(ierr(iibl,ichan1,isite2))) then
                     jchan = ichan1
                     jsite = isite2
                     ntemp = ntemp + 1
               endif
          endif
          if (isite2 .gt. 0 .and. ichan2 .gt. 0) then
              if (lbias(ierr(iibl,ichan2,isite2))) then
                     jchan = ichan2
                     jsite = isite2
                     ntemp = ntemp + 1
               endif
          endif

      return
      end
c************************************************************************
      subroutine check_gap (isite1,isite2,ichan1,ichan2,
     .                     jsite,jchan,ntemp,iibl)


      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         iibl,ntemp,jsite,jchan,igap
      integer         isite1,ichan1,isite2,ichan2
      logical         lgood


          igap = iibl - 1

          if (isite1 .gt. 0 .and. ichan1 .gt. 0) then
              if (.not.lgood(ierr(igap,ichan1,isite1))) then
                     jchan = ichan1
                     jsite = isite1
                     ntemp = ntemp + 1
               endif
          endif
          if (isite1 .gt. 0 .and. ichan2 .gt. 0) then
              if (.not.lgood(ierr(igap,ichan2,isite1))) then
                     jchan = ichan2
                     jsite = isite1
                     ntemp = ntemp + 1
               endif
          endif
          if (isite2 .gt. 0 .and. ichan1 .gt. 0) then
              if (.not.lgood(ierr(igap,ichan1,isite2))) then
                     jchan = ichan1
                     jsite = isite2
                     ntemp = ntemp + 1
               endif
          endif
          if (isite2 .gt. 0 .and. ichan2 .gt. 0) then
              if (.not.lgood(ierr(igap,ichan2,isite2))) then
                     jchan = ichan2
                     jsite = isite2
                     ntemp = ntemp + 1
               endif
          endif

      return
      end
c************************************************************************
      subroutine check_bgap (isite1,isite2,ichan1,ichan2,
     .                     jsite,jchan,ntemp,iibl)


      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         iibl,ntemp,jsite,jchan,igap
      integer         isite1,ichan1,isite2,ichan2
      logical         lgood,lbias


          igap = iibl - 1

          if (isite1 .gt. 0 .and. ichan1 .gt. 0) then
              if (.not.lgood(ierr(igap,ichan1,isite1)) .and.
     .                  lbias(ierr(iibl,ichan1,isite1))) then
                     jchan = ichan1
                     jsite = isite1
                     ntemp = ntemp + 1
               endif
          endif
          if (isite1 .gt. 0 .and. ichan2 .gt. 0) then
              if (.not.lgood(ierr(igap,ichan2,isite1)) .and.
     .                  lbias(ierr(iibl,ichan2,isite1))) then
                     jchan = ichan2
                     jsite = isite1
                     ntemp = ntemp + 1
               endif
          endif
          if (isite2 .gt. 0 .and. ichan1 .gt. 0) then
              if (.not.lgood(ierr(igap,ichan1,isite2)) .and.
     .                  lbias(ierr(iibl,ichan1,isite2))) then
                     jchan = ichan1
                     jsite = isite2
                     ntemp = ntemp + 1
               endif
          endif
          if (isite2 .gt. 0 .and. ichan2 .gt. 0) then
              if (.not.lgood(ierr(igap,ichan2,isite2)) .and.
     .                  lbias(ierr(iibl,ichan2,isite2))) then
                     jchan = ichan2
                     jsite = isite2
                     ntemp = ntemp + 1
               endif
          endif

      return
      end
c************************************************************************

