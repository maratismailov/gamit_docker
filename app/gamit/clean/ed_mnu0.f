c************************************************************************
      subroutine ed_menus (
     .         lmovie,stopmv,lstack,lspan,lfound,lcolor,quitit,allpts,
     .         seekob,newser,replot,redraw,penplt,newbox,domarg,hide_b,
     .         nbrack,nobs,nsave,nstack,ncfls,nsat,nmlfit,ngood,ntemp,
     .         npatch,nwins,iiend,iibeg,nelim,
     .         isite1,isite2,ichan1,ichan2,
     .         jsite1,jsite2,jchan1,jchan2,
     .         ksite1,ksite2,kchan1,kchan2,
     .         jsite,jchan,iwin,nl,luout,
     .         key,key0,pos,icomm,imenu,imenub,idump,ipoly,igraph,itemp,
     .         ymin2,ymax2,yadd,yend,ybeg,yslp,ymult,yold,yjump0,slip,
     .         nextslip,llist,ilist,iseek,
     .         ibmx0,ibmy0,nx,ny,l12v,l12e,l12p,ii0,ii1,iibl,iibr,iibl0,
     .         mplot,mrect,mtaxi,menus,mclip,csize,mwarn,mhelp,amsg,
     .         igerr,amenus,lbx,rbx,ix1,ix2,iy1,iy2,ix0,i0,dt,xs,ttt,
     .         plt,lslip,nslip,istart,ilast,epsilon,ytot,yflag0,yflag1,
     .         iprob,nprob,pbias,simple,imode,lcando)


      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      logical       lmovie,llmovie,stopmv,lstack,lspan,lfound,simple
      logical       lcolor,allpts,replot,hide_b
      logical       seekob,newser,redraw,penplt,newbox,domarg,quitit
      logical       slip,nextslip,lslip,mbias
      logical       lcando(-3:ntypes)
      integer       npatch,nslip,istart,ilast,iibl0,imode
      integer       llist,ilist,iseek
      integer       i,j,nl,luout,jsite,jchan
      integer       nbrack,nobs,nsave,nelim,nstack,ncfls,iwin
      integer       nsat,nmlfit,ngood,ntemp,nwins,iiend,iibeg
      integer       icomm,imenu,imenub,idump,ipoly,igraph,itemp
      integer       iibl,iibr,i0,ii0,ii1,l12v(maxwin),l12e,l12p
      integer       iprob,nprob,pbias(maxepc),isvmode
      integer       isite1,isite2,ichan1,ichan2
      integer       jsite1,jsite2,jchan1,jchan2
      integer       ksite1,ksite2,kchan1,kchan2
      integer       idd
      integer*2     ibmx0,ibmy0,nx,ny
      integer*2     lbx,rbx,ix1,iy1,ix2,iy2,ix0,pos(2)
      real*8        dt,xs,ttt,plt(maxepc),yflag0,yflag1,yjump0
      real*8        yadd,yend,ybeg,yslp,ymult,yold
      real*8        ymin2(maxwin),ymax2(maxwin),epsilon(4),ytot
      real*8        threshj,threshr
      character*1   key,key0
      integer*2     mplot(maxwin,4),mrect(4),mtaxi(4),menus(nmenus,4)
      integer*2     mwarn(4),mhelp(4),csize(2),mclip(4)
      integer*4     igerr
      character*16  amenus(nmenus)
      character*96  amsg

c     threshold for big jump in LC
      threshj = 0.9 * inter/60.d0
c     threshold for large RMS
      threshr = 1.d0


c**********************************************************************
c**********************************************************************
c
c     Now we are done plotting, so wait for input from mouse,
c     or go on to the next set of plots for a movie.
c
c**********************************************************************
c**********************************************************************


      if (lmovie) then
         call gmouse (key,pos,.false.)
         if (key.eq.MOUSE1 .or.
     .       key.eq.MOUSE2 .or.
     .       key.eq.MOUSE3 .or.
     .       (stopmv .and.
     .       (dabs(yjump0).gt.threshj.or.yflag1.gt.threshr)))
     .      then
             lmovie = .false.
             llmovie = .false.
             seekob = .false.
             newser = .true.
             redraw = .true.
             replot = .true.
            if (stopmv) then
               stopmv = .false.
               if (yflag1 .gt. threshr) then
                  call gmsg (mwarn,'HIGH RMS in LC!')
               else if (dabs(yjump0) .gt. threshj) then
                  call gmsg (mwarn,'BIG JUMP in LC!')
               else if (dabs(yjump0).gt.threshj
     .                  .and.yflag1 .gt.threshr) then
                  call gmsg (mwarn,'BIG JUMP AND HIGH RMS IN LC!')
               endif
c              restore channel and site indexes
c              to their values before after call to GETSER
               isite1 = jsite1
               isite2 = jsite2
               ichan1 = jchan1
               ichan2 = jchan2
            else
c              restore channel and site indexes
c              to their values before call to NEXTOB
               isite1 = ksite1
               isite2 = ksite2
               ichan1 = kchan1
               ichan2 = kchan2
            endif
            key = '!'

         else
            seekob = .true.
            newser = .true.
         endif
         redraw = .true.
      else if (nextslip) then
         key = '!'
         call gmouse (key,pos,.false.)
         if (key.eq.MOUSE1 .or. key.eq.MOUSE2 .or. key.eq.MOUSE3) then
            key     = '!'
            lslip    = .false.
            nextslip = .false.
            replot   = .true.
         else
            lslip    = .true.
         endif
         redraw = .true.
      else
         call gmouse (key,pos,.true.)
      endif

c**********************************************************************
c**********************************************************************
c
c     Decide what to do based on input from mouse.
c
c**********************************************************************
c**********************************************************************

c
      if (key .eq. MOUSE1 .or.
     .    key .eq. MOUSE2 .or.
     .    key .eq. MOUSE3) then


         if ((pos(1) .eq. 0) .and. (pos(2) .eq. 0)) then
c           check to see if we should refresh the screen
c           This works because GMOUSE returns 0,0 if screen has
c           be opened ("de-iconified") or uncovered.
            redraw = .true.
            replot = .true.
            newbox = .true.
         else
c           find out in which box the cursor lies and make it flash
            icomm = 0
            do 3020 i = 1, nmenus
               if (pos(1) .ge. menus(i,1) .and.
     .             pos(1) .le. menus(i,1)+menus(i,3) .and.
     .             pos(2) .ge. menus(i,2) .and.
     .             pos(2) .le. menus(i,2)+menus(i,4)) then
                  do j=1,4
                     mrect(j) = menus(i,j)
                  enddo
                  call gbox (mrect,igerr)
                  icomm = i
                  imenub = imenu
                  call gmsg (mwarn,' ')
               endif
 3020       continue
         endif

c        **************************************************************
         if (imenu .eq. -1) then
c        **************************************************************
c        **** OBSERVATION SELECTION MENU ****  MENU1
c                                   (see file ed_mnu1.ftn)
c        **************************************************************

c           **** SEEK ****
            if (icomm .eq. 1)  then
               call ed_seek (newser,seekob,nobs)

c           **** PLOT ****
            else if (icomm .eq. 2) then
               call ed_plot (newser,nobs)

c           **** LIST ****
            else if (icomm .eq. 3) then
                call ed_list (
     .           ichan1,ichan2,isite1,isite2,nsat,ncfls,
     .           llist,ilist,key,newser,redraw,mwarn,mhelp)

c           **** 1-WAY ***
            else if (icomm .eq. 4) then
                call ed_1way
     .                    (ichan1,ichan2,isite1,isite2,
     .                      mwarn,amsg,newser,redraw)

c           **** SAT1 ***
            else if (icomm .eq. 5) then
                call ed_sat1
     .          (MOUSE1,MOUSE2,MOUSE3,key,ichan1,nsat,redraw)

c           **** SAT2 ***
            else if (icomm .eq. 6) then
                call ed_sat2 (MOUSE1,MOUSE2,MOUSE3,
     .               key,ichan2,nsat,redraw)

c           **** SITE1 ****
            else if (icomm .eq. 7) then
                call ed_site1 (MOUSE1,MOUSE2,MOUSE3,
     .                         key,isite1,ncfls,redraw)

c           **** SITE2 ****
            else if (icomm .eq. 8) then
                call ed_site2 (MOUSE1,MOUSE2,MOUSE3,
     .                         key,isite2,ncfls,redraw)

c           **** FILE ***
            else if (icomm .eq. 9) then
                call ed_file (newser,penplt)

c           **** SEEK CHOICE ****
            else if (icomm .eq. 10) then
                call ed_seekch (MOUSE1,MOUSE2,MOUSE3,
     .                         key,iseek,newbox,redraw)

c           **** PRINT ***
            else if (icomm .eq. 110) then
                call ed_print(mwarn,ibmx0,ibmy0,
     .                        nx,ny,lcolor,idump)

c           **** L12V ****
c           select observables to edit and view
            else if (icomm .ge. nmenus/2+1 .and.
     .               icomm .le. nmenus/2+maxwin) then
                call ed_l12v (key,l12v,icomm,redraw)

c           **** IPOLY ****
c           choose the degree of the polynomial
c           -1 is first difference
            else if (icomm .eq. nmenus/2+maxwin+1) then
               call ed_ipoly (MOUSE1,MOUSE2,MOUSE3,
     .                        key,ipoly,maxply,newbox,redraw)

c           **** SKYMAP/TIME_SERIES ***
            else if (icomm .eq. nmenus/2+maxwin+2) then
               call ed_skymap (
     .         MOUSE1,MOUSE2,MOUSE3,key,igraph,maxgph,newbox,redraw)

c           **** STACK ***
            else if (icomm .eq. nmenus/2+maxwin+3) then
                call ed_stack (lstack,redraw,newbox,allpts,lspan,nstack)

c           **** MOVIE ***
            else if (icomm .eq. nmenus/2+maxwin+4) then
                call ed_movie (
     .          key,MOUSE3,lmovie,stopmv,redraw,seekob,newser)

c           **** STOP ****
            else if (icomm .eq. nmenus) then
                call ed_stop (imenu,quitit)

c           *** MOUSE2 key in data area should SEEK
            else if(pos(1) .ge. mclip(1) .and.
     .              pos(1) .le. mclip(1)+mclip(3) .and.
     .              pos(2) .ge. mclip(2) .and.
     .              pos(2) .le. mclip(2)+mclip(4) .and.
     .              key .eq. MOUSE2) then
                    call ed_seek (newser,seekob,nobs)
            endif

c           buffer channel and site indexes
ckf            jsite1 = isite1
ckf            jsite2 = isite2
ckf            jchan1 = ichan1
ckf            jchan2 = ichan2

c        **************************************************************
         else if (imenu .eq. 0) then
c        **************************************************************
c        **** MENU FOR PLOTTING ***** MENU2
c
c            (see files ed_mnu2.ftn ed_mnu3.ftn ed_save.ftn ed_bias ...)
c        **************************************************************

c           **** SAVE ****
            if ((icomm .eq. 1) .or. (icomm .eq. 15)) then
                  isvmode = 0
                  call ed_save (
     .                 nsave,nobs,imenu,mwarn,isite1,ichan1,
     .                 replot,newbox,redraw,lfound,
     .                 isite1,isite2,jsite1,jsite2,
     .                 ichan1,ichan2,jchan1,jchan2,isvmode)

c           **** ABORT ****
            else if (icomm .eq. 2) then
                  call ed_abort (
     .                 redraw,newbox,replot,
     .                 lfound,imenu,mwarn,
     .                 isite1,isite2,ichan1,ichan2,
     .                 jsite1,jsite2,jchan1,jchan2)

c           **** SPAN ***  look at part of series
            else if (icomm .eq. 3) then
                  call ed_span (
     .                 nbrack,lspan,allpts,redraw,replot,
     .                         nl,ii0,ii1,iibl,nobs,iibr)

c           **** ALL ****  Toggle showing all the data points
            else if (icomm .eq. 4) then
                  call ed_all (allpts,lspan,replot,redraw)

c           **** << T >> **** dilate the time series
            else if (icomm .eq. 5) then
                  call ed_gt (nl,ii1,ii0,replot)

c           **** >> T << ****
            else if (icomm .eq. 6) then
                  call ed_lt (nl,ii1,ii0,replot)

c           **** MARG ON/OFF ***
            else if (icomm .eq. 7) then
                  call ed_marg (domarg,redraw,replot)

c           **** HIDE ON/OFF ***
            else if (icomm .eq. 8) then
                  call ed_hide (hide_b,redraw,replot)

c           **** << **** Move left
            else if (icomm .eq. 9) then
                  call ed_mvleft (key,MOUSE2,nl,inter,ii0,
     .                                             ii1,replot)

c           **** >> ****  move right
            else if (icomm .eq. 10) then
                  call ed_mvright (key,MOUSE2,nl,inter,ii0,
     .                                              ii1,replot)

c           **** SLIP *** (search for a cycle slip)
            else if (icomm .eq. 11) then
                  call ed_slip(
     .                 replot,nextslip,lslip,npatch,nslip,
     .                                       simple,key,key0)

c           **** MOVE *** (set up menus for manual move)
            else if (icomm .eq. 12) then
                  call ed_move(
     .                 nbrack,nobs,iibr,iibl,l12e,imenu,
     .                 redraw,newbox,replot,yadd,mwarn)

c           **** PATCH **** (Automatic move)
            else if (icomm .eq. 13) then
                  idd = ichan1*ichan2*isite1*isite2
                  call ed_ptch(
     .                 nbrack,nobs,iibr,iibl,l12e,ii0,ii1,
     .                 ngood,yend,iiend,ybeg,iibeg,
     .                 yadd,yslp,itemp,isite1,ichan1,idd,
     .                 key,mwarn,mhelp,amsg,replot,domarg,nmlfit)

c           **** FIND *** (find where is the bias or a gap)
            else if (icomm .eq. 14) then
                 call ed_find(
     .                nbrack,ntemp,iibl,itemp,jsite,jchan,
     .                redraw,replot,newser,lfound,mwarn,
     .                amsg,isite1,isite2,ichan1,ichan2)

c           **** UNWEIGHT ***
            else if (icomm .eq. 16) then
                 call ed_unwe (
     .                replot,mwarn,nbrack,nobs,
     .                iibl,iibr,l12e)

c           **** REWEIGHT ****
            else if (icomm .eq. 17) then
                 call ed_rewe (
     .                replot,mwarn,nbrack,nobs,
     .                iibl,iibr,l12e)

c            **** BIAS ***
            else if (icomm .eq. 18) then
                  call ed_bias (
     .                 key,nbrack,nobs,iibr,iibl,ii0,ii1,l12e,
     .                 isite1,isite2,jsite,ichan1,ichan2,jchan,
     .                 ntemp,itemp,mwarn,amsg,replot,domarg)

c           **** UNDO ****
            else if (icomm .eq. 19) then
                 call ed_undo (replot,nobs,nelim,npatch,l12e)

c           **** CANCEL ****
            else if  (icomm .eq. 119) then
                  call ed_cncl (
     .                  redraw,replot,nbrack,nelim,npatch,nsave)

c           **** ELIM *** (remove all bias in span and save)
            else if (icomm .eq. 20) then
                  call ed_elim (
     .                 nelim,nbrack,nobs,imenu,mwarn,l12e,
     .                 iibl,iibr,ii0,ii1,iibeg,iiend,
     .                 replot,newbox,redraw,lfound,
     .                 isite1,isite2,jsite1,jsite2,
     .                 ichan1,ichan2,jchan1,jchan2)

c           **** SAVE1 *** (for automatic patch sequance)
            else if (icomm .eq. 21) then
                  isvmode = 1
                  call ed_save (
     .                 nsave,nobs,imenu,mwarn,jsite,jchan,
     .                 replot,newbox,redraw,lfound,
     .                 isite1,isite2,jsite1,jsite2,
     .                 ichan1,ichan2,jchan1,jchan2,isvmode)

c           **** Click in data area ***
            else if (pos(1) .ge. mclip(1) .and.
     .               pos(1) .le. mclip(1)+mclip(3) .and.
     .               pos(2) .ge. mclip(2) .and.
     .               pos(2) .le. mclip(2)+mclip(4)) then

c                Magic mouse keys only apply to time series
                 if (igraph .ne. 1) then
                    call ed_data(
     .                nbrack,nobs,igraph,mplot,mrect,mtaxi,mwarn,
     .                csize,lbx,rbx,ix1,ix2,iy1,iy2,ix0,iwin,nwins,
     .                pos,iibl,iibr,i0,ii0,ii1,l12p,l12e,l12v,dt,
     .                xs,ttt,yold,plt,amsg,igerr,replot,domarg)

                 else if (key .eq. MOUSE1) then
                    call ed_data(
     .                nbrack,nobs,igraph,mplot,mrect,mtaxi,mwarn,
     .                csize,lbx,rbx,ix1,ix2,iy1,iy2,ix0,iwin,nwins,
     .                pos,iibl,iibr,i0,ii0,ii1,l12p,l12e,l12v,dt,
     .                xs,ttt,yold,plt,amsg,igerr,replot,domarg)

                 else if (key .eq. MOUSE2) then
                    call ed_abort (
     .                 redraw,newbox,replot,
     .                 lfound,imenu,mwarn,
     .                 isite1,isite2,ichan1,ichan2,
     .                 jsite1,jsite2,jchan1,jchan2)

c                **************************************************
c                *** Right button to start  a sequence of functions

                 else if (key .eq. MOUSE3) then
                    call ed_data(
     .                nbrack,nobs,igraph,mplot,mrect,mtaxi,mwarn,
     .                csize,lbx,rbx,ix1,ix2,iy1,iy2,ix0,iwin,nwins,
     .                pos,iibl,iibr,i0,ii0,ii1,l12p,l12e,l12v,dt,
     .                xs,ttt,yold,plt,amsg,igerr,replot,domarg)

                    call ed_find1(
     .                  nbrack,ntemp,iibl,itemp,jsite,jchan,
     .                  redraw,replot,newser,lfound,mwarn,
     .                  npatch,amsg,isite1,isite2,ichan1,ichan2,
     .                  iprob,pbias,mbias)

                    if (lfound) then
                         idd = ichan1*ichan2*isite1*isite2
                         call ed_ptch(
     .                       nbrack,nobs,iibr,iibl,l12e,ii0,ii1,
     .                       ngood,yend,iiend,ybeg,iibeg,
     .                       yadd,yslp,itemp,jsite,jchan,idd,
     .                       key,mwarn,mhelp,amsg,replot,domarg,nmlfit)
                         imenu = 11
                         replot = .true.
                         write (amsg,'(a)')
     .                   'confirm PATCH: L-undo; R-save'
                         call gmsg(mwarn,amsg)
                    else if (mbias) then
                         call gmsg (mwarn,'More than 1 bias')
                         call multibias (
     .                       lfound,seekob,newser,replot,redraw,
     .                       newbox,domarg,lslip,nslip,nbrack,nobs,
     .                       ntemp,npatch,nwins,iiend,iibeg,itemp,
     .                       jsite1,jsite2,jchan1,jchan2,
     .                       isite1,isite2,ichan1,ichan2,
     .                       jsite,jchan,key,pos,l12e,l12v,l12p,
     .                       ngood,nsave,nmlfit,ymin2,ymax2,yadd,
     .                       yend,ybeg,yslp,ymult,yold,slip,nextslip,
     .                       ii0,ii1,iibl,iibr,iibl0,epsilon,istart,
     .                       ilast,imenu,mplot,mrect,mtaxi,menus,mclip,
     .                       csize,mwarn,amsg,igerr,amenus,lbx,rbx,
     .                       ix1,ix2,iy1,iy2,ix0,i0,dt,xs,ttt,ytot,
     .                       iprob,nprob,pbias,simple,nsat,ncfls,imode)
                         idd = ichan1*ichan2*isite1*isite2
                         call ed_ptch(
     .                       nbrack,nobs,iibr,iibl,l12e,ii0,ii1,
     .                       ngood,yend,iiend,ybeg,iibeg,
     .                       yadd,yslp,itemp,jsite,jchan,idd,
     .                       key,mwarn,mhelp,amsg,replot,domarg,nmlfit)
                         imenu = 11
                         replot = .true.
                         write (amsg,'(a)')
     .                   'confirm PATCH: L-undo; R-save'
                         call gmsg(mwarn,amsg)
                    else
                         call gmsg (mwarn,'No bias or gap')
                    endif
                 endif
            endif
c        **************************************************************
c        *** FAKE MENU LEVEL 11 to confirm automatic patch
c        **************************************************************
         else if (imenu .eq. 11) then
c            **** ABORT ****
             if (key .eq. MOUSE2 .or. icomm .eq. 2) then
                  call ed_abort (
     .                 redraw,newbox,replot,
     .                 lfound,imenu,mwarn,
     .                 isite1,isite2,ichan1,ichan2,
     .                 jsite1,jsite2,jchan1,jchan2)
c            **** UNDO ****
             else if (key .eq. MOUSE1 .or. icomm .eq. 18) then
                  call ed_undo (replot,nobs,nelim,npatch,l12e)
                  call gmsg (mwarn,'UNDID patch')
                  imenu = 0
c            **** SAVE ****
             else if (key .eq. MOUSE3 .or. icomm .eq. 1) then
                  isvmode = 1
                  call ed_save (
     .                 nsave,nobs,imenu,mwarn,jsite,jchan,
     .                 replot,newbox,redraw,lfound,
     .                 isite1,isite2,jsite1,jsite2,
     .                 ichan1,ichan2,jchan1,jchan2,isvmode)
                  call ed_plot (newser,nobs)
                  call gmsg (mwarn,'SAVED patch')
                  imenu = 0
             endif

c        **************************************************************
c        *** MENU LEVEL 14 (set up manual move)
c                                   (see file ed_move.ftn)
c        **************************************************************
         else if (imenu .eq. 12) then
                call ed_menu12(
     .               nobs,iibl,iibr,l12e,icomm,itemp,isite1,ichan1,
     .               imenu,yadd,ymult,key,replot,redraw,newbox,m
     .               warn,amsg)

         endif

c        **************************************************************
c        **** end of command switches ****
c        **************************************************************

c        put the cursor back to avoid confusing user
         call gmvcur (pos, igerr)

c        redraw command box
         if (icomm .ne. 0) then
            do 4300 i=1,4
               mrect(i) = menus(icomm,i)
 4300       continue
            call gmsg (mrect,amenus(icomm))
            icomm = 0
         endif

      endif

c     **************************************************************
c     *** Cycle-slip finder
c                           (see file ed_slip.ftn)
c     **************************************************************
      if (lslip) then

          key = key0
          call autoslip(
     .         lfound,seekob,newser,replot,redraw,newbox,domarg,lslip,
     .         nslip,nbrack,nobs,ntemp,npatch,nwins,iiend,iibeg,itemp,
     .         jsite1,jsite2,jchan1,jchan2,isite1,isite2,ichan1,ichan2,
     .         jsite,jchan,key,pos,l12e,l12v,l12p,ngood,nsave,nmlfit,
     .         ymin2,ymax2,yadd,yend,ybeg,yslp,ymult,yold,slip,nextslip,
     .         ii0,ii1,iibl,iibr,iibl0,epsilon,istart,ilast,imenu,
     .         mplot,mrect,mtaxi,menus,mclip,csize,mwarn,mhelp,
     .         amsg,igerr,
     .         amenus,lbx,rbx,ix1,ix2,iy1,iy2,ix0,i0,dt,xs,ttt,ytot,
     .         iprob,nprob,pbias,simple,nsat,ncfls,imode)


c         if (seekob) write(*,*)
c     .     'Done with ',isite2,isite1,ichan2,ichan1

      endif

      return
      end
c


c**********************************************************************

