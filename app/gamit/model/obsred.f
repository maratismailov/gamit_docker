Copyright (c) Massachusetts Institute of Technolog and the University of
California at San Diego, 1994.  All rights reserved.

      Subroutine obsred

c     Read one record from the input X- or C- file or get the epoch for simulated observations

c      Input
c            iobs  :  unit number of obs file (now in common /units/ of model.h
c            nchan :  number of channels on file  (model.h)
c            ischan:  PRN numbers of channels (dimension MAXSAT) (model.h)
c            ndat  :  number of data types (L1, L2, P1, P2, etc)  (model.h)
c            jd0   :  Julian Day of start epoch (GPST)  (model.h)
c            t0    :  seconds-of-day of start epoch (GPST) (model.h)
c            inter :  observation file interval (integer seconds) (model.h)
c            mtime :  type of time on obs file (1=UTC, 2=GPST)  (model.h)
c            iepoch:  epoch number expected to be read (check for consistency)
c       Output
c            msat  :  number of satellites this record  (model.h)
c            mprn  :  PRN numbers of satellite this record (dimension MAXSAT) (model.h)
c            jdobs :  Julian Day of this record (GPST) (model.h)
c            tobs  :  seconds-of-day this record (GPST) (model.h)
c            ier   :  data flags for channels of this record (dimension MAXSAT) (model.h)
c            isnr  :  signal-to-noise for each data type each channel (MAXDATxMAXSAT) (model.h)
c            obs   :  observations (MAXDATxMAXSAT)  (model.h)

C     Mod: MHM 870731  Error flagging from x-files
c     Mod: MHM 880330  Error flagging:  include low amplitude flags=3
c                         reset 98 and 4 to zero
c     Mod: RWK 881010  Error flagging:  apply different criteria to
c                         low amplitude acceptance based on receiver
c     Mod: RWK 891226  Return time of observation, either read from
c                         file or computed from start + increment
c     Mod: RWK 940712  Convert UTC to GPST

      implicit none

      include '../includes/dimpar.h'
      include '../includes/model.h'

      character*1  upperc,latflag,lonflag      
      character*256 message

      INTEGER*4 ndats,jsat,idoy,julday,nparts
     .        , iep,ihr,min,iyr,id,im,iprn,ierx
     .        , jdoy,jdread,i,j,m,is,kdat
     .        , iepochx,dlat,dlon,mlat,mlon,idcb,ioerr
c      this not used so removed from the x-file eventually
     .        , kflag 


      real*8  delt,sec,utcoff,taiutc
     .      , pres1,temp1,relhum1
     .      , tau(maxsat),tread,antaz,svcepc,svcl1,slat,slon
     .      , l1z,l1n,l1e,l2z,l2n,l2e

      character*128 xobs_line

C     Initialize variables
      
      jdoy= 0        
      svcepc = 0.d0
      svcl1 = 0.d0
      do j = 1, maxsat
         ier(j) = 1
         data_flag(j) = 0
         elev(j) = 0.d0
         elvdot(j) = 0.d0   
         atmdel(j) = 0.d0
         ampl1(j) = 0.d0
         ampl2(j) = 0.d0   
         do i = 1,maxdat
            isnr(i,j) = 0 
            obs(i,j) = 0.d0         
         enddo               
      enddo
                   
C     Read the Data from an X-File---

      if (upperc(obfiln(1:1)).eq.upperc('X')) then

        read ( iobs,'(/,2I4,I5,I4,2I3,F11.7)',iostat=ioerr)
     .         iepochx,msat,iyr,jdoy,ihr,min,sec  
        if( ioerr.gt.0 ) 
     .    call report_stat('WARNING','MODEL','obsred',' '
     .                 ,'Error reading X-file time line',ioerr)
        if( msat.gt.0 .or.iyr.ne.0 ) call fix_y2k(iyr)

c       Compare the epoch numbers input and read to detect a bogus x-file
        if( iepochx.ne.iepoch ) then  
           write(message,'(a,a4,a,i5,a,i5,a)') 
     .       'Epoch on X-file for site',sitecd,' (',iepochx
     .       ,') not expected (',iepoch,') --file corrupted?'
           call report_stat('FATAL','MODEL','obsred', ' ',message,0)
        endif

c       Compute the observation time from number of epochs
        delt= (iepoch-1)*dble(inter)
        jdobs= jd0
        tobs = t0
        call timinc (jdobs,tobs,delt)

c       Check this against the actual time if available (new-style X-files)
        if( jdoy.gt.0 ) then
         call monday(jdoy,im,id,iyr)
          jdread= julday(im,id,iyr)
          tread = ihr*3600.d0 + min*60.d0 + sec

c         if X-file is in UTC, convert to GPST
          if( mtime.eq.1 ) then
            utcoff = taiutc(jdread) - 19.d0
            call timinc(jdread,tread,utcoff)
          endif

C Assume the makex always gets the epoch correct. Skip this code since if the clock
C is drifting by mode than 1 sec/day it causes an unnecessary fatal stop in model. SCM 070402 !
Cc         let the time tag be up to 0.3 second off
Cc         currently, the only receiver having this problem is
Cc         the Trimble at SIG 3.25, which can drift up to 256 millisec.
Cc** rwk   some Trimble clocks are so bad that they can be 1 sec off!
Cc         For now allow greater slop.  In the long run, test on the
Cc         actual time (i.e., receiver time corrected by the known clock offset)
C          if (dabs((jdobs-jdread)*86400.d0+(tobs-tread)) .gt. 1.d0) then
C             write(message,65) sitecd,jdread,tread,jdobs,tobs,iepoch
C  65         format('Epoch mismatch: Time on input file for site ',a4
C     .        ,' (',I8,F10.2,') ','does not match expected time (',I8
C     .        , f10.2,') at epoch',I5)
C             call report_stat('FATAL','MODEL','obsred',' ',message,0)
C          endif

          jdobs= jdread
          tobs = tread
        endif

C        Now read the phase, amplitudes, and pseudo-range
         do m=1,msat
           read ( iobs,'(a)',iostat=ioerr ) xobs_line
           if( ioerr.eq.0 ) 
     .     read ( xobs_line,iostat=ioerr 
     .,           fmt ='(8x,i1,1x,2I2,2(1X,D22.15,1X,I3),2(2X,D22.15))')
     .            idcb,ierx,jsat
     .,           (obs(i,jsat),isnr(i,jsat),i=1,2)
     .,           (obs(i,jsat),i=3,ndat) 
cd              print *,'OBSRED ',iepochx,m,idcb,ierx,jsat
cd     .             ,(obs(i,jsat),i=1,ndat)
           if( ioerr.gt.0 ) then
              call report_stat('WARNING','MODEL','obsred',' '
     .                 ,'Error reading X-file data line',ioerr)
              ierx = 2   ! Mark as deleted observation
           endif
            
           ier(jsat)=ierx
           mprn(m) = ischan(jsat)
c          dummy the amplitudes for now
           ampl1(jsat) = 0.d0
           ampl2(jsat) = 0.d0  
c          set the data_flag bits for C1/P2 corrections (see obsmod)
           if( idcb.gt.0 ) call sbit(data_flag(jsat),28,1)
           if( idcb.eq.2 ) call sbit(data_flag(jsat),29,1)
             
c          convention for phase with increasing doppler:  Gamit currently
c          assumes negative, although RINEX assumes positive
           do i=1,2
               obs(i,jsat) = -obs(i,jsat)
            enddo  
         enddo   
                         
c      Read the Data from a C-File
      elseif( upperc(obfiln(1:1)).eq.upperc('C') ) then

         call readc4 (iobs
     .,               msat,mprn
     .,               iep,iyr,idoy,tobs
     .,               rclock,okmet,zendel,pres1,temp1,relhum1,atmlod
     .,               sitecd,latr_sph,lonr,radius
     .,               l1z,l1n,l1e,l2z,l2n,l2e,antaz
     .,               nsave,save )  
         call check_y2k(iyr)   
         call monday(idoy,im,id,iyr) 
         jdobs = julday(im,id,iyr)
      

c        convert from UTC to GPST if necessary
         if( mtime.eq.1 ) then
           utcoff = taiutc(jdobs) - 19.d0
           call timinc( jdobs,tobs,utcoff )
         endif
         if (iep .ne. iepoch) then
            write(message,'(a,2i6)')
     .      'C-file epoch counter mismatch.',iep,iepoch
            call report_stat('FATAL','MODEL','obsred',' ',message,0)
         endif

c        read the data, one sat at a time, put into ischan order
c        note that only ier,obs,isnr are passed back

         do 200 is = 1, msat
            i = 1
  210       if (mprn(is).eq.ischan(i)) then
               call readc5 (iobs
     .,                     iprn
     .,                     elev(i),elvdot(i)
     .,                     azim(i),azmdot(i),nadang(i)
     .,                     atmdel(i)
     .,                     svcepc,svcl1
     .,                     tau(i),drate(i)
     .,                     ier(i),data_flag(i)
     .,                     ndats
     .,                     obs(1,i)
     .,                     omc(1,i)
     .,                     isnr(1,i)
     .,                     ampl1(i),ampl2(i)
     .,                     nspare
     .,                     spare
     .,                     nparts
     .,                     tmpart(1,i))
c
c              if (ndats.ne.ndat .or. iprn.ne.ischan(i)) then
c                 write(message,'(a,4i5)')
c    .               'ndat or iprn missmatch ',
c    .                ndats,ndat,iprn,ischan(i)
c         call report_stat('FATAL','MODEL','obsred',' ',message,0)
c              endif
            else
               i=i+1
c              if (i.gt.nchan) then
c                 write(message,'(a,2i5)')
c    .               'msat missmatch ',mprn(is)
c         call report_stat('FATAL','MODEL','obsred',' ',message,0)
c              endif
               goto 210
            endif
  200    continue

      elseif(upperc(obfiln(1:1)).eq.upperc('S')) then

c       Compute the observation time from number of epochs
        delt= (iepoch-1)*dble(inter)
        jdobs= jd0
        tobs = t0
        call timinc (jdobs,tobs,delt)
        msat = nchan
        do i=1,msat
          mprn(i) = ischan(i)   
          ier(i) = 0 
          ampl1(i) = 1.0
          ampl2(i) = 1.0
        enddo    
c 

      else
         call report_stat('FATAL','MODEL','obsred',obfiln(1:1),
     .   'Input observation file type not X, C, or S ',0)
      endif

c**  With the move of the elevation cutoff code from MODEL to SOLVE
c    (also AUTCLN), there is no flagging of low-el points in MODEL
c    (code commented in model.f).  To preserve the fragile editing
c    of old rules, saved in 'clean' x-files, leave the low-el codes
c    as they were.   If the data are reedited by AUTCLN, the low-el
c    points might be recovered.   Removing this reset is particularly
c    useful for handling the problem of earlier data being cut using
c    an incorrect elevation angle, fixed in MODEL version 9.38 (28 Jun 95).
c    Also remove the reset of the 98 flag--no longer used.
c
c*   Change 98 flags to 0's (good observations)
C      A flag of 98 is temporarily assigned in CE (C-file editor)
C      at a possible cycle-slip to create a gap in LSQ

C    Change 4 flags to 0's as well.
C      A flag of 4 indicates an epoch limit or low elevation
C      These may have changed, and are reevaluated in the
C      the main program
c*     do i=1,nchan
c*       if( ier(i).eq.98 .or. ier(i).eq.4 ) ier(i)=0
c*     enddo
      
      return
      end
