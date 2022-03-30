Copyright (c) Massachusetts Institute of Technolog and the University of
California at San Diego, 1994.  All rights reserved.

      Subroutine obsred( obsfil,iobs,sitecd,xtype,skd
     .                 , nchan,ischan,ndat
     .                 , jd0,t0,inter,mtime
     .                 , iepoch,msat,mprn,jdobs,tobs,delt
     .                 , ier,data_flag,isnr,obs,obswgt
     .                 , ampl1,ampl2 )

c     Read one record from the input X- or C- file or get the epoch for simulated observations

c      Input
c            obsfil:  file type (X or C)
c            iobs  :  unit number of obs file
c            sitecd:  4-character site code (for error messages)
c            xtype :  X-file type: old (O) or extended (E)                     
c            skd   :  survey mode: static (S), kinematic (K), or dynamic (D)
c            nchan :  number of channels on file
c            ischan:  PRN numbers of channels (dimension MAXSAT)
c            ndat  :  number of data types (L1, L2, P1, P2, etc)
c            jd0   :  Julian Day of start epoch (GPST)
c            t0    :  seconds-of-day of start epoch (GPST)
c            inter :  observation file interval (integer seconds)
c            mtime :  type of time on obs file (1=UTC, 2=GPST) 
c            iepoch:  epoch number expected to be read (check for consistency)
c       Output
c            msat  :  number of satellites this record
c            mprn  :  PRN numbers of satellite this record (dimension MAXSAT)
c            jdobs :  Julian Day of this record (GPST)
c            tobs  :  seconds-of-day this record (GPST)
c            delt  :  time since first epoch (seconds)
c            ier   :  data flags for channels of this record (dimension MAXSAT)
c            isnr  :  signal-to-noise for each data type each channel (MAXDATxMAXSAT)
c            obs   :  observations (MAXDATxMAXSAT)
c            obswgt:  weights of observations (always 1.) (MAXDATxMAXSAT)

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

      character*1  obsfil,upperc,xtype,skd,latflag,lonflag      
      character*4  sitecd
      character*256 message

      INTEGER*4 iobs,jdobs,ndat,ndats,jsat,julday
     .        , inter,isnr(maxdat,maxsat),msat
     .        , mprn(maxsat),nchan,ischan(maxsat),ier(maxsat),jd0,nparts
     .        , iepoch,iep,nsave,nspare,ihr,min,iyr,id,im,iprn,ierx
     .        , jdoy,jdread,okmet,okwvr(maxsat),mtime,i,j,m,is,kdat
     .        , data_flag(maxsat),iepochx,dlat,dlon,mlat,mlon,idcb,ioerr

      real*4 ampl1(maxsat),ampl2(maxsat)

      real*8  tobs,delt,t0,sec,sod,rclock,utcoff,taiutc
     .      , azim(maxsat),azmdot(maxsat),elev(maxsat),elvdot(maxsat)
     .      , pres1,temp1,relhum1,wvrdel(maxsat)
     .      , atmdel(maxsat),tau(maxsat),drate(maxsat),tread
     .      , obs(maxdat,maxsat),omc(maxdat,maxsat)
     .      , obswgt(maxdat,maxsat)
     .      , tmpart(maxlab,maxsat),save(maxsav),spare(maxspr),antaz
     .      , svcepc,svcl1,latr,lonr,radius,slat,slon

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
         okwvr(j) = 0
         wvrdel(j) = 0.d0 
         atmdel(j) = 0.d0
         ampl1(j) = 0.d0
         ampl2(j) = 0.d0   
         do i = 1,maxdat
            isnr(i,j) = 0 
            obs(i,j) = 0.d0         
            obswgt(i,j) = 0.d0      
         enddo               
      enddo
                    

C     Read the Data from an X-File---

      if (upperc(obsfil).eq.upperc('X')) then

c       First read the time and (optionally) coordinates and antenna height
                                    
        if( upperc(xtype).eq.'E' .and.
     .     (upperc(skd).eq.'K' .or. upperc(skd).eq.'D') ) then
           read ( iobs, 10,iostat=ioerr,err=50)
     .          iepochx,msat,iyr,jdoy,ihr,min,sec
     .     , kflag,ksite,latflag,dlat,mlat,slat
     .     , lonflag,dlon,mlon,slon,radius
     .     , l1z,l1n,l1e,l2z,l2n,l2e
           if( msat.gt.0 .or.iyr.ne.0 ) call fix_y2k(iyr)
   10      format(/,2I4,I5,I4,2I3,F11.7,1x,i2,1x,a4,1x,
     .          a1,i2,1x,i2,1x,f8.5,1x,a1,i3,1x,i2,1x,f8.5,f13.4,
     .          1x,3F8.4,3X,3F8.4)
        elseif( upperc(xtype).eq.'O'.or.
     .          upperc(skd).eq.'S') then       
           read ( iobs, 20,iostat=ioerr,err=50 )
     .          iepochx,msat,iyr,jdoy,ihr,min,sec 
           if( msat.gt.0 .or.iyr.ne.0 ) call fix_y2k(iyr)
   20      format(/,2I4,I5,I4,2I3,F11.7)
        else
          write(message,30) sitecd,xtype,skd  
   30     format('Incompatible X-file type and survey mode for site '
     .     , a4,' X-file type: ',a2,' Survey mode: ',a2)
          call report_stat('FATAL','MODEL','obsred',' ',message,0)
        endif 

   50   if( ioerr.gt.0 ) then
          call report_stat('WARNING','MODEL','obsred',' '
     .                 ,'Error reading X-file time line',ioerr)
        endif
 

c       Compare the epoch numbers input and read to detect a bogus x-file
        if( iepochx.ne.iepoch ) then
           write(message,60) sitecd,iepochx,iepoch
   60      format('Epoch on X-file for site ',a4,' (',i5
     .         ,') not expected (',i5,' )--file corrupted?')
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
c          print *,iepochx,m,idcb,ierx,jsat ,obs(1,jsat),obs(2,jsat)
           if( ioerr.gt.0 ) then
              call report_stat('WARNING','MODEL','obsred',' '
     .                 ,'Error reading X-file data line',ioerr)
              ierx = 2   ! Mark as deleted observation
           endif
            
           ier(jsat)=ierx
           mprn(m) = ischan(jsat)
c          dummy the weights and amplitudes for now
           do kdat = 1,ndat
              obswgt(kdat,jsat)= 1.d0
           enddo                           
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
      elseif( upperc(obsfil).eq.upperc('C') ) then

         call readc4 (iobs
     .,               msat,mprn
     .,               iep,iyr,jdoy,sod
     .,               rclock,okmet,pres1,temp1,relhum1
     .,               nsave,save
     .,               kflag,ksite,latr,lonr,radius
     .,               l1z,l1n,l1e,l2z,l2n,l2e,antaz )  
         call check_y2k(iyr)
         call monday(jdoy,im,id,iyr)
         jdobs= julday(im,id,iyr)
         tobs = sod

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
     .,                     azim(i),azmdot(i)
     .,                     okwvr(i),wvrdel(i),atmdel(i)
     .,                     svcepc,svcl1
     .,                     tau(i),drate(i)
     .,                     ier(i),data_flag(i)
     .,                     ndats
     .,                     obs(1,i)
     .,                     omc(1,i)
     .,                     isnr(1,i)
     .,                     ampl1(i),ampl2(i)
     .,                     obswgt
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

      elseif(upperc(obsfil).eq.upperc('S')) then

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
          do j=1,ndat
            obswgt(i,j)= 1.d0
          enddo
        enddo    
c 

      else
         call report_stat('FATAL','MODEL','obsred',obsfil,
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
