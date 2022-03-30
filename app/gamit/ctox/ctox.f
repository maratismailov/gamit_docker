Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego, 1995. All rights reserved.

      Program CTOX
C
C     R.W. King  13 April 1987
C     New C-file format         :  12 Jul 1989    (R. King)
c     New C-file format         :  20 Dec 1991    (J.Genrich/Y. Bock)   
c     New C-file format (10.40) :  27 Aug 2010    (R. King)
c     New C-file format (10.60) :  28 Nov 2014    (R. King)
C
C     Write an X-File from a C-file.

c     Batch mode:        ctox <dfile> <c-letter> <x-letter>
c                        Translates all c-files to x-files
c                        No other options

c     Interactive mode:  ctox
c                        Dumping, decimating, normal-point options allowed

      implicit none

      include '../includes/dimpar.h' 
      include '../includes/units.h'
      include '../includes/errflg.h'
      include '../includes/makex.h'
      include '../includes/model.h'

      LOGICAL writex,owrite,bias(maxsat),lask,fcheck,wplot,wazel
      logical lgood,dump,lpart,allchan,allepochs

      CHARACTER*1 UPPERC,nbatch,cchar,xchar,fileid,latflag,lonflag
      character*3 daynum
      character*4 snames(maxsit)
      character*7 wildcard
      character*40 version
      CHARACTER*16 pickfn,fname,dfname,XFNAME,CFNAME,PFNAME 
      character*256 message

      integer*4 IY,IM,ID
     2      , IHR,MIN,IOERR,iux
     3      , i,icol,iarg1,iarg2
     4      , iarg3,istat,nstat,ndecim,idec,iepx,iep,newint,j
     5      , iidoy,idoy,jdoy,mchan,mstart,mstop
     6      , iclarg,mode,net,ishift,iyr,iflag,ichan,ityp,ipic
     7      , len,ians,nepochd,dlat,mlat,dlon,mlon
c           obsolete kinematic variable
     8      , kflag 
      integer lam(maxsit,maxsat,maxdat)

c     Function
      integer*4 julday 
      

      REAL*8    radkm,SEC,sod
     .        , pi,shift,relhum                  
     .        , slat,slon,l1z,l1n,l1e,l2z,l2n,l2e                            


      DATA ndecim/1/
      DATA owrite/.false./,writex/.false./,wazel/.false./
     1   , wplot/.false./,dump/.false./,lpart/.false./
     2   , allchan/.true./,allepochs/.true./
c  jfg change
      DATA PI/3.141592653589793D0/
c  jfg end change

c     Assign the unit numbers (in model.h commons)
      iterm = 5
      iscrn = 6
      iprnt = 8
      iud = 9  
c      input here (output in model)
      iuc = 20  
c      output here: rename  (input in model) 
      iobs = 40
      iux = iobs

C     Get and write the Program Name and Version Number

      call cversn(version)

c       Get the D-file name and input/output series letters from the
c       command-line argument if present, and translate C- to X-files

      iarg1 = iclarg(1,dfname)
      if( iarg1.gt.0 ) then
          iarg2 = iclarg(2,cchar)
          if ( iarg2.le.0 ) then
            call report_stat('WARNING','CTOX','ctox',' ',
     .'Incomplete runstring, use: ctox <dfile> <c-letter> <x-letter>',0)
          endif
          iarg3 = iclarg(3,xchar)
          if ( iarg3.le.0 ) then
            call report_stat('WARNING','CTOX','ctox',' ',
     .'Incomplete runstring, use: ctox <dfile> <c-letter> <x-letter>',0)
          endif
          open (unit=iud,file=dfname,status='old',err=900,iostat=ioerr)
          if (ioerr .ne. 0) then
            call report_stat('STATUS','CTOX','ctox',dfname,
     .      'Error opening D-File ',ioerr)
          endif
          call readdf( nbatch,iud,snames,daynum,nstat )
          fileid = 'D'
          writex = .true.
          wplot = .false.
c         Do not write everything to screen in copy mode
          dump = .false.
          iprnt = iscrn
          goto 50

      else

c       Interactive mode - check other options
              
   10   write(iscrn,11)
   11   format(/,1X,'Enter a C-file name, or',/
     .        ,1x,'      a D-file name to read a group of C-files :')
        wildcard = "[cd]*.*"    
        len = 7
        fname = pickfn(wildcard,len)
        fname = fname(1:len)
c        call pickfn('[cCdD]*.*',fname)   
        fileid = upperc(fname(1:1))
        if (fileid.eq.'C') then
          nstat=1
          cfname = fname
C         Find the period in the filename
          icol = index(cfname,'.')
          snames(1) = cfname(icol-5:icol-2)
          daynum    = cfname(icol+1:icol+3)
          cchar     = cfname(icol-1:icol-1)
        else if (fileid.eq.'D') then
          open (unit=iud,file=fname,status='old',err=10,iostat=ioerr)
          if (ioerr .ne. 0) then
            call report_stat('STATUS','CTOX','ctox',fname,
     .      'Error opening D-file ',ioerr)
          endif
          call readdf( nbatch,iud,snames,daynum,nstat )  
          write(iscrn,20)
   20     format(/,1X,
     .     'Series id (6th character) for input C-File names? ',$)
c           the following is a kluge to avoid reading a character-control
c           character from the format statement---Peng Fang / Bob King  990414
            read(iterm,'(a)') cchar,cchar
        else
          goto 10
        endif

c
c       Translate to X or dump to ASCII?
c
 25     continue              
        write(iscrn,30)
 30     format (1x,'Do you want to:',/,
     .  1x,'1 Translate the C-file(s) to X-file(s)',/,
     .  1x,'2 Dump the C-file(s) to an ASCII file',/,
     .  1x,'3 Make a plot of data availability',/,
c    .  1x,'4 Translate C-file(s) to X-file(s) with normal points',/,
     .  1x,'5 Write a file of az/el/SNR/ier values',//,
     .  1x,'Pick a number. ',$) 
        read(iterm,* ) ians 
        if ( ians .eq.1 .or. ians.eq.4 ) then
          write(iscrn,31)
 31       format(/,1X,
     .    'Series id (6th character) for output X-File names: ',$)
          read(iterm,'(a)') xchar
          writex = .true.
          wplot = .false.
c         Do not write everything to screen in copy mode
          dump = .false.
          iprnt = iscrn
          if( ians.eq.4 )
     .       call report_stat('STATUS','CTOX','ctox',' ',
     .   'Normal points no longer supported in CTOX--use AUTCLN',0)

        else if ( ians .eq.2 ) then
C         Open the output Print File
          write(iscrn,35)
 35       format (1x,'Enter file name for C-file dump.')
          read(iterm,'(a)') pfname
C         Open the output Print File 
c         check to make sure it's not the same as the input (typo blunder)
          if( pfname.eq.cfname ) call report_stat('FATAL','CTOX','ctox'
     .      ,' ','Output and input names the same',0)
          open(unit=iprnt,file=pfname,status='unknown',iostat=ioerr)
          if (ioerr .ne. 0) then
            call report_stat('STATUS','CTOX','ctox',pfname,
     .      'Error opening P-File ',ioerr)
          endif
          write(iscrn,36) pfname
 36       format(1x,'Opened: ',a)
          write(iprnt,37) version
 37       format(/,10x,'CTOX Version ',a40,//)
          dump = .true.
          writex = .false.
          wplot = .false.
          mstart = 1
          write(iscrn,'(a,$)') ' Dump all epochs (y/n)?'
          if (lask()) then
            allepochs = .true.
          else
            allepochs = .false.
            write(iscrn,'(a,$)') 'Enter start, stop epochs:'
            read(iterm,*)  mstart,mstop
          endif
          write(iscrn,'(a,$)') ' Dump all channels (y/n)?'
          if (lask()) then
            allchan = .true.
          else
            allchan = .false.
            write(iscrn,'(a,$)') 'Enter (1) channel to dump:'
            read(iterm,*)  mchan
          endif
          write(iscrn,'(a,$)') ' Dump partials (y/n)?'
          if (lask()) then
            lpart = .true.
          else
            lpart = .false.
          endif

        else if ( ians .eq.3 ) then
          dump = .false.
          writex = .false.
          wplot = .true.
          iprnt = iscrn

        else if ( ians.eq.5 ) then
          dump = .false.
          writex = .false.
          wplot = .false.
          iprnt = iscrn
          wazel = .true.
        else
          goto 25
        endif

c       See if decimation is desired
        write(iscrn,'(a,$)') ' Do you wish to decimate the data?'
        if (lask()) then
          write(iscrn,'(a,$)') 'New sampling interval in seconds? '
          read(iterm,'(i4)')  newint
c         set ndecim=0 as a flag until the old INTER read from the C-file
          ndecim = 0
          write(iscrn,'(a,$)') 'Choose an integer for shift(0,1,2...) '
          read(iterm,'(i4)')  ishift
        else
          ndecim = 1
          ishift = 0
          shift = 0.0d0
        endif

      endif

c     Begin loop over files
   50 do 300 istat=1,nstat
         if( fileid.eq.'D' ) then
            cfname= 'C'//snames(istat)//cchar//'.'//daynum
         endif
         call lowers(cfname)
         call copens(cfname,'OLD',IUC,IOERR)
         if (ioerr .ne. 0 ) then
            call report_stat('WARNING','CTOX','ctox',cfname,
     .      'Cannot open C-file, continuing with next C-file',ioerr)
             goto 300
         else
            call report_stat('STATUS','CTOX','ctox',cfname,
     .      'Opened C-file: ',ioerr)
         endif

                             
         if(writex) then
C           Open the Output X File
            xfname= 'X'//snames(istat)//xchar//'.'//daynum
            call lowers(xfname)
            if (fcheck(xfname) .and. .not. owrite) then
               WRITE(ISCRN,61) xfname
   61          FORMAT(/,1X,'Output X-File ',a16,' already exists.',/,
     .         1x,'Overwrite?',$)
               if (lask()) then
                  owrite = .true.
               else
                  owrite = .false.
               endif
            endif   
            if (owrite) then
               OPEN (UNIT=IUX,FILE=XFNAME,STATUS='OLD',iostat=ioerr)
            else
               OPEN (UNIT=IUX,FILE=xfname,STATUS='NEW',iostat=ioerr)
            endif
            owrite = .false.
            if (ioerr .eq. 0) then
               call report_stat('STATUS','CTOX','ctox',xfname,
     .        'Opened X-file: ',ioerr)
            else
               call report_stat('FATAL','CTOX','ctox',xfname,
     .        'Error opening X-file ',ioerr)
            endif
         endif

         if(wplot) then
C           Open the plot file
            xfname= 'a'//snames(istat)//cchar//'.'//daynum
            call lowers(xfname)
            OPEN (UNIT=IUX,FILE=XFNAME,STATUS='unknown',iostat=ioerr)
            if (ioerr .eq. 0) then
               call report_stat('STATUS','CTOX','ctox',xfname,
     .        'Opened X-file ',ioerr)
            else
               call report_stat('FATAL','CTOX','ctox',xfname,
     .        'Error opening X-file ',ioerr)
            endif
         endif

         if(wazel) then
C           Open the azel file
            xfname= 'z'//snames(istat)//cchar//'.'//daynum
            call lowers(xfname)
            OPEN (UNIT=IUX,FILE=XFNAME,STATUS='unknown',iostat=ioerr)
            if (ioerr .eq. 0) then
               call report_stat('STATUS','CTOX','ctox',xfname,
     .        'Opened X-file ',ioerr)
            else
               call report_stat('FATAL','CTOX','ctox',xfname,
     .        'Error opening X-file ',ioerr)
            endif
         endif
                 

C        Read and Print the C-File Header Records      

         call chdred ( dump,iy,im,id,ihr,min,sec )    
         latr_sph = preval(1)
         lonr     = preval(2)
         radius   = preval(3)*1.d3
         if( allepochs) mstop = nepoch
c        convert start time to GPST if in UTC
         if( mtime.ne.2 ) then
           jdoy = idoy( iy,im,id )
           call utc2gps( jdoy,iy,ihr,min,sec )
         endif 
c        compute (PEP) JD and seconds-of-day and store in model.h
         jd0 = julday(im,id,iy)
         t0  = ihr*3600.d0 + min*60.d0 + sec

c        Calculate the decimation factor

c        Check if an even multiple of inter
         if (ndecim .ne. 1) then
            ndecim = newint/inter
            if (newint .ne. ndecim*inter ) then
               write(message,70)newint,inter
   70          format('New sampling interval: ',i4,' is not an even'
     .                ,' multiple of the old interval: ',i4)
               call report_stat('FATAL','CTOX','ctox',' ',message,0)
            endif
            nepoch = nepoch/ndecim
            if (nepoch .gt. maxepc) nepoch = maxepc
c           shift the initial epoch to the middle of first new interval
            shift = dble(ishift*inter)
            inter = newint
            if( writex ) then 
               write (iscrn,*) 'Writing ',nepoch,' epochs to X-file'
            elseif ( dump ) then
              nepochd = (mstop-mstart)/ndecim + 1
              write (iscrn,*) 'Dumping ',nepochd,' epochs '
            endif
            write (iscrn,*) 'Decimating by a factor of ',ndecim
         endif

C        Write the X-File Header Records   
         if (writex) then
            call xhdrit( iux,cfname,shift )
         endif       

c        Write the header for the az/el file

         if( wazel) then
           write(iux,80) xfname,sitnam
     .                 , iy,im,id,ihr,min,sec,inter,nepoch
     .                 , nchan,(ischan(i),i=1,nchan)
* MOD TAH 200618: Updated 32I to 50I to allow for 35 Beidou satellites
   80      format(1x,'File: ',a16,'  Skymap values for site: ',a32,/
     .              ,1x,'Start: ',i4,2i3,1x,2i3,f5.2
     .              ,'  Interval: ',i3,'  # epochs:',i4,/
     .              ,1x,i2,' satellites: ',50i3)
           write(iux,'(/,a,/)')
     .      ' Epoch    For each SV:  PRN  IER  el  az  SNR(L1,L2)'
         endif


C        Loop over the data records
         ipic = ishift + 1   
         iepx = 0 
cd         print *,'CTOX nepoch ndecim shift ipic '
cd     .      ,nepoch,ndecim,shift,ipic
         DO 220 i=1,nepoch
           do 210 idec = 1, ndecim
             call cdtred( iidoy,lpart,allchan,mchan
     .                  , mstart,mstop,idec,ipic )
cd           print *,'CTOX nepoch i idec iepx ipic'
cd     .         ,nepoch,i,idec,iepx,ipic 
c          convert time to GPST if in UTC
           if( mtime.ne.2 ) then 
             call ds2hms(iy,jdobs,tobs,ihr,min,sec)
             call utc2gps( jdoy,iy,ihr,min,sec )
             sod = dble(ihr)*3600.d0 + dble(min)*60.d0 + sec
           endif

C         Station Coordinates       
         
          call raddms( latr,latflag,dlat,mlat,slat )
          if( latflag.eq.'-') then
            latflag = 'S'
          else
            latflag = 'N'
          endif
          call raddms( lonr,lonflag,dlon,mlon,slon )
          if( lonflag.eq.'-' ) then
           lonflag = 'W'
          else
           lonflag = 'E'
          endif

c         Set BIAS array to indicate that an extra bias flag
c         must be inserted at the next selected epoch.
c         BIAS flag set false after written in XDTRIT.
          do  j=1,nchan
             if (ier(j) .eq. igbias) bias(j)=.true.
          enddo

c         Write the record to the X-file
          if (idec .eq. ipic .and. writex) then    
             iepx = iepx + 1
              call xdtrit( iux,iepx,bias ) 
c RWK 160420
c     .            , latflag,dlat,mlat,slat,lonflag,dlon,mlon,slon,bias )
           endif          

           if (wplot) then
              do j = 1,nchan
                if (lgood(ier(j))) then
                  write (iux,'(f10.2,1x,i2)') sod,ischan(j)
                endif
             enddo
           endif
           if (idec.eq.ipic .and. wazel) then
              write(iux,fmt='(i4,32(2x,2i3,1x,f4.1,1x,f6.1,2i2))')
     .              iepx,(ischan(j),ier(j)
     .              , elev(j)*180./pi
     .              , azim(j)*180./pi
     .              , isnr(1,j),isnr(2,j), j=1,nchan)
c              add later L1/L2 RINEX SNR
           endif           



c--------end of loop over all output decimated epochs
  210       continue                           
c--------end of loop over all input epochs
  220    CONTINUE

c---- end loop over files---
  300 continue

c     END LOOP OVER FILES

 500  call report_stat('STATUS','CTOX','ctox',' ',
     .'Normal end of CTOX',0)
      stop

 900  call report_stat('FATAL','CTOX','ctox',' ',
     .'Error input D-file not found ',ioerr)

      END

