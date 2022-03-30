      Subroutine seschk( xorc, xfile, lxfil, lsesfo
     .                 , year, doy, isessn, nsod, nepoch, inter, ispan )
                            
c     Read the year, day-of-year, span, and SVs from session.info and check
c     these against the values in the headers of the X- or C-files.

c     R. King   6 April 1992             
c     Originally did not require session.info to be present and used X-, C-, or
c     I-file to get the session parameters.  Re-ordered by R. 15 August 2006 
c     to assume session.info is available.

C     Input:
c        XORC   : 1st character of X- or C-file name
C        XFILE  : X-file name
C        LXFIL  : unit number of X- or C-file  
c        IFILE  : I-file name
c        LIFIL  : unit number of I-file
c        IFLSW  : new I-file being written (T) or not (F)
c        LSESFO : unit number of session.info
c     Output     
c        year  : 4-digit integer year
c        doy   : integer day-of-year 
c        ISESSN : session number (assumed = 1 if session.info not available)
c        NSOD   : seconds-of-day at start rounded to the even minute   
c        ISPAN  : seconds in session
c
      implicit none

      include '../includes/dimpar.h'

      character* 1  xorc,skd,latflag,lonflag,gnss 
      character* 3  rcvrsw,rxobtyp(maxdat)
      character*4   antmod,svantmod(maxsat),drymap,wetmap,ionsrc
      character*5   frame,precmod,nutmod,gravmod,srpmod
     .,            eradmod,antradmod
      character* 6  antcod,magfield
      character*8   etidemod,atmlmod,otidemod,atmtide,hydrolmod
     .,             speopmod, cextra(maxext) 
      character*10  antmod_snx,svantmod_snx(maxsat)
      CHARACTER*16  XFILE,OBFILN,TFILN,jfiln,sitnam16,satnam(maxsat)
      character*20  rctype,rcvnum,anttyp,antnum,svantbody(maxsat)
      CHARACTER*32  SITNAM
      CHARACTER*80  TEXT(MAXTXT)
      character*256 message

      integer*4 lxfil,lsesfo,ioerr,iscrn,npart,nslip
     .        , idoy,ntext,ndat,lambda(maxsat,maxdat)
     .        , mtime,dattyp(maxdat),ischan(maxsat),avlmet,jde,jdr
     .        , ihr,min,nchan,nepoch,inter,ircint
     .        , nepoch1,inter1,doy1,ispan1,iy1,im1,id1
     .        , ihr1,min1,int1,iprnt,latd,latm,lond,lonm
     .        , isessn,nsod,nsod1,check_flg
     .        , doy,year,nyr,ndoy,norb,ispan
     .        , niextra,iextra(maxext),nrextra,ncextra
     .        , nchan1,ischan1(maxsat),ietide,isptide,nclock,i,j

      INTEGER*2    ISLIP(MAXCSB),ISLPST(MAXCSB)

      LOGICAL debug,fcheck,badsesfo

C
      REAL*8  offarp(3),offsl1(3),offsl2(3),te,tr,sec
     .     ,  ut1,ut1dot,xp,xpdot,yp,ypdot,psi,psidot,eps,epsdot
     .     ,  sec1,seclat,seclon,ht,elvcut,clock(4)
C MOD TAH 200126: Changed to svantdx(3,2,maxsat) (added L1/L2 index)
     .     ,  sod,sod1,atmlavg(3),hydrolavg(3),svantdx(3,2,maxsat)
     .     ,  rextra(maxext),fL1(maxsat),fL2(maxsat)
      real*8 antdaz  ! Antenna aligment from True N (deg).

      REAL*4  SWVER

      iscrn = 6

c     Initialize  seconds-of-day to 99999 as a flag for FIXDRV

      nsod = 999999  
      ispan = 86400
      
c------------------------------------------------------------------
                                                  
c     Read session.info and set the session ids 

      debug = .false.    
      check_flg = 0
      badsesfo = .false. 
c      print *,'SESCHK calling rsesfo yr doy ',year,doy
      call rsesfo( lsesfo,debug,check_flg,year,doy,isessn
     .           , ihr,min,inter,nepoch,nchan,ischan,badsesfo ) 
c      print *,'SESCHK ',lsesfo,debug,check_flg 
c      print *,'SESCHK ',year,doy,ihr,min,inter,nepoch,nchan 
      sod = dfloat(ihr*3600 + min*60)
      call even_minute(year,doy,sod,nyr,ndoy,nsod)           
    

c-------------------------------------------------------------------
c     If the X- or C-file available (not required), read the
c     year, day, and session information from the header

      if( fcheck(xfile) ) then
             
c       Read from the X-file

        if( xorc.eq.'x' ) then
C          open the X-file
           OPEN( LXFIL, FILE=XFILE, STATUS='OLD' )
           iprnt = 0
           call xhdred ( lxfil,iprnt,iscrn
     .         , nepoch1,inter1,ircint,isessn
     .         , mtime,iy1,im1,id1,ihr1,min1,sec1
     .         , nchan1,ischan1,satnam
     .         , ndat,dattyp,rxobtyp,lambda
     .         , offarp,sitnam16,rcvrsw,swver,antcod
     .         , rctype,rcvnum,anttyp,antnum
     .         , latflag,latd,latm,seclat,lonflag,lond,lonm,seclon,ht
     .         , ntext,text,gnss )
           doy1 = idoy(iy1,im1,id1)
           sod1 = ihr1*3600.d0 + min1*60 + sec1  
           ispan1 = (nepoch1-1)*inter1 
C          close the X-file
           close( lxfil )

c       Read from the C-file

        else if( xorc.eq.'c' ) then  
c          open the C-file
           CALL  COPENS( XFILE, 'OLD', LXFIL, IOERR )
C          read the headers
           CALL  READC1( LXFIL, NTEXT, TEXT )
           CALL  READC2( LXFIL, SITNAM 
     .,             rctype,rcvnum,rcvrsw,swver,anttyp,antnum 
     .,             NPART, norb, gnss, nchan1, ischan1, fL1, fL2
     .,             NDAT, DATTYP, LAMBDA, skd, nepoch1, inter1, ircint 
     .,             MTIME, isessn, iy1, im1, id1, ihr1,min1, sec1 
     .,             offarp, OFFSL1, OFFSL2, antdaz, svantdx 
     .,             OBFILN, TFILN, JFILN 
     .,             frame,precmod,nutmod,gravmod,srpmod  
     .,             eradmod,antradmod
     .,             ietide,isptide,speopmod
     .,             etidemod,otidemod,atmtide,atmlmod,hydrolmod  
     .,             atmlavg,hydrolavg,drymap,wetmap,ionsrc,magfield
     .,             antmod_snx,antmod,svantbody,svantmod_snx,svantmod
     .,             elvcut,nclock,clock
     .,             JDE, TE, JDR, TR 
     .,             UT1, UT1DOT, XP, XPDOT, YP, YPDOT 
     .,             PSI, PSIDOT, EPS, EPSDOT, AVLMET
     .,             NSLIP, ISLIP, ISLPST  
     .,             niextra,iextra,nrextra,rextra,ncextra,cextra)
           doy1 = idoy(iy1,im1,id1)
           sod1 = ihr1*3600.d0 + min1*60 + sec1    
           ispan1 = (nepoch1-1)*inter1
c          this probably not necessary--fixed prior to MODEL
           if( isessn.eq.0 ) isessn = 1
c          close the C-file
           CLOSE( LXFIL )

        else
          call report_stat('FATAL','FIXDRV','seschk',' '
     .          , 'Input observation file type neither X nor C',0)
        endif

c       if X or C-file is UTC, convert time to GPST
        if( mtime.eq.1 ) then
          call utc2gps(doy1,iy1,ihr1,min1,sec1)   
          sod1 = ihr1*3600.d0 + min1*60 + sec1
        elseif (mtime.ne.2 ) then
          call report_stat('WARNING','FIXDRV','seschk',' '
     .          , 'mtime neither 1 (UTC) nor 2 (GPST)',0)
        endif

c       Check the yr, day, and span  against session.info
 
        ispan1 = (nepoch1-1)*int1
        nsod1 = ihr1*3600.d0 + min1*60.d0 + nint(sec1)
        if( nchan.eq.nchan1 ) then 
c         sort X- or C-file SVs before checking against sesfo--ok since not passed out of this routine
          call sort1i( nchan,ischan )                
          do i=1,nchan
            if( ischan(i).eq.ischan1(i) ) then
              continue
            else
              call report_stat('WARNING','FIXDRV','seschk',' '
     .               ,'SVs in x-file do not match session.info:',0)
* MOD TAH 200618: Updated 32I to 50I to allow for 35 Beidou satellites
              write(message,'(a,a16,i2,3x,50i3)')
     .             '   x-file ',xfile,nchan1, (ischan1(j), j=1,nchan1)
              call report_stat('WARNING','FIXDRV','seschk',' '
     .                  ,message,0)
              write(message,'(a,a16,i2,3x,50i3)')
     .    '             session.info ',nchan, (ischan(j), j=1,nchan)
              call report_stat('WARNING','FIXDRV','seschk',' '
     .                    ,message,0)
                  call report_stat('FATAL','FIXDRV','seschk',' '
     . ,'SVs in X-file do not match session.info--see GAMIT.warning',0)
            endif
          enddo
        else
          write(message,'(a,i2,a,i2)')
     .       'Number of SVs in X- or C-file (=',nchan1
     .       ,') does not match number is session.info file (=',nchan
            call report_stat('FATAL','FIXDRV','seschk',' ',message,0)
        endif
                         
        if( inter.ne.inter1 .or. nepoch.ne.nepoch1 .or.
     .           iabs(nsod-nsod1).gt.10 ) then
          write(message,'(a,4i4,f7.2,i6,2i5)')
     . 'X- (or C-) file =',iy1,doy1,ihr1,min1,sec1,nsod1,inter1,nepoch1
          call report_stat('WARNING','FIXDRV','seschk',' ',message,0)
          write(message,'(a,4i4,f7.2,i6,2i5)')
     .     'session.info    =',nyr,ndoy,ihr,min,sec,nsod,inter,nepoch
          call report_stat('WARNING','FIXDRV','seschk',' ',message,0)
          write(message,'(a,a)')'Non-matching scenarios in '
     .         ,'session.info and input X- (or C-) file'
          call report_stat('WARNING','FIXDRV','seschk',' '
     .          , message,0)
        endif

c     endif if on existence of X-file
      endif

      RETURN
      END

Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      subroutine sort1i(ndat,iarray)
c
c     sort an integer array from the smallest to the bigest  
c     this version taken from /solve by rwk 010828; put both (and others)
c     in /lib eventually

      implicit none
                     
     
      integer ndat,iarray,iswtch,i,k
      dimension iarray(ndat)
c
      if(ndat.le.1) go to 100
      do 50 k=1,32000
         iswtch=0
         do 30  i=2,ndat
            if(iarray(i).ge.iarray(i-1)) go to  30
            iswtch=iarray(i)
            iarray(i)=iarray(i-1)
            iarray(i-1)=iswtch
 30      continue
         if(iswtch.eq.0) go to 100
 50   continue
c
 100  continue
c
      return
      end








