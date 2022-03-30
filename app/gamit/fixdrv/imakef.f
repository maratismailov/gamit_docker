      SUBROUTINE IMAKEF ( XFILE, XORC, KFILE, KPICK, CLKFIT,
     .                  LXFIL, LCFIL, formti, lprint )

C     This subroutine constructs an I-file for FIXDRV
C
C     Written  by R King 91/12/31  from SMAKEF (Bock, Shimada, Feigl)
c                                  and
C
C     arguments
C        XFILE : X-file name                         input
C        XORC  : first character of X,C-file         input
C        KFILE : clock file name                     input
C        KPICK : K-file flag                         input
C                  = Y Use K-file
C                  = N Do not use K-file
C        CLKFIT : clock fit polynomial function      input
C                  = 'C' Use cubic fitting
C                  = 'L' Use linear fitting
C                  = ' ' Use program default value
C        LXFIL  : unit number of X-file              input
C        LCFIL  : unit number of C-file              input
c        FORMTI : format statement for I-file        input
c        LPRINT : unit of print file                 input

C     files
C                   scratch file         LUN= 8      output
C                   S-file               LUN=16      output
C        XFILE    : X-file               LUN=LXFIL   input
C        CFILE    : C-file               LUN=LCFIL   input
C        KFILE    : K-file               LUN=20      input
C        leap.sec : leap second table    LUN=99      input

C     subroutines and functions
C        IMAKEF
C          +-CHDRED +-READC1 +-FERROR .
C          |        +-READC2 +-FERROR .
C          |        +-READC3 +-FERROR .
C          +-CLKPRM +-DAYNUM .
C          |        +-JULDAY
C          |        +-MONDAY .
C          +-COPENS --FERROR .
C          +-LFNTOP .
C          +-LOWERS .
C          +-NBLEN .
C          +-UPPERS .
C          +-XHDRED +-MONDAY .
C                   +-UPPERC .
C
      implicit none

      include '../includes/dimpar.h'

      logical pseudo_range

      real* 8 doy

      real*4 swver

      integer*4  ndat,dattyp(maxdat),lambda(maxsat,maxdat),ntext
     .         , lxfil,ioerr,inter,isprn,idyoyr
     .         , id,im,iy,idoy,ihr,imin,min,ircint,isessn
     .         , latd,latm,lond,lonm,nsat,nepoch,iclk
     .         , lcfil,nblen,idumb,mtime,it0,iscrn,iscratch,lprint
     .         , i0,i1,icol,lfntop,nslip,i

      integer*2  islip(maxcsb),islpst(maxcsb)

      real*8  twopi,t00,rate,sec,accel,cubic,epoch
     .      , offarp,offsl1,offsl2,height,seclat,seclon

      character* 1  unitcd, kpick, jumpfit, nojumpfit,clkfit, lowerc
     .            , xorc,latflag,lonflag,gnss
      character* 3  day,rcvrsw,rxobtyp(maxdat)
      character* 4  sname, clkmas
      character* 6  antcod
      character*16  xfile,kfile,kfile2,sitnam16,cfile,satnam(maxsat)
      character*20  rctype,rcvnum,anttyp,antnum
      character*32  sitnam
      character*48  formti
      character*64  coords
      character*80  text(maxtxt)

      dimension  it0(3), t00(3)
c     antenna offsets, up north east, ARP, L1 and L2
      dimension  offarp(3),offsl1(3), offsl2(3)
c     satellites
      dimension  isprn(maxsat)

      parameter( iscrn=6, iscratch=8 )

      twopi = 8.d0 * datan(1.d0)

c     initialize
      write (sitnam,'(32x)')
      do  10  i = 1, 63
         coords(i:i)=' '
   10 continue
      clkmas = '    '

c     determine day number
c     positions of first and last characters of x-file name
      i0 = lfntop( xfile )
      i1 = nblen( xfile )

c**     session number now read from X-file or set = 1 in XHDRED
C**     isessn = 1

c     determine series letter
      unitcd = xfile(i0+2:i0+2)
      call  uppers( unitcd )
      day(1:3) = xfile(i1-2:i1)

c     output site code (4 characters after first character)
      sname(1:4) = xfile(i0+1:i0+4)
c     render upper case
      call uppers (sname)

C     X- or C-file header read----------------------------------------

      open( 8, status='scratch' )
      if( lowerc(xorc) .eq. lowerc('x')) then
c        read x-file header
         open( lxfil, file=xfile, status='old',err=90 )
         call  xhdred ( lxfil, iscratch, iscrn, nepoch, inter, ircint
     .                , isessn, mtime, iy, im,id, ihr, min, sec
     .                , nsat, isprn, satnam
     .                , ndat, dattyp, rxobtyp,lambda
     .                , offarp, sitnam16,rcvrsw,swver,antcod
     .                , rctype,rcvnum,anttyp,antnum
     .                , latflag, latd, latm, seclat
     .                , lonflag, lond, lonm, seclon, height
     .                , ntext, text, gnss )
         close( lxfil )                 
         sitnam(1:16) = sitnam16
      else
c        read c file header
         cfile = xfile
         cfile(i0:i0) = 'c'
         call  lowers( cfile(i1-4:i1-4) )
         call  copens( cfile, 'old', lcfil, ioerr )
         call  chdred( lcfil, 8, nepoch, inter, mtime,
     .                 iy, im, id, ihr, min, sec,
     .                 nsat, isprn, ndat, dattyp, lambda,
     .                 offarp, offsl1, offsl2, sitnam, rcvrsw, swver,
     .                 rctype,rcvnum,anttyp,antnum,
     .                 ircint, isessn, ntext, text,
     .                 nslip, islip, islpst )
         close( lcfil )
      endif
      close(8)

c     if X or C times are UTC, convert to GPST
      idyoyr = idoy( iy,im,id )
      if( mtime.eq.1 ) then
        call utc2gps( idyoyr,iy,ihr,min,sec )
        call monday( idyoyr,im,id,iy )
      elseif( mtime.ne.2 ) then
        call report_stat('FATAL','FIXDRV','imakef',' '
     .     ,'mtime from X- or C-file neither 1 (UTC) nor 2 (GPST)',0)
      endif
c     put the times into the arrays
      it0(1) = im
      it0(2) = id
      it0(3) = iy
      t00(1) = dble(ihr)
      t00(2) = dble(min)
      t00(3) = sec


C     determine clock information-----------------------------------------------

c     check for pseudoranges: if they exist, assume k-files exist
      pseudo_range = .false.
      do i = 1, ndat
         if(  dattyp(i).eq.3
     .   .or. dattyp(i).eq.4
     .   .or. dattyp(i).eq.5) pseudo_range = .true.
      enddo

c     compute 3rd order polynomial from k-file
      epoch = 0.d0
      rate  = 0.d0
      accel = 0.d0
      cubic = 0.d0

c     see if k-file to be used    
      if( lowerc(kpick).eq.lowerc('y') )  then
         if( pseudo_range ) then
c          for instruments with pseudoranges
c          for now , set the new calling arguments to agree with old scheme:
c             linear no-jump fit, cubic jump fit; clock fit is input to 
c             imakef from the sittbl.  The current documentation describes
c             values of L (linear) for the nojump polynomial or C (cubic) for 
c             the with-jump polynomial, but we can recode the following for 
c             a variety of entries.  Clkera now recognizes only 'j' (with-jump),
c             'n' (no-jump), or blank (choose by residuals).
           nojumpfit = 'L'    
c rwk 971205: Cubic fits often fail to detect jumps; let's try linear for a while
c           jumpfit = 'C'
           jumpfit = 'L'
           if ( lowerc(clkfit).eq.'l' ) then
              clkfit = 'n'
           elseif ( lowerc(clkfit).eq.'c' ) then
              clkfit = 'j'
           else
              clkfit = ' '
           endif
           call clkera( kfile, sname, nojumpfit, jumpfit, clkfit
     .                , it0, t00, lprint
     .                , idumb, epoch, rate, accel, cubic )
         else
c          for instruments without pseudoranges
c          k-file name come from d-file or d-file+day no.
           kfile2 = kfile
           icol = index( kfile, '.' )
           kfile2(icol+1:icol+3) = day(1:3)
           call  clkprm( kfile, kfile2, it0, t00, unitcd, iclk,
     .                 epoch, rate )
c          first station is default master ?
           if( iclk .eq. 1 )  clkmas = sname
        endif
      endif

C     write I-file record for site -----------------------------------

C     start time
      IHR  = T00(1)
      IMIN = T00(2)
C     start year and doy
      call daynum(im,id,iy,doy)
      idyoyr=idint(doy)
c      write(buf3,'(a3)') day
c      read (buf3,'(i3)') idyoyr
      write(16,formti) sname,iy,idyoyr,isessn,ihr,imin,t00(3)
     .               , epoch,rate,accel,cubic

      RETURN

  90  call report_stat('FATAL','FIXDRV','imakef',xfile
     .   ,'In making I-file, error opening xfile:',0)

      END

