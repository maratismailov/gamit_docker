      subroutine wxhead (lu,stanam,coords,offarp,
c**     2           l1up,l1north,l1east,l2up,l2north,l2east,
     3           gnss,icount,satarray,ndat,dattyp,lambda, 
     .           rxobtyp,iobtypx,
     4           iwknstart,sowstart,interval,nepoch,isessn,
     5           version,irunt,xhdlin,ixhdln,
     6           uname,rcvrsw,swver,rcvcod,ircint,exptyp,antcod,
     7           dcb_override )
C
c     Subroutine to write the header portions of the Xfile
c
c     written by Peter Morgan for the Apollo  February 1987
C     modified by Kurt Feigl July 88
c     modified by Bob King   July 89
c     modified by Yehuda Bock April 25, 1990 (for coordinate conventions)
c
c
c     DATA DECLARATIONS

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'

      logical  dcb_override
      character      aline1*80,aline2*80,stanam*16
     .        ,      latflag*1,lonflag*1,version*40,yawbias*1
C     user name
      character*16   uname
c     receiver software
      character*3 rcvrsw
c     receiver software version
      real*4 swver
c     antenna code
      character*6 antcod,rcvcod
C     header lines from input file
      character*80   xhdlin(maxlin)     
c     GNSS code
      character*1 gnss
c     GNSS name
      character*16 satnam(maxsat)
c     SV antenna/body-type
      character*20 antbody     
c     Frequency channel for GLONASS
      integer*4 frqchn
c     svnav.dat start/stop times 
      integer*4 svnstart(5),svnstop(5)
c     experiment type
      character*9 exptyp
c     time-type for X-file
      character*4 x_time_flag
c     observable types and wavelength factors                        
      integer*4 ndat,dattyp(maxdat),iobtypx(6),lambda(maxsat,maxdat)
      character*3 rxobtyp(maxobt),obtypx(maxobt)
      integer*4      dlat,mlat,dlon,mlon,icount,satarray(maxsat)
     .        ,      itflag,isessn,jd,julday
     .        ,      iyear,iy2,imo,id,iday,ih,imin,ircint,interval
     .        ,      mchkey,istr,j
C     unit number for x file (opened)
     .        ,      lu
c     number of epochs
     .        ,      nepoch
C     number of lines in X file header
     .        ,      ixhdln
C     starting epoch week number
     .        ,      iwknstart
C     run time yr,mo,da,hr,min,sec
     .        ,      irunt(6)
     .        ,      i
c     sat IDs: SV and PRN
      integer nsn,prn

C     starting second of week
      real*8      sowstart
     .        ,   seclat,seclon,radius,seconds,utcoff,gpstutc
c**     .        ,   l1up,l1north,l1east,l2up,l2north,l2east
     .        ,   sbmass,taiutc,yawrate
     .        ,   coords(3),latr,lonr
     .        ,   offarp(3)  
      real*8 antpwr   ! Antenna transmitt power (W)       


cd      print *,'WXHEAD dattyp lambda(1); ',dattyp,(lambda(1,i),i=1,4)
                             
c        Write a dummy title line

      write(lu,10)
  10  format('GPS Phase and Pseudorange for GAMIT Processing ')


c        Write the MAKEX version number, user name and run time,
c        into X-file header

      write (lu, 11) version,uname,(irunt(i),i=1,6)
  11  format (/,1x,'MAKEX ver ',a40,/
     .   ,1x,'Run by ',a16,1x,'on',1x,
     .   i4,'-',i2.2,'-',i2.2,1x,i2.2,':',i2.2,':',i2.2,//)


c        Copy the input data file header lines to the X-file.
                                 
      dcb_override = .false.
      do i=1,ixhdln    
         write(lu,'(a80)') xhdlin(i)  
c        test each line for the use of a CC2NONCC'd RINEX file
         istr = mchkey(xhdlin(i),'CC C1, P2 converted',80,19)
         if( istr.gt.0 ) then
            dcb_override = .true.
         endif
      enddo    
      if( dcb_override) then
         call report_stat('WARNING','MAKEX','wxhead',' '
     .    ,'RINEX file has DCB corrections applied by CC2NONCC',0)
         call report_stat('WARNING','MAKEX','wxhead',' '
     .     ,' --set DCB flag on x-file to 0',0) 
         dcb_override = .true.
      endif

c       Write a blank character line

      write(lu,'(a)') ' '


c       Terminate the header block of the xfile with the END code.

      write(lu,'(a3)') 'END'


c       Write the coordinate file information. The coordinate file
c       contains the best apriori site information.

C     New conventions for coordinates are:
c        Right-handed coordinate system.
c        latitude can be  'N' or 'S'
C        longitude can be 'E' or 'W'
C        Yehuda Bock 4/25/90
c     For history's sake, the OLD convention was:
c        Left-handed coordinate system.
c        latitude could be  ' ' (positive) west  or '-' (negative) east
c        longitude could be ' ' (positive) north or '-' (negative) south
      write(lu,'(20x,a)') 'SPHERICAL COORDINATES'
      write(lu,13)
  13  format(1x,'STATION_NAME____     LATITUDE         LONGITUDE      RA
     .DIUS     RCVR SWVER  RCVCOD'
     .,/,1x,17x,'_DG MN SS.SSSSS __DG MN SS.SSSSS ____(M).____')

      call xyz2sph(coords,latr,lonr,radius)  
      call raddms( latr,latflag,dlat,mlat,seclat )
      if( latflag.eq.'-') then
        latflag = 'S'
      else
        latflag = 'N'
      endif
      call raddms( lonr,lonflag,dlon,mlon,seclon )
      if( lonflag.eq.'-' ) then
        lonflag = 'W'
      else
        lonflag = 'E'
      endif
      write(lu,14)
     *  stanam,latflag,dlat,mlat,seclat,lonflag,dlon,mlon,seclon,radius
     . ,rcvrsw,swver,rcvcod
 14   format(1x,a16,1x,a1,i2,1x,i2,1x,f8.5,1x,a1,i3,1x,i2,1x,f8.5,f13.4
     . ,2x,a3,1x,f5.2,1x,a6)


c       Write a null character line

      write(lu,'(80x)')


c       Write the antenna offset information
                          
c** rwk 050927: Write the ARP (mechanical) offsets, not the phase center offsets
      write(lu,60)
c*  60  format(' ANT OFFSETS (M) L1    UP     NORTH   EAST   L2   UP     N
c*     1ORTH   EAST') 
  60  format(' ANT ARP OFFSETS (M)   UP     NORTH   EAST ')              
c**      WRITE(lu,61) antcod,L1UP,L1NORTH,L1EAST,L2UP,L2NORTH,L2EAST
      write(lu,61) antcod,offarp
c**  61  FORMAT(19X,3F8.4,3X,3F8.4,/)
c**  61  FORMAT(1x,a6,12x,3F8.4,3X,3F8.4,/)
   61 format(1x,a6,12x,3f8.4,/) 


c       Convert from GPS week,seconds to date,time (now hardwired to GPST)

c**      itflag = +1  --old conversion to UTC
      itflag =  +4
      x_time_flag = 'GPST'
      call timcon (itflag,iwknstart,sowstart,iyear,
     1 iday,ih,imin,seconds,utcoff) 
c     utcoff will be returned zero for GPST/GPST case
      call monday(iday,imo,id,iyear)
      jd = julday(imo,id,iyear)
      gpstutc = taiutc(jd) - 19.d0    
c** rwk 150113: The loop here was on ndat, with format 7i3, but since we 
c        have no provision as yet to use more than 4 observables, and want
c        to add the observable codes as documentation (and for xtorx), 
c        hard-wire the loop and format to '4'.
      do i=1,4
        if(iobtypx(i).ne.0 ) then
          obtypx(i) = rxobtyp(iobtypx(i))
        else
          obtypx(i) = ''
        endif                                                         
      enddo
      write(lu,'(1x,i2,a,6x,a,4i3,2x,4(1x,a3))')
c** rwk 170601  icount,' SATELLITES',' 4 DATA TYPES:',(dattyp(j),j=1,4)
     .     icount,' SATELLITES   PRN',' 4 DATA TYPES:',(dattyp(j),j=1,4)
     .        ,(obtypx(i),i=1,4)
      do i=1,icount
         prn = satarray(i) 
* MOD TAH 190702: Added antpwr to snav_read call
         call svnav_read( -1,iyear,iday,ih,imin,gnss,prn,nsn
     .                  , frqchn,antbody,sbmass,yawbias,yawrate, antpwr
     .                  , svnstart,svnstop )
         write(satnam(i),'(a1,i2,1x,i3,1x,a8)') 
     .           gnss,prn,nsn,antbody(7:14)
c         write(lu,215)  i,prn,gnss_name,nsn
c     1               ,  (lambda(i,j),j=1,ndat)
c  215    format(' CHANNEL ',I2,'  PRN# ',I2,2X,a8,i2,2X,'LAMBDA'
c     1         ,  7I3 ) 
         write(lu,'(a,i2,a,a16,a,7i3)')  ' CHANNEL ',i,'  SV  '
     .         ,satnam(i),'LAMBDA ',(lambda(i,j),j=1,ndat)
       enddo

c       Write out the first epoch information
                                              
         iy2 = mod(iyear,100)
         write(lu,220) x_time_flag,gpstutc,iy2,iday,ih,imin,seconds
     .               , interval, ircint,isessn
  220    format(/,' FIRST EPOCH (',a4,')    GPST-UTC= ',f4.1,/
     . ,' YR DAY HR MN SECS   INTERVAL(SECS)   DATA INTERVAL  SESSION'
     . ,/,1X,I2,1X,I3,1X,I2,1X,I2,1X,F6.3,1X,I6,11x,i6,11x,i2/)   
       if( ircint.eq.0 ) call report_stat('WARNING','MAKEX','wxhead',' '
     .      ,'No sampling interval on RINEX header, possible problem in 
     .AUTCLN if different from X-file sampling',0)
  
C       Write out the number of epochs

      write(lu,240) nepoch
  240 format(i4,' EPOCHS' )


c       Write a blank line where the number of bands and type of
c       observvations used to be - now specified by the lambda array

      ALINE1 = '                                      '
      WRITE(lu,'(A80)') ALINE1


c       Write another blank line

      write(lu,'(a80)') aline1


c       Write the expanded data format indicator

      write(lu,250) exptyp
  250 format(' EXPANDED DATA FORMAT ',a9)


C       Write the general phase/pseuorange header line

      WRITE(lu,260)
 260  FORMAT(' EPOCH DCB IER CHN  L1 PHASE',11X,'AMP',7X,
     1'L2 PHASE',9X,'AMP',6X,'L1 PSEUDORANGE',11X,
     2'L2 PSEUDORANGE')

      return
      end

