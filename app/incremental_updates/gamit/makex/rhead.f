      Subroutine rhead( debug,gnss,rcvrsw,swver,rxver,rxpgm,rxtime
     .                , nobtyp,rxobtyp,nwave1,nwave2,ircint
     .                , nxhdlin,xhdlin )

c     Read the RINEX header to get information for the X-file header
c       R.W. King  20 Feb 1989

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'

      logical      debug, c2flag
      integer      nblen
                                   
c     Satellite system
      character*1 gnss

C     rcvr software name
      character*3  rcvrsw
      character*16 auser
C     user-input site code
     .,            asite
C     reciever serial number
     .,            arcvr
C     operator input antenna ht.
     .,            antht
C     receiver name
     .,            rcvrnm
C     number and strings for x-file header lines 
      integer*4 nxhdlin
      character*80 xhdlin(maxlin)

C     rcvr software version from input
      real*4       swver       
C     atmospheric pressure
      real*8       apress
C     atmospheric temperature
     .,            atemp
C     atmospheric pressure
     .,            ahumid
c     collection interval
     .,            sample_int
C     rcvr software version from rcvr
     .,            swveru
C     L1 phase center wrt mark
     .,            offsl1(3)
C     L2 phase center wrt mark
     .,            offsl2(3)
C     mult. factor for L1,L2 phases
      integer*4    l1fact,l2fact
c     receiver collection interval
      integer*4    ircint

c     RINEX defined items
      real*4 rxver
      character*20 rxpgm,rxusr,rxdat
c     comment
      integer irxcom
      character*60 rxcom(maxlin)
c     mark name
      character*60 rxmrk
c     observer
      character*20 rxobs
c     agency
      character*40 rxagy
c     receiver serial number, type and SW version
      character*20 rcvnum
      character*20 rctype
      character*20 rcvers
c     antenna serial number and type
      character*20 antnum,anttyp
c     aproximate coordinates
      real*8 apx,apy,apz
c     antenna offsets
      real*8 anth,ante,antn
c     wavelength factors
      integer*4 nwave1,nwave2
c     observation types
      integer nobtyp
      character*3 rxobtyp(maxobt)
c     data interval in seconds
      real*8 rxint
c     time type (always 'GPS' for mixed files)
      character*3 rxtime
c     data start time                                             
      integer irxyr0,irxmo0,irxdy0,irxhr0,irxmn0
      real*8 rxsec0
c     data stop time
      integer irxyr1,irxmo1,irxdy1,irxhr1,irxmn1
      real*8 rxsec1
c     report_stat routine message
      character*256 message  
c     other variables
      integer*4 i

c     Blank out the header to avoid writing nulls, interpreted by UNIX as binary

      do i=1,maxlin
        xhdlin(i) = ' '
      enddo        
      asite = ' '
      arcvr = ' '               
      c2flag = .false. 

c     Determine the receiver type from the input software type and version
                   
         call uppers(rcvrsw)
         rcvrnm = ' ' 
         if( rcvrsw.eq.'GES' .or. rcvrsw.eq.'COR' .or.
     .       rcvrsw.eq.'ROM' ) rcvrnm='TI 4100'
         if( rcvrsw.eq.'MAC' ) rcvrnm='MACROMETER II'
         if( rcvrsw.eq.'MIN' ) rcvrnm='MINI-MAC 2816' 
         if( rcvrsw.eq.'LEI' ) then
            if( swver.lt.5.0 ) then
               rcvrnm='Leica system 200' 
            else
               rcvrnm='Leica'
            endif 
         endif
         if( rcvrsw.eq.'TRM' ) then
            if (swver.lt.4.0)  rcvrnm='Trimble 4000 SDT'
            if (swver.ge.4.1)  rcvrnm='Trimble 4000 SST'
            if (swver.ge.5.5)  rcvrnm='Trimble 4000'
         endif
         write (uinfor,5) rcvrnm,rcvrsw,swver
         if( debug)  write (uscren,5) rcvrnm,rcvrsw,swver
   5     format (1x,'Receiver, Software from MAKEX Batch Controls: '
     .          , A16,3x,A6,f6.2)

c     Read the RINEX header

         call rrxhed ( debug,gnss
     .        , rxver,rxpgm,rxusr,rxdat,rxcom,irxcom,rxmrk,rxobs,rxagy
     .        , rcvnum,rctype,rcvers,antnum,anttyp,apx,apy,apz
     .        , anth,ante,antn,nwave1,nwave2,nobtyp,rxobtyp,rxint,rxtime
     .        , irxyr0,irxmo0,irxdy0,irxhr0,irxmn0,rxsec0
     .        , irxyr1,irxmo1,irxdy1,irxhr1,irxmn1,rxsec1)   
         if(debug)  
     .    print *,'RHEAD aft rrxhed rxver nwave1 nwave2 nobtyp rxobtyp '
     .          ,                  rxver,nwave1,nwave2,nobtyp,rxobtyp
                                    
c        For RINEX 2 with multiple GNSS, there is a problem for Glonass:
c        some stations label the low-frequency PR as P2 and omit C2 from
c        the header; others label it C2 but have P2 for GPS, so we have
c        to hide (by renaming) the P2 and select the C2.  This is kludge
c        and may not work for all cases.  We'll fix this with the new
c        x-file format, which will have the observables explicit at each epoch
* MOD TAH 200526: Block removed because not needed
C        if( rxver.lt.3.0 .and. gnss.eq.'R' ) then
C          do i=1,nobtyp
C            if(rxobtyp(i).eq.'C2 ') c2flag = .true.
C          enddo
C          if( c2flag ) then
C            do i=1,nobtyp
* MOD TAH 200526: Removed this "fix".  With sel_obtyp hierarchy C2 should 
*              get selected even if P2 present.
C              if(rxobtyp(i).eq.'P2 ' ) rxobtyp(i) = 'x2 '
C            enddo
C          endif 
C        endif  
         
c        get wavelength factors
         l1fact = nwave1
         l2fact = nwave2      
c** rwk 080522: If the wavelength factors is missing from the RINEX header
c               (not allowed by the standards, set them to 1)
         if( l1fact.eq.0 ) then
            call report_stat('WARNING','MAKEX','rhead',' '
     .        ,'Wavelength factors missing from RINEX header, set=1',0)
            l1fact = 1 
            l2fact = 1
         endif
     
c        Fill the xhdlin array with all the comment information in the X-file
c        First, a two-line header:
         write (xhdlin(1),'(1x,2a)') 'MADE FROM FILE: '
     .                             , frinex(1:nblen(frinex))
         nxhdlin = 2
         write (xhdlin(2),11)                      
         if (debug )  write (uscren,11)
         write (uinfor,11)
  11     format (1x,'HEADER FROM SOURCE RINEX FILE:')
c        Then the information from the RINEX header--first the comments
c        and then 7 lines of specific header information
         nxhdlin = nxhdlin + irxcom + 7
         if( nxhdlin.gt.maxlin ) then                               
           write(message,'(a,i3,a)') 'More than ',maxlin
     .      ,' comment lines for X-file header, omit some comments'
           call report_stat('WARNING','MAKEX','rhead',' ',message,0)
           irxcom = maxlin - 9
         endif
         if(debug) write(uscren,'(/,a,i3)') ' In RHEAD  irxcom=',irxcom
         do i=1,irxcom
           xhdlin(i+2) = rxcom(i)
         enddo
         nxhdlin = irxcom + 2

         write (uinfor,14) rxmrk
         write (uinfor,15) rxobs
         write (uinfor,16) rcvnum
         write (uinfor,20) anth,ante,antn
         write (uinfor,21) rctype,rcvers
         write (uinfor,22) anttyp,antnum
         write (uinfor,23) rxint
                
         if( debug ) then
           write (uscren,14) rxmrk
           write (uscren,15) rxobs
           write (uscren,16) rcvnum
           write (uscren,20) anth,ante,antn
           write (uscren,21) rctype,rcvers
           write (uscren,22) anttyp,antnum
           write (uscren,23) rxint
         endif

         write (xhdlin(nxhdlin+1),14) rxmrk
         write (xhdlin(nxhdlin+2),15) rxobs
         write (xhdlin(nxhdlin+3),16) rcvnum
         write (xhdlin(nxhdlin+4),20) anth,ante,antn
         write (xhdlin(nxhdlin+5),21) rctype,rcvers
         write (xhdlin(nxhdlin+6),22) anttyp,antnum
         write (xhdlin(nxhdlin+7),23) rxint
         nxhdlin = nxhdlin + 7

  14     format (1x,'Oper. site code  : ',a60)
  15     format (1x,'Operator name    : ',a20)
  16     format (1x,'Receiver serial #: ',a20 )
  20     format (1x,'Oper. ant. H,E,N : ',3f12.3)
  21     format (1x,'Receiver         : ',A20,'  Software: ',a20)
  22     format (1x,'Antenna          : ',A20,' Serial #: ',a20)
  23     format (1x,'Collection Interval: ',f10.3)
         ircint = int(rxint)

      return
      end
