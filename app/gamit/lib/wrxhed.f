      subroutine wrxhed ( lu,gnss,
     .   rxver,rxpgm,rxusr,rxdat,rxcom,irxcom,rxmrk,rxobs,rxagy,
     .   rcvnum,rctype,rcvers,antnum,anttyp,apx,apy,apz,
     .   anth,ante,antn,nwave1,nwave2,nobtyp,rxobtp,rxint,rxtime,
     .   irxyr0,irxmo0,irxdy0,irxhr0,irxmn0,rxsec0,
     .   irxyr1,irxmo1,irxdy1,irxhr1,irxmn1,rxsec1,
     .   numsv,totsv,numobs )

      implicit none
           
c     --dimensions (maxsat)
      include '../includes/dimpar.h'  
c     --dimensions (maxchn,maxobt,,maxlin)
      include '../includes/makex.h'
               
c     write RINEX header on logical unit lu, assumed open
      integer lu

c     strings for buffering data
      character*16  buff16

c     RINEX defined items    
      character*1 gnss 
      real*4 rxver
      character*20 rxpgm,rxusr,rxdat
c     comment
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
      character*20 antnum
      character*20 anttyp
c     aproximate coordinates
      real*8 apx,apy,apz
c     antenna offsets
      real*8 anth,ante,antn
c     wavelength factors
      integer nwave1,nwave2
c     observation types
      integer nobtyp
      character*3 rxobtp(maxobt)                          
c     time type 
      character*3 rxtime
c     data interval in seconds
      real*8 rxint
c     data start time
      integer irxyr0,irxmo0,irxdy0,irxhr0,irxmn0
      real*8 rxsec0
c     data stop time
      integer irxyr1,irxmo1,irxdy1,irxhr1,irxmn1
      real*8 rxsec1  
c     number of SVs 
      integer*4 numsv
c     satellite system (use as non-input scalar, set blank, for now) 
      character*1 satsys 
c     satellite numbers
      integer*4 totsv(maxsat)
c     number of obs of each type for each satellite
      integer*4 numobs(maxobt,maxsat)

c     --local variables
      integer irxcom,ioerr,len,rcpar,i,j  
      character*3 tmptyp(18)
      character*60 buff60
      character*80 prog_name
      character*256 message
         
 
c     get the calling module name for report_stat
      len = rcpar(0,prog_name)

      buff16 = 'OBSERVATION DATA'
      write (lu,'(f9.2,11x,a16,24x,a)',iostat=ioerr)
     .rxver,buff16,'RINEX VERSION / TYPE'

      write (lu,'(3a20,A)',iostat=ioerr)
     .rxpgm,rxusr,rxdat,'PGM / RUN BY / DATE'
                                               
      if( irxcom.gt.maxlin ) then
         write(message,'(a,i3,a)') 
     .        'Number of comment lines exceeds maxlin in makex.h (='
     .        , maxlin,')'
         call report_stat('WARNING',prog_name,'lib/wrxhed',' '
     .                   , message,0)
         irxcom = maxlin
      endif
      do  i = 1,irxcom   
         write (lu,'(a60,A)',iostat=ioerr)
     .   rxcom(i),'COMMENT'
      enddo
          
      write (lu,'(a60,a)',iostat=ioerr)
     .rxmrk,'MARKER NAME'
                             
      write (lu,'(a20,a40,A)',iostat=ioerr)
     .rxobs,rxagy,'OBSERVER / AGENCY'
                    
      write (lu,'(3a20,A)',iostat=ioerr)
     .rcvnum,rctype,rcvers,'REC # / TYPE / VERS'
       
      write (lu,'(2a20,20x,A)',iostat=ioerr)
     .antnum,anttyp,'ANT # / TYPE'

      write (lu,'(3f14.4,18X,A)',iostat=ioerr)
     .apx,apy,apz,'APPROX POSITION XYZ'

      write (lu,'(3f14.4,18X,A)',iostat=ioerr)
     .anth,ante,antn,'ANTENNA: DELTA H/E/N'

      write (lu,'(2i6,48X,A)',iostat=ioerr)
     .nwave1,nwave2,'WAVELENGTH FACT L1/2'
                        
      if( nobtyp.gt.maxobt ) then    
        write(message,'(a,i2,a,i2,a)') 'Number of observables = '
     .    ,nobtyp,' exceeds maxobt in makex.h (=',maxobt,')'
        call report_stat('FATAL',prog_name,'lib/wrxhed',' '
     .                   , message,0)
      endif
      if( nobtyp.gt.18 ) then  
        call report_stat('FATAL',prog_name,'lib/wrxhed',' '
     .      ,'Only 18 observables supported',0)
      endif
      do i = 1,nobtyp
        tmptyp(i) = rxobtp(i)
      enddo
      do i = nobtyp + 1,18
         tmptyp(i) = '   '
      enddo       
      if( rxver.lt.3. ) then
        write (lu,'(i6,9(4x,a2),A)',iostat=ioerr)
     .    nobtyp,(tmptyp(j)(1:2),j=1,9),'# / TYPES OF OBSERV'
        if( nobtyp.gt.9 ) then
          write (lu,'(6x,9(4x,a2),A)',iostat=ioerr)
     .    (tmptyp(j)(1:2),j=10,18),'# / TYPES OF OBSERV'
        endif
      else
        write(lu,'(a1,13(1x,a3),a)',iostat=ioerr) 
     .   gnss,(tmptyp(i),i=1,13),'SYS / # / OBS TYPES'
      
	if( nobtyp.gt.13 ) then
	  write(lu,'(a1,13(1x,a3),a)',iostat=ioerr) 
     .       gnss,(tmptyp(i),i=1,13),'SYS / # / OBS TYPES'
        endif
      endif
      
      write (lu,'(f10.3,50x,a)',iostat=ioerr)
     . rxint ,'INTERVAL'
         
      write (lu,'(5i6,f12.6,18X,A)',iostat=ioerr)
     .irxyr0,irxmo0,irxdy0,irxhr0,irxmn0,rxsec0,
     .'TIME OF FIRST OBS'
      
c     this entry is optional, use month to test for valid time            
      if( irxmo1.ne.0 ) then   
        write (lu,'(5i6,f12.6,18X,A)',iostat=ioerr)
     .        irxyr1,irxmo1,irxdy1,irxhr1,irxmn1,rxsec1,
     .        'TIME OF LAST OBS'
      endif  

c     this entry is also optional, test for validity
      if( numsv.ne.0 ) then 
         satsys = ' '
         write (lu,'(i6,54x,a)',iostat=ioerr)
     .          numsv,'# OF SATELLITES'
         do i = 1,numsv                
           buff60 = ' '
           write (buff60,'(3x,a,i2,9i6)') 
     .           satsys,totsv(i),(numobs(j,i),j=1,nobtyp) 
           write(lu,'(a60,a)') 
     .           buff60,'PRN / # OF OBS'
          enddo
      endif

      write (lu,'(60x,a19)')'END OF HEADER      '

      return
      end


