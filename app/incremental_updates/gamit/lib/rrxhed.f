      Subroutine rrxhed ( debug,gnss
     .   , rxver,rxpgm,rxusr,rxdat,rxcom,irxcom,rxmrk,rxobs,rxagy
     .   , rcvnum,rctype,rcvers,antnum,anttyp,apx,apy,apz
     .   , anth,ante,antn,nwave1,nwave2,nobtyp,rxobtyp,rxint,rxtime
     .   , irxyr0,irxmo0,irxdy0,irxhr0,irxmn0,rxsec0
     .   , irxyr1,irxmo1,irxdy1,irxhr1,irxmn1,rxsec1 )

c     Read a RINEX header on logical unit urinex (in includes/makex.h),
c     kurt feigl for RINEX Version 1
c     modified by r.king   19 Jan 91 for RINEX Version 2
c     modified by r.king    1 Mar 01 for RINEX Version 2.10 
c     modified by r.king   21 Dec 07 for RINEX Version 2.11 (no changes, just warning)
c     modified by r.king   22 Sep 15 for RINEX Version 3 (3.01 and 3.02 same)
c     modified by r.king   19 Feb 16 for RINEX Version 2.12 observables (lib/sel_obtyp.f only)
c     modified by M. Floyd 16 Oct 20 for RINEX Version 3.04 (no changes, just warning)
c     modified by M. Floyd 15 Apr 21 for RINEX Version 3.05 (no changes, just warning)

c     assumed open

      implicit none

      include '../includes/dimpar.h' 
      include '../includes/makex.h'

      logical debug 

c     variables for status reporting
      integer*4 ioerr,nerr,len,rcpar
      character*80 fname,prog_name
      character*256 message

c     strings for buffering data
      character*20 buff20
      character*80 buff80,blnk80

c     L1-only flag
      logical l1_only,ok_obtyp  

c     satellite system requested
      character*1 gnss 

c     RINEX defined items
c     file type and version number
      character*1  filtyp
      integer itype
      real*4 rxver
c     satellite system
      character*1 irxsvs
c     file identifiers
      character*20 rxpgm,rxusr,rxdat
c     comment
      character*60 rxcom(maxlin)
c     mark name
      character*60 rxmrk
c     mark number
      character*20 rxmrkn
c     observer
      character*20 rxobs
c     agency
      character*40 rxagy
c     receiver unit number,type, and SW version
      character*20 rcvnum,rctype,rcvers
c     antenna serial number and type
      character*20 antnum,anttyp
c     aproximate coordinates
      real*8 apx,apy,apz
c     antenna offsets
      real*8 anth,ante,antn
c     wavelength factors
      integer*4 nwave1,nwave2
c     different wavelength factors for satellites 
      integer nsvwave,isvwave(maxsat),nw,nw1,nw2
     .       ,iwave1(maxsat),iwave2(maxsat) 
c     flag for receiver clock offsets applied
      integer*4 irxrcvclk
c     SV system and observation types
      integer*4 nobtyp 
      character*1 asvid 
      character*3 rxobtyp(maxobt)
c     data interval in seconds
      real*8 rxint
c     data start time
      integer irxyr0,irxmo0,irxdy0,irxhr0,irxmn0
      real*8 rxsec0                             
      character*3 rxtime
c     data stop time
      integer irxyr1,irxmo1,irxdy1,irxhr1,irxmn1
      real*8 rxsec1
c     number of satellites
      integer irxnsv
c     satellite names (number designation)
      character*1 rxstyp(maxsat)
      integer irxsvn(maxsat)
c     observations for each satellite
      integer irxnob(maxsat,maxobt)

      integer irxcom,irxsat,n1,n2,i,j


c  Get calling program name and RINEX file name for report_stat

       len = rcpar(0,prog_name)
       inquire( unit=urinex, name=fname, iostat=ioerr )
       if(ioerr.ne.0) call report_stat('WARNING',prog_name,'lib/rrxhed'
     .   ,' ','Cannot get name of RINEX file for error reporting',ioerr)
     
c  Initialize buffers to avoid using info from previous file if line missing
      write (blnk80,'(80x)')
      write (rxpgm ,'(20x)')
      write (rxusr ,'(20x)')
      write (rxmrk ,'(20x)')
      write (rxobs ,'(20x)')
      write (rxagy ,'(40x)')
      write (rcvnum,'(20x)')
      write (rctype,'(20x)')
      write (rcvers,'(20x)')
      write (anttyp,'(20x)')
      write (antnum,'(20x)')
      apx = 0.d0
      apy = 0.d0
      apz = 0.d0
      anth = 0.d0
      ante = 0.d0
      antn = 0.d0
      nwave1 = 0
      nwave2 = 0
      irxyr0 = 0
      irxmo0 = 0
      irxdy0 = 0
      irxhr0 = 0
      irxmn0 = 0
      rxsec0 = 0.d0
      irxyr1 = 0
      irxmo1 = 0
      irxdy1 = 0
      irxhr1 = 0
      irxmn1 = 0
      rxsec1 = 0.d0
      nsvwave = 0 
      do i=1,maxobt
       rxobtyp(i) = '   '
      enddo    
      rxtime = '   '

c     set RINEX version to 1 until read, for initial check of end-of-header
      rxver = 1.0
      irxcom = 0
      irxsat = 1 
      rxint = 0.
      nerr = 0              
      

c rwk 060707: Read first line of RINEX header separately to guard against a bogus blank first line
      read( unit=urinex,iostat=ioerr,fmt='(a)') buff80
      if( buff80(61:73).eq.'rinex version' ) then
          call report_stat('WARNING',prog_name,'rrxhed',fname
     .                     ,'First line of RINEX lower case',ioerr) 
      endif
      call uppers(buff80)       
      if( buff80(61:73).ne.'RINEX VERSION' )  then 
           call report_stat('FATAL',prog_name,'rrxhed',fname
     .                     ,'Bogus first line of RINEX file',ioerr) 
      else
c          If the version number is an integer written under version 1 format
c          ( I6 ) it will not be read correctly by the version 2 format ( f9.2)
c          so read the version number free-format before reading the rest of
c          the line.  rwk 060526
          read( buff80(1:9),*,iostat=ioerr) rxver  
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'rrxhed'
     .      ,fname,'Error reading RINEX version number ',ioerr) 
          read (buff80,'(20x,a1,19x,a1,19x)',iostat=ioerr) filtyp,irxsvs
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'rrxhed'
     .      ,fname,'Error reading RINEX file type ',ioerr) 
      endif


 10   continue                
         buff80 = ' '
         read (unit   = urinex,
     .         iostat = ioerr,
     .         fmt    = '(a)') buff80
         BUFF20 = buff80(61:80)       
         call uppers(BUFF20)
            if (ioerr .ne. 0) then
                call report_stat('FATAL',prog_name,'rrxhed',fname
     .                     ,'Error reading header for RINEX file',ioerr)
            else if (
     .         (rxver.le.1.1.and. buff20.eq.'                    ') .or.
c              set to accept END OF HEADER with RINEX 1 (non-standard)
c     .         (irxver.eq.2 .and. buff20(:13).eq.'END OF HEADER') .or.  
     .         ( buff20(:13).eq.'END OF HEADER') .or.
     .         ioerr.eq.-1 ) then
                goto 20


            else if (BUFF20(1:19) .eq. 'PGM / RUN BY / DATE') then
               read (buff80,'(3a20)',iostat=ioerr)
     .         rxpgm,rxusr,rxdat

            else if (BUFF20(1:7) .eq. 'COMMENT') then
               read (buff80,'(a60)',iostat=ioerr)
     .         rxcom(irxcom+1)
*              MOD TAH 980115: +1 added to irxcom in if to ensure no bound
*                              violation.
               if( irxcom+1.lt.maxlin ) irxcom = irxcom + 1

            else if (BUFF20(1:11) .eq. 'MARKER NAME') then
               read (buff80,'(a60)',iostat=ioerr)
     .         rxmrk

            else if (BUFF20(1:13) .eq. 'MARKER NUMBER') then
               read (buff80,'(a20)',iostat=ioerr)
     .         rxmrkn

            else if (BUFF20(1:17) .eq. 'OBSERVER / AGENCY') then
               read (buff80,'(a20,a40)',iostat=ioerr)
     .         rxobs,rxagy

            else if (BUFF20(1:19) .eq. 'REC # / TYPE / VERS') then
               read (buff80,'(3a20)',iostat=ioerr)
     .         rcvnum,rctype,rcvers

            else if (BUFF20(1:12) .eq. 'ANT # / TYPE') then
               read (buff80,'(2a20)',iostat=ioerr)
     .         antnum,anttyp

            else if (BUFF20(1:19) .eq. 'APPROX POSITION XYZ') then
               read (buff80,'(3f14.4)',iostat=ioerr)
     .         apx,apy,apz

            else if (BUFF20 .eq. 'ANTENNA: DELTA H/E/N') then
               read (buff80,'(3f14.4)',iostat=ioerr)
     .         anth,ante,antn

            else if (BUFF20 .eq. 'WAVELENGTH FACT L1/2') then
c              check to see if different factors for different SVs  
               read(buff80(18:18),'(i1)') nw
               if( nw.ne.0 ) then  
c** rwk 190805: Since GAMIT won't support different wavelength factors for
c** different SVs, skip these lines
                 go to 10
                 nw1 = nsvwave + 1
                 nw2 = nw1 + nw - 1
                 if( nw2.gt.maxsat ) call report_stat('FATAL',prog_name
     .              ,'lib/rrxhed',fname,'Bad WAVELENGTH FACT',0)
                 read (buff80,'(2i6,i6,7(3x,a1,i2))',iostat=ioerr)  
     .              n1,n2,nw,(isvwave(i),i=nw1,nw2) 
                 do i=nw1,nw2 
                   iwave1(i) = n1
                   iwave2(i) = n2
                 enddo
                 nsvwave = nw2
               else
                read (buff80,'(2i6,i6)',iostat=ioerr) n1,n2
                nwave1 = n1
                nwave2 = n2
               endif
                      
c           RINEX 1/2 will have 'TYPES OF OBSERV', RINEX 3 will have 'OBS TYPES'
c           so we don't have to check which version here
            else if (buff20(1:19) .eq. '# / TYPES OF OBSERV') then
               read (buff80,'(i6)',iostat=ioerr) nobtyp     
               if( nobtyp.gt.maxobt ) then   
                 write(message,'(a,i2,a,i2)') '# obs types =',nobtyp
     .             ,' > maxobt = ',maxobt  
                 call report_stat('FATAL',prog_name,'lib/rrxhed',fname
     .                            ,message,0)
               elseif( nobtyp.le.9 ) then      
                 read (buff80,'(i6,9(4x,a2))',iostat=ioerr)
     .              nobtyp,(rxobtyp(j)(1:2),j=1,nobtyp)   
               else       
                 read (buff80,'(i6,9(4x,a2))',iostat=ioerr)
     .              nobtyp,(rxobtyp(j)(1:2),j=1,9)     
c                nobtyp > 9, continue to next line    
                 if( nobtyp.gt.9 ) then
                    read (unit=urinex,iostat=ioerr,fmt='(a)') buff80
cd                    print *,'buff80' ,buff80
                   if( ioerr.ne.0.or.buff80(65:72).ne.'TYPES OF' ) 
     .                call report_stat('FATAL',prog_name,'lib/rrxhed'
     .               ,fname,'Unexpected TYPES OF OBSERV line 2',ioerr)
* MOD TAH 120125: Read up to 18 first and then the remainder
                   read (buff80,'(6x,9(4x,a2))',iostat=ioerr)
     .                (rxobtyp(j)(1:2),j=10,min(18,nobtyp))       
                   if( ioerr.ne.0 ) 
     .                call report_stat('FATAL',prog_name,'lib/rrxhed'
     .                  ,fname,'Error reading TYPES OF OBSERV line 2'
     .                  ,ioerr)
                   endif  
* MOD TAH 120125   See if there is still more (19-nobtyp)
                   if( nobtyp.gt.18 ) then  ! Read next line
                      read (unit=urinex,iostat=ioerr,fmt='(a)') buff80
                      if( ioerr.ne.0.or.buff80(65:72).ne.'TYPES OF' ) 
     .                call report_stat('FATAL',prog_name,'lib/rrxhed'
     .                    ,fname,'Unexpected TYPES OF OBSERV line 2',
     .                     ioerr)                      
                      read (buff80,'(6x,9(4x,a2))',iostat=ioerr)
     .                   (rxobtyp(j)(1:2),j=19,min(27,nobtyp))       
                      if( ioerr.ne.0 ) 
     .                   call report_stat('FATAL',prog_name,'lib/rrxhed'
     .                     ,fname,'Error reading TYPES OF OBSERV line 3'
     .                     ,ioerr)
                    end if
               endif
c              observables are required to be uppercase, but at least one old version (1.7)
c              of TRRINEXO produced lowercase ids, so put this in for safety (Weiping Jiang/RWK 031016)
               do i=1,nobtyp
                 call uppers(rxobtyp(i))
               enddo                                         
               if(debug) then
                 write(*,*) 'RRXHED rxver nwave1 nwave2 nobtyp rxobtyp '
     .             , rxver,nwave1,nwave2,nobtyp,(rxobtyp(i),i=1,nobtyp)
               endif 
            
            elseif ( buff20(1:19).eq.'SYS / # / OBS TYPES' .and.
     .               buff80(1:1).eq.gnss  ) then
              read (buff80,'(a1,i5)',iostat=ioerr) asvid,nobtyp     
cd              print *,'RRXHED buff80 ',buff80
cd              print *,'gnss asvid nobtyp ',gnss,asvid,nobtyp
cd              print *,' decoding '
              if( nobtyp.gt.maxobt ) then   
                write(message,'(a,i2,a,i2)') '# obs types =',nobtyp
     .            ,' > maxobt = ',maxobt  
                call report_stat('FATAL',prog_name,'lib/rrxhed',fname
     .                          ,message,0)
              elseif( nobtyp.le.13 ) then
                read (buff80,'(1x,i5,13(1x,a3))',iostat=ioerr)
     .             nobtyp,(rxobtyp(j),j=1,nobtyp)   
              else 
                read (buff80,'(1x,i5,13(1x,a3))',iostat=ioerr)
     .            nobtyp,(rxobtyp(j),j=1,13)     
c               nobtyp > 13, continue to next line    
                if( nobtyp.gt.13 ) then
                  read (unit=urinex,iostat=ioerr,fmt='(a)') buff80
cd                  print *,'buff80' ,buff80
                  if( ioerr.ne.0.or.buff80(71:79).ne.'OBS TYPES' ) 
     .              call report_stat('FATAL',prog_name,'lib/rrxhed'
     .             ,fname,'Unexpected OBS TYPES line 2',ioerr)
c                   Read up to 26 first and then the remainder
                  read (buff80,'(6x,13(1x,a3))',iostat=ioerr)
     .                (rxobtyp(j),j=14,min(26,nobtyp))       
                  if( ioerr.ne.0 ) 
     .              call report_stat('FATAL',prog_name,'lib/rrxhed'
     .                 ,fname,'Error reading OBS TYPES line 2'
     .                 ,ioerr)
                endif
c               see if there are still more (27-nobtyp)
                if( nobtyp.gt.26 ) then  ! Read next line
                  read (unit=urinex,iostat=ioerr,fmt='(a)') buff80
                  if( ioerr.ne.0.or.buff80(71:79).ne.'OBS TYPES' ) 
     .               call report_stat('FATAL',prog_name,'lib/rrxhed'
     .                     ,fname,'Unexpected OBS TYPES line 2',
     .                     ioerr)                      
                   read (buff80,'(6x,13(1x,a3))',iostat=ioerr)
     .                  (rxobtyp(j),j=27,min(27,nobtyp))       
                   if( ioerr.ne.0 ) 
     .                call report_stat('FATAL',prog_name,'lib/rrxhed'
     .                  ,fname,'Error reading OBS TYPES line 3'
     .                  ,ioerr)
                endif  
              endif
              if(debug) then
                 write(*,*) 'RRXHED rxver nwave1 nwave2 nobtyp rxobtyp '
     .            , rxver,nwave1,nwave2,nobtyp,(rxobtyp(i),i=1,nobtyp)
              endif 

            else if (BUFF20(1:8) .eq. 'INTERVAL') then  
               read (buff80,'(f10.0)',iostat=ioerr)
     .           rxint   

            else if (BUFF20(1:17) .eq. 'TIME OF FIRST OBS') then
               read (buff80,'(5i6,f13.7,5x,a3)',iostat=ioerr)
     .         irxyr0,irxmo0,irxdy0,irxhr0,irxmn0,rxsec0,rxtime
               call fix_y2k(irxyr0)
                
            else if (BUFF20(1:16) .eq. 'TIME OF LAST OBS') then
               read (buff80,'(5i6,f12.6)',iostat=ioerr)
     .         irxyr1,irxmo1,irxdy1,irxhr1,irxmn1,rxsec1
               call fix_y2k(irxyr1)
  
            else if (BUFF20(1:16) .eq. '# OF SATELLITES') then
               read (buff80,'(i6)',iostat=ioerr)
     .         irxnsv

            else if (BUFF20(1:14) .eq. 'PRN / # OF OBS') then
               read (buff80,'(3x,a1,i2,9i6)',iostat=ioerr)
     .         rxstyp(irxsat),irxsvn(irxsat)
     .      , (irxnob(irxsat,i),i=1,nobtyp)
               if( irxsat.lt.maxsat ) irxsat = irxsat + 1
            
            elseif (BUFF20(1:19) .eq. 'RCV CLOCK OFFS APPL' ) then
c             warn but do nothing more for this entry
              read (buff80,'(2x,i4)') irxrcvclk
              if( irxrcvclk.ne.0 ) 
     .           call report_stat('WARNING',prog_name,'lib/rrxhed'
     .                ,fname,'Receiver clock offsets applied',0)
            endif
                            
            goto 10

  20       continue
           
c** debug
cd         do i=1,nobtyp
cd            print *,i,rxobtyp(i)
cd         enddo
     
c    RINEX 1 or 2, change special-case observation types
      if( rxver.lt.3.0 ) then
        l1_only = .true.        
        do itype=1,nobtyp
c         P2 called C2 for serial P-code Trimble  
c         rwk mod 060111: Do this only for old observations since there is now a true C2
c         if( rxobtyp(itype).eq.'C2 ' ) rxobtyp(itype) = 'P2 '
          if( rxobtyp(itype).eq.'C2 ' .and.irxyr0.lt.2004 )  
     .        rxobtyp(itype) = 'P2 '
          if( rxobtyp(itype).eq.'L2 ' ) l1_only = .false.
        enddo                          
        if( l1_only ) call report_stat('WARNING',prog_name
     .     ,'lib/rrxhed',' ','Single frequency observations',0)   
      endif   

c     RINEX 2 or 3 mislabeled Beidou second frequency (B2 but RINEX C7/L7/S7)
      if( gnss.eq.'C' ) then
        do itype=i,nobtyp
          if( rxobtyp(itype)(2:2).eq.'2' ) rxobtyp(itype)(2:2) = '7'
           call report_stat('WARNING',prog_name,'lib/rrxhed',fname
     .    ,'Beidou B2 obs changed from C2/L2 to RINEX standard C7/L7',0)
        enddo
      endif
          
c     if multiple wavelength factors input, set single values  (RINEX 2 only?) 
      if( nsvwave.gt.0 ) then  
         nwave1 = 1
         nwave2 = 1
         do i=1,nsvwave
           if( iwave1(i).ne.nwave1 ) then  
c             L1 should always be 1 (except Macrometer II) 
             call report_stat('WARNING',prog_name,'lib/rrxhed',fname
     .,'SV L1 wavelength factor NE 1-something wrong in RINEX header',0)
           endif
           if( iwave2(i).ne.nwave2 ) then
             call report_stat('WARNING',prog_name,'lib/rrxhed',fname
     .   ,'At least one SV L2 wavelength factor =2, set all =2',0)  
             iwave2(i) = 2
           endif
         enddo
         if( iwave2(1).eq.2 ) nwave2 = 2
      endif

c     check for reasonability
      if (rxver .gt. 3.05 ) then
         write (message,'(a,f9.2,a)') 'RINEX version = ',rxver
     .                    ,' greater than 3.05'
         call report_stat('WARNING',prog_name,'lib/rrxhed',fname
     .                   , message,0)

      endif

      if (filtyp .ne. 'O') then
         write (message,'(a,a1)')  'File type not O ',filtyp
         call report_stat('WARNING',prog_name,'lib/rrxhed',fname
     .          ,message,0)
      endif

      if ( anth.eq. 0.d0 ) then
         call report_stat('WARNING',prog_name,'lib/rrxhed',fname
     .                   ,'Antenna height is zero',0)

      endif

      if (anth .gt. 2.d0 ) then
         write (message,'(a,d12.3)') 'Antenna height > 2 meters: ',anth
         call report_stat('WARNING',prog_name,'lib/rrxhed',fname
     .                    ,message,0)

      endif
                                   
      if ( rxver.lt.3.0 .and. (nwave1.lt.1.or.nwave1.gt.2) )  then
         write (message,'(a,i2)') 'Invalid nwave1 ',nwave1
         call report_stat('WARNING',prog_name,'lib/rrxhed',fname
     .                   , message,0)   
c        add special fix for column misalignment
         if( nwave1.gt.0.and.nwave1.ne.2 ) then 
           call report_stat('WARNING',prog_name,'lib/rrxhed',fname
     .      ,'  Setting nwave1=1 ',0)  
           nwave1 = 1 
         endif
      endif

      if (.not.l1_only .and. rxver.lt.3.0 .and. 
     .    (nwave2.lt.1.or.nwave2.gt.2) ) then
         write (message,'(a,i2)') 'Invalid nwave2 ',nwave2
         call report_stat('WARNING',prog_name,'lib/rrxhed',fname
     .                   , message,0)  

c        add special fix for column misalignment
         if( rxver.lt.3.0 .and. nwave1.gt.0.and.nwave1.ne.2 ) then 
           call report_stat('WARNING',prog_name,'lib/rrxhed',fname
     .      ,'  Setting nwave1=1 ',0)  
           nwave1 = 1 
         endif
      endif

      if (dabs(ante) .gt. 0.002) then
         write (message,'(a,f10.1)') 'East antenna offset > 2 mm: '
     .                               ,ante*1000
        call report_stat('WARNING',prog_name,'lib/rrxhed',fname
     .                   , message,0)

      endif

      if (dabs(antn) .gt. 0.002) then
         write (message,'(a,f10.1)') 'North antenna offset > 2 mm: '
     .                                ,antn*1000
        call report_stat('WARNING',prog_name,'lib/rrxhed',fname
     .                   , message,0)
      endif      
                  
      if( rxver.lt.3.0 .and. rctype(1:15).eq.'ASHTECH GG-XXIV' ) then 
        ok_obtyp = .true.
        do itype=1,nobtyp
          if(rxobtyp(itype).eq.'P1 ') ok_obtyp = .false.
          if(rxobtyp(itype).eq.'P2 ') ok_obtyp = .false.
          if(rxobtyp(itype).eq.'L2 ') ok_obtyp = .false.
          if(rxobtyp(itype).eq.'D2 ') ok_obtyp = .false.
        enddo
        if( .not.ok_obtyp) call report_stat('WARNING',prog_name
     .    ,'lib/rrxhed',fname
     .  ,'Rcvr type indicates L1-only but obs type P or L2 in header',0)
      endif                                              

      if( rxtime.eq.'   ' ) rxtime = 'GPS'
                        
      return
      end

