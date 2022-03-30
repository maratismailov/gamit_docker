      Program ARGO2NAV

c     Translate NGS ARGO (CIGNET) orb format to RINEX nav format
                                                          
c     R. King    22 March 2000
                           
c     Runs from command-line:  argo2nav <argo-filename> 

      implicit none
       
c     --dimensions 
      include '../includes/argo2fic.h'

c     ARGO file items
         
      integer*4 iprn,iweek0,iyr,ihr,imin,ncnt
     .        , ittag(maxeph,7),locate(maxeph)
      real*8 sectag(maxeph),time(maxeph),ephr(maxeph,27)


c     RINEX defined items     

      real*4 rxver
      character*20 rxpgm,rxusr,rxdat
      real*8 bcephem(16),bcclock(6),subfr1(8)  
      real*8 sowclk,howsow


c     File names and units

      character*80 argofile,rinexfile
      integer uargo,urinex

      parameter (uargo=1)
      parameter (urinex=2)

c     Controls and other variables

      logical eof 
      integer*4 jyr,jmo,jdy,ios
      integer*4 iweek,iarg,iclarg,ioerr,idoy,jdoy,kdoy,i,j,k
      real*8 utcoff,sec
      character*15 buff15
      character*16 uname                   
      character*20 version  
      character*256 message

      data version/'ARGO2NAV version 1.0'/



***** Remove old versions of the status, warning, and error files

      call report_stat('CLEAR','ARG2NAV',' ',' ', ' ',0)
      call report_stat('CLEAR','LIB',' ',' ', ' ',0)
                             

***** Set the RINEX header info

c     -- RINEX version
      rxver = 2.0
c       identify ARGO2NAV version in 20 characters
      write(rxpgm,'(a20)') version
c     -- user name, agency, and run date
      call getusr(uname) 
      rxusr = '    '//uname    
      call getdat(jyr,jmo,jdy) 
      write (rxdat,'(i4.4,"/",i2.2,"/",i2.2)') jyr,jmo,jdy 

           
****** Read the name of the ARGO file and open it

      iarg = iclarg(1,argofile)
      if( iarg.le.0 ) then              
        write(6,'(1x)') 
        write(6,'(a)') 
     .     'ARGO2NAV requires ARGO file name as command-line input'
        stop
      endif
      open(unit=uargo,file=argofile,status='old'
     .    ,form='formatted',iostat=ioerr) 
      if (ioerr.ne.0) then 
        call report_stat('FATAL','ARGO2NAV','argo2mav',argofile
     .      ,'Cannot open input ARGO file',0) 
      else
        call report_stat('STATUS','ARGO2NAV','argo2mav',argofile
     .                       ,'Opening: ',0)
      endif

       
         
*****   Need to sort the ephemris records, so read all into storage

      eof = .false.
      ncnt = 1
      do while( .not.eof )
        call rdeph(uargo,ittag,sectag,ephr,ncnt,ios)
        if (ios .eq. -1) then
           ncnt = ncnt - 1
           close(uargo)   
           eof = .true.
        else
           if( ncnt.ge.maxeph ) then
             write(message,'(2a,i5)') 'Number of ephemeris records'
     .               ,' exceeds MAXEPH = ',maxeph
             call report_stat('FATAL','ARGO2NAV','argo2nav',' '
     .                       ,message,0)   
            else
              ncnt = ncnt + 1
            endif
        endif
      enddo  
c     set the initial week number from the orb file header
      iweek0 = ittag(1,1)
c     convert the time tags to week,sow for sorting
      do i=1,ncnt
        call fix_y2k(ittag(i,3))
        jdoy = idoy(ittag(i,3),ittag(i,4),ittag(i,5))
        call timcon(-4,iweek,sowclk
     .        , ittag(i,3),jdoy,ittag(i,6),ittag(i,7),sectag(i),utcoff )
c        print *,'prn y m d h m s w sow ',ittag(i,2),ittag(i,3)
c     . ,ittag(i,4),ittag(i,5),ittag(i,6),ittag(i,7),sectag(i),iweek,sowclk
c       if the data don't belong to this week flag time as negative
        if (iweek.eq.iweek0) then
          time(i) = sowclk
        else       
          write(message,'(a,i4,a,i4,a)') 'Epoch week (',iweek,
     .         ') not equal the initial week (',iweek0,')' 
          call report_stat('WARNING','ARGO2NAV','argo2nav',' '
     .                       ,message,0)   
          time(i) = -sowclk
        endif
      enddo     
      call indexx( ncnt,time,locate )
                

******  Begin loop over all ephemeris records
                                        
c     set 'last' doy to 0 to start
      kdoy = 0
      do i = 1,ncnt
        
        k = locate(i)                                                   
        sowclk = time(k)       
        if (sowclk.ge.0 ) then
          call timcon(4,iweek0,sowclk,iyr,jdoy,ihr,imin,sec,utcoff )
c
c         if new day, open a new RINEX file

          if( jdoy.ne.kdoy ) then   
            write (rinexfile,'(a4,i3.3,a1,".",i2.2,"n")')  
     .                       'brdc',jdoy,'0',mod(iyr,100)
            call lowers(rinexfile)  
            kdoy = jdoy
            open(unit =urinex,file=rinexfile,status='unknown'
     .          ,form='formatted',iostat = ioerr)
            if (ioerr.ne.0) then 
              call report_stat('FATAL','ARGO2NAV','argo2nav'
     .         ,rinexfile,'Cannot open output RINEX file',ioerr) 
            else   
              call report_stat('STATUS','ARGO2NAV','argo2nav'
     .                      ,rinexfile,'Opening: ',0)   
            endif
            kdoy = jdoy
            buff15= 'NAVIGATION DATA'
            write (urinex,'(f9.2,11x,a15,25x,a)',iostat=ioerr)
     .                  rxver,buff15,'RINEX VERSION / TYPE'
            write (urinex,'(3a20,A)',iostat=ioerr)
     .                   rxpgm,rxusr,rxdat,'PGM / RUN BY / DATE' 
            write (urinex,'(a,27x,a7,13x)') 
     .                'Broadcast elements from ARGO file','COMMENT'
            write (urinex,'(60x,a20)') 'END OF HEADER       '
          endif

c         map and write the navigation record
c         is howsow = sowclk?
          howsow = sowclk
          do j=1,16
            bcephem(j) = 0.d0
          enddo
          do j=1,6
            bcclock(j) = 0.d0
          enddo
          do j=1,6
            subfr1(j) = 0.d0      
          enddo          
c           xetoe             
          bcephem(1) = ephr(k,12)
c           xem0
          bcephem(2) = ephr(k,7)
c           xedn
          bcephem(3) = ephr(k,6)
c           xeart
          bcephem(4) = ephr(k,11)
c           xeecc
          bcephem(5) = ephr(k,9)
c           xei0
          bcephem(6) = ephr(k,16)
c           xeidt
          bcephem(7) = ephr(k,20)
c         xeom0
          bcephem(8) = ephr(k,14)
c           xeomd
          bcephem(9) = ephr(k,19)
c           xew
          bcephem(10) = ephr(k,18)
c           xecuc
          bcephem(11) = ephr(k,8)
c           xecus
          bcephem(12) = ephr(k,10)
c           xecrc   
          bcephem(13) = ephr(k,17)
c           xecrs
          bcephem(14) = ephr(k,5)
c         xecic
          bcephem(15) = ephr(k,13)
c           xecis
          bcephem(16) = ephr(k,15)
c           xetoc         
          bcclock(1) = sowclk
c           xeaf0
          bcclock(2) = ephr(k,1)
c           xeaf1
          bcclock(3) = ephr(k,2)
c           xeaf2
          bcclock(4) = ephr(k,3)
c           xeadc
          bcclock(5) = ephr(k,27)      
c         xetdg -- not used
c          bcclock(6) =   
c           cflgl2 
          subfr1(1) = ephr(k,21)
c           svaccr
          subfr1(2) = ephr(k,24)
c           svhlth
          subfr1(3) = ephr(k,25)
c           aodc
          subfr1(4) = ephr(k,27)
c           pflgl2
          subfr1(5) = ephr(k,23)
c           tgd
          subfr1(6) = ephr(k,26)
c           aode
          subfr1(7) = ephr(k,4)
c           howsow
          subfr1(8) = howsow   
          iprn = ittag(k,2)         
          call wrxnav ( urinex,iprn,iweek0,bcephem,bcclock,subfr1 )
        endif
c     end loop over all records
      enddo

      stop
      end


      subroutine rdeph(lu,ittag,sec,ephr,icnt,ios)
c     Read CIGNET ephemeris data file
c     Note that the CIGNET format changed on week 424
            
      include '../includes/argo2fic.h'

      integer*4 ittag(maxeph,7),icnt,igpswk,i,iwke
      integer lu,ios,nerr,idoy,idoyd,iflag
      real*8 sec(maxeph)
      real*8 ephr(maxeph,27),soweph,utcoff
      character*80 buff80
      character*256 message

      save nerr
      data nerr/0/
c
      ios = 0
c
c     Read the first line containing the GPS week number
c     Do this only once.
      if (icnt .eq. 1) then     
         read(lu,'(a80)',iostat=ios) buff80
         if( ios.ne.0 ) then
           call report_stat('FATAL','ARGO2NAV','rdeph'
     .     ,' ','Error reading first record of ARGO orb file',ios)
         else        
           call report_stat('STATUS','ARGO2NAV','rdeph',buff80
     .     ,'ARGO file header: ',0)   
            read(buff80,'(i4)',end=600,err=600,iostat=ios) ittag(1,1)
c           New format blank line
            if (ittag(1,1).ge.424) then
              read(lu,'(1x)',iostat=ios) 
              if( ios.ne.0 )
     .          call report_stat('FATAL','ARGO2NAV','rdeph'
     .       , ' ','Error reading second record of ARGO orb file',ios)
            endif
     
         endif
      endif

c     Read the date information
c     Different format begins during gps week 424 (M. Chin)
c
 150  continue
      if (ittag(1,1).lt.424) then
         read(lu,200,end=600,err=500,iostat=ios)(ittag(icnt,i),i=2,7),
     $       sec(icnt),(ephr(icnt,i),i=1,19)
 200     format(/,i2,5i3,f5.1,3d19.12,4(/,3x,4d19.12))
      else
         read(lu,201,end=600,err=500,iostat=ios)(ittag(icnt,i),i=2,7),
     $       sec(icnt),(ephr(icnt,i),i=1,27)
 201     format(i2,5i3,f5.1,3d19.12,6(/,3x,4d19.12))
c         print *,'RDEPH ',(ittag(icnt,i),i=1,7),sec(icnt)
      endif  
      call fix_y2k(ittag(icnt,3))

c     protect from bad week numbers
c     if bad, read another record, overwriting the bad one

      if (ittag(icnt,3) .gt. 0) then
         idoyd = idoy(ittag(icnt,3),ittag(icnt,4),ittag(icnt,5))
c        convert to GPST week number and second of week
         iflag = -4
         call timcon( iflag
     .              , igpswk,soweph
     .              , ittag(icnt,3),idoyd,ittag(icnt,6),ittag(icnt,7)
     .              , sec(icnt),utcoff)
         if (ittag(1,1).lt.424) then
            iwke = igpswk
         else
            iwke = ephr(icnt,22)
         endif
      else
         iwke = 0
      endif

      if (iwke .ne. ittag(1,1)) then
         if (mod(nerr,20) .eq. 0) then  
           write(message,'(a,2i5)') 'Week number mismatch: ',igpswk,iwke
           call report_stat('WARNING','ARGO2NAV','rdeph',' '
     .                       ,message,0)  
           write(message,'(a,6i5)') '   ',(ittag(icnt,i),i=2,7)
           call report_stat('WARNING','ARGO2NAV','rdeph',' '
     .                       ,message,0)   
         endif
         nerr = nerr + 1
         goto 150
      endif

500   if (ios .ne. 0) then   
         call report_stat('WARNING','ARGO2NAV','rdeph',' '
     .                   ,'File error',ios)   
        nerr = nerr + 1
        call ferror (ios,6)
        if (nerr .lt. 1000) then
           goto 150
        else          
         call report_stat('FATAL','ARGO2NAV','rdeph',' '
     .                   ,'Too many errors',ios)   
        endif
      end if
600   return
      end

