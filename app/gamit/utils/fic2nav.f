      Program FIC2NAV

c     Translate FICA Block 9 format (GAMIT e-files) to RINEX nav format
                                                          
c     R. King    11 June 2001
                           
c     Runs from command-line:  fic2nav <efile-name> 

      implicit none
       
c     --dimensions 
      include '../includes/argo2fic.h' 
      include '../includes/makex.h'

c     Ephemeris items  

      real*8 ephem(16),clock(6),subfr1(8),trans_sow
      integer*4 iwkn,nprn
         
c     File names and units

      character*80 efile,rxfile   
      integer*4 iue,iurx  
      parameter (iue=1)
      parameter (iurx=2)
 
c     RINEX header
      real*4 rxver  
      character*20 rxpgm,rxusr,rxdat

c     Controls and other variables

      logical eof 
      integer*4 icall,iflag,jyr,jmo,jdy,iargc,ioerr,nhead,ncnt,iyr,idoy   
      character*1  gnss
      character*3  adoy
      character*4  ayr
      character*15 buff15
      character*16 uname
      character*20 version 
      character*80 head(maxlin) 
      character*256 message

      data version/'FIC2NAV version 1.0'/



***** Remove old versions of the status, warning, and error files

      call report_stat('CLEAR','FIC2NAV',' ',' ', ' ',0)
      call report_stat('CLEAR','LIB',' ',' ', ' ',0)
      
***** Get the command-line arguments

      if( iargc().lt.2 ) then
         write(*,*)' '
         write(*,*)' '
         write(*,*)'Program FIC2NAV:'
         write(*,*)' '
         write(*,*)' Convert a FICA Blk 9 to RINEX navigation file'
         write(*,*)' '
         write(*,*)'Example:'  
         write(*,*) ' ' 
         write(*,*)'fic2nav <FICA-filename> <yr> <doy> '
         write(*,*)' '
         stop
      endif

      call rcpar(1,efile)  
      call rcpar(2,ayr)  
      call rcpar(3,adoy) 

***** Open the input and output files

      open(unit=iue,file=efile,status='old'
     .    ,form='formatted',iostat=ioerr) 
      if (ioerr.ne.0) then 
        call report_stat('FATAL','FIC2NAV','fic2nav',efile
     .      ,'Cannot open input e-file',0) 
      else
        call report_stat('STATUS','FIC2NAV','fic2nav',efile
     .                       ,'Opening: ',0)
      endif  
c     assume that the FICA file begins with the 4-character station name
      rxfile = ' '
      rxfile(1:4) = efile(1:4) 
      read(ayr,'(i4)') iyr 
      iyr = mod(iyr,100)
      write(ayr,'(i2)') iyr     
      read(adoy,'(i3)') idoy 
      write(adoy,'(i3)') idoy 
      if( adoy(1:1).eq.' ') adoy(1:1) = '0'
      if( adoy(2:2).eq.' ') adoy(2:2) = '0'    
      rxfile(5:7)= adoy
      rxfile(8:9) = '0.'
      rxfile(10:11) = ayr(1:2)
      rxfile(12:12) = 'n'
      open(unit=iurx,file=rxfile,status='unknown'
     .    ,form='formatted',iostat=ioerr) 
      if (ioerr.ne.0) then 
        call report_stat('FATAL','FIC2NAV','fic2nav',rxfile
     .      ,'Error opening output RINEX file',0) 
      else
        call report_stat('STATUS','FIC2NAV','fic2nav',rxfile
     .                       ,'Opening: ',0)
      endif  

  
***** Set the RINEX header info

c     -- RINEX version
      rxver = 2.0
c       identify FIC2NAV version in 20 characters
      write(rxpgm,'(a20)') version
c     -- user name, agency, and run date
      call getusr(uname) 
      rxusr = '    '//uname    
      call getdat(jyr,jmo,jdy) 
      write (rxdat,'(i4.4,"/",i2.2,"/",i2.2)') jyr,jmo,jdy 

           
****** Read and write the header and ephemeris records
       
      eof = .false.
      ncnt = 1  
c     read the first record to get the header   
      icall = 0
c  TEMPORARY: Force GPS
      gnss = 'G'
      call reade( iue,icall,gnss
     .          , iflag,trans_sow,nprn,iwkn,ephem,clock,subfr1 )  
      icall = 1  
      buff15= 'NAVIGATION DATA'
      write (iurx,'(f9.2,11x,a15,25x,a)',iostat=ioerr)
     .                rxver,buff15,'RINEX VERSION / TYPE'
      write (iurx,'(3a20,A)',iostat=ioerr)
     .                rxpgm,rxusr,rxdat,'PGM / RUN BY / DATE' 
      write (iurx,'(a,27x,a7,13x)') 
     .                'Broadcast elements from FICA file','COMMENT'
      write (iurx,'(60x,a20)') 'END OF HEADER       '   
c     map and write the first navigation record
      call wrxnav ( iurx,nprn,iwkn,ephem,clock,subfr1 )
c     now do the rest
      do while( .not.eof )  
        call reade( iue,icall,gnss
     .            , iflag,trans_sow,nprn,iwkn,ephem,clock,subfr1 )
        if (iflag.eq. -1) then
           ncnt = ncnt - 1
           close(iue)   
           eof = .true.
        else
           if( ncnt.ge.maxeph ) then
             write(message,'(2a,i5)') 'Number of ephemeris records'
     .               ,' exceeds MAXEPH = ',maxeph
             call report_stat('FATAL','FIC2NAV','fic2nav',' '
     .                       ,message,0)   
            else
              ncnt = ncnt + 1
            endif
        endif
        call wrxnav ( iurx,nprn,iwkn,ephem,clock,subfr1 )
      enddo

      stop
      end
