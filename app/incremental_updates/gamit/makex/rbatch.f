      subroutine rbatch
     .   ( debug,iyear,idoy,gnsslf,numsit,sites,rcvrs,vers )
* MOD TAH 200511: Return gnsslf not just gnss
c
c     Get parameters out of the batch file   
c
c     P. Morgan and K. Feigl   Jan-Jul 1987
c     R. King Jul 1999: Modified for 4-digit years and new format by R. King Jul 1999 
c     R. King Dec 2014: Modified for a single call, storing values in arrays
c                                            
c     The format of the file is a list of I/O files to be used,
c     followed by descriptor lines for each X-file to be generated.
c
c     Example
c
c
c infor 1 //bowie/data/gps/trex/makex/
c sceno 1 //bullen/data/tables/makex.sceno
c rinex 1 //bowie/data/gps/trex/makex/
c       0 no longer used
c coord 1 //bullen/data/gps/tables/ltrex2.all
c sited 1 station.info
c xfile 1 x
c svclk 1 //bowie/data/trex/makex/jtrex2.007
c clock 1 k
c sp3   1
c nav   1 //bowie/data/trex/makex/brdc007.92n
c snsys 1 G
c MOD TAH 200511: Mod line above to allow C7, G5 type entries.
c site        sw ver        
c (a4,1x,a4,1x,a3,1x,a1,2x,a3,1x,f4.2)
c FIBR 1986 365 1  GES 1.0

c In the I/O listing, all possible files must be listed in the order
c given.  If the variable name is followed by a '1', the file is to be used;
c if by a  '0' it is not used.  

c The next line is an unread header (e.g.): site  year  day session  sw  ver

c The next line is the format for the X-file descriptor lines to follow.
c Prior to release 9.90, the descriptor lines looked like
c cato871.278  GES 1.20 
c and were read by format (A4,A2,A1,1X,A3,2X,A3,1X,F4.2)  
c After 9.90 (July 1999), the descriptor lines look like
c cato 1987 278 1  GES 1.20
c and are read by format (a4,1x,a4,1x,a3,1x,a1,2x,a3,1x,f4.2)

      implicit none
                                                  
      include '../includes/dimpar.h' 

C     IO status
      integer    ioerr
c     0 for no file, 1 for file
      integer    iswitches(12)
c     every program needs it
      integer     i

c     number of sites
      integer*4 numsit  
c     sites
      character*4 sites(maxsit)
C     receiver software name  GES, COR, MAC, MIN, TRM, ASH, ROG
      character*3 rcvrs(maxsit)
c     firmware version
      real*4 vers(maxsit)
           
      character*1 gnss
      character*(*) gnsslf   ! GNSS with lower frequency added (optionally)
      character  fform1*70,line*80,year*4
      character  day*3,site*4,sesh*1
      character*120 fnames(12)
      character*256 message

c     return variables now integer
      character  buf2*2, buf3*3, buf4*4
      integer*4 iyear,idoy,isessn

      logical     debug,eof 

      include '../includes/makex.h'

c  Read 12 lines of switches and file names

         do i = 1,12
            read(ubatch,'(6x,i1,1x,a120)',iostat=ioerr)
     .        iswitches(i),fnames(i)                 
            if( ioerr.ne.0 ) then
              write(message,'(a,i2,a,i2,a,i2,a,a120)')
     .        'Error reading batch file--iswitches(',i,') ',iswitches(i)
     .                 ,', fnames(',i,') ',fnames(i)
              write (uinfor,'(a)') message       
              call report_stat('FATAL','MAKEX','rbatch',' ',message
     .                 ,ioerr)     
            endif
            if(debug)  write (uscren,*) i,iswitches(i),fnames(i)
         enddo
         if (iswitches(1) .eq. 1) then
            qinfor = .true.
            finfor = fnames(1)
         else
            qinfor = .false.
            finfor = 'junk'
         endif
         if (iswitches(2) .eq. 1) then
            qsceno = .true.
            fsceno = fnames(2)
         else
            qsceno = .false.
            fsceno = 'junk'
         endif
         if (iswitches(3) .eq. 1) then
            qrinex = .true.
            frinex = fnames(3)
         else
            qrinex = .false.
            frinex = 'junk'
         endif
c         if (iswitches(4) .eq. 1) then
c         else
c        endif
         if (iswitches(5) .eq. 1) then
            qcoord = .true.
            fcoord = fnames(5)
         else
            qcoord = .false.
            fcoord = 'junk'
         endif
         if (iswitches(6) .eq. 1) then
            qsited = .true.
            fsited = fnames(6)
         else
            qsited = .true.
            fsited = 'junk'
         endif
         if (iswitches(7) .eq. 1) then
            qxfile = .true.
            fxfile = fnames(7)
         else
            qxfile = .false.
            fxfile = 'junk'
         endif
         if (iswitches(8) .eq. 1) then
            qsvclk = .true.
            fsvclk = fnames(8)
         else
            qsvclk = .false.
            fsvclk = 'junk'
         endif
         if (iswitches(9) .eq. 1) then
            qclock = .true.
            fclock = fnames(9)
         else
            qclock = .false.
            fclock = 'junk'
         endif
         if (iswitches(10) .eq. 1) then
            qsp3 = .true.
            fsp3 = fnames(10)
         else
            qsp3 = .false.
            fsp3 = 'junk'
         endif
         if (iswitches(11) .eq. 1) then
            qnav   = .true.
            fnav   = fnames(11)   
         else
            qnav   = .false.
            fnav   = 'junk'
         endif
cd         print *,'RBATCH1 iswitches(12) fnames(12) '
cd     .        ,iswitches(12),fnames(12)(1:1)
         if (iswitches(12) .eq. 1 ) then
cd            print *,'fnames(12)(1:1) ',fnames(12)(1:1) 
* MOD TAH 200511: Get the gnsslg value from string. (2 characters)
            gnsslf = fnames(12)(1:2)
         else
            gnsslf = 'G'
         endif     
cd         print *,'gnss ',gnss
cd       stop 1
c        read the header line and the format line
         read(ubatch,'(a80)',iostat=ioerr) line   
cd         print *,'header line ',line
         if( ioerr.ne.0 ) then
           write(message,'(a)') 'Error reading batch file header line'
           write (uinfor,'(a)') message       
           call report_stat('FATAL','MAKEX','rbatch',' ',message,ioerr)
         endif
         read(ubatch,'(a70)',iostat=ioerr) fform1
cd         print *,'format line ',fform1
         if( ioerr.ne.0 ) then
           write(message,'(a)')'Error reading batch file format line'
           write (uinfor,'(a)') message       
           call report_stat('FATAL','MAKEX','rbatch',' ',message,ioerr)
         endif
         call lowers (finfor)
         call lowers (fsceno)
         call lowers (frinex)
         call lowers (fcoord)
         call lowers (fsited)
         call lowers (fxfile)
         call lowers (fsvclk)
         call lowers (fclock)
         call lowers (fsp3)
         call lowers (fnav  )
         call lowers (fextr3)
      
c  Loop over all site entries

c        Prior to release 9.90 (Y2K complient) we used a format of the form
c            cato871.278  GES 1.20 
c          read by the format (A4,A2,A1,1X,A3,2X,A3,1X,F4.2) 
c        Now we use a less artificial format of the form
c            cato 1987 278 1  GES 1.20
c          read by the format (a4,1x,a4,1x,a3,1x,a1,2x,a3,1x,f4.2)
c        Since the old format was never changed, we'll test on it to tell the difference.
c        The year and date must be the same for all sites (hence redundant), and 
c         session number is no longer used, so not returned.
                           
      eof = .false. 
      numsit = 0 
      do while (.not.eof ) 
        if( fform1(5:6).eq.'A2' .or. fform1(5:6).eq.'a2' ) then 
c         old format
          year = '    '                                      
          read(ubatch,'(a)',iostat= ioerr) line
cd          print *,'old format line ',line
* MOD TAH 200302: Check if site name is blank and if so treat this like EOF.
*         This was the old behavior of makex which changed for some 
*         reason.  Implemented with .or. statement below.
          if( ioerr.eq.-1 .or. line(1:4).eq.'    ' ) then 
            eof = .true.
          elseif( ioerr.ne.0 ) then
            call report_stat('FATAL','MAKEX','rbatch',' '
     .         ,'Error reading site line of batch file ',ioerr)
          else
            numsit = numsit + 1
            if(numsit.gt.maxsit) then 
                write(message,'(a,i3,a,i3,a)') 
     .         'Number of sites in makex batch file (',numsit
     .         ,')  greater than MAXSIT in dimpar.h (',maxsit,')'
              call report_stat('FATAL','MAKEX','rbatch',' ',message
     .             ,ioerr)  
            endif 
            read(line,fform1,iostat=ioerr) 
     .       sites(numsit),year(3:4),sesh,day,rcvrs(numsit),vers(numsit)
            if( ioerr.ne.0 ) then 
              write(message,'(a)') 'Error reading batch file site line'
              write (uinfor,'(a)') message       
              call report_stat('FATAL','MAKEX','rbatch',' ',message
     .             ,ioerr)
            endif
          endif
        elseif (fform1(5:6).eq.'1X'.or.fform1(5:6).eq.'1x') then
c         presumed new format
          read(ubatch,'(a)',iostat=ioerr) line
cd          print *,'new format line ',line
* MOD TAH 200624: Check if site name is blank and if so treat this like EOF.
*         This was the old behavior of makex which changed for some 
*         reason.  Implemented with .or. statement below. 
*         Added March update to the correct file format version
          if( ioerr.eq.-1 .or. line(1:4).eq.'    ' ) then 
            eof = .true.
          elseif( ioerr.ne.0 ) then    
            write(message,'(a)') 'Error reading batch file site line'
            write (uinfor,'(a)') message       
            call report_stat('FATAL','MAKEX','rbatch',' ',message
     .         ,ioerr)
          else
            numsit = numsit + 1  
            if(numsit.gt.maxsit) then 
                write(message,'(a,i3,a,i3,a)') 
     .         'Number of sites in makex batch file (',numsit
     .         ,')  greater than MAXSIT in dimpar.h (',maxsit,')'
              call report_stat('FATAL','MAKEX','rbatch',' ',message
     .             ,ioerr)                                         
            endif 
            read(line,fform1,iostat=ioerr)
     .         sites(numsit),year,day,sesh,rcvrs(numsit),vers(numsit)      
cd               print *,'numsit sites year day sesh rcvrs vers ',
cd     .    numsit,sites(numsit),year,day,sesh,rcvrs(numsit),vers(numsit)
            if( ioerr.ne.0 ) then 
              write(message,'(a)') 'Error reading batch file site line'
              write (uinfor,'(a)') message       
              call report_stat('FATAL','MAKEX','rbatch',' ',message
     .             ,ioerr)
            endif
          endif
        endif
        if( .not.eof ) then 
          if( numsit.eq.1 ) then   
            call lowers(year)
            call lowers(sesh)
            call lowers(day)   
            write(buf4,'(a4)') year
            read(buf4,'(i4)') iyear
            call fix_y2k(iyear)
            write(buf3,'(a3)') day
            read(buf3,'(i3)') idoy
          endif
          call lowers(sites(numsit))
          call lowers(rcvrs(numsit))
          if (debug) then
            print *,'RBATCH: '
            print *, fform1,numsit,sites(numsit),year,day
     .           ,rcvrs(numsit),vers(numsit)
          endif
        else
           write (uinfor,200) 
  200     format (//,' End of batch file ',//)
          call report_stat('STATUS','MAKEX','rbatch',' ',
     .       'End of batch file reached',ioerr)
          return
        endif 
c     ---end of loop on site lines
      enddo

      return
      end

