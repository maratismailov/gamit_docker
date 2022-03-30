      Subroutine rficahd( debug,found,rcvrsw,swveru,sample_int
     .                  , auser,arcvr,asite,antht
     .                  , ncomments,comments 
     .                  , nwave1,nwave2 )

c     Read header block types FICA files 
c     R.W. King  7 March 1997 from makex/rhead.f

      implicit none

      include '../includes/makex.h'

c ----Passed variables
        
c     Input
c       debug     logical   : extra print if true

c     Output
c       rcvrsw            c*3  : GAMIT-coded receiver firmware
c       swveru            r*8  : firmware version
c       sample_int        r*8  : receiver collection interval 
c       auser             c*8  : operator id   
c       arcvr             c*16 : receiver serial number
c       asite             c*16 : station name
c       antht             c*16 : operator antenna height
c       ncomments         i*4  : number of non-blank elements in comments array
c       comments(maxlin)  c*80 : array of comment lines
c       nwave1            i*4  : L1 frequency factor (RINEX definition)
c       nwave2            i*4  : L2 frequency factor (RINEX definition)

      character*3  rcvrsw 
      character*80 comments(maxlin)
      character*256 message

      integer*4 ncomments,nwave1,nwave2,nlines,nwords,i1,i2,i,j

      logical      debug,found,nerr

c     receiver serial number, station name, and operator antenna height
      character*16 arcvr,asite,antht

C     FICA block number
      integer*4 iblkid,ioerr
c     0=doc, 1=data, 2=ephem
      integer*4 iflag
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
c     user id
      character*8 auser

c     FICA ARRAYS
      real*8 ff(maxfic)
      integer*4 fi(maxfic)
      character*8 fc(maxfic)
C     number of elements in FICA arrays
      integer*4 nf,ni,nc
   
c     Blank out the character variables to avoid writing nulls, interpreted by UNIX as binary
     
      do i=1,maxlin
        comments(i) = ' '
      enddo    
      ncomments = 0 
      asite = ' '
      auser = ' '
      arcvr = ' ' 
   
      found = .false.
      nerr = .false. 
      ioerr = 0
      do while ( .not.found ) 
          
       call rfica (uficaf,uinfor,iblkid,ff,fi,fc,nf,ni,nc,ioerr)
       if( ioerr.eq.0 ) then

         if (iblkid .eq. 0) then 
             write (uinfor,'(a)')  'FICA Blk 0 comments:'
             nlines = nc/7 + 1 
             do j=1,nlines
                ncomments = ncomments+1
                if( ncomments.lt.maxlin ) then  
                  if( j.lt.nlines ) then
                    nwords = 7 
                  else 
                    nwords = mod(nc,7)
                  endif     
                  i1 =  (j-1)*7 + 1 
                  i2 = i1 + (nwords-1)
c                  print *,'BLK 0 nlines nwords j ncomments '
c     .                   , nlines,nwords,j,ncomments 
                  write(comments(ncomments)(1:60),'(7a8)') 
     .                          (fc(i),i=i1,i2)     
                  if( debug) write (*,*) (fc(i),i=1,nc)
                  write (uinfor,'(a)') comments(ncomments)(1:60)
                else          
                   write(message,'(a,i3,a,i3,a)') 
     .                'Number of comment lines from Blk 0 (=',ncomments
     .               ,') exceeds maxlin in makex.h (=',maxlin,')'
                   call report_stat('WARNING','FIC2RX','rficahd',' '
     .                             , message,0)
                endif
             enddo  
c            Try setting the receiver software from the BLK 0 comment lines.
c            For safety and efficiency, look for specific strings, which should be 
c            added to the code for each example encountered.
c            Case 1:  TI-ROM files with only BLK 0 and BLKs 401.  The BLK 0 entry reads
c                       BLK      0    0    0    8
c                       TI-NAV  88/02/12 ARL:UT 92/05/07 19:28:47 MS-DOS  MS-FOR  V 1.2 
               if( fc(1).eq.'TI-NAV  ' )  rcvrsw = 'ROM'   
c            Case 2: GESAR or ASDAP files with reduced obs blocks (BLK 80).  The BLK 0 entries read
c                       BLK      0    0    0    8                                                       
c                       FICALABL1.1     MIT GPS 89/02/18 16:01  moja8069 kurt    FICA                   
c                       BLK      0    0    0    8                                                       
c                       FICASLIM1.5     MIT GPS 88/11/06 22:10  moja8069 kurt    FICA                   
c                       BLK      0    0    0    8                                                       
c                       NGS2FIC 1.6     MIT GPS 1988072912:06:05moj426  mhm       TRADAS 
c            Case 3: GESAR files translated with UNPACK but with BLK 101 missing.  The BLK 0 entry reads
c                       BLK      0    0    0    8                                                       
c                       UNPACK  1.7     MIT GPS 1989062922:26:20   gavi9kurt       DAT80                
               if( fc(8).eq.'  TRADAS') rcvrsw = 'ROM'  
               if( fc(1).eq.'UNPACK  ') rcvrsw = 'GES'
         else if (iblkid .eq. 101) then
               call blk101 (debug,iflag,apress,atemp,ahumid,sample_int
     .                    ,swveru,auser,asite,arcvr,antht
     .                    ,ff,fi,fc,nf,ni,nc ) 
               found = .true.
               rcvrsw = 'GES'
c              set wavelength factors for TI 4100
               nwave1 = 1
               nwave2 = 1 
               write (uinfor,'(a)') 'FICA Blk 101 values:'
               write (uinfor,113) swveru
               write (uinfor,114) asite
               write (uinfor,115) auser
               write (uinfor,116) arcvr
               write (uinfor,120) antht
               write (uinfor,117) ahumid
               write (uinfor,118) atemp
               write (uinfor,119) apress
               write (uinfor,121) sample_int
               if( debug ) then
                 write (uscren,113) swveru
                 write (uscren,114) asite
                 write (uscren,115) auser
                 write (uscren,116) arcvr
                 write (uscren,120) antht
                 write (uscren,117) ahumid
                 write (uscren,118) atemp
                 write (uscren,119) apress
               endif
               if( (ncomments+9).le.maxlin ) then
                write (comments(ncomments+1),113) swveru
                write (comments(ncomments+2),114) asite
                write (comments(ncomments+3),115) auser
                write (comments(ncomments+4),116) arcvr
                write (comments(ncomments+5),120) antht
                write (comments(ncomments+6),117) ahumid
                write (comments(ncomments+7),118) atemp
                write (comments(ncomments+8),119) apress
                write (comments(ncomments+9),121) sample_int
                ncomments = ncomments + 9  
 113            format (1x,'Software version                    ',f4.2)
 114            format (1x,'Operator input site code            ',a16)
 115            format (1x,'Operator initials                   ',a16)
 116            format (1x,'Operator input receiver serial num  ',a16)
 120            format (1x,'Operator input antenna height (m)   ',a16)
 117            format (1x,'Operator input humidity (%)         ',f16.4)
 118            format (1x,'Operator input temperature (deg C)  ',f16.4)
 119            format (1x,'Operator input pressure (Mb)        ',f16.4)
 121            format (1x,'Operator input collection int (sec) ',f16.4)


               else        
                 write(message,'(a,i3,a,i3,a)') 
     .               'Number of comment lines from Blk 101 (=',ncomments
     .              ,') exceeds maxlin in makex.h (=',maxlin,')'
                   call report_stat('WARNING','FIC2RX','rficahd',' '
     .                             , message,0)  
               endif

         else if (iblkid .eq. 1001 .or. iblkid.eq.1101
     .          .or. iblkid.eq.1201 .or. iblkid.eq.1301 ) then
               if( iblkid.eq.1001) then
c                  MACROMETER II
                   call blk1001 ( debug,iflag,offsl1,offsl2,sample_int
     .                          , arcvr,asite,ff,fi,fc,nf,ni,nc )
                   rcvrsw = 'MAC'
                   found = .true.
                   nwave1 = 2
                   nwave2 = 2
               else if( iblkid.eq.1101) then 
c                  MiniMac
                   call blk1101 ( debug,iflag,offsl1,offsl2,sample_int
     .                          , arcvr,asite,ff,fi,fc,nf,ni,nc )   
                   rcvrsw = 'MIN'
                   found = .true.
                   nwave1 = 1
                   nwave2 = 2
               else if( iblkid.eq.1201) then
c                  Trimble SST (not tested)
                   call blk1201 ( debug,iflag,offsl1,offsl2,sample_int
     .                          , arcvr,asite,ff,fi,fc,nf,ni,nc ) 
                   call report_stat('WARNING','FIC2RX','rficahd',' '
     .                             ,'Trimble SST FICA not tested',0)
                   rcvrsw = 'TRM'
                   found = .true.
                   nwave1 = 1
                   nwave2 = 2
               else if( iblkid.eq.1301) then    
c                  Rogue
                   call blk1301 ( debug,iflag,offsl1,offsl2,sample_int
     .                          , arcvr,asite,ff,fi,fc,nf,ni,nc ) 
                   call report_stat('WARNING','FIC2RX','rficahd',' '
     .                             ,'Rogue FICA not tested',0) 
                   rcvrsw = 'ROG'
                   found = .true.
                   nwave1 = 1
                   nwave2 = 1
               endif
               write (uinfor,134) asite
               write (uinfor,136) arcvr
               write (uinfor,137) sample_int
               if( debug ) then
                 write (uscren,134) asite
                 write (uscren,136) arcvr
                 write (uscren,137) sample_int
               endif  
               if( (ncomments+3).le.maxlin) then
                 write (comments(ncomments+1),134) asite
                 write (comments(ncomments+2),136) arcvr
                 write (comments(ncomments+3),137) sample_int
                 ncomments = ncomments + 3   
 134             format (1x,'FICA file site code            ',a16)
 136             format (1x,'FICA file receiver serial num  ',a16)
 137             format (1x,'FICA file collection interval  ',f8.2)  
               else                 
                 write(message,'(a,i3,a,i3,a)') 
     .          'Number of comment lines from header block (=',ncomments
     .         ,') exceeds maxlin in makex.h (=',maxlin,')'
                 call report_stat('WARNING','FIC2RX','rficahd',' '
     .                           , message,0)  
               endif
         endif

       else                         
          call report_stat('WARNING','FIC2RX','rficahd',fficaf
     .                      ,'Missing header block',0)
          found = .true.
       endif
          
c     end of FICA read  
      enddo

      if( found ) then
c       check and convert sampling interval
        if( dmod(sample_int,1.d0).gt.1.d-06 ) then
           write(message,141) sample_int
  141      format('Receiver sampling interval (',f12.6
     .              , ' not an integer: ')
           call report_stat('FATAL','MAKEX','rhead',' ',message,0)
        endif
      endif   

      return
      end


