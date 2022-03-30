      Subroutine RIXHED ( luinx,version
     . , int,ixmaps,mapdim,iexp
     . , hgt1,hgt2,dhgt,lat1,lat2,dlat,lon1,lon2,dlon
     . , ixyr0,ixmo0,ixdy0,ixhr0,ixmn0,ixsec0
     . , ixyr1,ixmo1,ixdy1,ixhr1,ixmn1,ixsec1 )

c     Read an IONEX header on logical unit luinx, assumed open

c     Based on Feigl/King routine RRXHED using descriptors provided
c     by AIUB IONEX definition.   R. King 17 April 2007

c     assumed open

      implicit none

      include '../includes/dimpar.h'   
                 
c     unit number for IONEX file 
      integer*4 luinx
c     variables for status reporting
      integer*4 ioerr,nerr,len,rcpar
      character*12 fname
      character*80 prog_name
      character*256 message

c     strings for buffering data
      character*20 buff20
      character*80 buff80,blnk80

c     IONEX defined items
c       file type and version number
      character*1  filtyp
      real*4 version 
c       satellite system
      character*3 asvsys
c       file identifiers
      character*20 pgm,usr,date
c         comments and description: not read
c      character*60 comments(maxlin)
c      character*60 description(maxlin)
c       data interval in seconds  
      integer*4 iint
      real*4 int   
c       number of TEC maps in file
      integer*4 ixmaps
c       dimension of ionsophere maps
      integer*4 mapdim 
c       exponent for TEC values
      integer*4 iexp
c       spatial range and interval
      real*4 hgt1,hgt2,dhgt,lat1,lat2,dlat,lon1,lon2,dlon
c       data start time
      integer ixyr0,ixmo0,ixdy0,ixhr0,ixmn0,ixsec0
c       data stop time
      integer ixyr1,ixmo1,ixdy1,ixhr1,ixmn1,ixsec1
c       base radius
      real*8 brad
c  c    number of satellites -- not read
c      integer nsv 
c     satellite names (number designation) ---not read
c      character*3 asv(maxsat)
c       integer isv(maxsat)    


c  Get calling program name and IONEX file name for report_stat

       len = rcpar(0,prog_name)
       inquire( unit=luinx, name=fname, iostat=ioerr )
       if(ioerr.ne.0) call report_stat('WARNING',prog_name,'lib/rrxhed'
     .   ,' ','Cannot get name of IONEX file for error reporting',ioerr)
     
c  Initialize buffers to avoid using info from previous file if line missing
      write (blnk80,'(80x)')
      write (pgm ,'(20x)')
      write (usr ,'(20x)')
      ixyr0 = 0
      ixmo0 = 0
      ixdy0 = 0
      ixhr0 = 0
      ixmn0 = 0
      ixsec0 = 0.d0
      ixyr1 = 0
      ixmo1 = 0
      ixdy1 = 0
      ixhr1 = 0
      ixmn1 = 0
      ixsec1 = 0.d0

c     set IONEX version to 1 until read, for initial check of end-of-header
      version = 1.0
c      ixsat = 1 
      int = 0.    
      iexp = 0
      nerr = 0      

c Read first line of IONEX header separately to guard against a bogus blank first line
      read( unit=luinx,iostat=ioerr,fmt='(a)') buff80
      if( buff80(61:73).ne.'IONEX VERSION' ) then
           call report_stat('FATAL',prog_name,'rixhed',fname
     .                     ,'Bogus first line of IONEX file',ioerr) 
      else
c          If the version number is an integer written under version 1 format
c          ( I6 ) it will not be read correctly by the version 2 format ( f9.2)
c          so read the version number free-format before reading the rest of
c          the line.  rwk 060526
          read( buff80(1:9),*,iostat=ioerr) version  
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'rixhed'
     .      ,fname,'Error reading IONEX version number ',ioerr) 
          read (buff80,'(20x,a1,19x,a3)',iostat=ioerr) filtyp,asvsys
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'rixhed'
     .      ,fname,'Error reading IONEX file type ',ioerr) 
      endif


 10   continue 
      buff80 = ' '
      read (unit   = luinx,
     .     iostat = ioerr,
     .        fmt    = '(a)') buff80
      BUFF20 = buff80(61:80)
      call uppers(BUFF20)
      if (ioerr .ne. 0) then
         call report_stat('FATAL',prog_name,'rixhed',fname
     .                  ,'Error reading header for IONEX file',ioerr)
      elseif((version.le.1.1.and. buff20.eq.'                    ') .or.
c             set to accept END OF HEADER with IONEX 1 (non-standard)
c             (version.eq.2 .and. buff20(:13).eq.'END OF HEADER') .or.  
     .   ( buff20(:13).eq.'END OF HEADER') .or. ioerr.eq.-1 ) then
        goto 20

      elseif (BUFF20(1:19) .eq. 'PGM / RUN BY / DATE') then
         read (buff80,'(3a20)',iostat=ioerr)
     .         pgm,usr,date

      elseif (BUFF20(1:7) .eq. 'COMMENT') then 
c         for now
        continue
                
      elseif (BUFF20(1:12) .eq. 'DESCRIPTION') then
c         for now
        continue

      elseif (BUFF20(1:18) .eq. 'EPOCH OF FIRST MAP') then
        read (buff80,'(6i6)',iostat=ioerr)
     .       ixyr0,ixmo0,ixdy0,ixhr0,ixmn0,ixsec0  
        call fix_y2k(ixyr0)
               
      elseif (BUFF20(1:17) .eq. 'EPOCH OF LAST MAP') then
        read (buff80,'(5i6,f12.6)',iostat=ioerr)
     .         ixyr1,ixmo1,ixdy1,ixhr1,ixmn1,ixsec1
         call fix_y2k(ixyr1)
                                          
      elseif (BUFF20(1:8) .eq. 'INTERVAL' ) then 
        read (buff80,'(i6)',iostat=ioerr) iint
        int = real(iint)

      elseif (BUFF20(1:17) .eq. '# OF MAPS IN FILE') then
        read (buff80,'(i6)',iostat=ioerr)  ixmaps             
                                                        
      elseif (BUFF20(1:16) .eq. 'MAPPING FUNCTION' ) then 
c           for now
        continue         

      elseif (BUFF20(1:16) .eq. 'ELEVATION CUTOFF' ) then 
c              for now
        continue         
                               
      elseif (BUFF20(1:16) .eq. 'OBSERVABLES USED' ) then 
c            for now
        continue         

      elseif (BUFF20(1:13) .eq. '# OF STATIONS' ) then 
c            for now
        continue         

      elseif (BUFF20(1:15) .eq. '# OF SATELLITES' ) then 
c            for now
        continue         
                                
      elseif (BUFF20(1:11) .eq. 'BASE RADIUS' ) then 
        read( buff80,'(f8.1)',iostat=ioerr) brad

      elseif (BUFF20(1:13) .eq. 'MAP DIMENSION' ) then
        read( buff80,'(i6)') mapdim
        if( mapdim.ne.2 ) call report_stat('FATAL',prog_name
     .            ,'lib/rixhed',' ','Map dimension not 2',0)   
            
      elseif (BUFF20(1:19) .eq.'HGT1 / HGT2 / DHGT') then                 
        read( buff80,'(2x,3f6.1)',iostat=ioerr) hgt1,hgt2,dhgt

      elseif (BUFF20(1:19) .eq.'LAT1 / LAT2 / DLAT') then                 
          read( buff80,'(2x,3f6.1)',iostat=ioerr) lat1,lat2,dlat
                                
      elseif (BUFF20(1:19) .eq.'LON1 / LON2 / DLON') then                 
        read( buff80,'(2x,3f6.1)',iostat=ioerr) lon1,lon2,dlon
 
      elseif (BUFF20(1:8) .eq. 'EXPONENT' ) then 
           read( buff80,'(i6)',iostat=ioerr) iexp
                   
      elseif (BUFF20(1:  ) .eq. 'START OF AUX DATA' ) then
c          skip all of this for now
                 
      else
         continue
      endif

      if( ioerr.ne.0 )  then
         write(message,'(a,a20)') 
     .     'Bad IONEX header record, ID=',buff20
         call report_stat('FATAL',prog_name,'lib/rixhed'
     .                       ,fname,message,ioerr)
      endif

c     go read another record
      goto 10

  20  continue 
c     encountered: END OF HEADER

      return
      end

