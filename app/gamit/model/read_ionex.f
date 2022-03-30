      Subroutine read_ionex

c     Read all values from an IONEX file into storage in model.h
c     R. King 8 May 2007       

      implicit none
         
      include '../includes/dimpar.h'
      include '../includes/model.h'

c     unit number for IONEX file is iuf in common /units/ in model.h
        
c     variables for status reporting
      integer*4 ioerr,nerr
      character*12 fname
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
      real*4 dint   
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

c  Functions

      integer*4 idoy 

c  Local variables

      logical eof 
      integer*4 nlatcnt,nlonlines,ival1,ival2,lonval(maxilon)
     .        , imapx,myr,month,mday,mhr,min,msec 
     .        , i,j,il    
      real*4 ion_rad,lat
      character*80 line 


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
      dint = 0.    
      iexp = 0
      nerr = 0      

c.. Read first line of IONEX header separately to guard against a bogus blank first line
      read( unit=iuf,iostat=ioerr,fmt='(a)') buff80
      if( buff80(61:73).ne.'IONEX VERSION' ) then
           call report_stat('FATAL','MODEL','read_ionex',fname
     .                     ,'Bogus first line of IONEX file',ioerr) 
      else
c          If the version number is an integer written under version 1 format
c          ( I6 ) it will not be read correctly by the version 2 format ( f9.2)
c          so read the version number free-format before reading the rest of
c          the line.  rwk 060526
          read( buff80(1:9),*,iostat=ioerr) version  
          if( ioerr.ne.0 ) call report_stat('FATAL','MODEL','read_ionex'
     .      ,fname,'Error reading IONEX version number ',ioerr) 
          read (buff80,'(20x,a1,19x,a3)',iostat=ioerr) filtyp,asvsys
          if( ioerr.ne.0 ) call report_stat('FATAL','MODEL','read_ionex'
     .      ,fname,'Error reading IONEX file type ',ioerr) 
      endif

c.. Read all of the header records

 10   continue 
      buff80 = ' '
      read (unit   = iuf,
     .     iostat = ioerr,
     .        fmt    = '(a)') buff80
      BUFF20 = buff80(61:80)
      call uppers(BUFF20)
      if (ioerr .ne. 0) then
         call report_stat('FATAL','MODEL','read_ionex',fname
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
        read (buff80,'(6i6)',iostat=ioerr)
     .         ixyr1,ixmo1,ixdy1,ixhr1,ixmn1,ixsec1
         call fix_y2k(ixyr1)
                                          
      elseif (BUFF20(1:8) .eq. 'INTERVAL' ) then 
        read (buff80,'(i6)',iostat=ioerr) iint
        dint = real(iint)

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
        if( mapdim.ne.2 ) call report_stat('FATAL','MODEL'
     .            ,'read_ionex',' ','Map dimension not 2',0)   
            
      elseif (BUFF20(1:19) .eq.'HGT1 / HGT2 / DHGT') then                 
        read( buff80,'(2x,3f6.1)',iostat=ioerr) hgt1,hgt2,dhgt

      elseif (BUFF20(1:19) .eq.'LAT1 / LAT2 / DLAT') then                 
          read( buff80,'(2x,3f6.1)',iostat=ioerr) lat1,lat2,dlat  
          nilat = int((lat2-lat1)/dlat) + 1
cd          print *,'READ_IONEX header lat1 lat2 dlat nilat '
cd     .                             , lat1,lat2,dlat,nilat
          if( nilat.gt.maxilat ) then
            write(message,'(a,i3,a,i3,a)') 
     .       '# lat values in grid (',nilat,') gt maxilat (',maxilat,')'
            call report_stat('FATAL','MODEL','read_ionex',' ',message,0)
          endif
                                 
      elseif (BUFF20(1:19) .eq.'LON1 / LON2 / DLON') then                 
        read( buff80,'(2x,3f6.1)',iostat=ioerr) lon1,lon2,dlon  
        nilon = int((lon2-lon1)/dlon) + 1 
cd        print *,'READ_IONEX header lon1 lon2 dlon nilon '
cd     .                             , lon1,lon2,dlon,nilon
        if( nilon.gt.maxilon ) then
          write(message,'(a,i3,a,i3,a)') 
     .      '# lon values in grid (',nilon,') gt maxilon (',maxilon,')'
          call report_stat('FATAL','MODEL','read_ionex',' ',message,0)
        endif
 
      elseif (BUFF20(1:8) .eq. 'EXPONENT' ) then 
           read( buff80,'(i6)',iostat=ioerr) iexp 
cd           print *,'EXPONENT ',iexp
                   
      elseif (BUFF20(1:  ) .eq. 'START OF AUX DATA' ) then
c          skip all of this for now
                 
      else
         continue
      endif

      if( ioerr.ne.0 )  then
         write(message,'(a,a20)') 
     .     'Bad IONEX header record, ID=',buff20
         call report_stat('FATAL','MODEL','read_ionex'
     .                       ,fname,message,ioerr)
      endif

c     go read another record
      goto 10

  20  continue 
c     encountered: END OF HEADER

                        
c .. Now read all the values into storage

      eof = .false.
      ntion = 0  
      nlatcnt = 0
      do while ( .not.eof )

        read(iuf,'(a)',iostat=ioerr) line  
        if( ioerr.eq.-1 ) then
          eof = .true.
        elseif( ioerr.ne.0 ) then 
          call report_stat('FATAL','MODEL','read_ionex',' '
     .                    ,'Error reading IONEX file',ioerr)
        elseif( line(61:74).eq.'END OF TEC MAP' ) then   
          nlatcnt = 0
        elseif(line(61:76).eq.'START OF TEC MAP') then
          ntion = ntion + 1    
cd          print *,'read START OF TEC MAP ntion =',ntion
          if( ntion.gt.maxion ) then  
            write(message,'(a,i4,a)') 
     .        '# maps on IONEX file > maxion (',maxion,')'
            call report_stat('FATAL','GRDTAB','grdtab/get_ion_grid'
     .                     ,' ',message,0)
          endif
          read(line,'(i6)',iostat=ioerr) imapx  
          if( ioerr.ne.0 )  call report_stat('FATAL','GRDTAB'
     .         ,'get_ion_grid',' ','Error reading ion map #',ioerr)
cd          print *,'map # index map # from file ',ntion,imapx
          read(iuf,'(a)',iostat=ioerr) line 
          if( ioerr.ne.0 )  call report_stat('FATAL','GRDTAB'
     .           ,'get_ion_grid',' ','Error reading EPOCH line',ioerr) 
cd           print *,'DEBUG epoch ',line   
          if( line(61:80).ne.'EPOCH OF CURRENT MAP') then
             call report_stat('FATAL','GRDTAB','get_ion_grid',' ' 
     .       ,'Line after START OF TEC MAP not EPOCH OF CURRENT MAP',0)
          else 
c            get the time of the map
             read(line,'(6i6)',iostat=ioerr) myr,month,mday,mhr,min,msec
             ion_time(ntion) = real(idoy(myr,month,mday)) + mhr/24. 
     .                        + min/1440. + msec/86400.  
cd             print *,'read new epoch ',myr,month,mday,mhr,min,msec
cd     .          ,ion_time(ntion)
          endif     
c         read in the values    
          do j=1,nilat
            read(iuf,'(a)',iostat=ioerr) line  
            if( ioerr.ne.0.or.line(61:68).ne.'LAT/LON1') 
     .        call report_stat('FATAL','MODEL','read_ionex',' ' 
     .          ,'Line after EPOCH not LAT/LON1',ioerr)
            read(line,'(2x,5f6.1)',iostat=ioerr) 
     .            lat,lon1,lon2,dlon,ion_rad
cd            print *,'read lat lon1 lon2 dlon ',lat,lon1,lon2,dlon  
            if( ioerr.ne.0)  call report_stat('FATAL','MODEL'
     .           ,'read_ionex',' ','Error reading LAT/LON1/LON2 line'
     .           , ioerr) 
            nlatcnt = nlatcnt + 1   
            if( nlatcnt.gt.nilat ) then
               write(message,'(a,i3)') 
     .          '# lat values gt predicted #: ',nilat
               call report_stat('FATAL','MODEL','read_ionex',' '
     .                          ,message,0)
            endif
            nlonlines = (nilon-1)/16 + 1
cd            print *,'nlonlines ',nlonlines
            do il=1,nlonlines  
              ival1 = 16*il - 15
              ival2 = 16*il
              if( ival2.gt.nilon) ival2 = nilon
              read(iuf,'(a)',iostat=ioerr) line
              if( ioerr.ne.0 )  
     .          call report_stat('FATAL','GRDTAB','get_ion_grid'
     .                ,' ','Error reading longitude lines',ioerr)      
cd                print *,'il ival1 ival2  ',il,ival1,ival2
              read(line,'(16i5)',iostat=ioerr) (lonval(i),i=ival1,ival2)
              if( ioerr.ne.0)  call report_stat('FATAL','MODEL'
     .           ,'read_ionex',' ','Error reading TEC values',ioerr) 
            enddo    
            do i=1,nilon
              ion_val(i,nlatcnt,ntion) = float(lonval(i))*10.**iexp      
            enddo
cd            print *,'nilon,nlatcnt ntion ion_val '
cd     .          ,nilon,nlatcnt,(ion_val(i,nlatcnt,ntion),i=1,nilon)
          enddo
c       endif on reading start/end of TEC map 
        endif

c     enddo on all maps (over time)
      enddo
                     
c.. Save the grid coordinates in common (assume regular spacing)
      
      ilat1 = lat1
      dilat = dlat
      ilon1 = lon1  
      dilon = dlon
C      print *,'ilat1 dilat ilon1 dilon ',ilat1,dilat,ilon1,dilon                                     

C      print *,'End reading IONEX ntion ',ntion

      return
      end


