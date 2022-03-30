      Subroutine lread ( site, epoch )  

c     Read site coordinates from an l-file for apr file.  Adapted
c     from J. Ginrich routine of 911009 written for kinematic mode
c     with an l-file only.  Priority changed to take the last rather
c     than the first occurrence of a site.  R. King 020807  
c     Ability to read apr file and use velocities added, calling arguments
c     added, and common revised.  R. King 030215/030509/030605

c     Input and output stored in common'modkin.h'. 
c     Input is 4-character code 'site'; output is spherical coordinates
c     if site is not in the L-file, a warning is issued and the coordinates
c     are not updated.
                 
c     Input  
c       site :  4-character site id requested
c       epoch:  epoch of observations in decimal years
c       iul  :  unit number of the L-file (in commons of units.h)
c       kfflg:  type of L-file (spherical=0, Cartesian=1) (in model.h)
c       iueqrn: unit number of eq_rename (set in open.f and stored in model.h)

c     Output in commons of ../includes/model.h
                                             
c       sitecd  : 4-character site ID found (should be same as input)
c       kepoch0 : apr file epoch     
c       kpos    : Cartesian position from the L-file (m)
c       kvel    : Cartesian velocity from the (apr-style) L-file (m/yr)   
c       kvflg   : flag for velocities to be used (no=0, yes=1)          
c       kstartr(5): start epoch for validity of site coordinates (yr doy hr min sec)
c       kstop(5)  : stop  epoch for validity of site coordinates (yr doy hr min sec)

      implicit none
                                    
      include '../includes/dimpar.h'   
      include '../includes/units.h'
      include '../includes/model.h'
                     
      logical found, eof, first_call
              
      integer*4 maxnam
* MOD TAH 080626: Increased maxnam to 256 (from 10) becuase this dimensions not only
*     only site names but also numbes of earthquakes.
      parameter(maxnam=256)  


      integer*4 dlat, dlon, mlat, mlon, len,rcpar,ioerr, indx
     .        , mline, i 
            
      real*8 epoch,end_epoch,rad,slat,slon,rvalue

      character*1 upper1,cvalue,latflag,lonflag
      character*4 site,upperc
      character*8 sitematch
c *** I cannot get this variable passed through as an argument: 
c     temporarily put it in common here and in rename_match
      integer*4 nlfsite
      character*8 lfsites(maxnam)     
      common /lsite/ lfsites,nlfsite
c     keep the l-file site name local so that the station.info name takes precedence
      character*12 sname
      character*80 prog_name
      character*256 message,line    
      character*256  saveline(maxnam)  
                       
      data first_call/.true./
      save first_call

c     get calling module name for report_stat
      len = rcpar(0,prog_name)
         
      rewind iul
      eof = .false. 
      found = .false.    
      sitnam = ' '    
      sitecd = '    '     
      call uppers(site)
      do i=1,maxnam
        lfsites(i) = 'abcdefgh'
      enddo
                   
c     default is no velocities
      kvflg = 0               

c     initialize times of validity
      if( first_call ) then
        kstartr(1) = 1900
        kstopr(1)  = 2100
        kstopr(2) = 1
        do i=3,5
         kstartr(i) = 0
         kstopr(i) = 0   
        enddo  
cd        print *,'first_call kstopr ',kstopr
        first_call = .false.
      endif
             
c     kfflg set by program that opens the coordinate (L-) file
c       = 0  old-style (spherical) L-file
c       = 1  apr (Cartesian) L-file
 
cd      print *,' LREAD kfflag iul iueqrn ',kfflg,iul,iueqrn
                           
      if( kfflg.eq.0 ) then
c       read an L-file     
        do while (.not.eof ) 
          read(iul,'(a)',iostat=ioerr) line
          if( ioerr.eq.-1 ) then
            eof = .true.
          elseif ( ioerr.ne.0 ) then
            call report_stat('FATAL',prog_name,'lib/lread',' '
     .                      ,'Error reading l-file',ioerr) 
          else                    
            if(upperc(line(1:4)).eq.site) then 
              found = .true.                           
              sitecd = site   
              sitnam(1:12) = line(6:17)
              latflag = line(18:18)
              read(line(19:20),'(i2)',iostat=ioerr) dlat  
              read(line(22:23),'(i2)',iostat=ioerr) mlat
              read(line(25:32),'(f8.5)',iostat=ioerr) slat
              lonflag = line(34:34)
              read(line(35:37),'(i3)',iostat=ioerr) dlon
              read(line(39:40),'(i2)',iostat=ioerr) mlon
              read(line(42:49),'(f8.5)',iostat=ioerr) slon
              read(line(51:63),'(f13.4)',iostat=ioerr) rad  
              if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .          ,'lib/lread',' ','Error decoding l-file line',ioerr)
c             check radius to see if (illegal) geodetic entry
              if( rad.lt.6.d6 ) then
                write(message,'(a,a4,a)') 'Radius for site '
     .              ,sitecd,' < 6000 km --illegal L-file entry'
                call report_stat('FATAL',prog_name,'lib/lread',' '
     .             ,message,0)
              endif                                                        
c             Code for L-files created before coordinate convention change
              if(latflag.eq.' ') latflag=upper1('n')
              if(lonflag.eq.' ') lonflag=upper1('w')
              if(latflag.eq.'-' .or. dlat.lt.0) latflag=upper1('s')
              if(lonflag.eq.'-' .or. dlon.lt.0) latflag=upper1('w')    
c             convert to Cartesian coordinates for modkin.h commons
              call dmsrad( latflag,dlat,mlat,slat,latr )  
              call dmsrad( lonflag,dlon,mlon,slon,lonr )
              call sph2xyz( latr,lonr,rad,kpos ) 
              do i=1,3
                kvel(i) = 0.d0
              enddo  
c              write(*,'(a,i2,3f16.5,1x,3f9.4)'),
c     .              'LREAD kfflg kpos kvel ',kfflg,kpos,kvel
cd              print *,' rad lonr latr ',rad,lonr,latr
cd              print *,' dlon mlon slon,lonflag '
cd     .            ,dlon,mlon,slon,lonflag 
cd              print *,'  dlat mlat slat,latflag '
cd     .            ,dlat,mlat,slat,latflag
            endif  
          endif
        enddo

      elseif ( kfflg.eq.1 ) then
c       read and convert an apr file line   
        nlfsite = 0
        do while (.not.eof )  
          line = ' ' 
          read(iul,'(a)',iostat=ioerr) line    
cd          print *,'ioerr line ',ioerr,line
* MOD TAH 080626: Check to see if EXTENDED globk line.  These are not used
*         but gamit yet and so are commented out here
          call uppers(line)
          if( index(line,'EXTENDED').gt.0 ) then 
              line(1:1) = '*'
          endif
* END MOD TAH 080626:
          if( line(1:1).ne.' ' ) then
c            skip a comment line
             continue
          elseif( ioerr.eq.-1 ) then
            eof = .true.  
cd           print *,'eof set true by eof, ioerr ',ioerr
          elseif ( ioerr.ne.0 ) then
            call report_stat('FATAL',prog_name,'lib/lread',' '
     .                      ,'Error reading apr file',ioerr) 
          elseif( line(1:1).eq.' ' ) then 
            indx = 1  
            asite = ' ' 
            call read_line(line,indx,'CH',ioerr,rvalue,asite) 
            call uppers(asite)
            if( ioerr.eq. 0 ) then   
              if(  asite(1:4).eq.site ) then    
c               found a matching 4-character code; save and see if there are others
                nlfsite = nlfsite + 1 
                found = .true.   
                sitecd = site    
                sitnam(1:8) = asite  
                if( nlfsite.gt.maxnam ) 
     .            call report_stat('FATAL',prog_name,'lib/lread',' '
     .       ,'Too many sites with the same code in apr-style L-file',0)
                saveline(nlfsite)= line 
                lfsites(nlfsite) = asite
              endif
            endif  
          else   
cd              print *,'set eof true by blank, ioerr ',ioerr
            eof = .true.
          endif 
        enddo  
c       if more than one entry for the 4-character code, find the right one
        if( nlfsite.gt.1 ) then 
          if( iueqrn.eq.0 ) then
             write(message,'(a,a4,a)') 'Multiple entries for site '
     .        ,site,' but no eq_rename file; last entry used'
             call report_stat('WARNING',prog_name,'lib/lread',' ',
     .         message,0)
* MOD TAH 050506: Restore the coordinate line for the first occurence of
*            the coordinate.
             line = saveline(nlfsite)
          else                        
c           find line that matches the current epoch
cd            print *,'calling rename_match lfsites ',lfsites
            call rename_match( iueqrn,epoch,sitematch,end_epoch )
cd            print *,'return from rename_match',epoch,sitematch,end_epoch
cd            print *,  'lfsites ',lfsites
            if( sitematch.eq.'        ' ) then   
              mline = 1  
              line = saveline(mline)     
cd              print *,'no match mline line ',mline,line
              indx = 1
              call read_line( line,indx,'CH',ioerr,rvalue,asite )   
              call uppers(asite)
cd              print *,'read asite ',asite  
              write(message,'(a,a8)')  ' ----using ',asite
              call report_stat('WARNING',prog_name,'lib/lread',' '
     .                         ,message,0)   
            else     
cd              print *,'sitematch not blank '
cd              do i=1,nlfsite
cd                print *,i,lfsites(i)
cd              enddo
              do i=1,nlfsite
                if( sitematch.eq.lfsites(i) ) mline = i
              enddo
              line = saveline(mline) 
cd              print *,'entry found mline line ',mline,line
              indx = 1
              call read_line( line,indx,'CH',ioerr,rvalue,asite )   
              call uppers(asite)
cd              print *,'read asite ',asite    
              sitnam(1:8) = asite    
              call decyrs_to_ydhms( epoch,kstartr) 
              call fix_y2k(kstartr(1)) 
              call decyrs_to_ydhms( end_epoch,kstopr)   
cd              print *,'sitnam end_epoch kstopr ',sname,end_epoch,kstopr
              call fix_y2k(kstopr(1))   
cd              print *,'after fixy2k kstopr ',kstopr             
c**rwk 080513: Don't write this message anymore since it's confusing in MODEL
cd              write(message,'(a,a8)') 
cd     .              'Multiple entries in L-file, using ',asite
cd              call report_stat('STATUS',prog_name,'lib/lread',' ' 
cd     .                         ,message,0)
            endif
          endif   
        else
          line = saveline(1)
cd          print *,'saveline ',line(1:80)
        endif                                       
c       reset indx to avoid problem with extra blank lines 
        indx = 1                      
cd        print *,'line used ',line
        call read_line( line,indx,'CH',ioerr,rvalue,asite )   
        call uppers(asite)
cd        print *,'indx,ioerr,asite ',indx,ioerr,asite
        call read_line(line,indx,'R8',ioerr,kpos(1),cvalue)
cd        print *,'indx ioerr,kpos1 ',indx,ioerr,kpos(1)
        call read_line(line,indx,'R8',ioerr,kpos(2),cvalue) 
cd        print *,'indx ioerr,kpos2 ',indx,ioerr,kpos(2)
        call read_line(line,indx,'R8',ioerr,kpos(3),cvalue)   
cd        print *,'indx ioerr,kpos3 ',indx,ioerr,kpos(3)
        if( ioerr.ne.0 ) then  
          write(message,'(a,a4,a)') 'Error extracting position for '
     .        ,site,' from apr-style L-file'
          call report_stat('FATAL',prog_name,'lib/lread',' '
     .      ,message,ioerr)  
        endif
        call read_line(line,indx,'R8',ioerr,kvel(1),cvalue)
        call read_line(line,indx,'R8',ioerr,kvel(2),cvalue)
        call read_line(line,indx,'R8',ioerr,kvel(3),cvalue)  
        if( ioerr.ne.0 ) then
          write(message,'(a,a4,a)') 'Error extracting velocities for '
     .        ,site,' from apr-style L-file'
          call report_stat('WARNING',prog_name,'lib/lread',' '
     .       ,message,ioerr) 
          do i=1,3
             kvel(i) = 0.d0                   
          enddo             
          kepoch0 = epoch
          write(message,'(a,a,a4)') 'Epoch of coords set to '
     .        ,'current epoch and velocities to 0.',site
          call report_stat('WARNING',prog_name,'lib/lread',' '
     .   ,message,ioerr) 
        else
          call read_line(line,indx,'R8',ioerr,kepoch0,cvalue)  
          if( ioerr.ne.0 ) then
            write(message,'(a,a,a4)') 'Error extracting reference '
     .          ,'epoch from apr-style L-file for site ',site
            call report_stat('FATAL',prog_name,'lib/lread',' '
     .         ,message,ioerr)
          endif                     
        endif
      else
        call report_stat('FATAL',prog_name,'lib/lread',' '
     .                  ,'Invalid flag for coordinate-file type',0) 
      endif
      if( .not.found)  then
         write(message,'(a,a4,a)' ) 'Cannot find site code ',site
     .       ,' on L-file  '
         call report_stat('FATAL',prog_name,'lib/lread',' ',message,0)
      endif

      return
      end

c************************************************************************************

      Subroutine rename_match( iueqrn,epoch,sitematch,end_epoch)

c     Match a site read from the L- (apr-) file with the range of validity
c     implied by the eq_rename file

c     Input:  iueqrn           unit number of eq_rename file  
c             epoch            epoch of observations (decimal years)        
c             maxnam           dimension of lfsites
c    In common: 
c             nlfsite          number of L-file entries with matching 4-char codes to be checked
c             lfsites(maxnam)  8-character codes of L-file entries  (upper case)
c       
c     Output: sitematch        8-character code that matches the renames and earthquakes (blank if no match)
c             end_epoch        end time of validity of matching 8-charcter site name (decimal yrs)
                                          
      implicit none
                                      
      integer*4 maxnam
* MOD TAH 080626: Increased maxnam to 256 (from 10) becuase this dimensions not only
*     only site names but also numbes of earthquakes.
      parameter(maxnam=256)  

      integer*4 iueqrn,neq,nren,neqrn,len,rcpar,i,j

      real*8 epoch,end_epoch,eqepoch(maxnam)
     .     , rnstart(maxnam),rnstop(maxnam)
     .     , eqrnstart(maxnam),eqrnstop(maxnam)
               
      character*2 eqcodes(maxnam)      
      character*4 sitecd
      character*8 sitematch     
c *** I cannot get this variable passed through as an argument: 
c     temporarily put it in common here and in rename_match
      integer*4 nlfsite
      character*8 lfsites(maxnam)     
      common /lsite/ lfsites,nlfsite

      character*8 rnsites(maxnam),eqrnsites(maxnam)
      character*256 prog_name,message

      logical neweq,found
           
c     get calling module name for report_stat
      len = rcpar(0,prog_name)
       
cd      print *,'RENAME_MATCH iueqrn,epoch ',iueqrn,epoch 
cd      print *,' nlfsite lfsites ',(lfsites(i),i=1,nlfsite)  
        
   
c  Get a complete list of 8-character sites from the eq_file renames
           
      sitecd = lfsites(1)(1:4)                          
      call get_renames( iueqrn,sitecd,maxnam
     .                , nren,rnsites,rnstart,rnstop ) 
cd      print *,'returned from get_renames '
cd      do i=1,nren
cd        write(*,'(a8,1x,2f14.7)') 
cd     .       rnsites(i),rnstart(i),rnstop(i)
cd      enddo    
cd      print *,' lfsites ',lfsites
                      

c  Get from the l-file sites a list of EQs to check 

      neq = 0    
      do i=1,maxnam
        eqcodes(i) = ' '   
        eqepoch(i) = 2100.d0
      enddo
cd      print *,'checking for EQ names: nlfsite lfsites '
cd     .   ,nlfsite,(lfsites(i),i=1,nlfsite)
      do i=1,nlfsite
        if( lfsites(i)(7:8).ne.'PS' ) then  
          if( neq.eq.0 ) then
            neq = neq + 1  
            if( neq.gt.maxnam ) then 
              write(message,'(a,i3,a)') 
     .          'Number of earthquakes  > maxnam (',maxnam,')'
              call report_stat('FATAL',prog_name,'lib/lread',' '
     .              ,message,0)
            endif
            eqcodes(neq) = lfsites(i)(7:8) 
          else
c           check for duplicates  
            neweq = .true.       
cd            print *,'i lfsites ',i,lfsites(i)
            do j=1,neq  
              if( lfsites(i)(7:8).eq.eqcodes(j) ) neweq=.false.
            enddo    
cd            print *,'neweq ',neweq
            if( neweq ) then
               neq = neq+1
               eqcodes(neq) = lfsites(i)(7:8)    
cd               print *,'neq eqcodes ',neq,eqcodes(neq)
            endif
          endif   
        endif
      enddo         
cd      print *,'neq eqcodes ',neq,(eqcodes(i),i=1,neq)   
                            
c  Read the eq_rename file to get the epochs of the eathquakes
cd      print *,'calling get_equakes '
      if( neq.gt.0 ) call get_equakes(iueqrn,maxnam,neq,eqcodes,eqepoch)
cd      print *,'return from get_equakes '   

c  Merge the earthquakes into the rename site list if necessary      
            
cd      print *,'rename list nren ',nren                                 
cd    do i=1,nren
cd       write(*,'(a8,2f51.7)') rnsites(i),rnstart(i),rnstop(i)                                         
cd      enddo
cd      print *,'eq list neq ',neq
cd      do i=1,neq
cd        write(*,'(a2,f51.7)') eqcodes(i),eqepoch(i)
cd      enddo  
      do i=1,maxnam
        eqrnsites(i) = ' '
        eqrnstart(i) = 0.d0
        eqrnstop(i) = 0.d0
      enddo
      neqrn = 0
      do i=1,nren 
        neqrn = neqrn + 1 
        if( neqrn.gt.maxnam ) then
         write(message,'(a,i3,a)') 
     .      'Number of eqrenames > maxnam (',maxnam,')'
         call report_stat('FATAL',prog_name,'lib/lread',' ',message,0)
        endif
        eqrnsites(neqrn) = rnsites(i)
        eqrnstart(neqrn) = rnstart(i)
        eqrnstop(neqrn) = rnstop(i)  
        if( neq.gt.0 ) then
          do j=1,neq
            if( eqepoch(j).le.rnstart(i) ) then 
              eqrnsites(neqrn)(7:8) = eqcodes(j) 
            elseif( eqepoch(j).gt.rnstart(i) .and.
     .              eqepoch(j).le.rnstop(i) ) then 
c             change the stop time for the site and insert a new site  
              eqrnstop(neqrn) = eqepoch(j)
              neqrn = neqrn + 1         
              if( neqrn.gt.maxnam ) then
                write(message,'(a,i3,a)') 
     .           'Number of eqrenames > maxnam (',maxnam,')'
                 call report_stat('FATAL',prog_name,'lib/lread',' '
     .            ,message,0)
              endif
              eqrnsites(neqrn) = eqrnsites(neqrn-1)
              eqrnsites(neqrn)(7:8) = eqcodes(j) 
              eqrnstart(neqrn) = eqepoch(j) 
              eqrnstop(neqrn) = rnstop(i)     
            endif
          enddo 
        endif
      enddo
cd      print *,'epoch neqrn ',epoch,neqrn
cd      do i=1,neqrn
cd        write(*,'(i2,1x,a8,2f15.9)') 
cd     .     i,eqrnsites(i),eqrnstart(i),eqrnstop(i)
cd     enddo                   

c  Trap the unusual and un-coded case of two events in one day

      do i=1,neqrn-1 
        if( (dabs(eqrnstart(i)-epoch).lt..00274d0 ) .and.
     .      (dabs(eqrnstart(i)-eqrnstart(i+1)).lt..00274d0)  ) 
     .   call report_stat('WARNING',prog_name,'lib/lread',' '
     .    ,'Two eq-rename events within one day--coords may be wrong',0)
      enddo

c  Find the right name for the input epoch and return its end epoch
            
      sitematch = '        '    
      do i=1,neqrn 
        if( epoch.ge.eqrnstart(i).and.epoch.le.eqrnstop(i) ) then
         sitematch = eqrnsites(i)
          end_epoch = eqrnstop(i)   
        endif
      enddo   
cd      print *,' epoch sitematch end_epoch '
cd     .         ,epoch,sitematch,end_epoch 
      found =.false.     
      do i=1,nlfsite
       if( sitematch.eq.lfsites(i) ) found = .true.
      enddo
      if( .not.found ) then
        write(message,'(a,a8,a)') 
     .   'Multiple entries in L-file but no match to eq_file ('
     .      ,sitematch,')'
         call report_stat('WARNING',prog_name,'lib/lread',' ',message,0)
         sitematch = '        ' 
      endif  
cd      print *,'Final sitematch ',sitematch
cd      print *,'   lfsites ',lfsites
      return
      end    
 
c************************************************************************************

      Subroutine get_renames( iueqrn,sitecd,maxnam
     .                      , nren,sites,start,stop )

c     Read the eq_file and get a list of 8-character site names represented
c     by the explicit renames.  Assumes that the original names always end
c     in _GPS (i.e., a single rename is applied for a given interval)

c     Input:  iueqrn   unit number for eq_file
c             sitecd   4-charcter code for the site

c     Output: nren       number of 8-character sites after renames applied
c             sites    list if 8-character sites
c             start   start times (decimal years)
c             stop      stop times (decimal years)

      implicit none
                      
                                   
      integer*4 iueqrn,maxnam
     .        , nren,indx,ioerr,idate(5),doy,len,i
      
      character*4 sitecd 
      character*8 sites(maxnam),asite1,asite2       
      character*16 string
      character*80 line,prog_name

      real*8 start(maxnam),stop(maxnam),rvalue,sod

      logical eof          

c     function
      integer*4 idoy,rcpar
      real*8 decyrs

               
c     get calling module name for report_stat
      len = rcpar(0,prog_name)

c       Always create at least one entry

      nren = 1
      sites(1) = sitecd//'_GPS'
      start(1) = 1900.d0
      stop(1) = 2100.d0

c       Now add entries as renames are found

      rewind(iueqrn)       
      eof = .false.              
      do while( .not.eof )   
        read(iueqrn,'(a)',iostat=ioerr) line
cd        print *,'read iueqrn line ',line 
        if( ioerr.eq.-1 ) then 
          eof = .true.    
        elseif ( ioerr.ne.0 ) then
          call report_stat('FATAL',prog_name,'lib/lread',' '
     .                    ,'Error reading eq_rename file',ioerr) 
        elseif ( line(1:1).eq.' ' ) then
c         remove any comments
          indx = index(line,'!')
          if ( indx.gt.0 ) line(indx:) = ' '
          indx = 1  
          string = ' ' 
          call read_line( line,indx,'CH',ioerr,rvalue,string ) 
cd          print *,'read line string ',line,string
          call uppers(string)
          if( string(1:6).eq.'RENAME' ) then
            call read_line( line,indx,'CH',ioerr,rvalue,asite1 )  
cd            print *,'sitecd read rename asite1 ',sitecd,asite1
            call uppers(asite1)     
            if( asite1(1:4).eq.sitecd ) then 
              nren = nren + 1
* MOD TAH 080626: Added check on number of renames so that bounds are
*             not violated.  If too many simply stop adding new ones
*             and generate warning
              if( nren.gt.maxnam ) then
                  call report_stat('WARNING',prog_name,'lib/lread',' ',
     .               'Too many renames in eq_rename',maxnam)
                  nren = nren - 1
              endif 
              call read_line( line,indx,'CH',ioerr,rvalue,asite2 )
              call uppers(asite2)    
              sites(nren) = asite2  
c             next token can be an h-file string or an integer year
              call read_line( line,indx,'CH',ioerr,rvalue,string )
              read(string(1:4),'(i4)',iostat=ioerr) idate(1) 
              if( ioerr.ne.0) then
c               assume it's an h-file name and skip to the next token
                call read_line( line,indx,'CH',ioerr,rvalue,string )   
                read(string,'(i4)',iostat=ioerr) idate(1)    
c               check for bad read at this point
                if( idate(1).lt.1980.or.idate(1).gt.2100 ) 
     .             call report_stat('FATAL',prog_name,'lib/lread',' '
     .                    ,'Bad year in rename command',ioerr)
              endif    
              do i=2,5        
                call read_line( line,indx,'I4',ioerr,idate(i),string )  
                if(ioerr.ne.0) idate(i) = 1
              enddo  
              if( ioerr.ne.0 ) then
                 call report_stat('WARNING',prog_name,'lib/lread',' '
     .          ,'Error reading start time from eq_rename file',ioerr)    
c               if there is an error, default the start time to 1900
                if( ioerr.ne.0 ) idate(1) = 1900
              endif 
              call fix_y2k(idate(1))
              doy = idoy(idate(1),idate(2),idate(3))
              sod = idate(4)*3600.d0 + idate(5)*60.d0    
              start(nren) = decyrs(idate(1),doy,sod)   
              stop(nren-1) = start(nren)
cd              print *,'nren doy sod start '
cd     .                ,nren, doy,sod,start(nren) 
c             get the end time   
              do i=1,5
                call read_line( line,indx,'I4',ioerr,idate(i),string ) 
                if (ioerr.ne.0 ) idate(i) = 1
              enddo   
cd              print *,'stop idate ioerr ',idate,ioerr
              call fix_y2k(idate(1))
              if( ioerr.ne.0 ) then
                 call report_stat('WARNINGL',prog_name,'lib/lread',' '
     .          ,'Error reading stop time from eq_rename file',ioerr)    
               endif
c               if there is an error, default the end time to 2100 
              if( ioerr.ne.0 ) idate(1) = 2100
              doy = idoy(idate(1),idate(2),idate(3))
              sod = idate(4)*3600.d0 + idate(5)*60.d0 
              stop(nren) = decyrs(idate(1),doy,sod) 
cd              print *,'nren ioerr,doy sod stop '
cd     .               ,nren,ioerr,doy,sod,stop(nren) 
c           endif on matching site
            endif
c         end if on check for 'rename'
          endif
c       end if on error-free read of site name
        endif      
c     end do on records of file
      enddo   
cd      print *,'Finished renames nren ',nren                             
cd      do i=1,nren
cd       write(*,'(a8,1x,2f14.7)') 
cd     .       sites(i),start(i),stop(i)
cd      enddo  
      return
      end

c***********************************************************************************

      Subroutine get_equakes( iueqrn,maxnam,neq,eqcodes,epochs )

c     Read the eq_file and get the start and stop times of for the
c     earthquakes in the input list.
c    

c     Input:  iueqrn   unit number for eq_file
c             neq      number of eathquakes to be found
c             eqcodes  2-character codes for earthquakes requested

c     Output: epochs   times of  earthquakes (decimal years)

      implicit none
                      
                                   
      integer*4 iueqrn,maxnam,neq,indx,ioerr,doy
     .        , rcpar,yr,mon,day,hr,min,len,i
  
      character*2 eqcodes(maxnam),eqcode1
      character*16 string 
      character*80 prog_name
      character*256 line

      real*8 epochs(maxnam),rvalue,sod,lat,lon,dist,mag

      logical eof          

c     function
      integer*4 idoy,nblen
      real*8 decyrs

c     get calling module name for report_stat
      len = rcpar(0,prog_name)
           
cd      print *,'GET_EQUAKES maxnam ',maxnam
cd      print *,'  eqcodes ',eqcodes
      rewind(iueqrn)       
      eof = .false.     
cd      print *,'LREAD iueqrn ',iueqrn
      do while( .not.eof )  
        read(iueqrn,'(a)',iostat=ioerr) line   
cd        print *,'ioerr line ',ioerr,line
        if( ioerr.eq.-1 ) then
          eof = .true.    
        elseif( ioerr.ne.0 ) then
          call report_stat('FATAL',prog_name,'lib/lread',' '
     .                    ,'Error reading eq_rename file',ioerr) 
        elseif( line(1:1).eq.' ' ) then 
          string = ' '
          indx = 1
          call read_line( line,indx,'CH',ioerr,rvalue,string )    
cd          print *,'read_line ioerr string',ioerr,string
          call uppers(string) 
cd        print *,'string ',string    
          if( string(1:6).eq.'EQ_DEF' ) then
            call read_line( line,indx,'CH',ioerr,rvalue,eqcode1 )
cd            print *,'reading eqcode1 ioerr ',eqcode1,ioerr
            call uppers(eqcode1) 
cd            print *,'reading file neq eqcodes ',neq,(eqcodes(i),i=1,neq)
            do i=1,neq  
              if( eqcode1.eq.eqcodes(i) ) then 
cd                print *,'indx nblen ',indx,nblen(line)  
                read(line(indx:nblen(line)),*,iostat=ioerr) 
     .              lat,lon,dist,mag,yr,mon,day,hr,min  
                if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .          ,'lib/lread',' ','Error reading earthquake epoch',ioerr)
cd                print *,'eq line ',lat,lon,dist,mag,yr,mon,day,hr,min
                call fix_y2k(yr)
                doy = idoy(yr,mon,day)
                sod = hr*3600.d0 + min*60.d0 
                epochs(i) = decyrs(yr,doy,sod)  
cd                print *,'i eq epochs ',i,epochs(i)
              endif
            enddo
          endif 
        endif
      enddo 
c     put the list in time order   
cd      print *,'neq eqcodes eqepochs '
cd     .   ,neq,(eqcodes(i),epochs(i),i=1,neq)  
      call sorteq( maxnam,neq,epochs,eqcodes )
cd      print *,'neq eqcodes eq epochs '
cd     .    ,neq,(eqcodes(i),epochs(i),i=1,neq)
cd      print *,'ARGS ',iueqrn,maxnam,neq,epochs,eqcodes
      return              
      end

     
c************************************************************************************8

      subroutine sorteq( maxnam,neq,dates,names)  
                                         
c     Sort a real array with names into ascending order.
c     R. King, based on gamit/solve/sort1i.f
                 
       implicit none

       integer*4 maxnam,neq,i,k

       real*8 dates(maxnam),swap

       character*2 names(maxnam),nameswap

c
      if(neq.le.1) go to 100
      do 50 k=1,32000
         swap = 0.d0
         do 30  i=2,neq
            if(dates(i).ge.dates(i-1)) go to  30
            swap = dates(i)
            nameswap = names(i)
            dates(i)=dates(i-1) 
            names(i)=names(i-1)
            dates(i-1)=swap 
            names(i-1)=nameswap
 30      continue
         if(swap.eq.0.d0) go to 100
 50   continue
 100  continue
      return
      end
 

