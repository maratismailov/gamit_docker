c     Program to read an ensum SUM file and extract a station once, filtering on min vel sigma, 
c     and outputting three nrms values, the number of observations, and the epoch.   

c     R. King 15 February 2006

c The calling sequence is 
c
c    sortrms[SUM-file] [type] [sort] [min rms]  [maxsig]    
c
c     Required:
        
c     [SUM-file] is the name of the input file            
c     [type] is 'nrms' or 'wrms'                           
c     [sort] is 'a' for alphabetically 'n' for numerical 
c     [min rms]  is the minimum to include (dimensionaless or mm)
            
c     Optional:   

c     [maxsig] is maximum velocity sigma mm/yr  (default 10.) 
            

c The output file is sortrms.out, containing a list ordered alphabetically 
c or numerically by rms.  Sites removed are listed in sortrms.log

      implicit none
           
      integer*4 maxsit 
      parameter(maxsit=5000)
    
      character*1 sort,comp
      character*4 type               
      character*8 arg,site,sitea(maxsit)    
      character*12 sortmsg
      character*50 sumf,outfile,logfile
      character*256 line

      integer*4 iclarg,iarg,nblen,ns,nobs,icomp,indx
     .        , nobsa(maxsit),nsite,msite,nline,ioerr,i,j,k

      real*4 minrms,maxsig,vsig,wrms,nrms,dur,epoch
     .     , vsiga(3,maxsit),wrmsa(3,maxsit),nrmsa(3,maxsit)
     .     , dura(maxsit),epocha(maxsit)


      logical eof,oldsite,usesite(maxsit)
           

c Print help if no arguments

      iarg = iclarg(1,arg)
      if( iarg.eq.0 ) then
        print *,
     .  'sortrms needs at least 4 arguments'      
        stop
      endif
           
c Print some stuff to the screen

      print *,' '
      print *,' Starting sortrms '
      print *,' ' 

c Get the run-string arguments

      iarg = iclarg(1,sumf)
      open(unit=1,file=sumf,status='old',iostat=ioerr) 
      if( ioerr.ne.0 ) then
        print *,'Error opening sumf ',sumf,ioerr
        stop
      else
        print *,'Opened SUM file ',sumf
      endif

      iarg = iclarg(2,arg)   
      if( iarg.gt.0 ) then
        read(arg(1:4),'(a)') type
      else
        print *,'Missing type (nrms or wrms) in runstring'
        stop
      endif
              
      iarg = iclarg(3,arg)   
      if( iarg.gt.0 ) then
        read(arg(1:1),'(a)') sort
      else
        print *,'Missing sort (a or n) in runstring'
        stop
      endif

      iarg = iclarg(4,arg)
      read(arg,'(f8.0)') minrms 
      
      maxsig = 10.
      iarg = iclarg(5,arg)
      if( iarg.gt.0 ) then
        read(arg,'(f8.0)') maxsig 
      endif

c Print the input to the screen
               
      print *,' '                   
      if( sort.eq.'a') then
        sortmsg = 'alphabetical'
      else 
        sortmsg = 'numerical   '
      endif
      print *,'Sort on ',type,', order is ',sortmsg
      write(*,'(a,f4.2,a,f4.1)') ' Minimum rms ',minrms
     .   ,'   Maximum velocity sigma ',maxsig
      print *,' '
         

c Open the output files                           
      outfile = 'sortrms.out'
      open(unit=2,file=outfile,status='unknown',iostat=ioerr)
      if( ioerr.ne.0 ) then
        print*, 'Error opening outfile ',outfile,ioerr
      else
        print *,'Opened output file ',outfile
      endif 
      logfile = 'sortrms.log'
      open(unit=3,file=logfile,status='unknown',iostat=ioerr)
      if( ioerr.ne.0 ) then
        print*, 'Error opening the log file ',logfile,ioerr
      else
        print *,'Opened log file ',logfile
      endif 

    

c Initialize the arrays

      do i=1,maxsit
        sitea(i) = ' '
        epocha(i) = 0.
        dura(i) = 0.
        nobsa(i) = 0
        do j=1,3
          wrmsa(j,i) = 0.
          nrmsa(j,i) = 0.
          vsiga(j,i) = 0.
        enddo 
      enddo

c Read the values from the SUM file and store for sorting
                 
      eof = .false.
      nline = 0 
      nsite = 1
      do while (.not.eof )                                            
        line = '  ' 
        nline = nline + 1
        read(1,'(a)',iostat=ioerr) line     
c         print *,'Read LINE ',line 
        if( ioerr.eq.-1 ) then
          eof = .true.
          print *,'EOF, read ',nline,' lines and ',nsite,' sites'
        elseif ( ioerr.ne.0 ) then
          print *,'Error reading line ',nline,'  ioerr ',ioerr
        elseif ( line(1:1).eq.' ') then
c         all site names begin in column 1
          nline = nline + 1
          continue
        else
          read(line(1:8),'(a8)',iostat=ioerr) site   
          if( nline.eq.1 ) sitea(1) = site
          if( ioerr.ne.0 ) then
            print *,'Error decoding site name from line ',nline
            stop
          endif
          read(line(17:17),*,iostat=ioerr) comp 
          if( ioerr.ne.0 ) then
             print *,'Error reading component at line ',nline
             print *,'LINE: ',line 
             stop
          endif  
c         skip the position values to avoid **** lines
          read(line(18:25),*,iostat=ioerr) ns,nobs
          read(line(72:122),*,iostat=ioerr) vsig,wrms,nrms,dur,epoch   
          if( ioerr.ne.0 ) then
c            may have overflows in some fields
             print *,'Error reading values at line ',nline,' site ',site   
          else
c           see if we already have this site 
            oldsite =.false. 
c            print *,'DEBUG nsite,site ',nsite,site 
            do i=1,nsite
c              print *,'DEBUG checksite i,sitea ',i,site,sitea(i)
              if(sitea(i).eq.site) then
                indx = i
                oldsite = .true.
              endif
            enddo   
c            print *,'DEBUG oldsite ',oldsite
            if( .not.oldsite ) then 
              nsite = nsite + 1
              if(nsite.gt.maxsit) then
                print *,'nsite = ',nsite,' > maxsit'
                stop
              endif
              indx = nsite
            endif                
c            print *,'DEBUG checked nsite indx ',nsite,indx
            sitea(indx) = site
            nobsa(indx) = nobs
            dura(indx) = dur    
            epocha(indx) = epoch
            if( comp.eq.'N') then
              icomp = 1
            elseif( comp.eq.'E' ) then
              icomp = 2
            elseif( comp.eq.'U' ) then
              icomp = 3
            else
              print *,'Unidentified component at line ',nline
              stop
            endif
            wrmsa(icomp,indx) = wrms
            nrmsa(icomp,indx) = nrms
            vsiga(icomp,indx) = vsig
          endif
        endif
      enddo                              
      print *,'Read ',nsite,' sites from ',sumf

c  Array should now be full, remove sites if the smaller of the two horizontal
c  velocity sigmas is too high or if the largest of the rms values is too small
c  (check all if nrms, horizontal only if wrms)
                   
      print *,' '                   
      do i=1,nsite    
        usesite(i) = .true.        
        if( vsiga(1,i).gt.maxsig.and.vsiga(2,i).gt.maxsig) then  
          write(3,'(a,a8,a,2x,2f7.1)') 'Removed ',sitea(i)
     .      ,' vsig(N,E) = ',vsiga(1,i),vsiga(2,i)    
          usesite(i) = .false.
        endif
        if( type.eq.'nrms' ) then
           if( nrmsa(1,i).lt.minrms.and.nrmsa(2,i).lt.minrms.and.
     .         nrmsa(3,i).lt.minrms ) then
             write(3,'(a,a8,a,3f7.1)') 'Removed ',sitea(i)
     .          ,' nrms(N,E,U) = ',nrmsa(1,i),nrmsa(2,i),nrmsa(3,i)    
             usesite(i) = .false.   
           endif
         elseif( type.eq.'wrms') then
           if( wrmsa(1,i).lt.minrms.and.wrmsa(2,i).lt.minrms ) then
             write(3,'(a,a8,a,2f7.1)') 'Removed ',sitea(i)
     .         ,' wrms(N,E) = ',wrmsa(1,i),wrmsa(2,i)    
             usesite(i) = .false.   
           endif
         endif  
      enddo   
      msite = 0
      do i=1,nsite 
        if( usesite(i)) then
          msite = msite + 1
          sitea(msite) = sitea(i)
          nobsa(msite) = nobsa(i)
          dura(msite) =  dura(i)
          epocha(msite) =  epocha(i)
          do k=1,3
            wrmsa(k,msite) = wrmsa(k,i)
            nrmsa(k,msite) = nrmsa(k,i)
            vsiga(k,msite) = vsiga(k,i)
          enddo
        endif
      enddo  
      nsite = msite    
      print *,' ',nsite,' sites remaining'

c Sort alphabetically or by rms

      print *,' '
      print *,'Array filled for ',nsite,' sites '  
      call sortan( maxsit,nsite,sort,sitea,nobsa,dura,epocha,nrmsa
     .            , wrmsa,vsiga )

              
c  Write the output print file 

c     --write the header  
      write(2,'(a,a32,a,a12,a,a4)') '* Sites from ',sumf,'; order is '
     .   ,sortmsg 
      write(2,'(a,f5.1,a,f5.1)') '*  minimum rms ',minrms
     .  ,'  maximum velocity sigma ',maxsig
      write(2,'(a)') '*'
      write(2,'(3a)') ' Site     # Obs   Duration   Date  '
     .                 ,'       NRMS  N  E  U            WRMS  N  E U ',
     .                 '            VSig N  E  U '  
      write(2,'(2a)') '------------------------------------------------'
     .   ,'------------------------------------------------------------'
c     --write the values
      do i=1,nsite
        write(2,'(a8,i7,f11.1,f9.1,3(3x,3f7.1) )') sitea(i),nobsa(i)
     .    , dura(i),epocha(i),(nrmsa(j,i),j=1,3),(wrmsa(j,i),j=1,3)
     .    , (vsiga(j,i),j=1,3)   
      enddo

      end

c------------------------------------------------------------------------

      Subroutine sortan( ndim,n,sort,names,nobs,dur,epoch,nrms,wrms
     .                , vsig )

                               
c       Sort arrays, keying on either names (alphabetical) or nrms (numerical)
c        variable 'sort' is 'a' or 'n' to say which way to sort.

      implicit none
 
      integer*4 ndim,n,i,j,k
     .        , nobs(ndim),nbuf1,nbuf2 

      real*4 dur(ndim),epoch(ndim),nrms(3,ndim),wrms(3,ndim)
     .     , vsig(3,ndim),rbuf1(3),rbuf2(3),wbuf1(3),wbuf2(3)
     .     , ebuf1,ebuf2,dbuf1,dbuf2,vbuf1(3),vbuf2(3),magbuf1,magbuf2

      character*1 sort

      character*8 names(ndim),buf1,buf2      

c     Sort alphabetically
               
      do i = 1,n-1 
        do j = 1,n-i
           buf1 = names(j)
           buf2 = names(j+1)  
           dbuf1 = dur(j)
           dbuf2 = dur(j+1)       
           ebuf1 = epoch(j)
           ebuf2 = epoch(j+1) 
           nbuf1 = nobs(j)
           nbuf2 = nobs(j+1) 
           do k=1,3
             rbuf1(k) = nrms(k,j)
             rbuf2(k) = nrms(k,j+1)
             wbuf1(k) = wrms(k,j)
             wbuf2(k) = wrms(k,j+1)
             vbuf1(k) = vsig(k,j)
             vbuf2(k) = vsig(k,j+1)  
           enddo
           magbuf1 = sqrt(rbuf1(1)**2+rbuf1(2)**2+rbuf1(3)**2)
           magbuf2 = sqrt(rbuf2(1)**2+rbuf2(2)**2+rbuf2(3)**2)
           if( ( sort.eq.'a'. and. lle(buf1,buf2) ) .or.
     .         ( sort.eq.'n'. and. (magbuf1.le.magbuf2) ) ) then  
             names(j) = buf1
             names(j+1) = buf2  
             epoch(j) = ebuf1
             epoch(j+1) = ebuf2
             dur(j) = dbuf1
             dur(j+1) = dbuf2 
             nobs(j) = nbuf1
             nobs(j+1) = nbuf2    
             do k=1,3
               wrms(k,j) = wbuf1(k)
               wrms(k,j+1) = wbuf2(k)  
               nrms(k,j) = rbuf1(k)
               nrms(k,j+1) = rbuf2(k) 
               vsig(k,j) = vbuf1(k)  
               vsig(k,j+1) = vbuf2(k)     
             enddo
           else  
             names(j) = buf2
             names(j+1) = buf1  
             epoch(j) = ebuf2
             epoch(j+1) = ebuf1
             dur(j) = dbuf2
             dur(j+1) = dbuf1     
             nobs(j) = nbuf2
             nobs(j+1) = nbuf1   
             do k=1,3
               wrms(k,j) = wbuf2(k)
               wrms(k,j+1) = wbuf1(k)  
               nrms(k,j) = rbuf2(k)
               nrms(k,j+1) = rbuf1(k) 
               vsig(k,j) = vbuf2(k) 
               vsig(k,j+1) = vbuf1(k)  
             enddo
          endif 
        enddo
      enddo   

      return
      end


