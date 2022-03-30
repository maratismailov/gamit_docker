c     Program to extract stations from a glist output and form use_site, 
c     equate, and sig_neu lists for globk and glorg
c      
c     R. King 24 November 2001
c
      implicit none
       
      integer*4 maxsit
      parameter(maxsit=5000)

      logical autocol,contsig,found,eof,usesite(maxsit)
          
      integer*4 iun,iuse,iusig,minobs,mincobs,ioerr,is,iss
     .        , istart,istop,nusesites,ncontsites,numsit 
      integer*4 numobs(maxsit)


      character*1 ans  
      character*8 site(maxsit)
      character*16 glistf,usesitef,signeuf  
      character*64 line
                      
      real*4 span,lon,lat,ht,nsig,esig,usig
     .     , pos(3,maxsit),dur(maxsit),dist,start(maxsit),stop(maxsit)         

c     initialize 
      iun = 1
      iuse = 2
      iusig = 3
      autocol = .false. 
      contsig = .false. 
      found = .false. 
      do is = 1,maxsit
        usesite(is) = .false.
        numobs(is) = 0
        dur(maxsit) = 0.
      enddo  
      usesitef = 'usesite_file'  
      signeuf  = 'signeu_file'        

      write(*,*) 'Enter GLIST file name'
      read(*,*) glistf   
      open(unit=iun,file=glistf,status='old',iostat=ioerr) 
      if( ioerr.ne.0 ) then
        write(*,'(a,a,i5)') 'Error opening input file ',glistf,ioerr   
        stop
      endif  
      write(*,'(a,a)') ' Opened input file ',glistf 
      open(unit=iuse,file=usesitef,status='unknown',iostat=ioerr) 
      if( ioerr.ne.0 ) then
        write(*,'(a,a,i5)') 'Error opening usesite filefile '
     .        ,usesitef,ioerr   
        stop 
      else
         write(*,'(a,a)') ' Opened output file ',usesitef 
      endif

      write(*,*) 'Enter minimum observation span for station (yr)'
      read(*,*) span
      write(*,*) 'Enter minimum number of observations'
      read(*,*) minobs    
      write(iuse,'(a,i3,a,f5.2)') 
     .              '* Sites selected by GLIST2CMD   min obs =',minobs
     .             ,'  min span = ',span 
      write(*,*) 
     .    'Do you want to automatically include colocated sites (y/n)?'
      read(*,*) ans
      if( ans.eq.'y') then
        autocol = .true.  
        write(iuse,'(a)') 
     .     '*   --plus sites collocated with those meeting the criteria'  
      endif
      write(iuse,'(a)') '*'
      write(*,*) 'Do you want sig_neu for continuous stations?'
      read(*,*) ans
      if( ans.eq.'y') contsig = .true.    
      if( contsig ) then   
        open(unit=iusig,file=signeuf,status='unknown',iostat=ioerr) 
        if( ioerr.eq.0 ) then 
          write(*,'(a,a)') ' Opened output file ',signeuf  
        else
          write(*,'(a,a,i5)') 'Error opening usesite filefile '
     .        ,usesitef,ioerr   
          stop
        endif
        write(*,*) 'Enter minimum # of observations for continuous'
        read(*,*) mincobs
        write(*,*) 'Enter sig_neu values ( N E U milimeters)'
        read(*,*) nsig,esig,usig   
        write(iusig,'(a,i4)') 
     .       '* sig_neu list from GLIST2CMD  min obs =',mincobs   
        write(iusig,'(a)') '*'
      endif

c     read all the values into storage
       
c     find the first line of the section   
      eof = .false.
      do while (.not.found .and. .not.eof ) 
         read(iun,'(a)',iostat=ioerr) line 
         if( ioerr.ne.0 ) then
            write(*,*) 'Error looking for SUMMARY ', ioerr 
            if( ioerr.eq.-1 ) eof = .true.
         endif 
         if( line(3:9).eq.'SUMMARY' ) found = .true. 
      enddo  
c     skip the next two lines
      read(iun,'(1x)')
      read(iun,'(1x)')
      is = 0  
      eof = .false.
      do while (.not.eof)            
        line = ' ' 
        read(iun,'(a)',iostat=ioerr) line  
        if( ioerr.eq.-1 .or. line(3:5).eq.'  ' ) then
           eof = .true.
        elseif ( ioerr.ne.0 ) then 
           write(*,*) 'Error reading line, ioerr: ',ioerr
        else
          is = is + 1   
          if( is.le.maxsit) then
               read(line,'(1x,3f9.4,i5,2f8.2,f6.2,1x,a8)'
     .                           ,iostat=ioerr) 
     .        lon,lat,ht,numobs(is),start(is),stop(is),dur(is),site(is) 
               if( ioerr.ne.0 ) then
                   write(*,*) 'Error decoding line, ioerr: ',ioerr 
                   print *,'ErrLINE:',line 
               endif
              if( autocol )  call getpos(lon,lat,ht,pos(1,is))
          else
            write(*,*) '# sites > maxsit, STOP'
            stop 
          endif
        endif
      enddo  
      numsit = is
             
                      
c     Now loop through all the stations, writing use_site commands for all those 
c     that meet the duration and # obs criteria 
       
      do is = 1, numsit  
       if( numobs(is).ge.minobs .and. dur(is).ge.span ) 
     .     usesite(is) = .true. 
      enddo

c     Optionally, add all stations within 1 km of a used stations
       
      if( autocol ) then
        do is = 1,numsit 
          if( .not.usesite(is) ) then
c           since the sites are ordered by longitude, we need to search only a few
            istart = is - 10
            if( istart.le.0 ) istart = 1
            istop = is + 10
            if( istop.gt.numsit ) istop = numsit
            do iss = istart,istop 
              if( usesite(iss)) then
                dist = sqrt( (pos(1,is)-pos(1,iss))**2 +
     .                       (pos(2,is)-pos(2,iss))**2 +
     .                       (pos(3,is)-pos(3,iss))**2 )  
                if( dist.le.1e3  )  usesite(is) = .true.  
              endif   
            enddo
          endif
        enddo
       endif

             
c     Now write out the list

      write(*,*) 'Do you want to alphabetize the site list (y/n) ? '
      read(*,*) ans
      if( ans.eq.'y' ) then  
        call sorta(maxsit,numsit,site,usesite,start,stop,dur,numobs)  
      endif
      do is = 1,numsit 
         if( usesite(is) ) write(iuse,'(a,a8,a,2f8.1,f5.1,i4)') 
     . ' use_site ',site(is),' ! ',start(is),stop(is),dur(is),numobs(is)  
      enddo      


c     Finally, loop through again to write out sig_neu commands for continuous stations

      if( contsig ) then      
        ncontsites = 0
        do is = 1, numsit
          if( numobs(is).ge.mincobs .and. usesite(is) ) then
            write(iusig,'(a,a8,3f7.3)') 
     .          ' sig_neu ',site(is),nsig/1000.,esig/1000.,usig/1000. 
            ncontsites = ncontsites + 1  
          endif
        enddo  
      endif
        
c     Count the number of sites used
            
      nusesites = 0
      do is = 1, numsit
        if( usesite(is) ) nusesites = nusesites + 1
      enddo   

c     Write the summary
                
      write(*,*) ' '
      write(*,'(a,i4,a)') 'There are ',numsit,' sites in the glist file'
      write(*,'(1x,i4,a,a)')   nusesites,' sites written into ',usesitef
      if( contsig ) write(*,'(1x,i4,a,a)') 
     .     ncontsites,' continuous sites written into ',signeuf


      stop
      end    

      Subroutine GETPOS( lon,lat,ht,pos)
c     Convert lon,lat,ht to xyz 
c     **warning:  since we're only interested in small differences here, we're
c                 going to convert as though the coordinates are spherical, even
c                 though they are in fact geodetic (ellipsoidal).  Do not use
c                 this routine for anything else                  

      implicit none            

      real*4 lon,lat,ht,radius,pos(3),convd
        
      convd = (3.1415926/180)                      
      radius = 6378000. + ht  
      pos(1) = radius*cos(lon*convd)*cos(lat*convd)
      pos(2) = radius*sin(lon*convd)*cos(lat*convd)
      pos(3) = radius*sin(lat*convd)     
                                        
      return
      end
     
c
      Subroutine sorta( ndim,n,names,use,start,stop,dur,numobs )
                      
      implicit none
 
      integer*4 ndim,n,i,j
     .        , numobs(ndim),nbuf1,nbuf2 
      real*4 start(ndim),stop(ndim),dur(ndim)
     .      ,sbuf1,sbuf2,ebuf1,ebuf2,dbuf1,dbuf2
      character*8 names(ndim),buf1,buf2      
      logical use(ndim),lbuf1,lbuf2

c     Sort alphabetically
               
      do i = 1,n-1 
        do j = 1,n-i
           buf1 = names(j)
           buf2 = names(j+1)  
           sbuf1 = start(j)
           sbuf2 = start(j+1)
           ebuf1 = stop(j)
           ebuf2 = stop(j+1)
           dbuf1 = dur(j)
           dbuf2 = dur(j+1) 
           nbuf1 = numobs(j)
           nbuf2 = numobs(j+1)
           lbuf1 = use(j)
           lbuf2 = use(j+1)
           if( lle(buf1,buf2) ) then  
             names(j) = buf1
             names(j+1) = buf2  
             start(j) = sbuf1
             start(j+1) = sbuf2
             stop(j) = ebuf1
             stop(j+1) = ebuf2 
             dur(j) = dbuf1
             dur(j+1) = dbuf2 
             numobs(j) = nbuf1
             numobs(j+1) = nbuf2
             use(j) = lbuf1   
             use(j+1) = lbuf2
           else    
             names(j) = buf2
             names(j+1) = buf1  
             start(j) = sbuf2
             start(j+1) = sbuf1
             stop(j) = ebuf2
             stop(j+1) = ebuf1 
             dur(j) = dbuf2
             dur(j+1) = dbuf1     
             numobs(j) = nbuf2
             numobs(j+1) = nbuf1
             use(j) = lbuf2 
             use(j+1) = lbuf1
           endif 
         enddo
      enddo   

      return
      end
