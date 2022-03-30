c     Program to extract stations from a glist output and form use_site, 
c     unify_apr,  equate, and sig_neu lists for globk and glorg
c      
c     R. King 24 November 2001; last modified 31 October 2013
c                    

c     The program is command-line driven with the form:

c         glist2cmd [glist-file] [min-span] [min-num] [edit_file] colloc [use_file] [uni-file] [vel-file] [eqdef-file]

c           usually invoked by sh_glist2cmd
                   
c       Required:
c          glist-file   Name of input file from program glist

c       Optional:
c          min-span     Minumum number of span in years for site to be used (default 0)
c          min-num      Minumum number of observations (h-files) for site to be used (default 1)
c          edit_file    Input controls on editing of list (format below)
c          colloc       Include sites collocated with sites included by min-span and min-num   
c          use-file     Name of output use_site file for globk (default is 'glist.use_site')
c          unify        Create unify_pos_mode and unify_vel_mode files (default is to omit)
c          vel-file     Name of output vel file for plotting (default is to omit)
c          eqdef-file   Name of output file of site list for input to sh_make_eqdef
                     
c          Collocation distance hard-wired to 2 km

c          Format of edit_file (column 1 blank if not comment)
c           USE site1 site2 site3 ..
c           XCL site1 site2 site3 ..
c           BOX minlat minlon maxlat maxlon
c          (case insensitive, multiple USE/XCL commands allowed, 8 sites per line)
     
      implicit none                                                    
                                                                
      integer*4 maxsit,maxgrp
      parameter(maxsit=5000,maxgrp=30)

      logical autocol,contsig,found,eof,usesite(maxsit)
     .    ,skip,unify,endsite,startgroup, endgroup,vel,eqdef,sort
          
      integer*4 iun,iuedit,iuse,iusig,iunipos,iunivel,iuvel,iueqdef
     .        , mincobs,ioerr,is,iss,istart,istop,nusesites,ncontsites
     .        , numsit,iarg,js,j,unisite(maxsit),isu,lastuni,nuse,nxcl
     .        , jp1, jp2,numobs(maxsit)
c debug
      integer*4 i
                                                          
c  Functions
      integer*4 iclarg,nblen
      logical check_use 
      character*1 velmode

      character*1 ans,col1
      character*8 site(maxsit),grpsite(maxgrp),grppos1(maxgrp)
     .          , grppos2(maxgrp)
      character*32 glistf,editf,usesitef,signeuf,uniposf,univelf,velf
     .          , eqdeff,arg
      character*68 line
                      
      real*4 lon(maxsit),lat(maxsit),ht,nsig,esig,usig
     .     , pos(3,maxsit),dur(maxsit),dist,start(maxsit),stop(maxsit)         
     .     , evel,nvel,hvel,eres,nres,hres,rho,hsig,grpdur(maxgrp)
     ,     , grpdur1(maxgrp),grpdur2(maxgrp),colloc_dist
         
      include 'glist2cmd.h'

c     initialize 
      iun = 1
      iuse = 2
      iusig = 3 
      iuedit = 4
      iunipos = 7
      iunivel = 8
      iuvel = 9
      iueqdef = 10
      autocol = .false. 
      contsig = .false. 
      vel = .false. 
      eqdef = .false.
      found = .false.     
      do is = 1,maxsit
        usesite(is) = .false.
        numobs(is) = 0
        dur(maxsit) = 0.  
        unisite(is) = 0
      enddo             
      colloc_dist = 2000.
      minlat = -90.
      maxlat = 90. 
      minlon = 0.
      maxlon = 360.
      glistf = ' '
      editf = ' '  
      usesitef = 'usesite_file'  
      signeuf  = 'signeu_file'        
      uniposf = 'unify_pos_mode'
      univelf = 'unify_vel_mode'
      velf = 'vel_file'

c Get the run-string arguments and open the files

      iarg = iclarg(1,glistf)
      if( iarg.gt.0 ) then
        open(unit=iun,file=glistf,status='old',iostat=ioerr) 
        if( ioerr.ne.0 ) then
          write(*,'(a,a,i5)') 'Error opening input file ',glistf,ioerr   
          stop
        endif  
        write(*,'(a,a)') ' Opened input file ',glistf 
      else
        call proper_runstring('glist2cmd.hlp','glist2cmd',1)
      endif
      iarg = iclarg(2,arg)
      if( iarg.gt.0 ) then
        read(arg,*) span
      endif
      iarg = iclarg(3,arg)
      if( iarg.gt.0 ) then
        read(arg,*) minobs
      endif         
      write(*,'(a,f5.2,a,i3)') "Minimum span = "
     .     ,span,"  Minimum # = ",minobs

      iarg = iclarg(4,arg)      
      if( nblen(arg).gt.4 ) then
c      if( iarg.gt.0 ) then
          editf = arg(1:nblen(arg))
          open(unit=iuedit,file=editf,status='old',iostat=ioerr)
          if( ioerr.ne.0 ) then
            write(*,'(a,a,i5)') 'Error opening edit file ',editf,ioerr
            stop
          endif  
          write(*,'(a,a)') ' Opened edit file ',editf
c        endif  
      endif

      iarg = iclarg(5,arg) 
      if( arg(1:2).eq.'co' ) then
        autocol = .true.
        write(*,'(a)') "Including collocated sites"
      endif
            
      iarg = iclarg(6,arg) 
      if( nblen(arg).gt.4 ) then
        usesitef = arg(1:nblen(arg))
      else
        usesitef = 'glist.use_site'
      endif
      open(unit=iuse,file=usesitef,status='unknown',iostat=ioerr) 
      if( ioerr.ne.0 ) then
        write(*,'(a,a,i5)') 'Error opening usesite filefile '
     .        ,usesitef,ioerr   
        stop 
      else
         write(*,'(a,a)') ' Opened output file ',usesitef 
      endif 

      iarg = iclarg(7,arg)
      if( arg(1:2).eq.'un') then 
        unify = .true.      
        open(unit=iunipos,file=uniposf,status='unknown',iostat=ioerr)
        if( ioerr.eq.0 ) then 
          write(*,'(a,a)') ' Opened output unify pos_mode file '
     .                           ,uniposf       
          write(iunipos,'(a)') '* Entries for unify pos_mode file '
          write(iunipos,'(a)') '*  Caution: list may contain errors'
          write(iunipos,'(a)') ' POS_MODE'
        else
          write(*,'(a,a)') "Error opening unify pos_mode file: "
     .                        ,uniposf
        endif  
        open(unit=iunivel,file=univelf,status='unknown',iostat=ioerr)
        if( ioerr.eq.0 ) then 
          write(*,'(a,a)') ' Opened output unify vel_mode file '
     .                           ,univelf       
          write(iunivel,'(a)') '* Entries for unify vel_mode file '
          write(iunivel,'(a)')     '*  Caution: list may contain errors'
          write(iunivel,'(a)') ' VEL_MODE'
        else
          write(*,'(a,a)') 'Error opening unify vel_mode file: '
     .                   ,  univelf,ioerr
        endif  
      endif     
 

      iarg = iclarg(8,velf) 
      if( nblen(velf).gt.4) then
        vel = .true.      
        open(unit=iuvel,file=velf,status='unknown',iostat=ioerr)
        if( ioerr.eq.0 ) then 
          write(*,'(a,a)') ' Opened output vel file ',velf 
        else
          write(*,'(a,a)') "Error opening vel file: ",velf,ioerr
        endif
      endif

      iarg = iclarg(9,eqdeff)  
      if( nblen(eqdeff).gt.4 ) then
        eqdef = .true.      
        open(unit=iueqdef,file=eqdeff,status='unknown',iostat=ioerr)
        if( ioerr.eq.0 ) then 
          write(*,'(a,a)') ' Opened output eqdef file ',eqdeff
        else
          write(*,'(a,a)') "Error opening eqdef file: ",eqdeff,ioerr
        endif
      endif
                          

c** Code not implemented  to write sig_neu commands for continuous stations
c      write(iuse,'(a)') '*'
c      write(*,*) 'Do you want sig_neu for continuous stations?'
c      read(*,*) ans
c      if( ans.eq.'y') contsig = .true.    
c      if( contsig ) then   
c        open(unit=iusig,file=signeuf,status='unknown',iostat=ioerr) 
c        if( ioerr.eq.0 ) then 
c          write(*,'(a,a)') ' Opened output file ',signeuf  
c        else
c          write(*,'(a,a,i5)') 'Error opening usesite file '
c     .        ,signeuf,ioerr   
c          stop
c        endif                
c        write(*,*) 'Enter minimum # of observations for continuous'
c        read(*,*) mincobs
c        write(*,*) 'Enter sig_neu values ( N E U milimeters)'
c        read(*,*) nsig,esig,usig   
c        write(iusig,'(a,i4)') 
c     .       '* sig_neu list from GLIST2CMD  min obs =',mincobs   
c        write(iusig,'(a)') '*'
c      endif                
                         

c Write file headers

      write(iuse,'(3a,i3,a,f5.2)') 
     .              '* GLIST2CMD  input file: ',glistf(1:nblen(glistf)),
     .         '  min obs =',minobs,'   min span = ',span        
      if( eqdef ) then
         write(iueqdef,'(2a)') 
     .              '* GLIST2CMD  input file: ',glistf(1:nblen(glistf))
      endif
      if( autocol )   write(iuse,'(a)') 
     .     '*   --plus sites within 2 km of those meeting the criteria'    
      if( vel ) then
        write(iuvel,'(a)') ' SUMMARY VELOCITIES FROM GLIST2CMD'
        write(iuvel,'(3a)') '  Long.     Lat.        E & N Rate '
     .   ,'     E & N Adj.      E & N +-   RHO        H Rate  '
     .   ,' H adj.    +-  SITE'
        write(iuvel,'(2a)') '  (deg)    (deg)          (mm/yr)'
     .  ,'       (mm/yr)       (mm/yr)                 (mm/yr)'
        evel = 0.
        nvel = 0.
        hvel = 0
        eres = 0.
        nres = 0.
        hres = 0.
        hsig = 10.
        rho = 0.
      endif
 

c Read all the values into storage
       
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
c              Sometime between Nov 2004 and Sep 2005, the glist format was
c              changed to add a third decimal place to the start/stop times 
cd               read(line,'(1x,3f9.4,i5,2f9.3,f6.2,1x,a8)' 
cc** rwk 161029 fixed for current format
                 read(line,'(1x,2f9.4,1x,f9.4,i5,2f9.3,f7.3,1x,a8)'
     .                           ,iostat=ioerr) 
     .             lon(is),lat(is),ht,numobs(is),start(is),stop(is)
     .           , dur(is),site(is)      
               if( ioerr.ne.0 .or. stop(is).lt.1985.) then
                   write(*,*) 'Error decoding line, ioerr: ',ioerr 
                   write(*,*) 'Old-format glist file? , stop'
                   stop
               endif
               call getpos(lon(is),lat(is),ht,pos(1,is)) 
          else
            write(*,*) '# sites > maxsit, STOP'
            stop 
          endif
        endif
      enddo  
      numsit = is         
cd      print *,'start numsit usesite(1-3) ',numsit,(usesite(is),is=1,30)
                     

c  Read the input edit (use-list) into storage

      call read_edits( iuedit )
                      
c  Now loop through all the stations, marking the sites to use or unify
       
      do is = 1, numsit  
       usesite(is) = 
     .     check_use( site(is),lat(is),lon(is),numobs(is),dur(is) )   
      enddo
cd      print *,'aft checkr usesite(1-30) ',(usesite(is),is=1,30)
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
                if( dist.le.colloc_dist  )  then
                   usesite(is) = .true.  
                endif
              endif   
            enddo
          endif  
        enddo
       endif 
cd       print *,'after autocol usesite(1-30) ',(usesite(is),is=1,30)

      if( unify ) then    
        endsite = .false.
        is = 1
        isu = 1
        do while(.not.endsite)   
c         keep unisite = 0 for singular sites, then set a common numeral 
c         for sites to be grouped; after grouping, make the site with 
c         the longest span primary.   
          if( (is+1).le.numsit ) then        
            dist = sqrt( (pos(1,is)-pos(1,is+1))**2 +
     .                   (pos(2,is)-pos(2,is+1))**2 +
     .                   (pos(3,is)-pos(3,is+1))**2 )   
            if( dist.le.colloc_dist ) then
              unisite(is) = isu
              unisite(is+1) = isu 
            else
              isu = isu + 1
            endif
            is = is + 1   
          else
            endsite = .true. 
          endif  
        enddo               
c       now group and print the sites
        lastuni = 0   
        js = 1
        do is = 1,numsit +1   
cd          print *,'LOOP is ',is
cd          if( is.gt.30 ) stop 
cd         print *,'is usesite ',is,usesite(is)        
          if( usesite(is) ) then   
cd            print *,'is lastuni unisite js',is,lastuni,unisite(is),js
            if( unisite(is).ne.lastuni ) then 
               if( js.gt.1 ) then
c                sort each sub-group to make the site with the longest span primary and write it out
c                but only if at least one span is non-zero
                 sort = .false.
                 do i=1,js
                   if(grpdur(i).gt.0. ) sort = .true.
                 enddo
                 if( sort ) call sortgrp( maxgrp,js,grpsite,grpdur )
                 write(iunivel,'(1x,30(a8,1x))') (grpsite(j),j=1,js) 
c                for pos_mode, find the sites with matching 4-character codes and write them out
                 call write_pos( maxgrp,js,grpsite,grpdur,iunipos )
               endif
                js = 1  
               grpsite(1) = site(is) 
               grpdur(1) = dur(is)    
cd              print *,'Restarting js=1 ',grpsite(js)
            elseif( unisite(is).ne.0 ) then 
              js = js + 1
              grpsite(js) = site(is) 
              grpdur(js) = dur(is)
cd              print *,'Adding site is js ',is,js,grpsite(js)
            endif
            lastuni = unisite(is)
          endif
        enddo
      endif


c     Now write out the lists

c*      write(*,*) 'Do you want to alphabetize the site list?  (y/n) '
c*     read(*,*) ans
c*      print *,'Pause'    
c*       ans = ' ' 
c*       print *,' ans set blank' 
c*       read(*,'(a1)',iostat=ioerr) ans
c*       print *,'read ans ',ans
c      read(*,*,iostat=ioerr) ans
c*       if( ioerr.ne.0 ) then
c*          print *,'ioerr ',ioerr
c*          print *,'ans ',ans 
c*        else
c*          print *,'ioerr ok'    
c*          print *,'ans ',ans
c*       endif
c*      if( ans.eq.'y' ) then  
c        call sorta( maxsit,numsit,site,usesite,start,stop,dur,numobs
c     .            , lon,lat)  
c*      endif
      do is = 1,numsit 
        if( usesite(is) ) then
          write(iuse,'(a,a8,a,2f8.1,f5.1,i5,2f9.2)') 
     .       ' use_site ',site(is),' ! '
     .       ,start(is),stop(is),dur(is),numobs(is)  
     .       ,lon(is),lat(is)       

          if( eqdef ) then
            write(iueqdef,'(1x,2f9.4,1x,a8)') lon(is),lat(is),site(is)
          endif
     
          if (vel ) then
            col1 = ' ' 
c           set the sigma based on span length
c RWK 140626: For plotting with a velocity solution, it's better to 
c             have the added poorly observed have zero sigma (dots only),
c             so comment this out, at least temporarly
c           if( dur(is).lt.1.0 ) then
c             esig = 10.
c             nsig = 10.               
c             col1 = 'x'
c           elseif( dur(is).ge.5. ) then
c             esig = 1.
c             nsig = 1.
c           else
c             esig = 5./dur(is)
c             nsig = 5./dur(is)
c           endif                                    
            esig = .05
            nsig = .05
            write(iuvel,'
     .    (a1,f8.3,1x,f8.3,1x,4f8.2,2f8.2,f7.3,2x,f8.2,2f8.2,1x,a8)'
     .             ,iostat=ioerr)
     .              col1,lon(is),lat(is),evel,nvel,eres,nres,esig,nsig
     .             ,rho
     .             ,hvel,hres,hsig,site(is)
          endif    
          
        endif
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
                                                                             
c---------------------------------------------------------------------------

      Subroutine read_edits( iuedit )

      implicit none
                          
      integer*4 iuedit,indx,ioerr,idum,i
      logical eof       
      character*3 keywrd
      character*80 line


      include 'glist2cmd.h'
        
cd      print *,'in read_edits iuedit ',iuedit         
      eof = .false.      
      numu = 0 
      numx = 0     
      do while (.not.eof ) 
        keywrd = ' '      
        read(iuedit,'(a)',iostat=ioerr) line
        indx = 2
        if( ioerr.eq.-1 ) then
          eof = .true.
        elseif( ioerr.ne.0 ) then
          write(*,'(a,i5)') 'Error reading edit file ',ioerr                  
        elseif( line(1:1).ne.' ' ) then
          continue 
        else          
cd          print *,'calling read_line indx ',indx
          call read_line(line,indx,'CH',ioerr,idum,keywrd)
          call uppers(keywrd)  
cd          print *,'keywrd ',keywrd
          if( keywrd.eq.'USE' ) then  
            do while(ioerr.eq.0)
              numu = numu + 1  
cd              print *,'reading numu ',numu 
              if(numu.gt.maxedit ) then
                write(*,'(a,i4)') '# sites > maxedit ',maxedit
                stop
              endif   
              ulist(numu) = ' '
              call read_line(line,indx,'CH',ioerr,idum,ulist(numu))
              if( ioerr.eq.0 ) then
                call uppers(ulist(numu))                
c** leave blank for wildcard  if( ulist(numu)(5:8).eq.'    ' ) ulist(numu)(5:8)='_GPS'
              else
                numu = numu - 1  
cd                print *,'end numu ',numu
              endif
            enddo
          elseif( keywrd.eq.'XCL' ) then
            do while(ioerr.eq.0)
              numx = numx + 1           
cd              print *,'reading numx ',numx
              if(numx.gt.maxedit ) then
                write(*,'(a,i4)') '# sites > maxedit ',maxedit
                stop
              endif                   
              xlist(numx) = ' '
              call read_line(line,indx,'CH',ioerr,idum,xlist(numx))
              if( ioerr.eq.0 ) then
                call uppers(xlist(numx))                                
c** leave blank for wildcard  if( xlist(numx)(5:8).eq.'    ' ) xlist(numx)(5:8)='_GPS'
              else
                numx = numx - 1 
cd                print *,'end numx ',numx
              endif
            enddo
          elseif( keywrd.eq.'BOX') then
            read(line(6:80),*,iostat=ioerr) minlat,minlon,maxlat,maxlon
            if( ioerr.ne.0 ) then
              write(*,'(a)') 'Error reading bounding box: ',line
              stop
            endif   
          endif
        endif
      enddo           
      write(*,'(a,4f9.3)') 'Bounding box: ',minlat,minlon,maxlat,maxlon
      write(*,'(a,8(1x,a8))') 'Sites to always include: '
     .     ,(ulist(i),i=1,8)
      if(numu.gt.8) write(*,'((25x,8(1x,a8)))') (ulist(i),i=9,numu)
      write(*,'(a,8(1x,a8))') 'Sites to exclude: '
     .     ,(xlist(i),i=1,8)
      if(numx.gt.8) write(*,'((18x,8(1x,a8)))') (xlist(i),i=9,numx)
      
      return
      end

c---------------------------------------------------------------------------

      Function check_use( site,lat,lon,numobs,dur )  

               
      implicit none

      logical check_use  
      integer*4 numobs,i
      character*8 site
      real*4 lat,lon,dur

      include 'glist2cmd.h' 

      check_use = .false.
               
cd      print *,'CHECK_USE ',site,lat,lon,numobs,dur

c  Is the site within the bounding box?      
                                   
      if( lat.ge.minlat.and.lon.ge.minlon .and.
     .    lat.le.maxlat.and.lon.le.maxlon ) check_use = .true.
      

c  Are the number and span of observations sufficient?

      if( check_use ) then
         if( numobs.lt.minobs .or. dur.lt.span ) then
           check_use = .false.           
         endif
      endif
cd      print *,'aft obs dur box ',check_use          
    
c  Is the site on the explict use list?

      if( .not.check_use ) then
        do i=1,numu     
          if( site(1:4).eq.ulist(i)(1:4) ) then
            if( site(5:8).eq.ulist(i)(5:8) .or.
     .         ulist(i)(5:8).eq.'    ' )
     .      check_use = .true.
          endif
        enddo
      endif                           
cd      print *,'aft explicit use ',check_use

c  Is the site on the explicit exclude list?

      if( check_use ) then
        do i=1,numx   
          if( site(1:4).eq.xlist(i)(1:4) ) then
            if( site(5:8).eq.xlist(i)(5:8) .or.
     .         xlist(i)(5:8).eq.'    ')
     .        check_use = .false. 
          endif
        enddo
      endif           
cd      print *,'aft explicit unuse ',check_use

      return
      end

c----------------------------------------------------------------------------

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
     
c-----------------------------------------------------------------------------
      Subroutine sorta( ndim,n,names,use,start,stop,dur,numobs
     .                , lon,lat )
                      
      implicit none
 
      integer*4 ndim,n,i,j
     .        , numobs(ndim),nbuf1,nbuf2 
      real*4 start(ndim),stop(ndim),dur(ndim)   
     .      ,lon(ndim),lat(ndim)
     .      ,sbuf1,sbuf2,ebuf1,ebuf2,dbuf1,dbuf2
     .      ,lnbuf1,lnbuf2,ltbuf1,ltbuf2
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
           lnbuf1 = lon(j)
           lnbuf2 = lon(j+1)
           ltbuf1 = lat(j)
           ltbuf2 = lat(j+1)
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
             lon(j) = lnbuf1
             lon(j+1) = lnbuf2
             lat(j) = ltbuf1
             lat(j+1) = ltbuf2
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
             lon(j) = lnbuf2
             lon(j+1) = lnbuf1
             lat(j) = ltbuf2
             lat(j+1) = ltbuf1
           endif 
         enddo
      enddo   

      return
      end   

c---------------------------------------------------------

      Subroutine sortgrpa( ndim,n,names,dur )
                      
      implicit none
 
      integer*4 ndim,n,i,j   
      character*4 names(ndim),buf1,buf2      
      real*4 dur(ndim),dbuf1,dbuf2

c     Sort alphabetically
               
      do i = 1,n-1 
        do j = 1,n-i
           buf1 = names(j)
           buf2 = names(j+1) 
           dbuf1 = dur(j)
           dbuf2 = dur(j+1)
           if( lle(buf1,buf2) ) then  
             names(j) = buf1
             names(j+1) = buf2 
             dur(j) = dbuf1
             dur(j+1) = dbuf2 
           else    
             names(j) = buf2
             names(j+1) = buf1 
             dur(j) = dbuf2
             dur(j+1) = dbuf1 
           endif 
         enddo
      enddo   

      return
      end          

       

c---------------------------------------------------------
      Subroutine sortgrp( ndim,n,site,dur )
    
      implicit none

      integer*4 ndim,n,i,j
      character*8 site(ndim),site1,site2
      real*4 dur(ndim),dur1,dur2
                             
c     Sort by duration 

      do i = 1,n-1
        do j = 1,n-i
          dur1 = dur(j)
          dur2 = dur(j+1)
          site1 = site(j)
          site2 = site(j+1)
          if( dur1.ge.dur2 ) then
            dur(j) = dur1
            dur(j+1) = dur2
            site(j) = site1
            site(j+1) = site2      
          else
            dur(j) = dur2
            dur(j+1) = dur1
            site(j) = site2
            site(j+1) = site1    
          endif
        enddo
      enddo
      return
      end

c----------------------------------------------------------------
               

      Subroutine write_pos( maxgrp,numsit,site,dur,iunipos )
              
c     Write out pos_mode entries for matching 4-character names within the group

      integer*4 maxgrp,numsit,unisite(maxgrp),iunipos,is,isu
     .        , lastuni,js
      real*4 dur(maxgrp),grpdur(maxgrp)
      character*8 site(maxgrp),grpsite(maxgrp)
      character*4 site4(maxgrp) 
      logical endsite,sort

c      sites need to be alphabetized for this to work
c      (but resort below after selection to get longest span primary)
       do i=1,numsit
         site4(i) = site(i)(1:4)        
       enddo
       call sortgrpa(maxgrp,numsit,site4,dur)                             
       do i=1,numsit + 1
         unisite(i) = 0
       enddo 
       endsite = .false.
       is = 1
       isu = 1
       do while( .not.endsite )   
cd         print *,'SUB numsit is ',numsit,is
c        keep unisite = 0 for singular sites, then set a common numeral 
c        for sites to be grouped; after grouping, make the site with 
c        the longest span primary.   
         if( (is+1).le.numsit ) then        
           if( site(is+1)(1:4).eq.site(is)(1:4) ) then
             unisite(is) = isu
             unisite(is+1) = isu 
           else
             isu = isu + 1
           endif
           is = is + 1   
         else
           endsite = .true. 
         endif  
       enddo      
cd       print *,'SUB is unisite ',is,(unisite(i),i=1,is)         
c     now group and print the sites
      lastuni = 0   
      js = 1
      do is = 1,numsit +1   
cd        print *,'SUB loop is ',is
cd        print *,' is lastuni unisite js',is,lastuni,unisite(is),js
        if( unisite(is).ne.lastuni ) then 
           if( js.gt.1 ) then  
c            sort each sub-group to make the site with the longest span primary and write it out
c            but only if at least one span is non-zero
             sort = .false.
             do i=1,js
               if(grpdur(i).gt.0. ) sort = .true.
             enddo
              if( sort) call sortgrp( maxgrp,js,grpsite,grpdur )
              write(iunipos,'(1x,30(a8,1x))') (grpsite(j),j=1,js) 
            endif
             js = 1  
            grpsite(1) = site(is)  
            grpdur(1) = dur(is)
cd            print *,'SUB restarting js=1 ',grpsite(js)
        elseif( unisite(is).ne.0 ) then 
          js = js + 1
          grpsite(js) = site(is)  
          grpdur(js) = dur(is)
cd         print *,'SUB Adding site is js ',is,js,grpsite(js),grpdur(js)
        endif
          lastuni = unisite(is)
       enddo
                            
       return
       end

c---------------------------------------------------------------------------
         
CTITLE RCPAR
 
      integer*4 function rcpar( iel, arg )

 
*     Routine to emulate RCPAR using the igetarg UNIX subroutine
*     Modified to use getarg
 
*         iel       - Element of runstring to get
*       igetarg     - UNIX routine to read runstring
*       len_arg     - Length of arg.
*       trimlen     - Get length of string
*       offset      - Offset to be applied to the passsed element
*                    (0 is assumed to program name, 1 first argument)
 
      integer*4 iel, len_arg, trimlen, offset
 
*             arg   - Arg of runstring
 
      character*(*) arg
      character*4 test_arg
      
      data offset / -1 /
 
****  Get length of argument and runstring
* MOD TAH 010610: To see where the count starts for getarg
      if( offset.lt.0 ) then
          call getarg(0, test_arg)
	  len_arg = trimlen(test_arg)
	  if( len_arg.eq.0 ) then
	      offset = 1
	  else
	      offset = 0
	  end if
      end if
      
      len_arg = LEN(arg)
      call getarg( iel+offset, arg )
      rcpar = trimlen( arg )
 
***** Thats all
      return
      end

