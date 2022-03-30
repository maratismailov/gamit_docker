c     Program to extract stations from a GLOBK velocity and (optionally) position
c     print file and form a stab_site list  for repeatabilities, and optionally 
c     an unweight list for hfupd.  No longer accepts a prt or org file directly.
c      
c     R. King 1 February 2007; last modified 30 September 2009

c     The program is command-line driven with the form:
           
c         vel2stab [file] [hvsig] [vvsig] [hpsig] [vpsig] -staball -uw -uwdup 

c       Required:
c          file     Name of input vel file from GLOBK  (from sh_exglk or tsfit)
c                   Note: if positions also to be considered, append the position
c                   output of sh_exglk with its headers

c       Optional: 
c          hvsig        Max/min horizontal velocity sigma for stab/unweight (default 2 mm/yr)   
c          vvsig        Max/min vertical velocity sigma for stab/unweight (default 10 mm/yr)
c          hpsig        Max/min horizontal position sigma for stab/unweight (default 2 mm)
c          vpsig        Max/min vertical position sigma for stab/unweight   (default 10 mm)
c          -staball     Create the stab list with @ (necessary for tsfit) 
c          -uw          Create unweight commands for an hfupd edit file based on sigmas (default no)
c          -uwall       Include in the unweight list any sites that duplicate a 4-character id (default no)  
c
c       Note: Use zero for hvsig and vvsig to skip velocities (position-only file from glred/glorg)
c             Use zero for hpsig and vpsig to skip position (velocity-only file from tsfit)

      implicit none
       
      integer*4 maxsit
      parameter(maxsit=5000)

      logical eof,endvel,unweights,duplicates,uwflg(maxsit)
     .       , usevel,usepos,staball
         
      integer*4 iun,istab,ioerr,is,numsit,goodcount,xvcount,xpcount
     .        , numusit,iuwgt,iarg,iclarg,nscol,nr,len_run 
      integer*4 numobs(maxsit)

* Function
      integer*4 nblen,rcpar

      character*1 ans,col1(maxsit)  
      character*3 comtxt
      character*8 site(maxsit),sitep,usite(maxsit),word,siteout
      character*10 stabtxt,unwgttxt
      character*32 orgvelf,stabsitef,unwgtf,arg
      character*120 line,runstring

                      
      real*4 span,lon(maxsit),lat(maxsit),ht,nvsig(maxsit),evsig(maxsit)
     .     , hvsig(maxsit),npsig(maxsit),epsig(maxsit),hpsig(maxsit)
     .     , evel,nvel,hvel,hvmax,hpmax,vpmax,vvmax,devel,dnvel,dhvel
     .     , epos,npos,hpos,dnpos,depos,dhpos,rho,lonp,latp

* Initialize 
      iun = 1
      istab = 2    
      iuwgt = 3
      stabsitef = 'stabsite_file'  
      unwgtf = 'unweight_file'      
      goodcount = 0 
      hpmax = 2.
      hvmax = 2.
      vpmax = 10.
      vvmax = 10.     
      unweights = .false.                         
      duplicates = .false.       
      usevel = .true.
      usepos = .true.
      do is=1,maxsit
        col1(is) = ' '
      enddo

*  Get the run-string arguments and open the files

      iarg = iclarg(1,orgvelf)
      if( iarg.gt.0 ) then
        open(unit=iun,file=orgvelf,status='old',iostat=ioerr) 
        if( ioerr.ne.0 ) then
          write(*,'(a,a,i5)') 'Error opening input file ',orgvelf,ioerr
          stop
        endif  
        write(*,'(a,a)') ' Opened input file ',orgvelf 
      else             
        call proper_runstring('vel2stab.hlp','vel2stab',1)
      endif   
      open(unit=istab,file=stabsitef,status='unknown',iostat=ioerr) 
      if( ioerr.ne.0 ) then
        write(*,'(a,a,i5)') 'Error opening usesite filefile '
     .        ,stabsitef,ioerr   
        stop 
      else
         write(*,'(a,a)') ' Opened output file ',stabsitef 
      endif
      iarg = iclarg(2,arg)
      if( iarg.gt.0 ) then
        read(arg,*) hvmax
      endif
      iarg = iclarg(3,arg)
      if( iarg.gt.0 ) then
        read(arg,*) vvmax
      endif 
      iarg = iclarg(4,arg)
      if( iarg.gt.0 ) then
        read(arg,*) hpmax
      endif
      iarg = iclarg(5,arg)
      if( iarg.gt.0 ) then
        read(arg,*) vpmax
      endif
      write(*,'(a,2f10.1)') " Horizontal position,velocity cutoff "
     .                        ,hpmax,hvmax
      write(*,'(a,2f10.1)') " Vertical position,velocity cutoff   "
     .                       ,vpmax,vvmax

*     Loop decoding the rest of the runstring     
      nr = 5
      len_run = 1
      do while ( len_run.gt.0 )
         nr = nr + 1
         len_run = rcpar(nr, runstring )
         if( len_run.gt.0 ) then

*          See if we have new option
           word = runstring           
           call casefold(word)                             
cd           print *,'word ',word

*          See if we want to write an unweight file for hfupd        
           if( word(1:3).eq.'-UW ' ) then
             unweights = .true.
           elseif( word(1:6).eq.'-UWDUP' ) then
             duplicates = .true.

*          See if we want to use splat for stabilization sites
           elseif( word(1:8).eq.'-STABALL' ) then
             staball = .true.

           endif
         endif
      enddo  

      if( staball ) then
        write(*,'(a)') " Use splat with stab_site list"
      endif
      if( unweights ) then
        write(*,'(a)') " Writing an unweight file for hfupd "
      endif
      if( duplicates ) then    
        write(*,'(a)') " -- unweight duplicate 4-char IDs"
      else       
        write(*,'(a)') " -- retaining duplicate 4-char IDs"
      endif
          
* Check the inputs to see if velocities and/or positions to be included
                        
      if( hvmax.le.0. .or. vvmax.le.0. ) then 
        usevel = .false.              
      endif                                               
      if( hpmax.le.0 .or. vpmax.le.0 ) then
        usepos = .false.
      endif

* Write the headers for the stab_site file        

      write(istab,'(a,f5.1,1x,f5.1,a,f5.1,1x,f5.1,a)') 
     .     '* Sites selected by ORG2STAB  max vel sig (horiz vert) ='
     .    ,hvmax,vvmax,' mm/yr    max pos sig = ',hpmax,vpmax,' mm'  
      write(istab,'(2a)') '--1st column comment means omit for too-high'
     .                   ,' sigma: V velocity P position'
      write(istab,'(a)') '*'
      write(istab,'(a)') ' stab_site clear'
      write(istab,'(a,25x,2a)') '*',' Lon     Lat      Evsig    Nvsig  '
     .                  ,'  Uvsig    Epsig    Npsig   Upsig '


*  Read the velocity and optionally position values into storage,
*    --current code assumes a format with all header lines 
*      comments (non-blank column 1) and a (commented) position-header 
*      sequence beginning with '* Position' (sh_exglk output)
                                                   
      eof = .false.  
      endvel = .false.  
      is = 0  
cd      print *,'DEBUG on eof endvel ',eof,endvel
      do while (.not.eof ) 
        line = ' ' 
        read(iun,'(a)',iostat=ioerr) line 
        if( ioerr.eq.-1 .or. line(1:5).eq.'     ' ) then
          eof = .true.   
          numsit = is        
          if( endvel ) write(*,'(a,i5,a)') 'Read ',is,' positions'
        elseif( ioerr.ne.0 ) then  
          if(.not.endvel) write(*,*) 'Error reading velocity line',ioerr
          if( endvel )    write(*,*) 'Error reading position line',ioerr
          write(*,*) line
          stop
        elseif( line(1:5).eq.'* Pos') then   
          endvel = .true.     
          write(*,'(a,i5,a)') 'Read ',is,' velocities'
          numsit = is
          is = 0      
        elseif( line(1:1).ne.' ' ) then
          continue
        else  
          is = is + 1 
          if(is.gt.maxsit) then
            write(*,'(a,i6)') '# sites > maxsit ',maxsit
            stop    
          endif   
*         read a velocity or position (sh_exglk outputs the same number of columns but not the same precision or column width
*         -- first time through find the column length to get the site name
          if( is.eq.1 ) then
            if( line(nblen(line):nblen(line)).eq.'*') then
              nscol = nblen(line) - 8
            else
              nscol = nblen(line) - 7   
             endif  
          endif
          if( .not.endvel ) then   
            read(line(nscol:nscol+7),'(a8)',iostat=ioerr) site(is)
            if( ioerr.ne.0 ) then                      
              write(*,'(a)') 'Error reading site name for velocity'
              write(*,'(i3,1x,a)') 'nscol line ',nscol,line
              stop
            endif    
cd            print *,'nscol LINE',nscol,line
            read(line(1:nscol-1),*,iostat=ioerr) lon(is),lat(is)
     .         ,evel,nvel,devel,dnvel,evsig(is),nvsig(is)
     .         ,rho,hvel,dhvel,hvsig(is)
            if(ioerr.ne.0 ) then
              write(*,'(a,a8,a)') 'Error reading vel values for site '
     .             ,site(is),' set sigmas to 99'
              evsig(is) = 99.
              nvsig(is) = 99.
              hvsig(is) = 99.                       
            endif
          else
cd            print *,'reading sitep nscol  LINE ',nscol,line
            read(line(nscol:nscol+7),'(a8)',iostat=ioerr) sitep  
cd            print *,'read sitep ',sitep
            if( ioerr.ne.0 ) then                      
              write(*,'(a)') 'Error reading site name for position'
              write(*,'(i3,1x,a)') 'nscol line ',nscol,line
              stop
            endif  
            if( sitep.ne.site(is) ) then
              write(*,'(a,a4,1x,a4,a,i4)') 
     .           'Velocity, position site names (',site(is),sitep
     .           ,' ) do not match for site ',is
              stop
            endif    
cd            print *,'reading pos values nscol LINE ',nscol,line
            read(line(1:nscol-1),*,iostat=ioerr) lonp,latp,depos,dnpos
     ,           ,epsig(is),npsig(is),rho,dhpos,hpsig(is)
            if(ioerr.ne.0 ) then
              write(*,'(a,a8,a)') 'Error reading pos values for site '
     .             ,site(is),' set sigmas to 99'               
cd               print *,'ioerr lonp epos ',ioerr,lonp,epos 
              epsig(is) = 99.
              npsig(is) = 99.
              hpsig(is) = 99.                       
            endif
            if( abs(lon(is)-lonp).gt..0005 ) then 
              write(*,'(a,2f10.5,a,i4 )') 
     .             'Velocity, position longitudes ('
     .             ,lon(is),lonp,' ) do not match for site ',is
               stop
            endif    
*         endif for velocity or position
          endif
*       endif for valid read
        endif    
      enddo          

*  Check consistency of velocity and position counts

      if( usevel.and.usepos ) then
       if( is.ne.numsit ) then 
         write(*,'(a,2i6)') 'Something wrong # pos .ne. # vel',is,numsit
         stop   
       endif
      endif

            
*  Now loop through the entries, commenting with 'V' or 'P' if a site 
*  is excluded from stabilization for exceeding the maximum sigmas

      xvcount = 0
      xpcount = 0 
      do is = 1,numsit   
        if( usevel ) then
          if( evsig(is).gt.hvmax .or. nvsig(is).gt.hvmax .or.
     .        hvsig(is).gt.vvmax ) then
             col1(is)  = 'V'
             xvcount = xvcount + 1
          endif   
        endif
        if( usepos ) then  
          if( col1(is).ne.'v' ) then
             if( (epsig(is).gt.hpmax .or. npsig(is).gt.hpmax .or.
     .            hpsig(is).gt.vpmax ).and. col1(is).ne.'V' ) then
               col1(is)='P' 
               xpcount = xpcount + 1
             endif
          endif
        endif                             
      enddo
      


cd          print *
cd     .      ,is,site(is),lat(is),lon(is),evsig(is),nvsig(is),hvsig(is)
           
              
*  Now write out the stabilization list
                         
      stabtxt = 'stab_site ' 
      comtxt = ' ! '  
cd      print *,'writing out, numsit = ',numsit, stabtxt, comtxt
      do is=1,numsit
         siteout = site(is)
         if(staball ) siteout(5:8) = '@   '
        write(istab,'(a1,a10,a8,a3,8f9.3)') col1(is),stabtxt,siteout
     .     ,comtxt,lon(is),lat(is),evsig(is),nvsig(is),hvsig(is)
     .     ,epsig(is),npsig(is),hpsig(is)
      enddo         
         
*  Optionally create an unweight list for hfupd   

*     First go through the list and accumulate the sites to be unweighted
*     by the sigma criteria; then optionally add additional unweights 
*     for sites with dupllicate 4-character ids (the simple logic used
*     will generate duplicate unweight entries which are harmless but
*     could be removed with an editor
      if( unweights ) then   
        open(unit=iuwgt,file=unwgtf,status='unknown',iostat=ioerr)
        if( ioerr.ne.0 ) then
          write(*,'(a,a,i5)') 'Error opening unweight file '
     .        ,unwgtf,ioerr   
          stop 
        else
           write(*,'(a,a)') ' Opened output unweight file ',unwgtf
        endif    
        write(iuwgt,'(a,f5.1,1x,f5.1,a,f5.1,1x,f5.1,a)') 
     .     '* Sites selected by ORG2STAB  max vel sig (horiz vert) ='
     .    ,hvmax,vvmax,' mm/yr    max pos sig = ',hpmax,vpmax,' mm'  
        if( duplicates ) write(iuwgt,'(a)')
     .      '* --also unweight duplicate 4-character IDs'
        do is = 1,numsit
          uwflg(is) = .false.
          if( col1(is).ne.' ') uwflg(is) = .true.
        enddo              
        if( duplicates ) then
          do is = 2,numsit 
           if( site(is)(1:4).eq.site(is-1)(1:4).and.
     .        .not.uwflg(is-1) ) uwflg(is) =.true.
           if( is.gt.2.and.site(is)(1:4).eq.site(is-2)(1:4).and.
     .        .not.uwflg(is-2) ) uwflg(is) = .true. 
           if( is.gt.3.and.site(is)(1:4).eq.site(is-3)(1:4).and. 
     .        .not.uwflg(is-3) ) uwflg(is) = .true.
          enddo
        endif
        unwgttxt = ' unweight '  
        numusit = 0
        do is = 1, numsit
          if( uwflg(is) ) then
            write(iuwgt,'(a10,a8)') unwgttxt,site(is)
            numusit = numusit + 1
          endif
        enddo 
      endif

      write(*,'(a,i5)') 'Normal end of VEL2STAB numsit = ',numsit 
      write(*,'(a,i5)') 'Sites excluded for velocity sigma = ',xvcount
      write(*,'(a,i5)') 'Sites excluded for position sigma = ',xpcount
      goodcount = numsit - xvcount -xpcount
      write(*,'(a,i5)') 'Sites retained for stab_site list = ',goodcount
      if( unweights) write(*,'(a,i5)') 'Sites unweighted for hfupd = '
     .     ,numusit

      stop
      end
   
c-------------

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

