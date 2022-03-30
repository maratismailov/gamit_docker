c     Program to extract co-seismic offsets from one or more files of offsets 
c     obtained from a prior seismic or geodetic model, and apply these to
c     an apr file. The program also generates a set of 'constrain'  commands 
c     for glorg to be assocated with each set of offsets.
c
c     The program assumes that sites with the same 4-character names are 
c     contiguous in the apr file, but not necessarily time-ordered. The
c     coordinates must all be the same (i.e generated using unify_apr).
       
c     The input EQ files (nomailly named, e.g. landers.enu) must have a header 
c     line that gives their 2-character code and the date (decimal year) of 
c     the earthquake. The outputs are a new apr file with displacements 
c     applied and a file 'disp2apr_glorg' of 'constrain' commands for glorg.

c     The revised version of November 2016 allows for multiple displacements
c     between apr entries; i.e. for an earthquake which is not explicitly
c     recognized by the site extent but which occurs between two apr entries.
c     This is accomplished by using the dates of the earthquakes (as well
c     as the site's presence in the an earthquake's displacemtn file)
c     to determine if a displacement should be applied
c     

c     The program is command-line driven with the form:

c         disp2apr [in-apr] [out-apr] [eq-file-list] [start] [constraint] 
c                  [min-horiz] [min-vert] [max-horiz] [max-vert] 
                                                    
c         The first four entries are required. 

c            in-apr :  name of input apr file
c            out-apr:  name of output apr file
c            eq-file-list:  list of earthquake files, single column beginning in column 1
c            start:   date in decimal years to be associated with an initial 'GPS' 
c                     entry, usually corresponding to the first date of data 
c            constraint: Constraint applied as a percentage of the displacement (default 100)
c            min-horiz:  Apply a displacement if it exceeds this value (mm, default 1)
c            min-vert :  Apply a displacement if it exceeds this value (mm, default 3)
c            max-horiz:  No horziontal constraint generated if the displacement exceeds this value (mm, defauult 200)
c            max_vert :  No vertical constraint generated if the displacement exceeds this value (mm, defauult 200)  

c     Format of the displacement files:
c
c        Header is [EQ] [Date in decimal years]  (1s,a2,f10.0)
c        Entries are displacement in mm E N U SITE LAT LONG free-format
c        Example file nisqually.enu:
cNI  2001.1640
c    0.5     0.1     0.0 149F 47.7973  239.5993 Y  189.886 1998.729 2013.715
c    0.0     0.2     0.0 A545 45.4745  239.2565 Y  241.174 1992.467 2011.553
c    0.1    -0.1     0.1 ACME 48.7107  237.7964 Y  178.062 1996.596 2001.630
c    1.0     3.8    -3.0 APSA 46.6708  237.0147 N   56.749 1993.649 2000.568
c  (fields 7-10 are comments, in this case distance from the earthquake and data range)

c     R. King  12 September 2007; revised 5 February 2009; last revised 14 November 2016

      implicit none

                  
c Maximum number of earthquakes allowed  - need also in apply_disp
      integer*4 maxeq
      parameter(maxeq=10)
    
c Maximum number of affiliated sites allowed (same 4-char name) - need also in apply_disp
      integer*4 maxns
      parameter(maxns=30)
     
      integer*4 ioerr,ns,i,ibug,iarg,iueqlist
                       
      character*32 eqf(maxeq),cmdf,aprinf,aproutf,eqlistf,tmpfile,arg
      character*256 lines(maxns),tmpline
   
      logical newsite,eof,finished

c  Function
      integer*4 nblen,trimlen,iclarg
      
c  Unit numbers in common                             
      integer*4 iuaprin,iuaprout,iucmd,iueq
      common/units/iuaprin,iuaprout,iucmd,iueq(maxns)

c  Constraint criteria in common
      real*4 start,minhdisp,minvdisp,maxhdisp,maxvdisp,con_level
      common/crit/start,minhdisp,minvdisp,maxhdisp,maxvdisp,con_level

c  Earthquake names and dates from file headers
      character*2 eqext(maxeq)
      integer*4 neq
      real*8 eqdate(maxeq)  
      common/earthquakes/eqdate,neq,eqext 

c Print help if no arguments

      iarg = iclarg(1,arg)
      if( iarg.eq.0 ) then
        print *,    
     .   'Example:  disp2apr [in-apr] [out-apr] [eq-filelist] '
     .   ,' [max horiz] [max-vert] [% constraint]'
        stop
      endif

c Initialize units, file names, and defaults
 
      do i=1,maxeq
        iueq(i) = i + 10
      enddo  
      iuaprin = 51
      iuaprout =52  
      iueqlist = 53
      iucmd = 54
      aprinf = 'apr.in'
      aproutf = 'apr.out'        
      eqlistf = 'eqlistf' 
      cmdf = 'disp2apr_glorg'    
      minhdisp = 1
      minvdisp = 3
      maxhdisp = 200
      maxvdisp = 200
      con_level = 100
                                                                     
c Write a message to the user
 
      write(*,*) 'Update an apr file with earthquake displacements'
      write(*,*) ' ' 
      write(*,*) '**NOTES : 1. All position values for a site must'
     .          ,' be the same (unify_apr). ' 
      write(*,*) '          2. No blank lines allowed.'
      write(*,*) '          3. Comments are not copied since they '
     .          ,'cannot be kept in the right place.'
      write(*,*) ' '                   


c Get the run-string arguments and open the files
      
      iarg = iclarg(1,aprinf)
      open(unit=iuaprin,file=aprinf,status='old',iostat=ioerr)
      if( ioerr.ne.0 ) then
        write(*,'(a,a,i5)') 'Error opening input apr file '
     .     ,aprinf,ioerr
        stop
      else  
        write(*,'(a,a)') ' Opened input apr file ',aprinf  
      endif                                              
      iarg = iclarg(2,aproutf)
      open(unit=iuaprout,file=aproutf,status='unknown',iostat=ioerr)
      if( ioerr.ne.0 ) then
        write(*,'(a,a,i5)') 'Error opening output apr file '
     .     ,aproutf,ioerr
        stop
      else  
        write(*,'(a,a)') ' Opened output apr file ',aproutf  
      endif                                              
      iarg = iclarg(3,eqlistf)
      open(unit=iueqlist,file=eqlistf,status='old',iostat=ioerr)
      if( ioerr.ne.0 ) then
        write(*,'(a,a,i5)') 'Error opening EQ list file '
     .     ,eqlistf,ioerr
        stop
      else  
        write(*,'(a,a)') ' Opened EQ list file ',eqlistf  
      endif                                                            
      open(unit=iucmd,file=cmdf,status='unknown',iostat=ioerr)
      if( ioerr.ne.0 ) then
        write(*,'(a,a,i5)') 'Error opening ouyput equate file '
     .     ,cmdf,ioerr
        stop
      else  
        write(*,'(a,a)') ' Opened output equate file ',cmdf  
      endif  

c Read the earthquake list file, open the displacement files, and read the headers
                    
      eof = .false.          
      neq = 0
      do while(.not.eof)    
        read(iueqlist,'(a)',iostat=ioerr) tmpfile 
        if( tmpfile(1:1).eq.' '.or. ioerr.ne.0 ) then
          eof =.true.
        else
          neq = neq + 1   
          if( neq.gt.maxeq ) then
            write(*,'(a,i3)') '# EQ files > maxeq ',maxeq
            stop 
          endif
          eqf(neq) = tmpfile    
        endif
      enddo 
      do i=1,neq 
        open(unit=iueq(i),file=eqf(i),status='old',iostat=ioerr)   
cd        print *,i,iueq(i),eqf(i)
        if( ioerr.ne.0 ) then
          write(*,'(a,a,i5)') 'Error opening input file ',eqf(i),ioerr   
          stop
        else  
          write(*,'(a,a)') ' Opened input file ',eqf(i)    
c         read the header
          read(iueq(i),'(a2,1x,f11.0)',iostat=ioerr) eqext(i),eqdate(i)
        endif
      enddo   
cd      print *,'eqf eqext  ',(eqf(i),eqext(i),i=1,neq)                       
c     put the earthquakes in chronological order
      call sorteq(maxeq,neq,eqdate,eqext,iueq)       
      write(*,*) 'Earthquakes with displacements provided: '
      do i=1,neq
        write(*,'(1x,a2,f10.4)') eqext(i),eqdate(i)
      enddo

c  Read the start date to be used if the first entry is _GPS

      iarg = iclarg(4,arg)
      if( iarg.le.0 ) then
        write(*,*) 'Start date in decimal years  now required'
      else
        read(arg,*) start
      endif

c  Read the percentage constraint and the displacement limits (optional)
                     
      iarg = iclarg(5,arg)
      if( iarg.gt.0 ) then
        read(arg,*) con_level
      endif
      iarg = iclarg(6,arg)   
      if( iarg.gt.0 ) then
        read(arg,*) minhdisp
      endif       
      iarg = iclarg(7,arg) 
      if( iarg.gt.0 ) then
        read(arg,*) minvdisp
      endif
      iarg = iclarg(8,arg)   
      if( iarg.gt.0 ) then
        read(arg,*) maxhdisp
      endif                
      iarg = iclarg(9,arg)     
      if( iarg.gt.0 ) then
        read(arg,*) maxvdisp
      endif       

c  Write the headers of the output apr and constraint files
                       
      write(iuaprout,'(a)') '*'
      write(iuaprout,'(a)') 
     .    '* Apriori coordinate file with EQ displacements added'         
      write(iuaprout,'(a,a16)') '* Original file : ',aprinf  
      write(iuaprout,'(a,10a16)') '* Earthquake files : ',
     .       (eqf(i),i=1,neq)         
      write(iuaprout,'(a,f8.3)') 
     .    '* Earliest date considered ',start
      write(iuaprout,'(a,f7.0,a,f7.0,a)') 
     .   '* Minimum displacements applied(mm) ',minhdisp
     .     ,' Hor  ',minvdisp,' Vert'
      write(iuaprout,'(a)') '* '                           
      write(iucmd,'(a)') '* Earthquake constraints for glorg'   
      write(iucmd,'(a,f4.0)') '* Level of constraint (%) ',con_level
      write(iucmd,'(a,f7.0,a,f7.0,a)') 
     .   '* No costraint if displacements (mm) greater than ',maxhdisp
     .     ,' Hor  ',maxvdisp,' Vert' 
       write(iucmd,'(a)') '*'  

c  Convert thresholds to internal units (m, fraction)
      maxhdisp = maxhdisp/1000.
      maxvdisp = maxvdisp/1000.
      minhdisp = minhdisp/1000.
      minvdisp = minvdisp/1000.
      con_level = con_level/100.

c  Read the apr file one group at a time, look for displacements, 
c  and write out the new lines (done in apply_disp even if no change)
                            
      ibug = 1              
      finished = .false.
      eof = .false.         
      newsite = .false.      
      do while (.not.finished )         
        ns = 1
        do while(.not.newsite.and. .not.finished )
          if( ns.gt.maxns ) then
            write(*,'(a,i3)') '# renamed sites > maxns ',maxns   
            stop
          endif         
          tmpline = ' '
          read(iuaprin,'(a)',iostat=ioerr) tmpline  
cd          print *,'ns ioerr  Read tmpline ',ns,ioerr,tmpline
          if( ioerr.eq.-1 ) then
            finished = .true.  
            newsite = .true.     
            ns = ns-1
cd            print *,'ioerr -1 ns newsite finished ',ns,newsite,finished
          elseif( ioerr.ne.0 ) then
            write(*,*) 'Error reading input apr file ',ioerr
            write(*,'(a)') tmpline
            stop   
          elseif( tmpline(1:1).ne.' '.or.tmpline(2:5).eq.'    ') then
cd            print *,'comment or null tmpline(1:5)= ',tmpline(1:5)  
c           Don't write comments since they will appear in the
c           wrong place due to the need to detect new sites'
c            write(iuaprout,'(a)') tmpline(1:trimlen(tmpline))  
             continue
          elseif( tmpline(2:9).eq.'EXTENDED') then
            write(iuaprout,'(a)') tmpline(1:trimlen(tmpline))  
          else  
cd            print *,'potential new site ns ',ns  
            if(ns.eq.1) then  
               lines(ns) = tmpline
cd              print *,'ns = 1 lines(1): ',lines(ns)(1:150) 
               ns = ns + 1
            elseif( tmpline(2:5).eq.lines(ns-1)(2:5)) then  
cd              print *,'sites match '
               lines(ns) = tmpline
               newsite = .false. 
               ns = ns + 1  
            else                     
cd              print *,'sites do not match newsite = T'  
              newsite = .true.                         
              ns = ns - 1
              backspace(iuaprin)
cd              print *,'backspaced iuaprin '
            endif
          endif
c---      read another line, potentially of the same group
        enddo 
c---    apply the diplacements (or not) and write out this group
          call apply_disp( ns,lines )         
c---      go read another group  
          newsite = .false.  
      enddo

c     EOF encountered, done.
      write(*,'(a)') 'Normal end'

      stop
      end

c-----------------------------------------------------------------------------

      Subroutine apply_disp(ns,lines)
                                      
c       Read the earthquake files and apply appropriate displacements
c       to the site coordinates for this group of sites

c         ns         # of sites in group  (dimension maxns)
c         line(1-ns) lines of the apr file for this group of sites

c       In common/earthquakes/
c          neq       # of earthquake files (dimenion maxeq)
c          eqext      2-character extents of the earthquakes available (maxeq)
c          eqdate    date of earthquakes in decimal years (maxeq)

c ** NB: Currently assumes that all coordinates are equal on the input apr file.

      implicit none
                 
c Maximum number of earthquakes allowed  - copied from main
      integer*4 maxeq
      parameter(maxeq=10)
    
c Maximum number of affiliated sites allowed (same 4-char name) - copied from main
      integer*4 maxns
      parameter(maxns=30)

      integer*4 ns,is,ieq,i,ioerr
  
      real*8 firstpos(3),readpos(3),sitpos(3),lastpos(3)
     .     , lastdate,newdate
     .     , sitvel(3),date,dispneu(3),dispxyz(3),disptotneu(3)
     .     , constr_neu(3),constr_inc_neu(3)
     .     , lon,lat,loc_coord(3),rot_matrix(3,3)
     .     , Econ,Ncon,Ucon
                 
      character*2 readeq,lasteq
      character*8 site,lastsite    
      character*256 lines(maxns),line_rem   

      logical found 

                               
                   
c  Unit numbers in common                             
      integer*4 iuaprin,iuaprout,iucmd,iueq
      common/units/iuaprin,iuaprout,iucmd,iueq(maxns)  

c  Constraint criteria in common
      real*4 start,minhdisp,minvdisp,maxhdisp,maxvdisp,con_level
      common/crit/start,minhdisp,minvdisp,maxhdisp,maxvdisp,con_level
                   
c  Earthquake names and dates from file headers
      character*2 eqext(maxeq)
      integer*4 neq
      real*8 eqdate(maxeq)  
      common/earthquakes/eqdate,neq,eqext 

                                               
c** DEBUG
cd        print *,'debug lines'
cd        do i=1,ns
cd          print *,lines(i)
cd        enddo

c  Order the site lines chronologically by earthquake
           
      call sortsites(maxeq,neq,eqext,maxns,ns,lines)
                                                       
c** DEBUG
cd         print *,'after sort'
cd         do i=1,ns
cd           print *,lines(i)
cd         enddo
                                                                       

c    Read and write the first site of the group
                                           
cd      do i=1,ns
cd        print *,'i lines',i,lines(i)(1:150)
cd      enddo
      call read_apr_line(lines(1),lastsite,firstpos,sitvel,date )
cc     .                  ,len_rem,line_rem,ioerr)   
      if( ioerr.eq.0 ) then
cc        write(iuaprout,10) lastsite,firstpos,sitvel,date
cc     .      ,line_rem(1:len_rem)
cc 10     format(1x,a8,3(1x,f14.4),1x,3(f9.5,1x),1x,f9.4,1x,a) 
        write(iuaprout,'(1x,a8,3(1x,f14.4),1x,3(f9.5,1x),1x,f9.4)')
     .         lastsite,firstpos,sitvel,date
cd        print *,'wrote out lastsite ',lastsite
        lasteq = lastsite(7:8)       
        call uppers(lasteq) 
        if(lasteq.eq.'PS') then
          lastdate = start
        else
          do i=1,neq
            if(eqext(i).eq.lasteq) lastdate = eqdate(i)
          enddo
        endif
cd         print *,'1st site, lasteq lastdate lastpos '
cd   .        ,lastsite,lasteq,lastdate,lastpos
      else
        write(*,'(2a)') 'Error reading apr line 1: ',lines(1)
        stop
      endif          
c ---   if only one site in the group, skip the rest
      if( ns.eq.1.or.ioerr.ne.0 ) return        
      do i=1,3
        lastpos(i) = firstpos(i)
      enddo

c    Read subsequent lines for the site and apply displacements 
                       
      do i=1,3
        sitpos(i) = lastpos(i)
      enddo
      do is=2,ns
        call read_apr_line(lines(is),site,readpos,sitvel,date )
cc  .                    ,len_rem,line_rem,ioerr )      
        if( ioerr.ne.0 ) then
          write(*,'(2a)') 'Error reading apr line ',is,' : ',lines(is)
          stop
        endif
        readeq = site(7:8)
        call uppers(readeq)
c       **temporary: make sure coordinates haven't changed
        if( readpos(1).ne.firstpos(1) ) then
          write(*,'(a,1x,a8,1x,a8,a)') '**Coordinates of common sites '
     .          ,site,lastsite,' do not match, skip'
          do i=1,3
            sitpos(i) = lastpos(i)
          enddo
        elseif( readeq.eq.lasteq ) then 
c         no new EQ: use coordinates from last set of displacements
cd        print *,'site ',is,site,' readeq lasteq ',readeq,lasteq
          do i=1,3
            sitpos(i) = lastpos(i)
          enddo     
        else                
c         get date of new EQ
          do i=1,neq
            if(eqext(i).eq.readeq )  newdate = eqdate(i)
          enddo 
c         get all displacements greater than the minimum since the last EQ
          do i=1,3
            dispneu(i) = 0.d0     
            disptotneu(i) = 0.d0
            constr_neu(i) = 0.d0
          enddo
          do ieq=1,neq                                               
            found = .false.
            if(eqdate(ieq).gt.lastdate.and.eqdate(ieq).le.newdate) then
              call get_disp( iueq(ieq),site(1:4),eqext(ieq)
     .                     , dispneu,found )
              if( found .and. 
     .           (dabs(dispneu(1)).gt.minhdisp.or.
     .            dabs(dispneu(2)).gt.minhdisp.or.
     .            dabs(dispneu(3)).gt.minvdisp) )  then
                call rotate_geod( dispneu,dispxyz,'NEU','XYZ',readpos
     .                        , loc_coord,rot_matrix)                      
                do i=1,3
                  sitpos(i) = sitpos(i) + dispxyz(i)  
                enddo                       
                write(iuaprout,'(a,a4,a,a2,a,3f9.4,a,3f9.4,a)') 
     .           '#  Displacement to ',site(1:4),' applied for '
     .            ,eqext(ieq)
     .            ,' : NEU ',(dispneu(i),i=1,3)
     .            ,'   XYZ ',(dispxyz(i),i=1,3),' m' 
c               accumulate the NEU displacements for the constraint comment
                do i=1,3
                  disptotneu(i) = disptotneu(i) + dispneu(i)
                enddo
c               get the constraint for this EQ 
                call get_constraint(dispneu,constr_inc_neu)
c               add the constraint for this EQ to any earlier ones between the renames
c               (should be a quadratic add, but don't bother)
                do i=1,3
                  constr_neu(i) = dsqrt( constr_neu(i)**2 
     .                                  + constr_inc_neu(i)**2 )
                enddo                                  
              endif
            endif  
          enddo   
c         write the line into the apr file
          write(iuaprout,'(1x,a8,3(1x,f14.4),1x,3(f9.5,1x),1x,f9.4)')
     .           site,sitpos,sitvel,date
c         write the constraints into the file for glorg 
          if( constr_neu(1).gt.0.d0 .and. constr_neu(1).lt.99.d0 )  
     .      write(iucmd,'(a,f7.4,1x,2(a8,a,a8,2a,f8.4))') 
     .          ' constrai ',constr_neu(1),site
     .            ,' npos ',lastsite,' npos', '  ! EQ dN ',disptotneu(1)     
          if( constr_neu(2).gt.0.d0 .and. constr_neu(2).lt.99.d0 ) 
     .      write(iucmd,'(a,f7.4,1x,2(a8,a,a8,2a,f8.4))') 
     .           ' constrai ',constr_neu(2),site 
     .            ,' epos ',lastsite,' epos', '  ! EQ dE ',disptotneu(2)
          if( constr_neu(3).gt.0.d0 .and. constr_neu(3).lt.99.d0 )     
     .       write(iucmd,'(a,f7.4,1x,2(a8,a,a8,2a,f8.4))') 
     .           ' constrai ',constr_neu(3),site
     .            ,' upos ',lastsite,' upos' ,'  ! EQ dU ',disptotneu(3)
c       endif on matching coordinates 
        endif
c       update eq name
        lasteq = readeq  
        lastdate = newdate 
c     enddo on reading apr lines 
      enddo
        
      return
      end
     
c-------------------------------------------------------------------------------

      Subroutine read_apr_line(line,site,sitpos,sitvel,date )
cd     ,                         ,len_rem,line_rem,ioerr )

c     Read a line of an apr file.
c        site, sitpos,sitvel,date are the required entires
c        line_rem of length len_rem is everything else
c     R. King from TAH kf/utils/unify_apr.f   070914

      character*8 site,cd
      character*256 line
cc   .   line_rem
                                 
      integer*4 indx,ioerr
      

      real*8 sitpos(3),sitvel(3),date,rd     

cd       print *,'READ_APR_LINE ',line 

      indx = 1       
      call read_line(line, indx, 'CH', ioerr, rd , site )   
      if( ioerr.ne.0.or.site(1:4).eq.'    ' ) then
        write(*,'(a,i4)') 'Error reading site from apr file ',ioerr
        return
      endif  
      call multiread(line, indx, 'R8', ioerr, sitpos, cd, 3 )    
      if( ioerr.ne.0 ) then
        write(*,'(a,a8,1x,i4)') 
     .     'Error reading position from apr file ',site,ioerr
        return
      endif                     
      call multiread(line, indx, 'R8', ioerr, sitvel, cd, 3 )  
      if( ioerr.ne.0 ) then
        write(*,'(a,a8,1x,i4)') 
     .     'Error reading velocity from apr file ',site,ioerr
        return
      endif                          
      call read_line(line, indx, 'R8', ioerr, date , cd )     
      if( ioerr.ne.0 ) then
        write(*,'(a,a8,1x,i4)') 
     .     'Error reading date from apr file ',site,ioerr
        return
      endif                          
cc      line_rem = line(indx:trimlen(line)) 
cc      len_rem = max(1,trimlen(line_rem))

      return
      end
   
c--------------------------------------------------------------------------

      Subroutine get_disp(iunit,site,eqext,disp,found)
                 

c     Read the earthquake displacement file and return values in meters

      implicit none    

      integer*4 iunit,ioerr,i
      character*2 eqext
      character*4 site,site1
      real*8 disp(3),N,E,U,lat,lon 
      logical found,eof         
      character*256 tmpline

      rewind(iunit)    
c     skip the header
      read(iunit,'(a)') 
      found = .false.  
      eof = .false.   
      do i=1,3
        disp(i) = 0.d0
      enddo
      do while(.not.found .and. .not.eof )       
c        read(iunit,'(f7.1,2f8.1,1x,a4,f8.4,f10.4)',iostat=ioerr)
c     .        E,N,U,site1,lat,lon  
        read(iunit,*,iostat=ioerr)  E,N,U,site1,lat,lon  
        disp(1)=N
        disp(2)=E
        disp(3)=U   
        call uppers(site1)
        if( ioerr.eq.-1) then
          write(*,'(a,a4,a,a2,a)') 'Site ',site,' not found on '
     .         ,eqext,' earthquake file, setting disp = 0.'  
          do i=1,3
            disp(i) = 0.
          enddo
          eof = .true.
        elseif( ioerr.ne.0 ) then
          write(*,'(a,1x,a2,1x,a,a4)') 'Error reading ',eqext
     .     ,' earthquake file for site ',site   
          backspace(iunit)
          read(iunit,'(a)') tmpline
          print *,'tmpline ',tmpline
          stop
        elseif(site1.eq.site) then
          found = .true.
        endif
      enddo                   
      do i=1,3
        disp(i) = disp(i)/1.d3
      enddo
      return
      end

c------------------------------------------------------------------------

      Subroutine get_constraint( dispneu,constr_neu )

c     Get a constraint value for selected displacements if criteria are met
c     If not, return 99. (m)  signalling no constraint 

      implicit none

      real*8 dispneu(3),constr_neu(3)

c  Constraint criteria in common
      real*4 start,minhdisp,minvdisp,maxhdisp,maxvdisp,con_level
      common/crit/start,minhdisp,minvdisp,maxhdisp,maxvdisp,con_level

c   Use Zheng-kang Shen's formula, applying half of the other component displacements
c        in addition to the main component (e.g. if input level is 60%,
c        then the constraint for E is sqrt(.6*dE**2 + .3*dN**2 + .3*dU**2)
       
      constr_neu(1) =  dsqrt( (      con_level*dispneu(1))**2 +
     .                        (0.5d0*con_level*dispneu(2))**2 +
     .                        (0.5d0*con_level*dispneu(3))**2 )
      constr_neu(2) = dsqrt( (      con_level*dispneu(2))**2 +
     .                        (0.5d0*con_level*dispneu(1))**2 +
     .                        (0.5d0*con_level*dispneu(3))**2 )
 
      constr_neu(3) =  dsqrt( (      con_level*dispneu(3))**2 +
     .                        (0.5d0*con_level*dispneu(1))**2 +
     .                        (0.5d0*con_level*dispneu(2))**2 )
    
c   Don't allow a constraint to be less than the displacement minimum 
      if( constr_neu(1).lt.minhdisp ) constr_neu(1) = minhdisp
      if( constr_neu(2).lt.minhdisp ) constr_neu(2) = minhdisp
      if( constr_neu(3).lt.minhdisp ) constr_neu(3) = minhdisp

c   Don't apply a constraint if the displacement exceeds the maximum                          
      if( dabs(dispneu(1)).gt.maxhdisp .or. 
     .   dabs(dispneu(2)).gt.maxhdisp ) then
        constr_neu(1) = 99.
        constr_neu(2) = 99.
      endif
      if( dabs(dispneu(3)).gt.maxvdisp ) constr_neu(3) = 99.

      return
      end

c-------------------------------------------------------------------------

   
      Subroutine sorteq( maxnam,neq,dates,names,units)  
                                         
c     Sort a real array with names into ascending order.
c     R. King, identical to sorteq in gamit/lib/lread.f
                 
       implicit none

       integer*4 maxnam,neq,units(maxnam),unitswap,i,k

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
            unitswap = units(i)
            dates(i)=dates(i-1) 
            names(i)=names(i-1)  
            units(i)=units(i-1)
            dates(i-1)=swap 
            names(i-1)=nameswap  
            units(i-1)=unitswap
 30      continue
         if(swap.eq.0.d0) go to 100
 50   continue
 100  continue
      return
      end   

c------------------------------------------------------------------------------------

      Subroutine sortsites(maxeq,neq,eqext,maxns,ns,lines)  
      
c     Sort the site names first by 6th character (normally ordered
c     G, 1, 2, 3...), then by same order as the earthquakes 
c     (chronological by a previous sort, though this doesn't matter here)
c     Assumes that the site names start in column 2

       implicit none

       integer*4 maxeq,neq,maxns,ns,nl,isitext(maxns),iswap,ioerr,i,j
                     
       character*1 sitext(maxns)                           
       character*2 eqext(maxeq),code
       character*256 lines(maxns),templines(maxns),lineswap
 
       logical found

c     Sort numerically but with G first (assign all alphamerics to 0 )
           
      do i=1,ns     
        call uppers( lines(i) )   
        read(lines(i)(7:7),'(i1)',iostat=ioerr) isitext(i)
        if( ioerr.ne.0 ) isitext(i) = 0
      enddo
      if(ns.le.1) go to 100
      do  j=1,32000
         iswap = -1
         do 30  i=2,ns
            if(isitext(i).ge.isitext(i-1)) go to  30
            iswap = isitext(i)
            lineswap = lines(i) 
            isitext(i)=isitext(i-1) 
            lines(i)=lines(i-1)  
            isitext(i-1)=iswap
            lines(i-1)=lineswap
 30      continue
         if(iswap.eq.-1) go to 100
      enddo  
 100  continue

c     Now sort by earthquake order

c     start the line counter
      nl = 0 
c     put any lines not corresponding to an EQ first
      do i=1,ns   
        code = lines(i)(8:9)  
        found = .false.     
        do j=1,neq
          if(code.eq.eqext(j)) found = .true.   
        enddo
        if( .not.found ) then  
           nl = nl + 1
           templines(nl) = lines(i)
        endif
      enddo

c     now add the EQ lines in the order of the EQs
      do i=1,neq 
        do j=1,ns
          if( lines(j)(8:9).eq.eqext(i) ) then
            nl = nl + 1  
            templines(nl) = lines(j)
          endif
        enddo
      enddo
  
c     check for missing lines
      if( nl.ne.ns ) then
         write(*,'(a,i2,a,i2)') 
     .     'Lines missed in resort; neq=',neq,' nl=',nl
         do i=1,neq
           write(*,'(a)') lines(i)
         enddo
      endif
         
      do i=1,nl
        lines(i) = templines(i)
      enddo

      return
      end


