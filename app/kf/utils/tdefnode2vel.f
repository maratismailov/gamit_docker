c Convert a TDEFNODE vsum file to a GLOBK vel file, selecting a subset of stations
c if requested.

c vsum format 
c format(a8,3(1x,a4),   3i2,       4f10.4, 3(3f10.2,2f8.2), f8.3, 10f8.2, 2(4f8.2,f8.4), 3f8.2 )
CHUR_FPS PNW1 NoAm P001 1 1 0  265.9113   58.7591 1990.0000 2020.0000      -1.0       0.0      -1.0    0.30   -3.43      -1.3       0.0      -1.3    0.30   -4.27       0.0       0.0       0.0    0.00    0.00  -0.008     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0  0.0000     0.0     0.0     0.0     0.0  0.0000    0.00    0.00    0.00
c
c  1 Site name 
c  2 Data code 
c  3 Block code     
c  4 Pole code     
c  5 E flag (1=used, 0=not)
c  6 N flag (1=used, 0=not)
c  7 U flag (1=used, 0=not)
c  8 Longitude
c  9 latitude
c 10 Year 1
c 11 Year 2
c 12 E observed
c 13 E calculated
c 14 E residual
c 15 E sigma
c 16 E res/sigma
c 17 N observed
c 18 N calculated
c 19 N residual
c 20 N sigma
c 21 N res/sigma
c 22 U observed
c 23 U calculated
c 24 U residual
c 25 U sigma
c 26 U res/sigma
c 27 rho observed (NE)
c 28 E calculated elastic velocity
c 29 N calculated elastic velocity
c 30 U calculated elastic velocity
c 31 E calculated elastic velocity sigma (not yet calculated)
c 32 N calculated elastic velocity sigma (not yet calculated)
c 33 U calculated elastic velocity sigma (not yet calculated)
c 34 E calculated strain rate velocity
c 35 N calculated strain rate velocity
c 36 E calculated strain rate velocity sigma (not yet calculated)
c 37 N calculated strain rate velocity sigma (not yet calculated)
c 38 E calculated rotation velocity
c 39 N calculated rotation velocity
c 40 E calculated rotation velocity sigma
c 41 N calculated rotation velocity sigma
c 42 NE calculated rotation correlation
c 43 E calculated network velocity
c 44 N calculated network velocity
c 45 E calculated network velocity sigma
c 46 N calculated network velocity sigma
c 47 NE calculated network correlation 
c 48 E component of MR velocity
c 49 N component of MR velocity
c 50 U component of MR velocity

      implicit none

      integer eflag,nflag,uflag,ioerr,iudef,iuvel,iusites
     .      , maxsit,nssites,i
      parameter(maxsit=2000)

      real*4 lat,long,evel,nvel,hvel,eres,nres,eres1,nres1
     .     , vres,vsig,rat,hres,esig,nsig,hsig,dum,rho
                            
      character*1 col1,col77,ans
      character*4 soln,block,pole,frame,blocksel,polesel,solsel
      character*7 select
      character*8 site,selsites(maxsit)
      character*30 deffile,velfile,selsitefile
      character*256 line              

                
      logical endfile,matchsite,match,found
                          
c     initialize units and strings
      iudef = 1
      iuvel = 2 
      iusites = 3
      deffile = ' '
      velfile = ' '
      line = ' '   
      selsitefile = ' ' 
      select = ' ' 

c     set the height values 
      hvel = 0.
      hsig = 20.
      hres = 0.

c     get and open the input files
 
      write(*,*) 'Enter name of the input TDEFNODE vsum file:'
      read(*,*) deffile
      open(unit=1,status='old',file=deffile,iostat=ioerr) 
      if( ioerr.ne.0 ) then
        write(*,'(a,a30)')  'Error opening input TDEFNODE file ',deffile
        stop
      endif
      write(*,'(a,a30)') 'Input TDEFNODE vsum file: ',deffile     
      write(*,'(a,a4)')  "Choose 4-character solution (or 'all')"
      read(*,*,iostat=ioerr) solsel
      if(ioerr.ne.0) solsel = 'all '  
      if( solsel.eq.'all ') then
        write(*,'(a)') 'All solutions used'
      else                                      
        call uppers(solsel)
        write(*,'(a,a4)') 'Selected solution: ',solsel
      endif                                     
      write(*,'(a,a4)')  "Choose 4-character pole (or 'all')"
      read(*,*,iostat=ioerr) polesel
      if(ioerr.ne.0) polesel = 'all '  
      if( polesel.eq.'all ') then
        write(*,'(a)') 'All poles/blocks used'
      else        
        call uppers(polesel)
        write(*,'(a,a4)') 'Selected pole ',polesel
      endif                                     
      write(*,'(a,a4)')  "Read a file of sites for selection? (y/n)"
      read(*,'(a)') ans
      if( ans.eq.'y' ) then
        write(*,'(a)') "Enter file name:"
        read(*,*,iostat=ioerr) selsitefile
        open(unit=iusites,status='old',file=selsitefile,iostat=ioerr)   
        if( ioerr.ne.0 ) then    
          write(*,'(a,a30)')  'Error opening select-site file '
     .       ,selsitefile
          stop
        endif                        
        nssites = 1      
        endfile = .false.
        do while(.not.endfile)  
          selsites(nssites)(5:8) = '    '  
          read(iusites,'(a1,a8)',iostat=ioerr) col1,selsites(nssites)
          if( ioerr.eq.-1) then
            endfile = .true.      
          elseif( col1.eq. ' ' ) then   
            nssites = nssites + 1
            if(nssites.gt.maxsit) then
              write(*,'(a,i4,a)') '# sites > maxsit (',maxsit,')'
              stop
             endif
          endif
        enddo   
        nssites = nssites -1        
        write(*,'(i4,a)') nssites," sites read from the select-file"
        write(*,'(a)') 
     .   "Do you want to include (i) or exclude (e) these sites?"
        read(*,'(a)') ans
        if( ans.eq.'i') then
          select = "include"
        else        
          select = "exclude"
        endif
      endif
      write(*,*) 'Enter name of the output GLOBK vel file:'
      read(*,*) velfile
      write(*,*) 'Opening output file ',velfile
      open(unit=iuvel,status='unknown',file=velfile,iostat=ioerr)
      if( ioerr.ne.0 ) then
        write(*,'(a,a30,i5)')  
     .     'Error opening output vel file ',velfile,ioerr
        stop
      endif
c     write the header lines
      write(iuvel,'(a)') ' SUMMARY VELOCITiES & RESIDUALS FROM DEF2VEL'
      write(iuvel,'(3a)') '  Long.     Lat.        E & N Rate '
     .   ,'     E & N Adj.      E & N +-   RHO        H Rate  '
     .   ,' H adj.    +-  SITE'
      write(iuvel,'(2a)') '  (deg)    (deg)          (mm/yr)'
     .  ,'       (mm/yr)       (mm/yr)                 (mm/yr)'


c     read and write each line of the files

c     no header on the TDEFNODE vsum file      
      endfile = .false.
      do while (.not.endfile ) 
        read(iudef,'(a)',iostat=ioerr) line 
        if( ioerr.eq.-1 ) then
          endfile = .true.
          print *,'Normal stop at EOF '
          stop
        elseif( ioerr.ne.0 ) then  
          print *,'Error reading line  ioerr ',ioerr     
          print *,line
          stop   
        else          
          read(line,'(a8,3(1x,a4),3i2,2f10.4,20x
     .              ,3(f10.1,10x,f10.1,f8.2,8x),f8.3)')
     .       site,soln,block,pole,eflag,nflag,uflag,long,lat
     .     , evel,eres,esig,nvel,nres,nsig,hvel,hres,hsig,rho 
cd          print *,'site pole ',site,pole
          if( ioerr.ne.0 ) then
            print *,'Error decoding line ',ioerr
            print *,line  
            stop
          endif       
          call uppers(pole)
          call uppers(block)
          call uppers(soln)
cd         print *,evel,nvel,esig,nsig,rho
cd         print *,site,soln
cd         print *,eres1,nres1,vres,vsig,rat
c** Don't replace the extent character: need to be able to pick out unique site
c**        site(5:8) = '_DEF'  
cd           print *,'solsel soln',solsel,soln 
c** but fill in something if it's incomplete in the residual file
          if(site(6:6).eq.' ') site(5:8) = '_DEF'   
c         see if a subset is to be selected  
          col1 = ' ' 
c         --  by names from a file  
cd          print *,'DEBUG select ',select
          if( select.ne.' ') then   
            i=1           
            found = .false.
            do while(.not.found.and.i.le.nssites)  
cd              print *,'match ',i,selsites(i),site
              match = matchsite(selsites(i),site)
              if( match ) then
                found = .true.
                if( select.eq.'include' ) then
                   col1 = ' '
                elseif( select.eq.'exclude') then
                   col1 = 'x'
                endif                              
              endif
              i=i+1
            enddo     
            if( .not.found ) then
              if( select.eq.'include' ) col1 = 'x'
              if( select.eq.'exclude' ) col1 = ' ' 
            endif 
cd            print *,'list col1 ',col1
          endif
c         -- by solution  
cd          print *,'soln solsel col1 ',soln,solsel,col1
          if( (soln.eq.solsel.or.solsel(1:3).eq.'all') .and.
     .        col1.ne.'x' ) then     
            col1 = ' '
          else
            col1 = 'x'
          endif                       
cd          print *,'soln col1 ',col1
c         -- by pole   
cd          print *,'pole polesel col1 ',pole,polesel,col1
          if( (pole.eq.polesel.or.polesel(1:3).eq.'all') .and.
     .        col1.ne.'x' ) then
           col1 = ' '
          else
            col1 = 'x'
          endif                    
cd          print *,'pole col1 ',col1 
c         finally, 'x' out the site if horizontal not used by TDEFNODE
          if( eflag.eq.0 .or. nflag.eq.0 ) then
            col1 = 'x'
          endif                      
cd          print *,'flags col1 ',col1
          write(iuvel,'
     .      (a1,f8.3,1x,f8.3,1x,4f8.2,2f8.2,f7.3,2x,f8.2,2f8.2,1x,a8)'
     .           ,iostat=ioerr)
     .         col1,long,lat,evel,nvel,eres,nres,esig,nsig,rho
     .       , hvel,hres,hsig,site  
          if(ioerr.ne.0 ) stop 4    
        endif
      enddo

      stop
      end                       
ccc
      Logical function matchsite(site1,site2)
     
      character*8 site1,site2
      
      matchsite = .false.   
cd      print *,'MATCHSITE site1 site2 ',site1,site2,matchsite   

      if( site1.eq.site2 ) matchsite = .true. 
cd      print *,'all char test ',matchsite 
      if( site1(5:8).eq.'    '.or. site1(5:8).eq.'_DEF'.or.
     .    site1(5:5).eq.'@' .or.site1(5:6).eq.'_@') then
        if(site1(1:4).eq.site2(1:4)) matchsite = .true.
      endif
cd      print *,'site1 test ',matchsite
      if( site2(5:8).eq.'    '.or. site2(5:8).eq.'_DEF'.or.
     .    site2(5:5).eq.'@' .or.site2(5:6).eq.'_@') then
        if(site1(1:4).eq.site2(1:4)) matchsite = .true.
cd      print *,'site2 test ',matchsite
      endif
      return
      end

      



