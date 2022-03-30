c     Create Nutation, Sun, and Moon Tables for GAMIT processing
c     from a PEP N-body Ephemeris or a JPL Planetary Ephemeris
c     R. King December 2007, from old emntab (1984-2005)+ JPL routines.
c     R. King December 2009, change to read the ascii (export) version of the nbody file.

      implicit none                       

 
c        Maximum number of entries for Sun and Moon
      integer*4 maxsun,maxmoon
      parameter(maxmoon=1500)
      parameter(maxsun=150)
           
c        Input file name and source
      character*80 nbodyf
      character*3  source  

c        Unit numbers for input ephemeris and output tables
      integer*4 ibody,isun,imoon,inut,iout
      data ibody/10/inut/11/,isun/12/,imoon/13/,iout/8/

c        Span for output files
      integer*4 istart(3),istop(3),jdstart,jdstop

c        Flags for writing output files
      logical lsun,lmoon,lnut
      data lsun/.false./,lmoon/.false./,lnut/.false./
      integer*4 nvel
      data nvel/0/

c        Values for output files 
      character*80 titsun,titmon,titnut
      real*8 dints,dintm,earth(6,maxsun),moon(6,maxmoon)
      real*4 ts(maxsun),tm(maxmoon),tn(maxmoon)
     .     , dpsi(maxmoon),deps(maxmoon)
      integer*4 nsun,nmoon,nnut,npr
                               
c        Local
      character*1 ans
      integer*4 ioerr

c        Function
      integer*4 julday
      character*1 lowerc
          
       real*4 testvar            
       character*128 title
            
       testvar = .4728
c        Read the name of the input ephemeris file 
       
      write(6,'(a)')'Program to write GAMIT tables from PEP N-body file'
     .             ,' or JPL Planetary Ephemeris'
      write(*,'(a)') 'Enter the name of the input ephemeris file:'
      read(*,'(a)') nbodyf 
      write(*,'(a)') 'What is the type (JPL or PEP)?'
      read(*,'(a3)') source  
      call uppers(source) 
      if( source.ne.'PEP'.and.source.ne.'JPL') then
        write(*,'(a)') 'Type of file not PEP or JPL '
        stop
      endif
         
c       Requested span 
                                
      write(*,'(a)') 'Enter start, stop dates (YYYY MM DD  YYYY MM DD )'
      read(*,*,iostat=ioerr) istart,istop
c     convert to PEP JD
      jdstart = julday(istart(2),istart(3),istart(1))
      write(*,'(a,3i4,a,i7)') 'Start date ',istart,'  PEP JD ',jdstart
      jdstop = julday(istop(2),istop(3),istop(1))
      write(*,'(a,3i4,a,i7)') 'Start date ',istop,'  PEP JD ',jdstop


c       What output files do we want?

      write(*,'(a)') 'Do you want to write a Sun table (y/n)?'
      read(5,'(a)') ans
      if( lowerc(ans).eq.'y') then
        lsun = .true.
        open( unit=isun,file='soltab',form='formatted',status='unknown'
     .    , iostat=ioerr ) 
        if( ioerr.ne.0 ) then
          write(*,'(a,i4)') 'Error opening soltab ',ioerr
          stop
        endif                         
        write(*,'(a)') 'Enter title of Sun table'
        read(*,'(a80)') titsun
      endif 
      write(*,'(a)') 'Do you want to write a Moon table (y/n)?'
      read(*,'(a)') ans
      if( lowerc(ans).eq.'y') then   
        lmoon = .true.
        open( unit=imoon,file='luntab',form='formatted',status='unknown'
     .    , iostat=ioerr ) 
        if( ioerr.ne.0 ) then
          write(*,'(a,i4)') 'Error opening luntab ',ioerr
          stop
        endif
        write(*,'(a)') 'Enter title of Moon table'
        read(*,'(a80)') titmon
      endif 
      write(*,'(a)') 'Do you want to write a nutation table (y/n)?'
      read(*,'(a)') ans
      if( lowerc(ans).eq.'y') then   
        lnut = .true.
        open( unit=inut,file='nutabl',form='formatted',status='unknown'
     .      , iostat=ioerr ) 
        if( ioerr.ne.0 ) then
          write(*,'(a,i4)') 'Error opening nutabl ',ioerr
          stop
        endif
        write(*,'(a)') 'Enter title of Nutation table'
        read(*,'(a80)') titnut
      endif 
               

c       Set the tabular interval (days) and number of values (3 for position-only)
      
      dints = 4.d0
      dintm = 0.5d0
      npr = 3
         

c       Read the values into arrays  (units are km, arcsec)

      if( source.eq.'PEP' ) then 
c** rwk 091211: assume that we're now reading the ascii (export) version, but keep 
c               the  old subroutine, renamed to pepephb (binary)
c        call pepephb( ibody,jdstart,jdstop,lsun,lmoon,lnut,dints,dintm
c     .  , maxsun,maxmoon,nsun,nmoon,nnut,ts,tm,tn,earth,moon,dpsi,deps
c     .  , testvar ) 
        call pepeph( ibody,nbodyf,jdstart,jdstop,lsun,lmoon,lnut
     .             , dints,dintm,maxsun,maxmoon,nsun,nmoon,nnut
     .             , ts,tm,tn,earth,moon,dpsi,deps,testvar ) 
      elseif( source.eq.'JPL') then    
        call jpleph( nbodyf,ibody,jdstart,jdstop,lsun,lmoon,lnut
     .      , dints,dintm,maxsun,maxmoon,nsun,nmoon,nnut
     .      , ts,tm,tn,earth,moon,dpsi,deps )
      endif



