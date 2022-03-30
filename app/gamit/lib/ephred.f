Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.

      Subroutine EPHRED( lun,pjd0,lbody,ivel )    

c     Read a PEP or JPL n-body ephermeris file to populate arrays of tabular
c     point for the Earth, Moon, Venus, and Jupiter

c     R. King 9 January 2018 from older routines emntab, soltab, luntab
   
c     Input:

c     lun     i*4  logical unit for the ephemeris file (different in ARC and MODEL)
                           
c     pjd0    r*8  PEP JD in coordinate time

c     lbody   i*4  0  Earth and Moon only
c                  1  Also Venus and Jupiter

c     ivel    i*4  0 position only
c                  1 position and velocity

c     Output: contents of common /ephtab.
  
c      In order to use the Everett interpolator, we need 8 values centered 
c      about the requested epoch.  For the Earth, Venus, and Jupiter, each 
c      at 4-day intervals on the PEP ephemeeris, this implies reading in a 
c      40-day span, or two records of 20 days each.  Since we'd like to 
c      accommodate GNSS integrations of at least 3 days, we'll read in a 
c      60-day span (3 records), giving us giving us 15 points, 10 of which 
c      will be used for interpolation at any given  epoch.  For the Moon, the 
c      tabular interval is 0.5 days, so only 4 days would be  needed, but to 
c      simplify the logic in reading the PEP ephemeris, which stores 40 lunar 
c      points per record, we'll read in 120 points for the Moon at the same 
c      time. Thus the start time in the arrays will  be the same for all four 
c      bodies.  For the JPl ephemeris, we use JPL-created routines that read 
c      interpolate the values from a direct-acces file; thus our scheme will 
c      be to call these routines to create tabular values analogous to those 
c      from the PEP ephemeris. 
            
c      NB: The time arguments passed from here and used throughout GAMIT
c      are PEP JDs (= true JD +0.5), either as an integer plus a UTC or
c      TDT seconds or fraction of a day, or as a decimal - 24400000 in
c      or to retain precision.  When this last value is used, is will be
c      one day greater than the conventional Modified Julian date. 

       implicit none
         
c  Calling arguments
      integer*4 lun,lbody,ivel
      real*8 pjd0,rvec(6,4)

c  Ephemeris file name
      character*5 nbodyf 
                    
c  Program name and function used to get it
      character*80 prog_name
      integer*4 rcpar 

c  File type (PEP ascii, PEP binary, JPL binary)
      character*4 etype
           
c  Coordinate arrays and tabular intervals 
      common/ephtab/earth(6,15),moon(6,120),venus(6,15),jupiter(6,15)
     .  ,tsun(15),tmoon(120),jdstart,intsun,intmoon,intvenus,intjupiter
      integer*4 jdstart,intsun,intmoon,intvenus,intjupiter
      real*8 earth,moon,venus,jupiter,tsun,tmoon

c     earth, moon, venus, jupiter:  Tablular values over 60 days
c     tsun, tmoon                :  PEP JD of the epochs of the tabular points
c                                   (tsun covers all three planets)
c     jdstart                    :  PEP JD of first tabular point in the arrays
c     intsun, intmoon,intvenus, intjupiter: Tabular interval (4 for planets, -1 implying 0.5 for Moon)

c Other variables
      integer*4 jd0int,len,ioerr,i,j
c      character*4 byte
      character*16 cform,cassess
      character*80 message
      logical lexist 

      character*3 char3 

      logical debug/.false./
  
c  Get the calling program name for report_stat
      len = rcpar(0,prog_name) 
                                  
c  Determine the start and stop times for the values to be read 
c    PEP Earth/Venus/Jupiter  ephemerides are tabulated at 4-day intervals 
c    beginning at 24xxxx01 so we need to align with these

      jd0int = int(pjd0)
      jdstart = jd0int - mod((jd0int),20) +1 - 20 
      if(debug) print *,'EPHRED pjd0 jdstart ',pjd0,jdstart
                           
cc      jdstart = int(pjd0) - 5*intsun 
cc     this value will be reset in PEPRED to correspond to a tabular interval 
cc      if(debug) print *,'EPHRED jdstart ',jdstart

c  Determine the type of ephemeris, open it, and the read values into storage             
                      
c     file name currently hardwired
      nbodyf = 'nbody' 
      open(unit=lun,file=nbodyf,iostat=ioerr)
      if( ioerr.eq.0 ) then
        call report_stat('STATUS',prog_name,'lib/ephred'
     .   ,nbodyf,'Opened planetary ephemeris file ',0)
      else
        call report_stat('FATAL',prog_name,'lib/ephred'
     .   ,nbodyf,'Error opening planetary ephemeris file',ioerr)
      endif 
              
c  On Unix, there seems to be no clear delineation of sequential vs 
c  direct-access that can be recognized by the 'inquire' command, so
c  determine whether this is a JPL direct-access binary file or
c  a PEP sequential or binary  ascii file by reading the first three characters

      read(lun,'(a3)',iostat=ioerr) char3
      if( ioerr.ne.0 ) then
        write(message,'(a,a3)') 
     .    'Error reading 1st 3 characters of ephemeris file ',char3
        call report_stat('FATAL',prog_name,'lib/ephred',nbodyf,message
     .                  ,ioerr)
      endif
      close(lun)
      if( char3.eq.'JPL') then
        etype = 'JPL ' 
      elseif( char3.eq.' N-' ) then
        etype = 'PEPa'
      elseif( char3.eq.'200' ) then 
c       not yet sure this works 
        etype = 'PEPb'
      else 
        write(message,'(a,a3)') 'Cannot identify ephemeris type; char3='
     .     ,char3
        call report_stat('FATAL',prog_name,'lib/ephred','nbody',message
     .     ,0)
      endif
      if(debug) print *,'EPHRED etype ',etype                         

c Read the tabular values of the ephemeris into the arrays
         
      if( etype.eq.'PEPa'.or.etype.eq.'PEPb' ) then
        call pepred( lun,etype,pjd0,lbody )  
                       
      elseif( etype.eq.'JPL ' ) then 
        call jplred( lun,pjd0,lbody )
      endif

      if(debug) then
        write(*,'(a)') 'EPHRED Earth'
        do j=1,15
          write(*,'(f12.3,3f16.1,3f6.2)') tsun(j),(earth(i,j),i=1,6)
        enddo     
        write(*,'(a)') 'EPHRED Venus'
        do j=1,15
          write(*,'(f12.3,3f16.1,3f6.2)') tsun(j),(venus(i,j),i=1,6)
        enddo

        write(*,'(a)') 'EPHRED Jupiter'
        do j=1,15
          write(*,'(f12.3,3f16.1,3f6.2)') tsun(j),(jupiter(i,j),i=1,6)
        enddo
        write(*,'(a)') '       Moon'
        do j=1,120
          write(*,'(f12.3,3f16.1,3f6.2)') tmoon(j),(moon(i,j),i=1,6)
        enddo  
      endif
                  
      return
      end 

c***************************************************************************************888

      Subroutine PEPRED( lun,etype,pjd0,lbody )
                   
c     Read a binary or ascii (export) PEP N-body ephemeris to get 
c     15 tabular points for the Earth, Moon, Venus, and Jupiter
c     symmetric about PEP JD pjd0 

      implicit none

c       Unit number for n-body file
      integer*4 lun
                     
c       Type of file ('PEPa' or 'PEPb')         
      character*4 etype

c       Midpoint of requested  span
      real*8 pjd0 

c       Coodinates requested (lbody = 0 for Earth and Moon, =1 to add Venus and Jupiter)
      integer*4 lbody
                
c  Ephemeris file name
      character*5 nbodyf         

c  Program name and function used to get it
      character*80 prog_name
      integer*4 rcpar,len 

c  Coordinate arrays and tabular intervals 
      common/ephtab/earth(6,15),moon(6,120),venus(6,15),jupiter(6,15)
     .  ,tsun(15),tmoon(120),jdstart,intsun,intmoon,intvenus,intjupiter
      integer*4 jdstart,intsun,intmoon,intvenus,intjupiter
      real*8 earth,moon,venus,jupiter,tsun,tmoon

c  Rotation matrix from B1950 to J2000
      real*8 a2000(3,3)
      data a2000/.9999257079523629D0, .0111789381264276D0,
     .           .0048590038414544D0, -.0111789381377700D0,
     .           .9999375133499888D0, -.0000271579262585D0,
     .           -.0048590038153592D0, -.0000271625947142D0,
     .           .9999881946023742D0/    
                                        
c  Local timing quantities 
      real*8 dints,dintm

c        Header values      
      integer*4 jd1,jd2,jdbd0(10),intrd
     .        , nmoon1,itrt
     .        , jdbd1,jdbd2,nbdy2  
      integer*2 npl(10),ncp(10),intb(10),nbdy,nbdy1,intbd
     .        , kbdy(40),kkbdy(80)  
      character*28 name(12)
c     on the binary file this read as integers
      integer*4 iname(6,10)
      real*8 mass1(10),relftb(10),beta(6,10),meqinc
      real*4 epsbd
      character*128 title
      character*4 level,levp(10) 

c        Data record values    
      integer*4 ivl,mvl,jdbd
      real*8   frct,merc(6,10),body(6,5,8),mon(6,40)  
      real*4  psid(40),epsd(40),librt(40,3)

c        Local                  
      integer*4 jds1,jds2,idir,jm,nskip,kskip,ioerr,iout,line,npg
     .  ,jvlbd,nbdyp2,kk,i,j,k,nsun,nmoon,msun,mmoon,nrec,mrec
      data  intrd/0/,iout/8/   
      data      nskip/0/,kskip/9999/   
           real*8 mass(10),tmp3(3)
      real*8  convds,ltvel,tt,aultsc,fract,aukm,aukmvl
      data   convds/4.8481368110953599d-6/,ltvel/2.99792458d5/
     .     , aultsc/499.00478370d0/,fract/0.d0/,aukm/0.d0/,aukmvl/0.d0/    
      character*5 inframe  
      character*80 message 

      logical first,finished

      logical debug/.false./     
                
c  Get the calling program name for report_stat

      len = rcpar(0,prog_name) 

c  Open the requested file type
                                  
c     File name currently hardwired to 'nbody'
      nbodyf = 'nbody'
      if( etype.eq.'PEPa') then        
        open( unit=lun,file=nbodyf,form='formatted',status='old'   
     .      , iostat= ioerr)   
        if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/ephred'
     .         ,'nbody','Error opening PEP ascii file',ioerr)           
      elseif( etype.eq.'PEPb') then
        open( unit=lun,file=nbodyf,form='unformatted',status='old'   
     .      , iostat= ioerr)   
        if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/ephred'
     .         ,'nbody','Error opening PEP binary file',ioerr)    
      endif 

c  Read the header records from the N-body file 
                                   
      title = ' '
      if( etype.eq.'PEPa') then
        read(lun,'(a80)',iostat=ioerr) title(1:80)
        read(lun,'(a48,a4)',iostat=ioerr) title(81:128),level   
        if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/ephred'
     .         ,'nbody','Error reading title on PEP ascii file',ioerr)    
        if(debug ) then                                                       
          write(6,'(a,a128)') 'Title on ASCII  N-body file: ',title
          write(6,'(a,a4)') '  Version:',level
        endif
        read(lun,'(i2,3x,10i3)',iostat=ioerr) nbdy,(npl(i),i=1,nbdy)
        read(lun,'(10i3)',iostat=ioerr) (ncp(i),i=1,nbdy) 
        read(lun,'(10i3)',iostat=ioerr) (intb(i),i=1,nbdy)  
        read(lun,'(2i10/10i8)',iostat=ioerr)  
     .     jdbd1,jdbd2,(jdbd0(i),i=1,nbdy)
        read(lun,'(3d25.18)',iostat=ioerr) ((beta(i,j),i=1,6),j=1,nbdy)
        read(lun,'(3d25.18)',iostat=ioerr) meqinc  
        nbdyp2=nbdy+2        
cc RWK 18718: Correct bad read
cc        read(lun,'(a24,a4)',iostat=ioerr) (name(j),levp(j),j=1,nbdyp2)
        read(lun,'(a28)',iostat=ioerr) (name(j),j=1,nbdyp2)
        read(lun,'(4i5)',iostat=ioerr) nmoon,nbdy1,intbd,jvlbd
        read(lun,'((d25.18))',iostat=ioerr) (mass1(i),i=1,nbdy)
        read(lun,'(e15.8,2i6)',iostat=ioerr) epsbd,itrt,npg
        read(lun,'(2(16i5/),8i5)',iostat=ioerr) kbdy
        read(lun,'(4(16i5/),16i5)',iostat=ioerr) kkbdy        
        read(lun,'((3d25.18))',iostat=ioerr) (relftb(i),i=1,nbdy)
        if(ioerr.ne.0) call report_stat('FATAL',prog_name,'lib/ephred'
     .   ,'nbody','Error reading header lines on PEP ascii file',ioerr)
      elseif( etype.eq.'PEPb' ) then
        read(lun,iostat=ioerr) title
        if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/ephred'
     .         ,'nbody','Error reading title on PEP binary file',ioerr)    
        if(debug)  write(*,'(/,1x,a,/,1x,a128)') 
     .    'Title on binary N-body file:',title
       read(lun,iostat=ioerr) nbdy,(npl(i),i=1,nbdy),(ncp(i),i=1,nbdy)
     .         , (intb(i),i=1,nbdy),jdbd1,jdbd2,(jdbd0(i),i=1,nbdy)
     .         , ((beta(i,j),i=1,6),j=1,nbdy)
     .         , ((iname(i,j),i=1,6),j=1,nbdy)
     .         , nmoon,nbdy1,intbd,jvlbd,epsbd,(kbdy(i),i=1,40),itrt,npg
     .         , (mass1(i),i=1,nbdy),(relftb(i),i=1,nbdy)
        if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,'lib/ephred'
     .   ,'nbody','Error reading header lines on PEP binary file',ioerr)
      endif


c Compute masses and constants

      if( mass1(10).eq.0.d0) mass1(10) = 82.3005883d0
      do i=1,nbdy
        mass(i) = 0.d0
        if( mass1(i).gt.1.d0 ) mass(i) = 1.d0/mass1(i)
      enddo  
      aukm = ltvel*aultsc
      aukmvl = aukm/86400.d0


c If debug print the N-body header
                  
      if( debug ) then 
        write(*,'(/,a,6x,2i10,i5)') 'JDBD1,JDBD2,NBDY=',jdbd1,jdbd2,nbdy
        idir = 1
        if(jdbd1.gt.jdbd2) idir = -1
        write(*,'(/,2a)')'   NPL    NPC       NAME          '
     .     ,'PEP Ver      MASS1'
        print *,'nbdy npl ',nbdy,(npl(k),k=1,nbdy)
        write(*,'(/,(2x,2i5,2x,a24,1x,a4,4x,1pd22.15))')  
     .      ( npl(k),ncp(k),name(k),levp(k),mass1(k),k=1,nbdy )
        write(*,'(/,a,10i5,/,3(6x,10i5,/))') 'kbdy= ',kbdy
        write(*,'(/,a,10i5,/,7(6x,10i5,/))') 'kkbdy=',kkbdy       
        write(*,'(/,a,3d25.18,/,3(8x,3d25.18,/))')  
     .          'relftb= ',(relftb(i),i=1,nbdy)    
        if( etype.eq.'PEPb') then
          write(*,'(a,3x,6i3,1x,1pd10.3)')
     .         'nmoon1,nbdy1,intbd,jvlbd,intrd,npg,epsbd'
     .      ,   nmoon1,nbdy1,intbd,jvlbd,intrd,npg,epsbd
        endif 
      endif    
        
c  Trap a file backwards-in-time or a reference frame not B1950 or J2000
    
      if( jdbd1.gt.jdbd2) call report_stat('FATAL',prog_name
     .  ,  'lib/ephred','nbody'
     .  ,  'Logic fails for a file backwards-in-time',0)
      if( kkbdy(70).eq.0 ) then 
        inframe = 'B1950' 
      elseif( kkbdy(70).eq.1 ) then
        inframe = 'j2000'
      else 
         call report_stat('FATAL',prog_name,'lib/ephred','nbody'
     .         ,'kkbdy(70) indicates an invalid frame',0)
      endif        
                                      
c       Set the start and stop times: assume that the interval for 
c       Earth, Venus, and Jupiter is the same (nominally 4 days)
c       and that the start for the Moon will be the same as for
c       the other bodies 

      intsun = intb(3)
      intvenus = intb(2)
      intjupiter = intb(5)
      intmoon = intb(10)
c      jdstart = int(pjd0) - 7*intsun
      dints = intb(3)
      dintm = 2.d0**intb(10)
c      if(debug) write(*,'(a,i9 )') 'Requested startjd ',jdstart


c  Read through the file and fill the arrays
c       nsun/nmoon count the values in the arrays to be returned
c       msun/mmoon count the values within each nbody record 
c       msun/mmoon counts the number of records read for the requested values     
c       nrec counts the total records on the nbody file (not used)

      first = .true.            
      finished = .false.                                 
      nrec = 0 
      nbdy2 = nbdy1-1
      nsun = 0 
      nmoon= 0   
      do while(.not.finished)          
        nrec = nrec + 1
        if( etype.eq.'PEPa' ) then 
          read(lun,'(i10,d25.18,2i5)',iostat=ioerr,end=999) 
     .       jdbd,frct,ivl,mvl  
          do i=1,10
            read(lun,'((3d25.18/3d25.18))',iostat=ioerr) 
     .        (merc(j,i),j=1,ivl)
          enddo  
          if( debug.and.nrec.le.2) then      
            write(*,'(/,a,i8,f6.3,i3)') 'jdbd fract ivl Mercury'
     .        ,jdbd,frct,ivl
            do i=1,10
              write(*,'(6d10.3)') (merc(j,i),j=1,6) 
            enddo 
          endif
          do i=1,nbdy2
           read(lun,'((3d25.18/3d25.18))',iostat=ioerr) 
     ,          ((body(k,j,i),k=1,ivl),j=1,5) 
          enddo   
          read(lun,'((3d25.18/3d25.18))',iostat=ioerr) 
     .      ((mon(i,j),i=1,mvl),j=1,40)
          read(lun,'((2e15.8))',iostat=ioerr) (psid(i),epsd(i),i=1,40)
          read(lun,'((3e15.8))',iostat=ioerr) 
     .       ((librt(j,i),i=1,3),j=1,40) 
          if( ioerr.ne.0 ) then
            write(message,'(a,i4,a)') 'Read error at data record ',nrec
     .      ,' on PEP ASCII file '
            call report_stat('FATAL',prog_name,'lib/ephred','nbody'
     .                     ,message,ioerr)
          endif
        elseif( etype.eq.'PEPb' ) then
          read(lun,iostat=ioerr,end=999)jdbd,frct,ivl,mvl
     .        , ((merc(i,j),i=1,ivl),j=1,10)
     .        , (((body(i,j,k),i=1,ivl),j=1,5),k=1,nbdy2)
     .        , ((mon(i,j),i=1,mvl),j=1,40)
     .        , (psid(j),epsd(j),j=1,40), ((librt(j,i),i=1,3),j=1,40)
         if( ioerr.ne.0 ) then
            write(message,'(a,i4,a)') 'Read error at data record ',nrec
     .             ,' on PEP binary file '
            call report_stat('FATAL',prog_name,'lib/ephred','nbody'
     .                     ,message,ioerr)
          endif
        endif
        if( jdbd.ge.jdstart ) then  
          mrec = 0 
          msun = 0 
          mmoon = 0         
          do j=1,5
            jm= 8*j - 7
            nsun = nsun + 1
            msun = msun + 1 
            if(debug) print *,'nsun msun ',nsun,msun
            if(nsun.gt.15) call report_stat('FATAL',prog_name
     .          ,'lib/ephred','nbody','# Earth values > 15 ',0)
            tsun(nsun)= dfloat(jdbd-2400000) + (dfloat(msun-1))*dints
            if( debug ) write(*,'(a,i3,f13.5,f6.1,6f16.8)') 
     .              'EMbary nsun tsun dints coords '
     .           ,nsun,tsun(nsun),dints,(body(i,j,2),i=1,6)
            do i=1,3
              earth(i,nsun) = (body(i,j,2)-mass(10)*mon(i,jm))*aukm
              earth(i+3,nsun) =
     .                  (body(i+3,j,2)-mass(10)*mon(i+3,jm))*aukmvl
            enddo   
            if( inframe.eq.'B1950' ) then 
              call matmpy( a2000,earth(1,nsun),tmp3,3,3,1)
              do i=1,3
                earth(i,nsun) = tmp3(i)
               enddo
            endif
            if(debug) 
     .          write(*,'(a,i3,f13.5,f6.1,6f16.3)') 
     .              'Earth nsun tsun dints coords  '
     .           ,nsun,tsun(nsun),dints,(earth(i,nsun),i=1,6)
            if( lbody.gt.0 ) then 
              do i=1,3
                venus(i,nsun) = body(i,j,1)*aukm
                venus(i+3,nsun) = body(i+3,j,1)*aukmvl
              enddo   
              if( inframe.eq.'B1950' ) then 
                call matmpy( a2000,venus(1,nsun),tmp3,3,3,1)
                do i=1,3
                  venus(i,nsun) = tmp3(i)
                 enddo
              endif
              if(debug) 
     .          write(*,'(a,i3,f13.5,f6.1,6f16.3)') 
     .              'Venus nsun tsun dints coords  '
     .           ,nsun,tsun(nsun),dints,(venus(i,nsun),i=1,6)

              do i=1,3
                jupiter(i,nsun) = body(i,j,4)*aukm
                jupiter(i+3,nsun) = body(i+3,j,4)*aukmvl
              enddo   
              if( inframe.eq.'B1950' ) then 
                call matmpy( a2000,jupiter(1,nsun),tmp3,3,3,1)
                do i=1,3
                  jupiter(i,nsun) = tmp3(i)
                 enddo
              endif
              if(debug) 
     .          write(*,'(a,i3,f13.5,f6.1,6f16.3)') 
     .              'Jupiter nsun tsun dints coords  '
     .           ,nsun,tsun(nsun),dints,(jupiter(i,nsun),i=1,6)
            endif  
          enddo  
          mmoon = 0 
          do j=1,40
            nmoon = nmoon + 1        
            mmoon = mmoon + 1
            if(nmoon.gt.120) call report_stat('FATAL',prog_name
     .          ,'lib/ephred','nbody','# Moon values > 120 ',0)
            tmoon(nmoon) = dfloat(jdbd-2400000) + dfloat(mmoon-1)*dintm 
            if(debug)  write(*,'(a,i3,f13.5,f6.1,6f16.10)') 
     .              'Moon nmoon tmoon dintm coords '
     .           ,nmoon,tmoon(nmoon),dintm,(mon(i,j),i=1,6)
            do i=1,3
              moon(i,nmoon) = mon(i,j)*aukm
              moon(i+3,nmoon)= mon(i+3,j)*aukmvl
            enddo      
            if( inframe.eq.'B1950' ) then 
              call matmpy( a2000,moon(1,nmoon),tmp3,3,3,1)
              do i=1,3
                moon(i,nmoon) = tmp3(i)
               enddo
            endif
            if(debug) write(*,'(a,i3,f13.5,f6.1,6f16.2)') 
     .              'Moon nmoon tmoon dintm coords '
     .           ,nmoon,tmoon(nmoon),dintm,(moon(i,nmoon),i=1,6)
          enddo
        endif
        if( nsun.ge.15..or.nmoon.ge.120 ) finished = .true. 
      enddo

 999  return
      end

c**********************************************************************************

      Subroutine JPLRED(lun,pjd0,lbody)
 
c     Read a binary direct-access JPL n-body ephemeris to get 15 tabular points 
c     at 4-day intervals for the Earth, Venus, and Jupiter, 120 at half-day 
c     intervals for the Moon, all with starting epoch on the nearest 20-day
c     tabular point at least 20 days prior to the epoch of the data or integration
c     (pjd0). This scheme chosen to match the one used for the PEP ephermeris
c     and therefore able to share the interpolation code.  See the documentation
c     at the top of PEPRED. 

c     R. King 30 January 2018 from the routines of gg/utils/emntab.

      implicit none

      integer*4 lun,lbody
  
      real*8 pjd0
     
c   Input:
                 
c     lun     i*4  logical unit for the ephemeris file (different in ARC and MODEL)
c                  file name hard-wired to 'nbody'
                           
c     pjd0    r*8  PEP JD in coordinate time

c     lbody   i*4  0  Earth and Moon only
c                  1  Also Venus and Jupiter

c  Output:  contents of common /ephtab/ - see documentation at top of ephred.

c  Coordinate arrays and tabular intervals 
      common/ephtab/earth(6,15),moon(6,120),venus(6,15),jupiter(6,15)
     .  ,tsun(15),tmoon(120),jdstart,intsun,intmoon,intvenus,intjupiter
      integer*4 jdstart,intsun,intmoon,intvenus,intjupiter
      real*8 earth,moon,venus,jupiter,tsun,tmoon

c  Values from the JPL ephemeris
      integer*4 ntarg,ncent
      real*8 et,rrd(6)

c  Local                      
      real*8 xjd,dintm
      integer*4 jdstop,jd,nsun,nmoon,rcpar,len,i  
      character*80 nbodyf
      character*16 prog_name
                           
      logical debug/.false./
            
                
c Get the calling program name for report_stat

      len = rcpar(0,prog_name) 


c Set the file name for calls to JPL routines

      nbodyf = 'nbody'
                    
c Get the start time and tabular intervals to match a PEP ephemeris
            
      intsun = 4
      intvenus = 4
      intjupiter = 4
      intmoon = -1 
      jdstart = int(pjd0) - 7*intsun 
c      if(debug) write(*,'(a,f16.2,i8)') ' pjd0, jdstart ',pjd0,jdstart
c      jdstart = jdstart - mod(jdstart,20) + 1     
       jdstop  = jdstart + 59
c      if(debug) write(*,'(a,i8)') 'modified jdstart ',jdstart

c Get the tabular values
                              
      jd = jdstart    
      nsun = 0
      do while( jd.le.jdstop)         
        et = dfloat(jd) - 0.5d0      
c       Earth 
        ntarg = 3 
        ncent = 11    
        nsun = nsun + 1
        if(nsun.gt.15) call report_stat('FATAL',prog_name
     .         ,'lib/ephred','nbody','# Earth values > 15 ',0)
        tsun(nsun) = float(jd) - 2400000.
        call pleph( nbodyf,lun,et,ntarg,ncent,rrd ) 
        do i=1,6   
          earth(i,nsun) = rrd(i)
        enddo
        if( lbody.gt.0 ) then
c         Venus
          ntarg = 2 
          ncent = 11    
          call pleph( nbodyf,lun,et,ntarg,ncent,rrd ) 
          do i=1,6   
            venus(i,nsun) = rrd(i)
          enddo
c         Jupiter
          ntarg = 5 
          ncent = 11    
          call pleph( nbodyf,lun,et,ntarg,ncent,rrd ) 
          do i=1,6   
            jupiter(i,nsun) = rrd(i)
          enddo
        endif             
        if(debug) then
          write(*,'(a,i8,f16.2,i4,f10.2)') 'jd et nsun tsun '
     .        ,jd,et,nsun,tsun(nsun)
          write(*,'(a,6f16.2)') 'Earth  ',(earth(i,nsun),i=1,6)
          write(*,'(a,6f16.2)') 'Venus  ',(venus(i,nsun),i=1,6)
          write(*,'(a,6f16.2)') 'Jupiter',(jupiter(i,nsun),i=1,6)
        endif
        jd = jd + intsun
      enddo            
c     Moon
      ntarg = 10 
      ncent = 3    
      nmoon = 0  
      jd = jdstart 
      xjd = dfloat(jd)   
      dintm = 2.0**intmoon
      if(debug) print *,'Calling PLEPH ntarg ncent jd jdstop xjd dintm '
     .                               , ntarg,ncent,jd,jdstop,xjd,dintm
      do while( jd.le.jdstop) 
        nmoon = nmoon + 1 
        et = xjd - 0.5d0      
        if(nmoon.gt.120) call report_stat('FATAL',prog_name
     .         ,'lib/ephred','nbody','# Moon values > 120 ',0)
        tmoon(nmoon) = xjd - 2400000.d0
        call pleph( nbodyf,lun,et,ntarg,ncent,rrd ) 
        if(debug ) print *,'et nmoon tmoon rrd '
     .                    , et,nmoon,tmoon(nmoon),rrd
        do i=1,6   
          moon(i,nmoon) = rrd(i)
        enddo                     
        if(debug) then
          write(*,'(a,i8,f16.2,i4,f10.2)') 'jd et nmoon,tmoon '
     .      ,jd,et,nmoon,tmoon(nmoon)
          write(*,'(a,6f16.2)') 'Moon  ',(moon(i,nmoon),i=1,6)
        endif
        xjd = xjd + dintm
        jd = idint(xjd)
      enddo              
                       
      return
      end
        
        
C++++++++++++++++++++++++++
C
      SUBROUTINE PLEPH ( namfil,nrfile,ET, NTARG, NCENT, RRD )
C
C++++++++++++++++++++++++++
C  NOTE : Over the years, different versions of PLEPH have had a fifth argument:
C  sometimes, an error return statement number; sometimes, a logical denoting
C  whether or not the requested date is covered by the ephemeris.  We apologize
C  for this inconsistency; in this present version, we use only the four necessary 
C  arguments and do the testing outside of the subroutine.
C
C
C
C     THIS SUBROUTINE READS THE JPL PLANETARY EPHEMERIS
C     AND GIVES THE POSITION AND VELOCITY OF THE POINT 'NTARG'
C     WITH RESPECT TO 'NCENT'.
C
C     CALLING SEQUENCE PARAMETERS:
C
C       ET = D.P. JULIAN EPHEMERIS DATE AT WHICH INTERPOLATION
C            IS WANTED.
C
C       ** NOTE THE ENTRY DPLEPH FOR A DOUBLY-DIMENSIONED TIME **
C          THE REASON FOR THIS OPTION IS DISCUSSED IN THE 
C          SUBROUTINE STATE
C
C     NTARG = INTEGER NUMBER OF 'TARGET' POINT.
C
C     NCENT = INTEGER NUMBER OF CENTER POINT.
C
C            THE NUMBERING CONVENTION FOR 'NTARG' AND 'NCENT' IS:
C
C                1 = MERCURY           8 = NEPTUNE
C                2 = VENUS             9 = PLUTO
C                3 = EARTH            10 = MOON
C                4 = MARS             11 = SUN
C                5 = JUPITER          12 = SOLAR-SYSTEM BARYCENTER
C                6 = SATURN           13 = EARTH-MOON BARYCENTER
C                7 = URANUS           14 = NUTATIONS (LONGITUDE AND OBLIQ)
C                            15 = LIBRATIONS, IF ON EPH FILE
C
C             (IF NUTATIONS ARE WANTED, SET NTARG = 14. FOR LIBRATIONS,
C              SET NTARG = 15. SET NCENT=0.)
C
C      RRD = OUTPUT 6-WORD D.P. ARRAY CONTAINING POSITION AND VELOCITY
C            OF POINT 'NTARG' RELATIVE TO 'NCENT'. THE UNITS ARE AU AND
C            AU/DAY. FOR LIBRATIONS THE UNITS ARE RADIANS AND RADIANS
C            PER DAY. IN THE CASE OF NUTATIONS THE FIRST FOUR WORDS OF
C            RRD WILL BE SET TO NUTATIONS AND RATES, HAVING UNITS OF
C            RADIANS AND RADIANS/DAY.
C
C            The option is available to have the units in km and km/sec.
C            For this, set km=.true. in the STCOMX common block.
C
                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      integer*4 ntarg,ncent,ncon,i,k


      DIMENSION RRD(6),ET2Z(2),ET2(2),PV(6,13)
      DIMENSION SS(3),CVAL(400),PVSUN(6),zips(2)
      data zips/2*0.d0/

      LOGICAL BSAVE,KM,BARY
      LOGICAL FIRST
      DATA FIRST/.TRUE./

      INTEGER LIST(12),IPT(39),DENUM
                                       
      character*80 namfil
      integer*4 nrfile

      COMMON/EPHHDR/CVAL,SS,AU,EMRAT,DENUM,NCON,IPT

      COMMON/STCOMX/KM,BARY,PVSUN  

      dimension pnut(4)

c** rwk 071231: Set logical 'km' to true to return values in km, not au
      km = .true.

C     INITIALIZE ET2 FOR 'STATE' AND SET UP COMPONENT COUNT
C
      ET2(1)=ET
      ET2(2)=0.D0
      GO TO 11

C     ENTRY POINT 'DPLEPH' FOR DOUBLY-DIMENSIONED TIME ARGUMENT 
C          (SEE THE DISCUSSION IN THE SUBROUTINE STATE)

      ENTRY DPLEPH(ET2Z,NTARG,NCENT,RRD)

      ET2(1)=ET2Z(1)
      ET2(2)=ET2Z(2)

  11  IF(FIRST) CALL STATE(namfil,nrfile,zips,list,pv,pnut)
      FIRST=.FALSE.

  96  IF(NTARG .EQ. NCENT) RETURN

      DO I=1,12
      LIST(I)=0
      ENDDO

C     CHECK FOR NUTATION CALL

      IF(NTARG.NE.14) GO TO 97
        IF(IPT(35).GT.0) THEN
          LIST(11)=2
          CALL STATE(namfil,nrfile,ET2,LIST,PV,pnut)
c          print *,'after nutation STATE ET2 LIST ',et2,list
c          print *,'PV '
c          print *,' pnut ',pnut
c          stop    
           do i=1,4
            rrd(i) = pnut(i)
           enddo
          RETURN
        ELSE
          do i=1,4
          rrd(i)=0.d0
          enddo
          WRITE(6,297)
  297     FORMAT(' *****  NO NUTATIONS ON THE EPHEMERIS FILE  *****')
          STOP
        ENDIF

C     CHECK FOR LIBRATIONS

  97  do i=1,6
      rrd(i)=0.d0
      enddo

      IF(NTARG.NE.15) GO TO 98
        IF(IPT(38).GT.0) THEN
          LIST(12)=2
          CALL STATE(namfil,nrfile,ET2,LIST,PV,pnut)
          DO I=1,6
          RRD(I)=PV(I,11)
          ENDDO
          RETURN
        ELSE
          WRITE(6,298)
  298     FORMAT(' *****  NO LIBRATIONS ON THE EPHEMERIS FILE  *****')
          STOP
        ENDIF

C       FORCE BARYCENTRIC OUTPUT BY 'STATE'

  98  BSAVE=BARY
      BARY=.TRUE.

C       SET UP PROPER ENTRIES IN 'LIST' ARRAY FOR STATE CALL

      DO I=1,2
      K=NTARG
      IF(I .EQ. 2) K=NCENT
      IF(K .LE. 10) LIST(K)=2
      IF(K .EQ. 10) LIST(3)=2
      IF(K .EQ. 3) LIST(10)=2
      IF(K .EQ. 13) LIST(3)=2
      ENDDO
                
C       MAKE CALL TO STATE

      CALL STATE(namfil,nrfile,ET2,LIST,PV,pnut)

      IF(NTARG .EQ. 11 .OR. NCENT .EQ. 11) THEN
      DO I=1,6
      PV(I,11)=PVSUN(I)
      ENDDO
      ENDIF

      IF(NTARG .EQ. 12 .OR. NCENT .EQ. 12) THEN
      DO I=1,6
      PV(I,12)=0.D0
      ENDDO
      ENDIF

      IF(NTARG .EQ. 13 .OR. NCENT .EQ. 13) THEN
      DO I=1,6
      PV(I,13)=PV(I,3)
      ENDDO
      ENDIF

      IF(NTARG*NCENT .EQ. 30 .AND. NTARG+NCENT .EQ. 13) THEN
      DO I=1,6
      PV(I,3)=0.D0
      ENDDO
      GO TO 99
      ENDIF

      IF(LIST(3) .EQ. 2) THEN
      DO I=1,6
      PV(I,3)=PV(I,3)-PV(I,10)/(1.D0+EMRAT)
      ENDDO
      ENDIF

      IF(LIST(10) .EQ. 2) THEN
      DO I=1,6
      PV(I,10)=PV(I,3)+PV(I,10)
      ENDDO
      ENDIF

  99  DO I=1,6
      RRD(I)=PV(I,NTARG)-PV(I,NCENT)
      ENDDO                

      BARY=BSAVE

      RETURN
      END

      
C++++++++++++++++++++++++++++++++
C
      SUBROUTINE STATE(namfil,nrfile,ET2,LIST,PV,PNUT)
C
C++++++++++++++++++++++++++++++++
C
C THIS SUBROUTINE READS AND INTERPOLATES THE JPL PLANETARY EPHEMERIS FILE
C
C     CALLING SEQUENCE PARAMETERS:
C
C     INPUT:
C
C         ET2   DP 2-WORD JULIAN EPHEMERIS EPOCH AT WHICH INTERPOLATION
C               IS WANTED.  ANY COMBINATION OF ET2(1)+ET2(2) WHICH FALLS
C               WITHIN THE TIME SPAN ON THE FILE IS A PERMISSIBLE EPOCH.
C
C                A. FOR EASE IN PROGRAMMING, THE USER MAY PUT THE
C                   ENTIRE EPOCH IN ET2(1) AND SET ET2(2)=0.
C
C                B. FOR MAXIMUM INTERPOLATION ACCURACY, SET ET2(1) =
C                   THE MOST RECENT MIDNIGHT AT OR BEFORE INTERPOLATION
C                   EPOCH AND SET ET2(2) = FRACTIONAL PART OF A DAY
C                   ELAPSED BETWEEN ET2(1) AND EPOCH.
C
C                C. AS AN ALTERNATIVE, IT MAY PROVE CONVENIENT TO SET
C                   ET2(1) = SOME FIXED EPOCH, SUCH AS START OF INTEGRATION,
C                   AND ET2(2) = ELAPSED INTERVAL BETWEEN THEN AND EPOCH.
C
C        LIST   12-WORD INTEGER ARRAY SPECIFYING WHAT INTERPOLATION
C               IS WANTED FOR EACH OF THE BODIES ON THE FILE.
C
C                         LIST(I)=0, NO INTERPOLATION FOR BODY I
C              `                  =1, POSITION ONLY
C                                =2, POSITION AND VELOCITY
C
C               THE DESIGNATION OF THE ASTRONOMICAL BODIES BY I IS:
C
C                         I = 1: MERCURY
C                           = 2: VENUS
C                           = 3: EARTH-MOON BARYCENTER
C                           = 4: MARS
C                           = 5: JUPITER
C                           = 6: SATURN
C                           = 7: URANUS
C                           = 8: NEPTUNE
C                           = 9: PLUTO
C                           =10: GEOCENTRIC MOON
C                           =11: NUTATIONS IN LONGITUDE AND OBLIQUITY
C                           =12: LUNAR LIBRATIONS (IF ON FILE)
C
C
C     OUTPUT:
C
C          PV   DP 6 X 11 ARRAY THAT WILL CONTAIN REQUESTED INTERPOLATED
C               QUANTITIES.  THE BODY SPECIFIED BY LIST(I) WILL HAVE ITS
C               STATE IN THE ARRAY STARTING AT PV(1,I).  (ON ANY GIVEN
C               CALL, ONLY THOSE WORDS IN 'PV' WHICH ARE AFFECTED BY THE
C               FIRST 10 'LIST' ENTRIES (AND BY LIST(12) IF LIBRATIONS ARE
C               ON THE FILE) ARE SET.  THE REST OF THE 'PV' ARRAY
C               IS UNTOUCHED.)  THE ORDER OF COMPONENTS STARTING IN
C               PV(1,I) IS: X,Y,Z,DX,DY,DZ.
C
C               ALL OUTPUT VECTORS ARE REFERENCED TO THE EARTH MEAN
C               EQUATOR AND EQUINOX OF J2000 IF THE DE NUMBER IS 200 OR
C               GREATER; OF B1950 IF THE DE NUMBER IS LESS THAN 200. 
C
C               THE MOON STATE IS ALWAYS GEOCENTRIC; THE OTHER NINE STATES 
C               ARE EITHER HELIOCENTRIC OR SOLAR-SYSTEM BARYCENTRIC, 
C               DEPENDING ON THE SETTING OF COMMON FLAGS (SEE BELOW).
C
C               LUNAR LIBRATIONS, IF ON FILE, ARE PUT INTO PV(K,11) IF
C               LIST(12) IS 1 OR 2.
C
C         NUT   DP 4-WORD ARRAY THAT WILL CONTAIN NUTATIONS AND RATES,
C               DEPENDING ON THE SETTING OF LIST(11).  THE ORDER OF
C               QUANTITIES IN NUT IS:
C
C                        D PSI  (NUTATION IN LONGITUDE)
C                        D EPSILON (NUTATION IN OBLIQUITY)
C                        D PSI DOT
C                        D EPSILON DOT
C
C           *   STATEMENT # FOR ERROR RETURN, IN CASE OF EPOCH OUT OF
C               RANGE OR I/O ERRORS.
C
C
C     COMMON AREA STCOMX:
C
C          KM   LOGICAL FLAG DEFINING PHYSICAL UNITS OF THE OUTPUT
C               STATES. KM = .TRUE., KM AND KM/SEC
C                          = .FALSE., AU AND AU/DAY
C               DEFAULT VALUE = .FALSE.  (KM DETERMINES TIME UNIT
C               FOR NUTATIONS AND LIBRATIONS.  ANGLE UNIT IS ALWAYS RADIANS.)
C
C        BARY   LOGICAL FLAG DEFINING OUTPUT CENTER.
C               ONLY THE 9 PLANETS ARE AFFECTED.
C                        BARY = .TRUE. =\ CENTER IS SOLAR-SYSTEM BARYCENTER
C                             = .FALSE. =\ CENTER IS SUN
C               DEFAULT VALUE = .FALSE.
C
C       PVSUN   DP 6-WORD ARRAY CONTAINING THE BARYCENTRIC POSITION AND
C               VELOCITY OF THE SUN.
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      SAVE

      DIMENSION ET2(2),PV(6,12),PNUT(4),T(2),PJD(4),BUF(1500),
     . SS(3),CVAL(400),PVSUN(6)

      INTEGER LIST(12),IPT(3,13)
      integer*4 nrfile,numde,ncon,nrecl,ksize,irecsz,ncoeffs,nrl,nr
     .        , i,j,k

      LOGICAL FIRST
      DATA FIRST/.TRUE./

      CHARACTER*6 TTL(14,3),CNAM(400)
      CHARACTER*80 NAMFIL

      LOGICAL KM,BARY,swap

      COMMON/EPHHDR/CVAL,SS,AU,EMRAT,NUMDE,NCON,IPT
      COMMON/CHRHDR/CNAM,TTL
      COMMON/STCOMX/KM,BARY,PVSUN

C
C       ENTRY POINT - 1ST TIME IN, GET POINTER DATA, ETC., FROM EPH FILE
C
      IF(FIRST) THEN
        FIRST=.FALSE.

C ************************************************************************
C ************************************************************************

C THE USER MUST SELECT ONE OF THE FOLLOWING BY DELETING THE 'C' IN COLUMN 1

C ************************************************************************

C        CALL FSIZER1(NRECL,KSIZE,NRFILE,NAMFIL)  
        CALL FSIZER2(NRECL,KSIZE,NRFILE,NAMFIL,swap)
C        CALL FSIZER3(NRECL,KSIZE,NRFILE,NAMFIL)

      IF(NRECL .EQ. 0) WRITE(*,*)'  ***** FSIZER IS NOT WORKING *****'

C ************************************************************************
C ************************************************************************

      IRECSZ=NRECL*KSIZE
      NCOEFFS=KSIZE/2

        OPEN(NRFILE,
     *       FILE=NAMFIL,
     *       ACCESS='DIRECT',
     *       FORM='UNFORMATTED',
     *       RECL=IRECSZ,
     *       STATUS='OLD')
c      print *,'opened in STATE',nrfile,namfil,irecsz,swap

      READ(NRFILE,REC=1)TTL,CNAM,SS,NCON,AU,EMRAT,
     . ((IPT(I,J),I=1,3),J=1,12),NUMDE,(IPT(I,13),I=1,3)    

      if( swap ) then   
        call swap_bytes(4,ncon,1)   
        call swap_bytes(4,numde,1)
        call swap_bytes(8,au,1)
        call swap_bytes(8,emrat,1) 
        do i=1,3
           call swap_bytes(8,ss,3)  
        enddo
c        print *,'swapped ncon,numde,ss,au,emrat',ncon,numde,ss,au,emrat
        do j=1,13   
          do i=1,3
            call swap_bytes(4,ipt(i,j),1)
          enddo
        enddo
      endif
c      do j=1,13
c        print *,(ipt(i,j),i=1,3)
c      enddo  
 

      READ(NRFILE,REC=2)CVAL   
      if( swap ) then
        call swap_bytes(8,cval,400)
      endif
c      print *,'CVAL 1-6 ',(cval(i),i=1,6) 


      NRL=0

      ENDIF


C       ********** MAIN ENTRY POINT **********


      IF(ET2(1) .EQ. 0.D0) RETURN

      S=ET2(1)-.5D0
      CALL SPLIT(S,PJD(1))
      CALL SPLIT(ET2(2),PJD(3))
      PJD(1)=PJD(1)+PJD(3)+.5D0
      PJD(2)=PJD(2)+PJD(4)
      CALL SPLIT(PJD(2),PJD(3))
      PJD(1)=PJD(1)+PJD(3)

C       ERROR RETURN FOR EPOCH OUT OF RANGE

      IF(PJD(1)+PJD(4).LT.SS(1) .OR. PJD(1)+PJD(4).GT.SS(2)) GO TO 98

C       CALCULATE RECORD # AND RELATIVE TIME IN INTERVAL

      NR=IDINT((PJD(1)-SS(1))/SS(3))+3
      IF(PJD(1).EQ.SS(2)) NR=NR-1

        tmp1 = DBLE(NR-3)*SS(3) + SS(1)
        tmp2 = PJD(1) - tmp1
        T(1) = (tmp2 + PJD(4))/SS(3)

C       READ CORRECT RECORD IF NOT IN CORE

      IF(NR.NE.NRL) THEN
        NRL=NR
        READ(NRFILE,REC=NR,ERR=99)(BUF(K),K=1,NCOEFFS)
        if( swap ) then
           call swap_bytes(8,buf,ncoeffs)
        endif
      ENDIF

      IF(KM) THEN
      T(2)=SS(3)*86400.D0
      AUFAC=1.D0
      ELSE
      T(2)=SS(3)
      AUFAC=1.D0/AU
      ENDIF

C   INTERPOLATE SSBARY SUN

      CALL INTERP(BUF(IPT(1,11)),T,IPT(2,11),3,IPT(3,11),2,PVSUN)

      DO I=1,6
      PVSUN(I)=PVSUN(I)*AUFAC
      ENDDO

C   CHECK AND INTERPOLATE WHICHEVER BODIES ARE REQUESTED

      DO 4 I=1,10
      IF(LIST(I).EQ.0) GO TO 4

      CALL INTERP(BUF(IPT(1,I)),T,IPT(2,I),3,IPT(3,I),
     & LIST(I),PV(1,I))

      DO J=1,6
       IF(I.LE.9 .AND. .NOT.BARY) THEN
       PV(J,I)=PV(J,I)*AUFAC-PVSUN(J)
       ELSE
       PV(J,I)=PV(J,I)*AUFAC
       ENDIF
      ENDDO

   4  CONTINUE

C       DO NUTATIONS IF REQUESTED (AND IF ON FILE)

      IF(LIST(11).GT.0 .AND. IPT(2,12).GT.0)
     * CALL INTERP(BUF(IPT(1,12)),T,IPT(2,12),2,IPT(3,12),
     * LIST(11),PNUT)

C       GET LIBRATIONS IF REQUESTED (AND IF ON FILE)

      IF(LIST(12).GT.0 .AND. IPT(2,13).GT.0)
     * CALL INTERP(BUF(IPT(1,13)),T,IPT(2,13),3,IPT(3,13),
     * LIST(12),PV(1,11))

      RETURN

  98  WRITE(*,198)ET2(1)+ET2(2),SS(1),SS(2)
 198  format(' ***  Requested JED,',f12.2,
     * ' not within ephemeris limits,',2f12.2,'  ***')

      stop

   99 WRITE(*,'(2F12.2,A80)')ET2,'ERROR RETURN IN STATE'

      STOP

      END

C++++++++++++++++++++++++
C
      SUBROUTINE FSIZER2(NRECL,KSIZE,NRFILE,NAMFIL,swap)
C
C++++++++++++++++++++++++
C  THIS SUBROUTINE OPENS THE FILE, 'NAMFIL', WITH A PHONY RECORD LENGTH, READS 
C  THE FIRST RECORD, AND USES THE INFO TO COMPUTE KSIZE, THE NUMBER OF SINGLE 
C  PRECISION WORDS IN A RECORD.  
C
C  THE SUBROUTINE ALSO SETS THE VALUES OF  NRECL, NRFILE, AND NAMFIL.

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      SAVE

      CHARACTER*6 TTL(14,3),CNAM(400)
      CHARACTER*80 NAMFIL 

      DIMENSION SS(3)

      INTEGER*4 IPT(3,13),nrecl,ksize,nrfile,mrecl,nbodyf,ncon
     .        , numde,kmx,khi,nd,itest,ioerr,i,j

      character*16 prog_name
      integer*4 rcpar,len


      logical swap,debug/.false./
                   
      len = rcpar(0,prog_name) 

C  *****************************************************************
C  *****************************************************************
C
C  THE PARAMETERS NRECL, NRFILE, AND NAMFIL ARE TO BE SET BY THE USER
C
C  *****************************************************************

C  NRECL=1 IF "RECL" IN THE OPEN STATEMENT IS THE RECORD LENGTH IN S.P. WORDS
C  NRECL=4 IF "RECL" IN THE OPEN STATEMENT IS THE RECORD LENGTH IN BYTES
C  (for UNIX, it is probably 4)
C
      NRECL= 4

C  NRFILE IS THE INTERNAL UNIT NUMBER USED FOR THE EPHEMERIS FILE

cxx      NRFILE=12

C  NAMFIL IS THE EXTERNAL NAME OF THE BINARY EPHEMERIS FILE


css      NAMFIL= 'unxp2000.405'

C  *****************************************************************
C  *****************************************************************

C  **  OPEN THE DIRECT-ACCESS FILE AND GET THE POINTERS IN ORDER TO 
C  **  DETERMINE THE SIZE OF THE EPHEMERIS RECORD

      MRECL=NRECL*1000

        OPEN(NRFILE,
     *       FILE=NAMFIL,
     *       ACCESS='DIRECT',
     *       FORM='UNFORMATTED',
     *       RECL=MRECL,
     *       STATUS='OLD',iostat=ioerr )    
      if( ioerr.ne.0 ) then    
        call report_stat('FATAL',prog_name,'lib/ephred',namfil
     .        ,'Error opening JPL ephemeris file ',ioerr)
      else               
        call report_stat('STATUS',prog_name,'lib/ephred',namfil
     .                  ,'Opened JPL ephemeris',0)
      endif  

      READ(NRFILE,REC=1)TTL,CNAM,SS,NCON,AU,EMRAT
     * ,((IPT(I,J),I=1,3),J=1,12),NUMDE,(IPT(I,13),I=1,3)   
c     see if byte-swapped    
      if( ncon.lt.0.or.ncon.gt.1000) then
        swap= .true.     
        if(debug) write(*,'(/,a,/)') '**Byte-swapping invoked '
      else
        swap = .false.
      endif       
      if( swap ) then   
        call swap_bytes(4,ncon,1)   
        call swap_bytes(4,numde,1)
        call swap_bytes(8,au,1)
        call swap_bytes(8,emrat,1) 
        do i=1,3
           call swap_bytes(8,ss,3)
        enddo
        do j=1,13   
          do i=1,3
            call swap_bytes(4,ipt(i,j),1)
          enddo
        enddo
      endif
 
      CLOSE(NRFILE)



C  FIND THE NUMBER OF EPHEMERIS COEFFICIENTS FROM THE POINTERS

      KMX = 0
      KHI = 0

      DO I = 1,13
         IF (IPT(1,I) .GT. KMX) THEN
            KMX = IPT(1,I)
            KHI = I
         ENDIF
      ENDDO

      ND = 3
      IF (KHI .EQ. 12) ND=2

      KSIZE = 2*(IPT(1,KHI)+ND*IPT(2,KHI)*IPT(3,KHI)-1)     
c      print *,'ksize ',ksize 

      RETURN

      END
C+++++++++++++++++++++++++++++++++
C
      SUBROUTINE INTERP(BUF,T,NCF,NCM,NA,IFL,PV)
C
C+++++++++++++++++++++++++++++++++
C
C     THIS SUBROUTINE DIFFERENTIATES AND INTERPOLATES A
C     SET OF CHEBYSHEV COEFFICIENTS TO GIVE POSITION AND VELOCITY
C
C     CALLING SEQUENCE PARAMETERS:
C
C       INPUT:
C
C         BUF   1ST LOCATION OF ARRAY OF D.P. CHEBYSHEV COEFFICIENTS OF POSITION
C
C           T   T(1) IS DP FRACTIONAL TIME IN INTERVAL COVERED BY
C               COEFFICIENTS AT WHICH INTERPOLATION IS WANTED
C               (0 .LE. T(1) .LE. 1).  T(2) IS DP LENGTH OF WHOLE
C               INTERVAL IN INPUT TIME UNITS.
C
C         NCF   # OF COEFFICIENTS PER COMPONENT
C
C         NCM   # OF COMPONENTS PER SET OF COEFFICIENTS
C
C          NA   # OF SETS OF COEFFICIENTS IN FULL ARRAY
C               (I.E., # OF SUB-INTERVALS IN FULL INTERVAL)
C
C          IFL  INTEGER FLAG: =1 FOR POSITIONS ONLY
C                             =2 FOR POS AND VEL
C
C
C       OUTPUT:
C
C         PV   INTERPOLATED QUANTITIES REQUESTED.  DIMENSION
C               EXPECTED IS PV(NCM,IFL), DP.
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      SAVE
C
      DOUBLE PRECISION BUF(NCF,NCM,*),T(2),PV(NCM,*),PC(18),VC(18)
      
      integer*4 ncf,ncm,na,ifl,np,nv,i,j,l

      DATA NP/2/
      DATA NV/3/
      DATA TWOT/0.D0/
      DATA PC(1),PC(2)/1.D0,0.D0/
      DATA VC(2)/1.D0/
C
C       ENTRY POINT. GET CORRECT SUB-INTERVAL NUMBER FOR THIS SET
C       OF COEFFICIENTS AND THEN GET NORMALIZED CHEBYSHEV TIME
C       WITHIN THAT SUBINTERVAL.
C
      DNA=DBLE(NA)
      DT1=DINT(T(1))
      TEMP=DNA*T(1)
      L=IDINT(TEMP-DT1)+1

C         TC IS THE NORMALIZED CHEBYSHEV TIME (-1 .LE. TC .LE. 1)

      TC=2.D0*(DMOD(TEMP,1.D0)+DT1)-1.D0

C       CHECK TO SEE WHETHER CHEBYSHEV TIME HAS CHANGED,
C       AND COMPUTE NEW POLYNOMIAL VALUES IF IT HAS.
C       (THE ELEMENT PC(2) IS THE VALUE OF T1(TC) AND HENCE
C       CONTAINS THE VALUE OF TC ON THE PREVIOUS CALL.)

      IF(TC.NE.PC(2)) THEN
        NP=2
        NV=3
        PC(2)=TC
        TWOT=TC+TC
      ENDIF
C
C       BE SURE THAT AT LEAST 'NCF' POLYNOMIALS HAVE BEEN EVALUATED
C       AND ARE STORED IN THE ARRAY 'PC'.
C
      IF(NP.LT.NCF) THEN
        DO 1 I=NP+1,NCF
        PC(I)=TWOT*PC(I-1)-PC(I-2)
    1   CONTINUE
        NP=NCF
      ENDIF
C
C       INTERPOLATE TO GET POSITION FOR EACH COMPONENT
C
      DO 2 I=1,NCM
      PV(I,1)=0.D0
      DO 3 J=NCF,1,-1
      PV(I,1)=PV(I,1)+PC(J)*BUF(J,I,L)
    3 CONTINUE
    2 CONTINUE
      IF(IFL.LE.1) RETURN
C
C       IF VELOCITY INTERPOLATION IS WANTED, BE SURE ENOUGH
C       DERIVATIVE POLYNOMIALS HAVE BEEN GENERATED AND STORED.
C
      VFAC=(DNA+DNA)/T(2)
      VC(3)=TWOT+TWOT
      IF(NV.LT.NCF) THEN
        DO 4 I=NV+1,NCF
        VC(I)=TWOT*VC(I-1)+PC(I-1)+PC(I-1)-VC(I-2)
    4   CONTINUE
        NV=NCF
      ENDIF
C
C       INTERPOLATE TO GET VELOCITY FOR EACH COMPONENT
C
      DO 5 I=1,NCM
      PV(I,2)=0.D0
      DO 6 J=NCF,2,-1
      PV(I,2)=PV(I,2)+VC(J)*BUF(J,I,L)
    6 CONTINUE
      PV(I,2)=PV(I,2)*VFAC
    5 CONTINUE
C
      RETURN
C
      END

C+++++++++++++++++++++++++
C
      SUBROUTINE SPLIT(TT,FR)
C
C+++++++++++++++++++++++++
C
C     THIS SUBROUTINE BREAKS A D.P. NUMBER INTO A D.P. INTEGER
C     AND A D.P. FRACTIONAL PART.
C
C     CALLING SEQUENCE PARAMETERS:
C
C       TT = D.P. INPUT NUMBER
C
C       FR = D.P. 2-WORD OUTPUT ARRAY.
C            FR(1) CONTAINS INTEGER PART
C            FR(2) CONTAINS FRACTIONAL PART
C
C            FOR NEGATIVE INPUT NUMBERS, FR(1) CONTAINS THE NEXT
C            MORE NEGATIVE INTEGER; FR(2) CONTAINS A POSITIVE FRACTION.
C
C       CALLING SEQUENCE DECLARATIONS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION FR(2)

C       MAIN ENTRY -- GET INTEGER AND FRACTIONAL PARTS

      FR(1)=DINT(TT)
      FR(2)=TT-FR(1)

      IF(TT.GE.0.D0 .OR. FR(2).EQ.0.D0) RETURN

C       MAKE ADJUSTMENTS FOR NEGATIVE INPUT NUMBER

      FR(1)=FR(1)-1.D0
      FR(2)=FR(2)+1.D0

      RETURN

      END

c*******************************************************************************

      Subroutine TABSUN( iout,title,npr,dint,maxval,n,t,x )

c     Write a GAMIT-format output Sun table from an array of values.
         
c     Input values
c       iout      output unit number
c       title     title for output file
c       npr       number values (3 for postion-only, 6 with velocities)
c       dint       interval (days) of values in output table
c       maxval    dimension of output arrays
c       n         number of values in input arrays
c       t         array of input dates (NB: PEP JD, not MJD)
c       x         array of input coordinates
              
c     Table values set here
c        outunt   output units (input assumed to be km)
c        frame    set to J2000

      integer*4 maxval,ix(6),npjd,ifix,jd1,jd2,nvel,iout,int,npr,i,j,n
      real*8 x(6,maxval),outunt,xx,dint
      real*4 t(maxval)
      character*80 title
      character*24 varfmt
      character*5  frame
      data frame/'J2000'/,varfmt/'(1x,i5,6i11)            '/ 
     .    ,outunt/1.0d0/

      jd1 = ifix(t(1)) + 2400000
      jd2 = ifix(t(n)) + 2400000
                                            
c      print *,'In TABSUN n t x ',n,(t(i),x(1,i),i=1,n)      
      

c     tabular interval is 4 days
      if( dint.eq.4.d0 ) then
        int = 4
      else
        write(*,'(a)') 'Output interval for soltab must be 4 days'
        stop
      endif

      write(iout,10) title,varfmt,nvel,jd1,jd2,npr,int,outunt,frame
   10 format(a80,/,a24,4x,i2,1x,i7,1x,i7,1x,i2,1x,i2,1x,1pe15.0,1x,a5)

      do  j = 1,n
        do i=1,npr
          xx = x(i,j)/outunt
          ix(i) = idint(dsign(dabs(xx)+0.5d0,xx))
        enddo
        npjd = ifix(t(j))
        write(iout,varfmt) npjd,(ix(i),i=1,npr)
      enddo

      return
      end

                       
c******************************************************************************

      Subroutine TABMON( iout,title,npr,dint,maxval,n,t,x )

c     Write a GAMIT-format output Moon table from an array of values.
             
c     Input values
c       iout      output unit number
c       title     title for output file
c       npr       number values (3 for postion-only, 6 with velocities)
c       dint       interval (days) of values in output table
c       maxval    dimension of output arrays
c       n         number of values in input arrays
c       t         array of input dates (NB: PEP JD, not MJD)
c       x         array of input coordinates
    
c     Table values set here
c        outunt   output units (input assumed to be km)
c        frame    set to J2000
          
      integer*4 maxval,ix(6),npjd,ifix,jd1,jd2,nvel,iout,int,npr,i,j,n
      real*8 x(6,maxval),outunt,xx,dint
      real*4 t(maxval)
      character*80 title
      character*24 varfmt
      character*5  frame
      data frame/'J2000'/,varfmt/'(1x,i5,6i11)            '/
     .    ,outunt/1.d-3/

      jd1 = ifix(t(1)) + 2400000
      jd2 = ifix(t(n)) + 2400000
                   

c     tabular interval is integer power-of-two days
      if( dint.eq.0.5d0 ) then
        int = -1
      else
        write(*,'(a)') 'Output interval for lutab must be 0.5 days'
        stop
      endif

      write(iout,10) title,varfmt,nvel,jd1,jd2,npr,int,outunt,frame
   10 format(a80,/,a24,4x,i2,1x,i7,1x,i7,1x,i2,1x,i2,1x,1pe15.0,1x,a5)

      do  j = 1,n
        do i=1,npr
          xx = x(i,j)/outunt
          ix(i) = idint(dsign(dabs(xx)+0.5d0,xx))
        enddo
        npjd = ifix(t(j))
        write(iout,varfmt) npjd,(ix(i),i=1,npr)
      enddo

      return
      end

       
c************************************************************************************

      Subroutine TABNUT( iout,title,npr,dint,maxval,n,t,x,y )

c     Input values
c       iout      output unit number
c       title     title for output file
c       npr       number values (3 for postion-only, 6 with velocities)
c       dint       interval (days) of values in output table
c       maxval    dimension of output arrays
c       n         number of values in input arrays
c       t         array of input dates (NB: PEP JD, not MJD)
c       x         array of input coordinates
    
c     Table values set here
c        outunt   output units (input assumed to be km)
c        frame    set to J2000

      integer*4 maxval,ix(4),iy(4),npjd(4),ifix,npr,jd1,jd2,int,iout
     .        , i,j,k,n

      real*4 t(maxval),x(maxval),y(maxval),units,outunt
      data outunt/1.e-4/

      real*8 xx,yy,dint

      character*80 title
      character*24  varfmt
      data          varfmt/'(1x,i5,8i8,8x,i2)       '/
                 
c      print *,'TABNUT npr n '
c      do i=1,n
c        print *,'t x y ',t(i),x(i),y(i)
c      enddo

      jd1 = ifix(t(1)) + 2400000
      jd2 = ifix(t(n)) + 2400000
      
c     tabular interval is integer power-of-two days
      if( dint.eq.0.5d0 ) then
        int = -1
      else
        write(*,'(a)') 'Output interval for nutabl must be 0.5 days'
        stop
      endif

      write(iout,10) title,varfmt,jd1,jd2,npr,int,outunt
   10 format(a80,/,a24,11x,i7,1x,i7,1x,i2,1x,i2,1x,1pe15.0)

      i = 0
   20 k = 0
   30 i = i + 1
      if( i.gt.n ) goto 50
      k = k + 1
      xx = x(i)/outunt
      yy = y(i)/outunt
      ix(k) = idint(dsign(dabs(xx)+0.5,xx))
      iy(k) = idint(dsign(dabs(yy)+0.5,yy))
      npjd(k) = ifix(t(i))
      if( k.ne.4 ) goto 30
      write(iout,varfmt) npjd(1),(ix(j),iy(j),j=1,4)
      goto 20
   50 continue
      if( k.eq.0 ) goto 70
      write(iout,varfmt) npjd(1),(ix(j),iy(j),j=1,k)
      write(6,60) k
   60 format(/,1x,'Need to edit nutation table to add ',i1, 'in col 80')

   70 return
      end
















  










