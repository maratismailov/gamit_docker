Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego. 2003.   All rights reserved.

      Subroutine get_svantpcv( jd,ichan,nadangd,yatt,first
     .                       , found_ant,svoffl1,svoffl2,corrl1,corrl2 )
c
c Purpose    :  To compute phase centre offsets of a satellite antennas as 
c               a function of nadir angle.  Calls read_antex to read an 
c               IGS-standard ANTEX table.
c
c Parameters :
c         In:  ichan       : channel number for requested SVes            i*4
c              gnss        : system  G R C E J I (global.h)               c*1   
* MOD TAH 200512: Added primary/secondary F2 choice. 
c              atxfrq(3)   : ANTEX frequency codes to be used (model.h)   c*4
c              svantbodyx  : SV body type BLOCK IIR-M, BEIDOU-2I etc      c*20
c                            (from svantbody(ichan) in model.h)
c              iprn        : prn number for SV  (from ischan(ichan))      i*4      
c              jd          : Julian day of observations                   i*4
c              nadangd     : nadir angle of SV antenna (deg)              r*8
c              yatt        : yaw angle, a proxy for azimuth  (deg)        r*8
c              first       : .true. if call from SETUP, read in values    logical
c                          .false. if call from MODEL, just interpolate
c                In common /lunits/ in model.h
c              ipcv        : unit number of antmod.dat or ANTEX file      i*4   
c                              (model.h)
c
c        Out : found_ant   : true if found requested svantbody model      logical    
c              svoffl1(3)  : L1 phase center offset                       r*8
c              svoffl2(3)  : L2 phase center offset                       r*8 
c              svelvtabl1(maxel,maxsat), svelvtabl2(maxel,maxsat) : 
c                               Azimuth-avergaged corrections (mm)        r*8
c              svtabl1(maxel,maxaz,maxsat), svtabl2(maxel,maxaz,maxsat)
c                               Corrections by azimuth & elevation (mm)   r*8
c              corrl1      : variable phase centre (vpc) corection l1     r*8
c              corrl2      : variable phase centre (vpc) corection l2     r*8
c                These all in model.h:
c              svantmod_in : requested model (NONE, ELEV, AZEL)           c*4
c              svantmod(maxsat) : model found NONE, ELEV, AZEL            c*4
c                          Note change from prior scheme of putting a  
c                          2-character code into the last four positions
c                          (e.g, I1 for IGS 01, N3 for NGS 03)   
c                (note change from old scheme of ELnn where nn is name of model)   
c              svantmod_snx(maxsat): Full name of model from ANTEX file   c*10
c                          (replaces old 'amodel' c*16)     

      Implicit none

      include '../includes/dimpar.h'  
      include '../includes/global.h'
      include '../includes/model.h'
      
      logical      first,found_ant,found_f1,found_f2
                                  
      character*1   pcvtype,lastchr,gnss_req 
      character*5   antftype
      character*10  header_snx
      character*20  refant,svantsn,svantbodyx
      character*60  comments(maxtxt)
      character*80  prog_name
      character*256 message   

      integer*4 ichan,iprn,len,rcpar,ncomments,ioerr,nblen,jd
     .        , nchr,i,j

      real*4 antex_vers               

      real*8 nadangd,yatt,az,corrl1,corrl2,svoffl1(3),svoffl2(3) 

      logical nadir_warning/.false./,pco_warning/.false./,debug/.false./

*     Make the azimuth positive for interpolating from the ANTEX file
      az = yatt        
      if(yatt.lt.0.d0) az = 360.d0 + yatt

      if( debug ) then
        print *,'Entering GET_SVANTPCV ipcv gnss svantbody '
     .        ,ipcv,gnss,svantbody(ichan)
        print *,'jd nadangd ',jdobs,nadangd
        print *,' yatt az ',yatt,az 
        print *,' svantmod_in first svantmod '
     .           ,svantmod_in,first,svantmod(ichan)
      endif
       
       iprn = ischan(ichan)                 
       svantbodyx = svantbody(ichan)
       svantsn = ' '                                         
c      serial number seems to be the system+PRN                    
c      
       write(svantsn(1:3),'(a1,i2)') gnss,iprn
       if(svantsn(2:2).eq.' ') svantsn(2:2) = '0' 

c**   If the first call for this satellite, read all the array values into storage
c        (this conditional ends near the bottom of the subroutine)
       
      if (first) then 
                      
        rewind(ipcv) 
 
c**     Initialize the offset arrays
                  
        do  i=1,maxel   
          svelvtabl1(i,ichan) = 0.d0
          svelvtabl2(i,ichan) = 0.d0
          do  j=1,maxaz
            svtabl1(i,j,ichan)=0.0d0
            svtabl2(i,j,ichan)=0.0d0
          enddo  
        enddo 
        svantmod(ichan) = '    '

c**     Determine the format of the table

        read(ipcv,'(60x,a5)',iostat=ioerr) antftype
        if( ioerr.ne.0) then
          call report_stat('FATAL','MODEL','get_svantpcv',' '
     .       ,'Error reading first line of antenna PCV table',ioerr)  
        else
          if( antftype.ne.'ANTEX') antftype = 'GAMIT' 
        endif 
        rewind(ipcv)  


c**     Read the tabular values into storage
                         
        if( antftype.ne.'ANTEX' )                
     .    call report_stat('FATAL','MODEL','get_svantpcv'
     .             ,' ','Old antmod.dat format no longer supported',0) 
        if( debug ) print *,'GET_SVANTPCV svantsn ',svantsn
        call read_antex_head( ipcv,antex_vers,pcvtype,refant
     .                      , header_snx,ncomments,comments )
cd        print *,'GET_SVANTPCV bef read_antex jd ',jd
       gnss_req = gnss
       call read_antex(ipcv,antex_vers,svantbodyx,svantsn
     .                , gnss_req,atxfrq,jd
     .                , found_ant,found_f1,found_f2
     ,                , svantmod_snx(ichan),svdazi(ichan)
     .                , svzen1(ichan),svzen2(ichan),svdzen(ichan)
     .                , svoffl1,svoffl2
     .                , svelvtabl1(1,ichan),svelvtabl2(1,ichan)
     .                , svtabl1(1,1,ichan),svtabl2(1,1,ichan) )
       if( debug ) print *
     .   ,'GET_SVANTPCV gnss_req atxfrq found_ant found_f1 found_f2 '
     .                , gnss_req,atxfrq,found_ant,found_f1,found_f2 

c       no swapping of order for offsets needed for SV (both ANTEX and GAMIT are X Y Z)
        if( .not.found_ant.and.pcvtype.eq.'A' ) then
          write(message,'(a,a20,a,i2,a)') 'Input SV antenna '
     .         ,svantbodyx,' for PRN ',iprn
     .         ,' not in ANTEX PCV file'
          call report_stat( 'WARNING','MODEL','get_svantpcv',' '
     .                    , message,0 )        
        else   
         svnel(ichan)=int((svzen2(ichan)-svzen1(ichan))/svdzen(ichan))+1
        endif                                               
        if( svnel(ichan).gt.maxel ) then
          write(message,'(a,i3,a,i3,a)') 
     .      'Number of nadir-angle values (',svnel(ichan)
     .       ,') on ANTEX file exceeds maxel (',maxel,')'
          call report_stat('FATAL','MODEL','get_svantpcv',' '
     .                    ,message,0)    
        endif  


c**     See if we got the antenna and model we requested 

        if( found_ant ) then 

          if( svantmod_in.eq.'NONE' ) then
             svantmod(ichan) = 'NONE'    

          elseif( svantmod_in(1:2).eq.'EL' ) then
           if( svdzen(ichan).gt.0.d0 ) then
             svantmod(ichan) = 'ELEV'            
             svnel(ichan) = 
     .         int((svzen2(ichan)-svzen1(ichan))/svdzen(ichan)) + 1
           else
             write(message,'(a,a20,a)') 
     .        'Requested nadir-angle-dependent phase center model for '
     .        ,svantbodyx,' but only offset available' 
              call report_stat('WARNING','MODEL','get_svantpcv'
     .                         ,' ',message,0)
             svantmod(ichan) = 'NONE'
           endif

          elseif( svantmod_in.eq.'AZEL' ) then
            if( svdazi(ichan).gt.0.d0 ) then
             svantmod(ichan) = 'AZEL'
             svnel(ichan) = 
     .         int((svzen2(ichan)-svzen1(ichan))/svdzen(ichan)) + 1
             svnaz(ichan) = int(360/svdazi(ichan)) + 1
            elseif( svdzen(ichan).gt.0.d0) then
             write(message,'(a,a20,a)') 
     .         'Requested AZ-EL phase center model for ',svantbodyx
     .           ,' but only nadir-angle corrections available'    
             call report_stat('WARNING','MODEL','get_svantpcv'
     .                     ,' ',message,0)
             svantmod(ichan) = 'ELEV'    
             svnel(ichan) = 
     .          int((svzen2(ichan)-svzen1(ichan))/svdzen(ichan)) + 1
            else
              write(message,'(a,a20,a)') 
     .       'Requested AZ-EL phase center model for ',svantbodyx
     .           ,' but only offset available' 
             svantmod(ichan) = 'NONE' 
             call report_stat('WARNING','MODEL','get_svantpcv'
     .                   ,' ',message,0)
            endif
          endif 
 
        else
          svantmod(ichan) = 'NONE' 
          svantmod_snx(ichan) = ' '  
          do i=1,3
            svoffl1(i) = 0.d0
            svoffl2(i) = 0.d0
          enddo  
          return
        endif     
               
c     End if for first call of routine call of subroutine (SETUP)
      endif

c**   Interpolate the values (start here if table previously read)     

c     check for reasonableness and warn if beyond the table limits
      if ( nadangd.lt.svzen1(ichan) ) then
         nadangd = svzen1(ichan)
         if( .not.nadir_warning ) then
           write(message,'(2a,f7.2,a,f7.2)') 
     .      'Observed nadir angle < table minimum '
     .       ,'Use antmod.dat value for ',svzen1(ichan)
     .       ,' degrees. NOT ',nadangd
           len = rcpar(0,prog_name)
           call report_stat('WARNING','MODEL','get_svantpcv',' '
     .        ,message,0)
           nadir_warning = .true.
         endif
       elseif ( nadangd.gt.svzen2(ichan) ) then
         nadangd  = svzen2(ichan)
         if( .not.nadir_warning ) then
           write(message,'(2a,f7.2,a,f7.2)') 
     .       'Observed zenith angle > tables maxium '
     .       ,'Use antmod.dat value for ',svzen2(ichan)
     .       ,' degrees. NOT ',nadangd
           len = rcpar(0,prog_name)
           call report_stat('WARNING','MODEL','get_svantpcv',' '
     .       ,message,0)
           nadir_warning = .true.
         endif
       endif

      if( debug ) then
         print *,'GET_SVANTPCV ichan  maxel svnel svelvtabl1(1-5) '
     .      ,ichan,maxel,svnel(ichan),(svelvtabl1(i,ichan),i=1,1) 
         print *,'GET_SVANTPCV ichan  maxel svnel svelvtabl2(1-5) '
     .      ,ichan,maxel,svnel(ichan),(svelvtabl2(i,ichan),i=1,5)    
        print *,'  nadangd svantmod_in svantmod '
     .      ,nadangd,svantmod_in,svantmod(ichan)
      endif

c     If azimuth as well as elevation, use a 2-d interpolation
      if ( svantmod(ichan).eq.'AZEL' ) then
c       L1   
        call interp_azel( nadangd,az,svnel(ichan),svnaz(ichan)
     .                 , svdzen(ichan),svdazi(ichan),svzen2(ichan)
     .                 , svtabl1(1,1,ichan),corrl1)
        if( debug ) 
     .   print *,'GET_SVANTPCV interpolating AZEL L1 nel naz dzen dazi '
     .         , ' nadangd az svtabl1(1,1,ichan),svtabl1(2,3,ichan) '
     .    ,svnel(ichan),svnaz(ichan),svdzen(ichan),svdazi(ichan)
     .    ,nadangd,az,svtabl1(1,1,ichan),svtabl1(2,3,ichan)
c       L2
        call interp_azel( nadangd,az,svnel(ichan),svnaz(ichan)
     .                 , svdzen(ichan),svdazi(ichan),svzen2(ichan)
     .                 , svtabl2(1,1,ichan),corrl2)
        if(debug) 
     .   print *,'GET_SVANTPCV interpolating AZEL L2 nel naz dzen dazi '
     .         , ' nadangd az svtabl2(1,1,ichan) svtabl2(2,3,ichan) '
     .    ,svnel(ichan),svnaz(ichan),svdzen(ichan),svdazi(ichan)
     .    ,nadangd,az,svtabl2(1,1,ichan),svtabl2(2,3,ichan)

      elseif(svantmod_in(1:2).eq.'EL'.and.
     .       svantmod(ichan)(1:2).eq.'EL') then
c       L1
        call linear( nadangd,svnel(ichan),svdzen(ichan)
     .             , svzen1(ichan),svzen2(ichan),svelvtabl1(1,ichan)
     .             , corrl1)
        if(debug)
     .   print *,'GET_SVANTPCV interpolating ELEV L1 nel dzen  '
     .        , ' nadangd svelvtabl1(1,ichan) svelvtabl1(3,ichan) '
     .    , svnel(ichan),svnaz(ichan),svdzen(ichan)
     .    , nadangd,svelvtabl1(1,ichan),svelvtabl1(3,ichan)

c       L2
        call linear( nadangd,svnel(ichan),svdzen(ichan)
     .            , svzen1(ichan),svzen2(ichan),svelvtabl2(1,ichan)
     .            , corrl2 )
        if(debug)
     .   print *,'GET_SVANTPCV interpolating ELEV L2 nel naz dzen  '
     .        , ' nadangd svelvtabl2(1,ichan) svelvtabl2(3,ichan) '
     .    , svnel(ichan),svnaz(ichan),svdzen(ichan)
     .    , nadangd,svelvtabl2(1,ichan),svelvtabl2(3,ichan)
 
      else
        corrl1 = 0.d0
        corrl2 = 0.d0
      endif
                 
      if(debug)  print *, 'GET_SVANTPCV ichan corrl1 corrl2 '
     .          ,  corrl1,corrl2

      return
      end

