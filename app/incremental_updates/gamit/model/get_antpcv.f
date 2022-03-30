Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego. 2003.   All rights reserved.

      Subroutine GET_ANTPCV( jd,el,az,debug
     .                     , pcoffl1,pcoffl2,corrl1,corrl2 )

c Purpose    :  To compute phase centre offsets, as a function of elevation
c               angle and (optionially) azimuth) for a receiver or satellite 
c               antenna.  Reads only an ANTEX file, no longer old GAMIT format.
c
c Parameters :
c         In:  jd        : PEP Julian day of observations            i*4
c              el        : elevation angle of sat degrees            r*8
c              az        : azmith of satellite degrees               r*8
c              antdaz    : Deviation of antenna from true north.
c                          Saved in model.h (TAH 200205)             r*8
c              newant    : true if tables not read earlier           l*4
c              debug     : true if debug printout requested          l*4
c                  In common /lunits/ in model.h       
c              ipcv      : unit number of antmod.dat (ANTEX) file    i*4
c              iepcv     : unit number of site-dependent ANTEX file  i*4
c                  In common /systyp/ in global.h           
c              gnss      : GNSS system calibrations requested        c*1 
c                  In common /rcvantcom/ in model.h                     
c              atxfrq(2) : ANTEX frequency code (in model.h)         c*3
c              anttyp    : 20-character antenna name                 c*20
c              radome_in : requested radome                          c*5
c              radome    : radome available in ANTEX file            c*5  
c              antmod_in : requested model (NONE, ELEV, AZEL)        c*4
c              antmod    : model available (NONE, ELEV, AZEL)        c*4
c                  In common /lfcom/ in model.h
c              sitecd    : site code                                 c*4
c                  In common /antcom/ in model.h
c              newant    : true if tables not read earlier           l*4

c        Out :
c              pcoffl1(3): L1 phase center wrt ARP, N E U (mm)       r*8
c              pcoffl2(3): L2 phase center wrt ARP, N E U (mm)       r*8
c              corrl1    : variable phase centre  corection l1 (mm)  r*8
c              corrl2    : variable phase centre  corection l2  (mm) r*8
c                 In common /rcvantcom/ in model.h
c              antmod    : model found ( NONE ELEV AZEL )            c*4
c                (note change from old scheme of ELnn where nn is name of model)
c              radome    : NONE if radome_in not found               c*5
c              antmod_snx: full name of model from table             ch*10
c              pcvminelev: minimum elevation angle for PCV table     r*8 
c              
c
c subroutines called : ant_alias, bilin8, linear, report_stat
c
c created    :  2003/4/14 by R. King, based on read_antmod.f written
c               by S. McClusky (1995/10/17 and modified 1997/05/09
c               Last modified by R. King 2010/9/6
c
3      Implicit none

      include '../includes/dimpar.h'   
      include '../../libraries/includes/const_param.h'
      include '../includes/units.h'      
      include '../includes/global.h'
      include '../includes/model.h'
      
      logical      found_ant,found_f1,found_f2,debug
                                  
      character*1   pcvtype,lastchr,buf1,gnss_req
      character*3   rxobtyp(maxdat)
* MOD TAH 200520: Changed associated with primary F2 choice
      character*4   atxfrq_req(3)  ! Requested values (not found)
      character*5   antftype
      character*6   buf6 
      character*10  header_snx
      character*20  refant,buf20
      character*60  comments(maxtxt)     
      character*80 prog_name
      character*256 message,line

      integer*4 len,rcpar,jd,ncomments,ioerr,nchr,nblen,i,j
                              
      real*4 antex_vers

      real*8 el,az,pcoffl1(3),pcoffl2(3),corrl1,corrl2,zen,offtmp

* Added TAH 200205: Rotation matrix to rotate to orientation
*     of the antenna
      real*8 rotmat(3,3)  ! Rotation by antsaz (for NEU)
      real*8 basoff1(3), basoff2(3)   ! Basevector of PVO offset before rotation
                          ! before rotation at L1 and L2 (NEU)

      real*8 eaz          ! Effective azimith (az-antdaz) (deg)

      integer*4 cnt
       
      logical warning
      data warning/.false./ 
      data cnt / 0 /             
      save warning

c         Initialize the SINEX codes
       header_snx = ' ' 
       antmod_snx = ' ' 

c**   If the first call for this station, read all the array values into storage
c        (this conditional ends near the bottom of the subroutine)
       if( debug ) print *,'GET_ANTPCV newant anttyp radome_in '
     .   ,newant,anttyp,radome_in

      if (newant) then                                      

c**     Initialize the offset arrays
          
        do i=1,3
          offl1(i) = 0.d0
          offl2(i) = 0.d0
        enddo
        do  i=1,maxel   
          elvtabl1(i) = 0.d0
          elvtabl2(i) = 0.d0
          do  j=1,maxaz
            tablel1(i,j)=0.0d0
            tablel2(i,j)=0.0d0
          enddo  
        enddo
                              
c**     Read a site-specific model if available
               
        found_ant = .true.
        if( iepcv.ne.0 ) then   
           call report_stat('STATUS','MODEL','get_antpcv',epcvfiln
     .        ,'Reading site-specific ANTEX file',ioerr)
          rewind( iepcv )
          read(iepcv,'(60x,a5)',iostat=ioerr) antftype
          if( ioerr.ne.0) then
            call report_stat('FATAL','MODEL','get_antpcv',' '
     .        ,'Error reading first line of site-specific ANTEX file'
     .        ,ioerr)
          endif            
          rewind ( iepcv )
          call read_antex_head( iepcv,antex_vers,pcvtype,refant
     .                        , header_snx,ncomments,comments )
          if( debug ) then
           print *,'ANTEX head: vers pcvtype refant header_snx '
     .          ,              antex_vers,pcvtype,refant,header_snx
           print *,'before READ_ANTEX anttyp ',anttyp   
          endif   
          call read_antex( iepcv,antex_vers,anttyp,antsn,gnss,atxfrq,jd
     .                   , found_ant,found_f1,found_f2
     .                   , antmod_snx,dazi,zen1,zen2,dzen
     .                   , pcoffl1,pcoffl2,elvtabl1,elvtabl2
     .                   , tablel1,tablel2 )         
          if( .not.found_ant ) then
            call report_stat('WARNING','MODEL','get_antpcv',epcvfiln
     .       ,'Site-specific ANTEX file requested but antenna not found'
     .       , 0)  
           else
             write(iprnt,'(a)') 'Using an empirical PCV model'
           endif
        endif
                                                 
             
c**     Read the usual ANTEX file
           
        if( iepcv.eq.0 .or. .not.found_ant ) then   
          rewind( ipcv ) 
c         check the format
          read(ipcv,'(60x,a5)',iostat=ioerr) antftype
          if( ioerr.ne.0) then
            call report_stat('FATAL','MODEL','get_antpcv',' '
     .        ,'Error reading first line of antenna PCV table',ioerr) 
          endif
          if( antftype.ne.'ANTEX' ) then
             call report_stat('FATAL','MODEL','get_antpcv'
     .             ,' ','Old antmod.dat format no longer supported',0)
          endif               
          rewind( ipcv )    
          call read_antex_head( ipcv,antex_vers,pcvtype,refant
     .                      , header_snx,ncomments,comments )
          if( debug ) then
            print *,'ANTEX head: vers pcvtype refant header_snx '
     .          ,              antex_vers,pcvtype,refant,header_snx
            print *,'before READ_ANTEX anttyp ',anttyp   
          endif  
          call read_antex( ipcv,antex_vers,anttyp,antsn
     .                   , gnss,atxfrq,jd
     .                   , found_ant,found_f1,found_f2
     .                   , antmod_snx,dazi,zen1,zen2,dzen
     .                   , pcoffl1,pcoffl2,elvtabl1,elvtabl2
     .                   , tablel1,tablel2 )   
          if( debug) print *
     .       ,'GET_ANTPCV gnss atxfrq found_ant found_f1 found_f2 '
     .                   ,gnss,atxfrq,found_ant,found_f1,found_f2

c         If no match of antenna+radome, try antenna alone   
          if( .not.found_ant ) then
            call report_stat('WARNING','MODEL','get_antpcv','antmod.dat'
     .      ,'No PCV model match of antenna+radome, try antenna alone'
     .      , 0) 
            anttyp(17:20) = 'NONE' 
            rewind(ipcv) 
            call read_antex_head( ipcv,antex_vers,pcvtype,refant
     .                           , header_snx,ncomments,comments )
            call read_antex( ipcv,antex_vers,anttyp,antsn
     .                 , gnss,atxfrq,jd
     .                 , found_ant,found_f1,found_f2 
     .                 ,antmod_snx,dazi,zen1,zen2,dzen
     .                 , pcoffl1,pcoffl2,elvtabl1,elvtabl2
     .                 , tablel1,tablel2 )      
            if( debug) print *
     .       ,'GET_ANTPCV gnss atxfrq found_ant found_f1 found_f2 '
     .                   ,gnss,atxfrq,found_ant,found_f1,found_f2
          else
            anttyp(17:20) = radome_in(1:4)
          endif
          if( .not.found_ant ) then
            write(message,'(a,a20,a,a20,a,a4)') 
     .             'Input antenna type ',anttyp, ' (',anttyp,
     ,          ') not in ANTEX PCV file; site ',sitecd
             call report_stat( 'FATAL','MODEL','get_antpcv'
     .           ,'antmod.dat',message,0 )
          endif              

c         if antenna found but no frequency match, use GPS L1 and L2 calibrations
          if( .not.found_f1 .or. .not.found_f2 ) then   
            gnss_req = 'G'    
            atxfrq_req(1) = atxfrq(1)
            atxfrq_req(2) = atxfrq(2)
            atxfrq_req(3) = atxfrq(3)
            atxfrq(1) = 'G01'
            atxfrq(2) = 'G02' 
            rewind(ipcv) 
            call read_antex_head( ipcv,antex_vers,pcvtype,refant
     .                          , header_snx,ncomments,comments )
            call read_antex( ipcv,antex_vers,anttyp,antsn
     .                   , gnss_req,atxfrq,jd
     .                   , found_ant,found_f1,found_f2
     .                   , antmod_snx,dazi,zen1,zen2,dzen
     .                   , pcoffl1,pcoffl2,elvtabl1,elvtabl2
     .                   , tablel1,tablel2 )  
            if( debug) print *
     .        ,'GET_ANTPCV gnss atxfrq found_ant found_f1 found_f2 '
     .                    ,gnss,atxfrq,found_ant,found_f1,found_f2
          endif
          if( found_ant.and.(.not.found_f1.or..not.found_f2) ) then 
            write(message,'(a,3(1x,a3),a)') 
     .                'No ANTEX calibrations for frequency ',atxfrq_req
     .               ,' use G01 and G02'
            write(iprnt,'(a)') message 
            call report_stat('WARNING','MODEL','get_antpcv'
     .                      ,'antmod.dat',message,0)
          endif 
        endif   
c       end of test for site-specific ANTEX file
        
         
c**     Save the minimum elevation angle to pass back to the calling routine

        pcvminelev = 90.d0 - zen2    
        nel = int((zen2-zen1)/dzen) + 1  

* MOD TAH 200205: Rotate the pcoffl1,pcoffl2 if the antdaz is not zero.
        if( antdaz.ne.0 ) then
*           Rotate the PCO using R = [cos(az) sin(az) 0 ; -sin(az) cos(az) 0 ; 0 0 ];
            call report_stat('STATUS','MODEL','get_antpcv',' ',
     .                       'Rotating PCO',int(antdaz) )
            rotmat(1,1) =  cos(antdaz*pi/180)
            rotmat(1,2) = -sin(antdaz*pi/180)
            rotmat(1,3) =  0
            rotmat(2,1) = +sin(antdaz*pi/180)
            rotmat(2,2) =  cos(antdaz*pi/180)
            rotmat(2,3) =  0
            rotmat(3,1) =  0
            rotmat(3,2) =  0
            rotmat(3,3) =  1.d0
*           Rotate L1 and L2; make copies of original vectors
            basoff1 = pcoffl1
            basoff2 = pcoffl2
            pcoffl1 = matmul(rotmat,basoff1)
            pcoffl2 = matmul(rotmat,basoff2)
!           write(*,200) cnt, pcoffl1, basoff1, pcoffl2, basoff2
!200        format('PCOFF ',i4,' L1 ',2(3F12.5,2x),' L2 ', 2(3F12.5,2x))
        endif
 
        if( debug ) then  
          print *,'after READ_ANTEX jd found_ant antmod_snx '
     .                            , jd,found_ant,antmod_snx
          print *,'dazi zen1 zen2 dzen '
     .           , dazi,zen1,zen2,dzen
          print *,'pcoffl1 pcoffl2 ',pcoffl1,pcoffl2   
          print *,'nel pcvminelev ',nel,pcvminelev
        endif
        if( nel.gt.maxel ) then
           write(message,'(a,i3,a,i3,a)') 'Number of zenith values ('
     .           ,nel,') on ANTEX file exceeds maxel (',maxel,')'   
          call report_stat('FATAL','MODEL','get_antpcv',' '
     .                    ,message,0)    
        endif  
   

c**     See if we got the model we requested 

        if( antmod_in.eq.'NONE' ) then
           antmod = 'NONE'

        elseif( antmod_in(1:2).eq.'EL' ) then
          if( dzen.gt.0.d0 ) then
            antmod = 'ELEV' 
            nel = int((zen2-zen1)/dzen) + 1
          else    
            write(message,'(a,a20,a)') 
     .         'Requested elevation-dependent phase center model for '
     .         ,anttyp,' but only offset available' 
               call report_stat('WARNING','MODEL','get_antpcv'
     .                         ,' ',message,0)    
            antmod = 'NONE'
          endif

        elseif( antmod_in(1:2).eq.'AZ' ) then    
           if( dazi.gt.0.d0 ) then
             antmod = 'AZEL'
             nel = int((zen2-zen1)/dzen) + 1
c             MOD SCM and MJMi 120104: to allow 0.5 degree antex files,
c             move int()
             naz = int(360/dazi) + 1 
c             naz = 360/int(dazi) + 1 
           elseif( dzen.gt.0.d0) then     
             write(message,'(a,a20,a)') 
     .         'Requested AZ-EL phase center model for ',anttyp
     .         ,' but only elevation corrections available'    
             call report_stat('WARNING','MODEL','get_antpcv'
     .                       ,' ',message,0)
             antmod = 'ELEV'    
             nel = int((zen2-zen1)/dzen) + 1
           else            
             write(message,'(a,a20,a)') 
     .         'Requested AZ-EL phase center  model for ',anttyp
     .             ,' but only offset available' 
             call report_stat('WARNING','MODEL','get_antpcv'
     .                        ,' ',message,0)  
             antmod = 'NONE' 
           endif
        endif 
         
cd        debug = .true.
        if( debug ) then
          print *,'GET_ANTPCV end of newant call '
          write(*,'(19f6.2)') (elvtabl1(i),i=1,nel)
          print *,'tablel1 : '
          do j=1,naz
            write(*,'(19f6.2)')(tablel1(i,j),i=1,nel)
          enddo      
          write(*,'(19f6.2)') (elvtabl2(i),i=1,nel)
          print *,'tablel2 : '
          do j=1,naz   
            write(*,'(19f6.2)')(tablel2(i,j),i=1,nel)
          enddo
        endif         
cd        stop
     
c     End if for call of subroutine for new antenna 
      newant = .false.
      endif
               
      if( debug ) print *,
     .  'GET_ANTPCV antmod_in antmod antmod_snx '
     .     , antmod_in,antmod,antmod_snx

c**   Interpolate the values (start here if table previously read)     

c     NOTE: Convention of tables changes with release 10.08 to
c           be zenith angle rather than elevation angle

      zen = 90.d0 -el
c     check for reasonableness and warn if beyond the table limits
      if ( zen.lt.zen1 ) then
         if( .not.warning ) then
           write(message,'(2a,f7.2,a,f7.2)') 
     .      'Observed zenith angle < table minimum '
     .       ,'Use antmod.dat value for ',zen1,' degrees. NOT ',zen
           len = rcpar(0,prog_name)
           call report_stat('WARNING',prog_name,'lib/get_antpsv',' '
     .        ,message,0)
           warning = .true.
         endif
       elseif ( zen.gt.zen2 ) then
         zen  = 90.d0               
         if( .not.warning ) then
           write(message,'(2a,f7.2,a,f7.2)') 
     .       'Observed zenith angle > tables maxium '
     .       ,'Use antmod.dat value for ',zen2,' degrees. NOT ',zen
           len = rcpar(0,prog_name)
           call report_stat('WARNING',prog_name,'MODEL','get_antpcv',' '
     .       ,message,0)
           warning = .true.
         endif
       endif

c     If azimuth as well as elevation, use a 2-d interpolation
      if ( antmod.eq.'AZEL' ) then
C MOD TAH 200205: Update azimuth for miss-aligned antennas.
C       (Adding 360 ensures eff_az is 0-360)
        eaz = mod(az - antdaz+360,360.d0)
        cnt = cnt + 1
c       L1 
        call interp_azel(zen,eaz,nel,naz,dzen,dazi,zen2,tablel1,corrl1)
        if( debug ) 
     .     print *,'GET_ANTPCV inter  AZEL nel naz dzen dazi '
     .         , ' zen az tables(1&6,1) '
     .    ,nel,naz,dzen,dazi,zen,az,eaz
     .    ,tablel1(1,1),tablel1(6,1),tablel2(1,1),tablel2(6,1)
c       L2
        call interp_azel(zen,eaz,nel,naz,dzen,dazi,zen2,tablel2,corrl2)  
!       write(*,210) cnt,zen,az,eaz,naz, corrl1,corrl2
!210    format('ANTPCV ',i6,3F9.2,1x,I5,2F10.4)   
c     If elevation only, use a 1-d interpolation
      elseif (antmod.eq.'ELEV') then
c       L1  
        if( debug ) 
     .    print *,'GET_ANTPCV interg EL-only zen az tables(1&6,1) ' 
     .    ,zen,az,tablel1(1,1),tablel1(6,1),tablel2(1,1),tablel2(6,1)
* MOD RWK 180720: Replace this with a new more general version using floating 
*       call linear(zen,nel,incelv,elvtabl1,corrl1)
        call linear(zen,nel,dzen,zen1,zen2,elvtabl1,corrl1)
c       L2
*      call linear(zen,nel,dzen,elvtabl2,corrl2)
       call linear(zen,nel,dzen,zen1,zen2,elvtabl2,corrl2)
      else
        corrl1 = 0.d0
        corrl2 = 0.d0
      endif 
      if(debug) print *,'GET_ANTPCV zen az corrl1 corrl2 '
     .                  ,zen,az,corrl1,corrl2
              
      return
      end
              





