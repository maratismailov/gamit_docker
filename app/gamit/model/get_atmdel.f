Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      Subroutine GET_ATMDEL( jd,t,prn,az,el,atdel,part )
C
C       Calculate the atmospheric delay using one of several models.
C       MODIFIED 27-MAY-82 BY RIA FOR GPS
C         R. Abbot - 27 May 82
C         J. Davis - 15 Feb 87
C         R. King  - 4 June 87, 17 May 96    
c         P. Tregoning/R. King - Mar-Jul 2006
C
C      Input
C
c       IUZ  :  Unit number for met print (z-) file  (in common /units/ of model.h)
c       JD   :  PEP Julian Day of observation    
c       T    :  Seconds of day observation
C       HEIGHT: Site geodetic height (km) (model.h)
C       LATR:   Site geodetic latitude (radians)  (model.h)
c       LONR:   Site longitude (radians)   
C       PRN :   Satellite PRN number
c       AZ  :   Satellite azimuth (radians)
C       EL  :   Satellite elevation (radians)
c         Met values/flags from from common /metcom/ in ../includes/model.h 
c       METSRC C*3  Source of met values (RNX UFL GPT STP PTH)
C       DRYZEN C*4  Model to be used for the dry zenith delay
C       WETZEN C*4  Model to be used for the wet zenith delay
C       DRYMAP C*4  Mapping function to be used for the dry delay
C       WETMAP C*4  Mapping function to be used for the wet delay
C       PRES0   R*8 Default (constant) pressure (sealevel or local) (millibars) 
c                    from STP or GPT
C       TEMP0   R*8 Default (constant) temperature (sealevel or local) (degrees C)    
c       WETVAR0 R*8 Relative humidity (0-1) from batch file humidity*100 or GPT e0      
c       **----local variables PRES, TEMP, WETVAR may come from time-dependent measurements
c             local variables SEALOC_P and SEALOC_T indicates whether sea-level (S) or 
c             local (L) values of P and T
c       SENSOR_HT R*8 Height of sensor on RINEX met file (m)
c
c         U-file variables available from common/ufcom/ in ../includes/model.h     
c       METMOD C*8  Source of P, T, H values
c       MAPMOD C*8  Source of mapping function values   
c       NTMET  I*4  Number of epochs of met values (RINEX of u-file) 
c       MET_VAL(ntmet,3): P,T,H  at times MET_TIME(ntmet)        
c       NTMAP  I*4  Number of epochs of mapping function values (u-file)
c       NMAP   I*4  Number of mapping function values each epoch
c       MAP_VAL(ntmap,nmap): Mapping function values at times MET_TIME 
c                            May also include ZHD if no met-file values on u-file

C      Output
C
c       ATDEL: Total atmospheric delay (seconds)
c       PART : Atmospheric delay partial (seconds / second)
c       OKMET: Bits set if values are availble without extrapolation from a
c              RINEX met file or from a VMF1 mode (0 if from GPT or STP)
c               Binary coded  Bits: 1 P   2 T   4 wet   8 wvr
c               (model.h)

c     References;
c       Saastamoinen, J., Atmospheric correction for the troposphere and 
c         statosphere in radio ranging of satellites, in "The Use of Artificial
c         Satellites for Geodesy", AGU Monograph Ser. 15, eds. S.W. Henrikson et al,
c         247-251, AGU, Washington, D.C., 1972. 

      implicit none         

      include '../includes/dimpar.h'
      include '../../libraries/includes/const_param.h'
      include '../includes/model.h'

      integer*4 jd,idoy,iyr,ihr,min,prn,imap,ngood
            
      logical grdvmf1            
c       Flags indicating whether pressure, temperature, water-vapor, and/or ZHD available
      logical pflg,tflg,wflg,zflg
          
      real*4 metval(3),mapval(3),exttim_met,exttim_map
      real*8 tkelv,tkelvs,e,gm,gm0,R,dmap,wmap
     .     , el,delev,hgt,dmjd
     .     , part,latd
     .     , zhd,zwd,ztd,atdel,radian
     .     , hmf(2),wmf(2),doy
     .     , zd,shd,swd,sec,t,az                         
     .     , bh,c0h,phh,c11h,c10h,ch,bw,cw
     .     , sine,beta,gamma,topcon
     .     , a_ht,b_ht,c_ht,ht_corr_coef,ht_corr
   
c       Functions                    
      integer*4 match_name
      logical fcheck 

 
c* WVR not used
c     real*4 awvr,bwvr,cwvr  
c       T1,T2 : WVR brightness temperatures
c       AWVR, BWVR, CWVR:  WVR algorithm coefficients

c       mapval : 3-element array of interpolated VMF1 coefficients 
c                Ah Aw ZHD
                               
      character*1 sealoc_p,sealoc_t
      character*8 metmodz 
      character*256 message

      logical first,early_met_warning,late_met_warning

c     debug only
c      integer*4 i

      data early_met_warning/.false./,late_met_warning/.false./

c     Initialize 'first' for setup call
      data first /.true./

cd      print *,'ATMDEL DEBUG prn el dryzen wetzen drymap wetmap '
cd    .       ,             prn,el,dryzen,wetzen,drymap,wetmap
cd      print *,'  pres0 temp0 wetvar0 ',pres0,temp0,wetvar0
cd      print *,' iuz ',iuz    
cd      print *,'jd t ',jd,t
                            
c.. Initialize return variables to zero
      atdel = 0.0D+00
      part  = 0.0D+00

c.. Degrees to radians
      radian = datan(1.d0)/45.d0

c..  Calculate date and time for mapping functions and z-file
      call dayjul( jd,iyr,idoy)         
      doy = dfloat(idoy) + t/86400.d0 
                           
c.. Local validity flags for P,T, Rel. Hum, ZHD
      pflg = .false.
      tflg = .false.  
      wflg = .false.
      zflg = .false.

   
c.. If values from RINEX file or u-file, interpolate to the current epoch
      if( metsrc.eq.'RNX'.or. lmet ) then  
cd       print *,'ATMDEL idoy t doy ',idoy,t,doy
cd       print *, 'ntmet met_time ',ntmet,(met_time(i),i=1,ntmet)
        call lininterp( maxmet,3,met_time,met_val,doy,ntmet,3,metval
     .                , exttim_met)    
cd       if( doy.gt.200.082 ) print *,'ATMDEL doy exttim '
cd    .          ,doy,exttim
      endif
      if( lmap ) then   
cd       print *,'map_name ntmap nmap ',map_name,ntmap,nmap
       call lininterp(maxmap,nmap,map_time,map_val,doy,ntmap,nmap,mapval
     .               , exttim_map )   
cd      print *,'doy mapval ',doy,mapval   
      endif

c.. Get P T H and/or ZHD from available values  
c         (if we combine 'map' and 'met' files, this code can be shortened considerably)

c....RINEX met file
                      
      if( metsrc.eq.'RNX') then   
c       For RINEX (already extrapolated 3 hours, use only if within range)
cd       print *,'TEST exttim ',exttim
        if( exttim_met.eq.0. ) then
c         check pressure validity
          call check_met(1,1,'P',metval(1),ngood)
          if( ngood.gt.0 ) then
            pres = metval(1)  
            sealoc_p = 'L'
            pflg = .true.
          else      
            write(message,'(a,f6.1,a,f8.3,a,f6.1)') 
     .        'Bad RINEX pressure value (',metval(1),') at ',doy,'
     .         , using ',pres0
            call report_stat('WARNING','MODEL','atmdel',' ' ,message,0) 
            pres = pres0   
c           sealoc assumes pres0 is from RINEX or GPT
            sealoc_p = 'L'  
          endif  
c         check temperature validity
          call check_met(1,1,'T',metval(2),ngood)
          if( ngood.gt.0 ) then
            temp = metval(2)  
            sealoc_t = 'L'
            tflg = .true.
          else      
            write(message,'(a,f6.1,a,f8.3,a,f6.1)') 
     .        'Bad RINEX temperature value (',metval(2),') at ',doy,'
     .         , using ',temp0
            call report_stat('WARNING','MODEL','atmdel',' ' ,message,0) 
            temp = temp0   
c           sealoc assumes temp0 is from RINEX or GPT
            sealoc_t = 'L'  
          endif  
c         check humidity validity
          call check_met(1,1,'H',metval(3),ngood)
          if( ngood.gt.0 ) then
            wetvar = metval(3)/100.      
            wflg = .true.
          else      
            write(message,'(a,f6.1,a,f8.3,a,f5.2)') 
     .        'Bad RINEX humidity value (',metval(3),') at ',doy,'
     .         , using ',wetvar0*100.
            call report_stat('WARNING','MODEL','atmdel',' ' ,message,0) 
            wetvar = wetvar0   
          endif
        else  
          if( exttim_met.lt.0.and..not.early_met_warning ) then
              write(message,'(a,f7.3,a,f7.3,a,2f6.1,f5.2)') 'Obs time ('
     .        ,doy,') too early for  RINEX met file start (',met_time(1)
     .        ,') setting values from defaults (',pres0,temp0,wetvar0
             call report_stat('WARNING','MODEL','atmdel',' ' ,message,0)
             early_met_warning = .true.
          endif
          if( exttim_met.gt.0.and..not.late_met_warning ) then
            write(message,'(a,f7.3,a,f7.3,a,2f6.1,f5.2)') 'Obs time ('
     .       ,doy,') too late for  RINEX met file end (',met_time(ntmet)
     .       ,') setting values from defaults (',pres0,temp0,wetvar0
            call report_stat('WARNING','MODEL','atmdel',' ' ,message,0)
            late_met_warning = .true.
          endif
          pres = pres0    
          temp = temp0
          wetvar = wetvar0
          sealoc_p = 'L'  
        endif    

c.... U-file (VMF1 for now)

      elseif( metsrc.eq.'UFL' ) then
        if( lmet ) then
          call report_stat('FATAL','MODEL','atmdel',' ' 
     .       ,'u-file P T H from met grid or list not yet coded',0)  
        endif
        if( lmap ) then           
c         For u-file, allow the observation to be 3 hrs before or after the entries
          if( abs(exttim_map)-0.125.le.0. ) then
            imap = match_name(nmap,2,map_name,'ZH' )
            if( imap.gt.0 ) then  
c             check reasonableness of zenith delay 
              call check_met(1,1,'Z',mapval(imap),ngood)
              if( ngood.gt.0 ) then
                zhd = mapval(imap)  
                zflg = .true.
              else      
                write(message,'(a,f6.1,a,f8.3,a,f6.1)') 
     .            'Bad u-file ZHD value (',mapval(imap),') at ',doy,'
     .             , using ',zhd0
                call report_stat('WARNING','MODEL','atmdel',' ' 
     .                           ,message,0) 
                zhd = zhd0   
              endif    
            endif
            imap = match_name(nmap,2,map_name,'SP' )
            if( imap.gt.0 ) then  
c             check pressure validity
              call check_met(1,1,'P',mapval(imap),ngood)
              if( ngood.gt.0 ) then
                pres = mapval(1)  
                sealoc_p = 'S'
                pflg = .true.
              else      
                write(message,'(a,f6.1,a,f8.3,a,f6.1)') 
     .            'Bad u-file pressure value (',mapval(imap),') at '
     .               ,doy,' using ',pres0
                call report_stat('WARNING','MODEL','atmdel',' ' 
     .               ,message,0) 
                pres = pres0   
c               sealoc assumes pres0 is from u-file
                sealoc_p = 'S'  
              endif    
            endif  
c           * add tests here for temperature and water vapor variables when available *
          else
            write(message,'(a,f8.3,a,2f8.3,a,2f6.1,f5.2)') 'Obs time ('
     .          ,doy,' more than 3 hrs outside range of u-file ('
     .          ,map_time(1),map_time(ntmap)
     .          ,') setting values from defaults (',pres0,temp0,wetvar0
            call report_stat('WARNING','MODEL','atmdel',' ' ,message,0) 
            pres = pres0    
            temp = temp0
            wetvar = wetvar0
            sealoc_p = 'L'  
          endif  
          if( .not.pflg .and. .not.zflg ) then   
            write(message,'(a,f8.1,a)') 
     .         'Missing U-file ZHD and pressure at ',doy,'
     .         , using initial value of ZHD '
             call report_stat('WARNING','MODEL','atmdel',' ' ,message,0)  
             zhd = zhd0  
             pres = pres0
             temp = temp0
             wetvar = wetvar0 
             sealoc_t = 'L'   
          endif  
          if( .not.tflg ) then
            temp = temp0
          endif
          if( .not.wflg ) then
             wetvar = wetvar0
          endif
        endif

c....Constant values 

      elseif( metsrc.eq.'PTH' .or. metsrc.eq.'GPT' .or.
     .         metsrc.eq.'STP' ) then
c        values assigned in setup.f 
         pres = pres0   
         temp = temp0   
         wetvar = wetvar0
         if( metsrc.eq.'GPT') then
           sealoc_t = 'L'
           sealoc_p = 'L'
         else
           sealoc_t = 'S'
           sealoc_p = 'S'
         endif
      else
        call report_stat('FATAL','MODEL','atmdel',' ' 
     .                  , 'Invalid met source',0) 
      endif
     

                                        
c.. Set flag for c-file, indicating whether the requested local values were available
      okmet = 0
      if( pflg .or. zflg ) okmet = okmet + 1
      if( tflg ) okmet = okmet + 2
      if( wflg ) okmet = okmet + 4  
                                      

c.. Get the acceleration of gravity for pressure mapping 
c       Saastamoinen [1972], Eqn. (22):  gm is at the centroid
c       of the atmsphere; gm0 is adjusted so that  height in the 
c       equation is station height  (km/sec**2)
      gm0 = 9.784d-3
      gm = gm0*(1.d0 - 0.266D-02 * cos(2.0D0 * latr ) 
     .                    - 0.28D-03 * height )     
c     R  is the specific gas constant for dry air (km**2/sec**2/degC = 28.9644d9 kg/kmol)
      R = 2.89644d-4
                     
c.. Assume a standard lapse rate for now (later read from a met file?)
c      old atmdel used 4.5 C/km (45N July); revert to standard value used for GPT
c**      alpha = 6.5d0 -- this now set (as 'lapse' in sb setup.f and stored in model.h)
                
c.. Get T in Kelvin for pressure mapping and z-file print
      tkelv = temp + 273.16   
      if( sealoc_t.eq.'S' ) then
        tkelvs = tkelv  
cd       print *,'height sensor_ht ',height,sensor_ht
        if( sensor_ht.gt.1.d5 ) then
c         sensor height set to large number if not available
          tkelv = tkelvs - lapse*height
        else
          tkelv = tkelvs - lapse*(height-sensor_ht/1.d3)
         endif
      endif                           
 
c.. Compute site pressure (mb) and temperature from sea-level or
c   sensor values if necessary; formula from Hopfield [1972]
c     
      if( sealoc_p.eq.'S' ) then 
cd       print *,'pres tkelv tkelvs gm R lapse '
cd    .     ,pres,tkelv,tkelvs,gm,r,lapse   
        pres = pres*(tkelv/tkelvs)**(gm/R/lapse)   
        temp = temp - lapse*height
      endif    
cd      print *,'sealoc_p,sealoc_t,height lapse tkelv sensor_ht pres  '
cd    .         ,sealoc_p,sealoc_t,height,lapse,tkelv,sensor_ht,pres     
cd    print *,'GET_ATMDEL zflg pflg wflg pres wetvar '
cd   .                   ,  zflg, pflg, wflg,pres,wetvar

c.. Compute the zenith hydrostatic delay 
          
      if( metsrc.ne.'RNX'.and. zflg ) then     
c       convert RINEX value from mm to m
        zhd = zhd/1.d3/vel_light    
      else     
c       Saastamoinen [1972]
        zhd = 0.2277d-2 * pres * (gm0/gm) / vel_light
      endif 
c     if zhd from u-file used, derive a pressure value for the c-file (kept 0 for the z-file)
      if( zsrc.ne.'    '.and.psrc.eq.'   ' ) then
        pres = zhd * vel_light / (.2277d-2 * (gm0/gm) )
      endif
                               
c.. Compute the zenith wet delay. 
                                  
c     Saastamoinen [1972]              
c     get partial pressure of water vapor, in mbar
c**      e = wpress(wetvar,tkelv-273.16)         
      call wpress(1,wetvar,e,tkelv-273.16)      
      zwd = 0.2277d-2 * (0.1255d4/tkelv+.05d0) * e * (gm0/gm) /vel_light
cd      print *,'wetvar e tkelv gm/gm0 zwd '
cd     .       , wetvar,e,tkelv,gm/gm0,zwd 

c.. WVR delay : not currently supported.  If values available on the u-file,
c                 read them and restore this code.
c      okwvr = -1   ** hardwired now in /cfout 
c      if (wetzen(1:3) .eq. 'WVR' ) then
c....   Calculate zenith delay in seconds
c        zwd = (awvr + bwvr * t1 + cwvr * t2) * 1.0D-02 / vel_light
c        uninitialized variable eliminated:  wetwvr
c        zwd = wetwvr
c        okwvr = 0


c...Get the hydrostatic (dry) mapping function    

       zd = radian*90.d0 - el  
       dmjd = jd + t/86400.d0 - 2400000.5d0
       hgt = height*1.d3

c.....Niell?
      if( drymap .eq. 'NMFH' ) then
        doy = idoy
        latd = latr/radian
        delev = el/radian
        call nmfh2p1( doy,latd,hgt,delev,hmf )
        dmap = hmf(1)

c.... GMF?  (Global Mapping Function by Boehm, Niell, Schuh, GRL 2005)
      elseif( drymap.eq.'GMFH'.or.drymap(1:2).eq.'GP' ) then       
c rwk 121212: if GPT2 used, no need to call GMF, but must apply the height correction
        if( fcheck('gpt.grid')) then  
          bh = 0.0029
          c0h = 0.062
          if (latd.lt.0.d0) then   ! southern hemisphere
            phh  = pi
            c11h = 0.007
            c10h = 0.002
          else                     ! northern hemisphere
            phh  = 0.d0
            c11h = 0.005
            c10h = 0.001
          endif
          ch = c0h + ((dcos(doy/365.25d0*2.d0*pi + phh)+1.d0)*c11h/2.d0 
     .             + c10h)*(1.d0-dcos(latd))  
          sine   = dsin(pi/2.d0 - zd)
          beta   = bh/( sine + ch  )
          gamma  = ah/( sine + beta)
          topcon = (1.d0 + ah/(1.d0 + bh/(1.d0 + ch)))
          dmap   = topcon/(sine+gamma)
cd          print *,'GMF  ah sine beta gamma topcon dmap '
cd     .                , ah,sine,beta,gamma,topcon,dmap
c         height correction for hydrostatic mapping function from Niell (1996)
          a_ht = 2.53d-5
          b_ht = 5.49d-3
          c_ht = 1.14d-3
          beta   = b_ht/( sine + c_ht )
          gamma  = a_ht/( sine + beta)
          topcon = (1.d0 + a_ht/(1.d0 + b_ht/(1.d0 + c_ht)))
          ht_corr_coef = 1/sine - topcon/(sine + gamma)
          ht_corr      = ht_corr_coef * height
          dmap = dmap + ht_corr
cd          print *,'ATMDEL ht_corr dmap for GPT2: ',ht_corr,dmap
        else
          call gmf( dmjd,latr,lonr,hgt,zd,hmf(1),wmf(1) )
          dmap = hmf(1)      
cd          print *,'ATMDEL dmap for GPT: ',ht_corr,dmap
        endif

c.... IMF?
      elseif( drymap.eq. 'IMFH' ) then   
        call report_stat('FATAL','MODEL','atmdel',' '
     .        ,'IMF mapping function no longer supported',0)
  
c.... VMF1?    
      elseif (drymap.eq. 'VMF1'.or.drymap.eq.'VMFH' )then   
        if( lmap ) then  
c         set the variables using the 2-character tokens from the u-file 
          ah = 0.d0
          imap = match_name(nmap,2,map_name,'AH' ) 
          if( imap.gt.0 ) ah = mapval(imap)   
cd         print *,'AH imap ah ',imap,ah
          aw = 0.d0
          imap = match_name(nmap,2,map_name,'AW' )  
          if( imap.gt.0 ) aw = mapval(imap) 
cd         print *,'AW imap aw ',imap,aw
c         if the values are from a grid, they need a height correction                          
          grdvmf1 = .false.    
          if( mapmod(5:5).eq.'G' ) grdvmf1 = .true.  
cd         if( first ) print *,'ah aw dmjd latr hgt zd ',
cd    .         ah,aw,dmjd,latr,hgt,zd
          call vmf1(ah,aw,dmjd,latr,hgt,zd,hmf(1),wmf(1),grdvmf1)    
c          if( first ) print *,'hmf1 wmf1 grdvmf1 ',hmf(1),wmf(1),grdvmf1
          dmap = hmf(1)     
        else
          call report_stat('FATAL','MODEL','atmdel',' '
     .       ,'VMF1 mapping function requested but not available',0)
        endif      
      else
        write(message,'(a,a4,a)') 'Dry mapping function ',drymap
     .    ,' not identified'
        call report_stat('FATAL','MODEL','atmdel',drymap,message,0) 
      endif

c.. Wet mapping function

c.... Neill?
      if ( wetmap. eq. 'NMFW' ) then
        latd = latr/radian
        delev = el/radian
        call nmfw2( latd, delev, wmf )
        wmap = wmf(1)

c.... GMF?  (Global Mapping Function by Boehm, Niell, Schuh, GRL 2005)
      elseif(wetmap.eq. 'GMFW'.or.wetmap(1:2).eq.'GP' ) then
        if( fcheck('gpt.grid') ) then 
          bw = 0.00146
          cw = 0.04391
          beta   = bw/( sine + cw )
          gamma  = aw/( sine + beta)
          topcon = (1.d0 + aw/(1.d0 + bw/(1.d0 + cw)))
          wmap  = topcon/(sine+gamma)
cd          print *,'ATMDEL wmap for GPT2 ',wmap
        else
          call gmf(dmjd,latr,lonr,hgt,zd,hmf(1),wmf(1) )
          wmap = wmf(1)                       
cd          print *,'ATMDEL wmap for GPT ',wmap
        endif

c.... IMF?
      elseif( drymap.eq.'IMFH' ) then 
           call report_stat('FATAL','MODEL','atmdel',' '
     .         ,'IMF mapping function no longer supported',0)

c.... VMF1?    

      elseif (wetmap.eq. 'VMF1'.or.wetmap.eq.'VMFW' ) then   
        if( drymap.eq.'VMF1'.or.drymap.eq.'VMFH') then
c          wet values already interpolated and evaluated
           wmap = wmf(1)
        else
          call report_stat('FATAL','MODEL','atmdel',' '
     .     ,'Problem using VMF1 for wet but not dry mapping function',0)
        endif

      else 
        write(message,'(a,a4,a)') 'Wet mapping function ',wetmap
     .    ,' not identified'
        call report_stat('FATAL','MODEL','atmdel',wetmap,message,0)
      endif


c.. Form the atmospheric delay
      atdel = zhd * dmap + zwd * wmap

cd          print *,'latr el height drymap dmap wetmap wmap zhd zwd ',
cd    .              latr,el,height,drymap,dmap,wetmap,wmap,zhd,zwd
cd          if( el.lt.1.57d0 ) stop

c...Write values into the weather print (z-) file 
          
      if( iuz.ne.0 ) then
        if( first ) then
          write(iuz,'(a,a4,a,f8.3)') 
     .        '* A priori atmospheric values for ',sitecd
     .       ,'   Geodetic height(m) = ',height*1.d3   
          if( metsrc.eq.'UFL'.and. lmet ) then 
             metmodz = metmod
          elseif( metsrc.eq.'UFL'.and.lmap ) then
             metmodz = mapmod   
          else
             metmodz = ' '
          endif        
          write(iuz,'(a,a3,1x,a8,1x,4(a,a4))')  
     .       '* Met source ',metsrc,metmodz
     .      ,'  Dry Zen ',dryzen,'  Wet Zen ',wetzen,'  Dry Map ',drymap
     .      ,'  Wet Map ',wetmap   
c          write(iuz,'(a)') '* Units: Az/el in deg, T in Kelvin, P in mb, delay in m'
          write(iuz,'(a)') '*'
          write(iuz,'(4a)') '* Yr  Doy Hr Mn Sec  PRN  Azimuth '
     .     ,' Elevation  Pres   Temp    WV Pres  Dry Zen (m) '
     .     ,'Wet Zen   Total Zen  Dry Map   Wet Map '
     .     ,' Dry Slant  Wet Slant Total Slant'
         first = .false.
        else           
         call ds2hms(iyr,idoy,t,ihr,min,sec)   
         ztd = (zhd + zwd)  * vel_light
         shd = (zhd * dmap) * vel_light
         swd = (zwd * wmap) * vel_light   
         write(iuz,'(1x,2i4,2i3,f4.0,i4,1x,2f9.4,2x,f8.1,f7.1,f7.1
     .        ,8f11.4)') 
     .       iyr,idoy,ihr,min,sec,prn,az/radian,el/radian
     .     , pres,tkelv,e,zhd*vel_light,zwd*vel_light,ztd,dmap,wmap
     .     , shd,swd,atdel*vel_light 
        endif  
      endif


c.... Partial derivative is wet mapping function
      part = wmap
      
c.... Debug:
cd    print *, ' PRES, TEMP, WETVAR, gd;at, height, el '
cd    print *, pres,temp,wetvar,latr,height,el
cd    print *, ' ZHD, ZWD, DMAP, WMAP',zhd,zwd,dmap,wmap
cd    print *, ' '
      return
      end
