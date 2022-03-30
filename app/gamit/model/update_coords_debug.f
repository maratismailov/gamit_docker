      Subroutine update_coords( sitecd, skd, iepoch, jdobs, tobs
     .                        , iuh,iprnt,iant
     .                        , semi,finv,shft,antcod,newant,evec0)

c     Update the Earth-fixed coordinates with new position (eq/rename)
c     velocity,  or antenna offsets (station.info)

c     R. King 8 May 2003 from (now-obsolete) routine update_offs (RWK 14 Feb 97)

c     Input
          
c         sitecd  C*4  : 4-character station name    
c         skd     C*1  : 'S', 'K', or 'D' for static, kinematic, or dynamic
c         iepoch  I*4  : epoch of observation
c         jdobs   I*4  : PEP JD of observation time
c         tobs    R*8  : GPST seconds-of-day of observation time  
c         iuh     I*4  : unit number for station.info  
c         iprnt   I*4  : unit number for P-file  
c         iant    I*4  : unit number for antmod.dat file  
c         semi    R*8  : semi-major axis of geodetic datum (m)
c         finv    R*8  : inverse flattening of datum
c         shft(3) R*8  : center-of-figure offset of datum (m)

        
c      Output 
                 
c         antcod     C*6 : Updated antenna type    
c         newant     L*4 : True if new station.info entry found 
c         evec0(3,2) R*8 : Updated cartesian coordinates for L1 and L2

c     Output new L-file coordinates through common /kinpar1/ in modkin.h

      implicit none

      include '../includes/dimpar.h'
      include '../includes/modkin.h'


      logical old_stinf
           
      integer*4 icall,iy,iday,isessn
     .        , iuh,iprnt,iepoch,jdobs,iant
     .        , nyr,ndoy,nsod,itimdif,itime(5),i
                                                                      
      real*4 swver
      real*8 tobs,anthgt,offstn,offste,dhpab,offsl1(3),offsl2(3)
     .     , evec0(3,2),sitepart0(3,3),gdlat,gdhgt,sec,semi,finv,shft(3)
     .     , pos(3),decyrs,tobs_yrs
                      
      character*1 skd,upper1
      character*4 sitecd,trkcod,prject,orbit,sitcod
      character*5 hgtcod,radome
      character*6 rcvcod,antcod
      character*16 stanam,amodel 
      character*256 message 
     
      logical newcrd,newant              

c  Variables and logic for stations coordinates:
c
c     L-file values read as spherical (old-style) or Cartesian (GLOBK-style) 
c         but now stored only as Cartesian in the commons of modkin.h.
c
c     Earth-fixed values corrected for antenna offsets and velocity (kvel)
c         are evec0(3,2) (km)
c
c     Inertial values are evec(6,2) (km)
c
c     The L-file values are read initially in SETUP and written into the
c     C-file header.  They can be updated during a session (UPDATE_COORDS
c     called by MODEL), however, if there is an earthquake or other rename 
c     (eq_rename.dat entry) or change in antenna type or offset (new station.info 
c     entry), or if velocities are included (GLOBK apr-file) and exceed 1 mm/day.  
c     For kinematic or dynamic measurements, they are updated at each epoch.
c     Whenever the L-file values are updated, the Earth-fixed and inertial
c     values are recomputed.  
c      
c  Keys for updates in modkin.h commons:
c      kstarts, kstops : current stnfo entry     defaults 
c      kstartr, kstopr:  current eq/rename entry   default 1900, 2100
c      kvflg : true if velocities to be used 
  
c  Initialize flags for new initial coordinates, offsets, or position 

      newcrd = .false.

C  Convert JD, sod to yr,doy,hr,min,sec [itime(1-5)] for station.info/eq-rename checks
c  Get time in years for coordinate updates
                              
        print *,'UPDATE_COORS jdobs ',jdobs
        call dayjul(jdobs,itime(1),itime(2))   
        print *, 'yr day ',itime(1),itime(2)
        call ds2hms( itime(1),itime(2),tobs,itime(3),itime(4),sec) 
        itime(5) = int(sec)   
        print *,  'itime kstart kstopr ',itime,kstartr,kstopr 
        tobs_yrs = decyrs( itime(1),itime(2),tobs )

              
c  If kinematic or earthquake/rename then read new values from the L-file
       
      if( upper1(skd).eq.'K' .or.upper1(skd).eq.'D' .or.
     .    itimdif(itime,kstopr).gt.0 ) then 
        newcrd = .true.                
        call lread( sitecd,tobs_yrs )  
c       note in the p-file that we've done an update
        write(message,'(a,i5,i4,2i3,a )') 'Coords updated at '
     .                                ,(itime(i),i=1,4),' --see P-file'
        call report_stat('STATUS','MODEL','update_coords',' ',message
     .                  ,0)
        write(iprnt,'(a,i5,2x,i4,1x,i3,1x,3i3)')
     .      'Coordinates updated at epoch ',iepoch,(itime(i),i=1,5)     
        write(iprnt,'(3f9.4,2x,3f9.4,f11.4)')
     .        ' New values (m  m/yr) :    '
     .        ,kpos,kvel,kepoch0  
      endif   
                 

c  If outside the bounds of the current station.info entry, read another one

      print *,'itime kstarts kstops ',itime,kstarts,kstops
      if( itimdif(itime,kstarts).lt.-60 ) then
c       obs time more than 1 min earlier than station.info time, entry is missing
        write(message,'(a,i5,i4,3i3,a,i5,i4,3i3,a)') 'Obs time ('
     .         ,itime,') too early for station.info start (',kstarts,')'
        call report_stat('FATAL','MODEL','model',' ',message,0) 
      elseif( itimdif(itime,kstops).gt.60 ) then
        newant = .true.
c       obs time more than 1 min later than station.info stop time, 
c       read another entry and update coordinates  
       print *,'updating station.info entry sitecd itime,kstops '
c     .        ,sitecd,itime,kstops
        call dayjul(jdobs,iy,iday)
        call check_oldstnfo( iuh, old_stinf )
        if( old_stinf ) then 
c         get the most recent entry prior to the current time
          icall = 3        
c         don't check project and orbit codes
          prject = '    '
          orbit  = '    '  
c         set sitcod = blank to check only on trkcod  (=sitcod for static)
          sitcod = '    '
          trkcod = sitecd
c         don't check session number
          isessn = 99    
          call even_minute(iy,iday,tobs,nyr,ndoy,nsod)
          call rstnfo( iuh, icall, prject, orbit
     .              , trkcod, sitcod, nyr, ndoy, isessn, nsod
     .              , stanam,  anthgt, offstn, offste, rcvcod, antcod
     .              , hgtcod, swver, kstarts, kstops )
c        print *,'iuh,icall,prject,orbit,trkcod,sitcod,nyr,ndoy,isessn,
c    .   nsod ',iuh,icall,prject,orbit,trkcod,sitcod,nyr,ndoy,isessn,nsod
c        print *,'stanam,anthgt,offstn,offste,rcvcod,antcod,hgtcod,swver'
c    .       ,stanam,anthgt,offstn,offste,rcvcod,antcod,hgtcod,swver
c        print *,'kstarts kstops ',kstarts,kstops
        else          
         print *,'calling even_minute iy iday tobs ', iy,iday,tobs
          call even_minute(iy,iday,tobs,nyr,ndoy,nsod)   
         print *,'called even_minute ',nyr,ndoy,nsod
          call rstnfo2( iuh, kinflg, trkcod, sitecd, nyr, ndoy, nsod, 0
     .              , stanam, anthgt, offstn, offste, rcvcod, antcod
     .              , hgtcod, radome, swver, kstarts, kstops )
        endif
       print *,'new anthgt kstarts kstops ',anthgt,kstarts,kstops
c        Convert the raw antenna offset position to phase-center-above-mark
c        The quantities offsl1 and offsl2 are offsets of L1 and L2 phase centres
c        from the monument, computed by hisub from the offsets measured in the
c        field (anthgt,offstn,offste) and the L1/L2 APR to phase centre offsets
c        read from antmod.dat.
        call hisub(iant,anthgt,offstn,offste,antcod
     .        , hgtcod,dhpab,trkcod,iy,iday,isessn,offsl1,offsl2,amodel)
c       Note in the P-file that we've done an update
        write(message,'(a,i5,i4,2i3,a )') 'Antenna offsets updated at '
     .                                ,(itime(i),i=1,4),' --see P-file'
        call report_stat('STATUS','MODEL','update_offs',' ',message,0)
        write(iprnt,'(a,i5,2x,i4,1x,i3,1x,3i3)')
     .       'Offset from monument updated at epoch '
     .       ,iepoch,(itime(i),i=1,5)
        write(iprnt,'(a,1x,a6,1x,a5,3f8.4)') 
     .    '  New values (antcod hgtcod U N E): ',antcod,hgtcod,anthgt
        write(iprnt,'(a,6f8.4)') 
     .'             (L1/L2 U N E): ',offsl1,offsl2
        l1z = offsl1(1)
        l1n = offsl1(2)
        l1e = offsl1(3)
        l2z = offsl2(1)
        l2n = offsl2(2)
        l2e = offsl2(3)    
      endif
                     

c  If velocities to be used, update the position  
              
      if( kvflg.gt.0 ) then 
        do i=1,3
          pos(i) = kpos(i) + kvel(i)*(tobs_yrs-kepoch0)
        enddo   
      else
        do i=1,3
          pos(i) = kpos(i)
        enddo
      endif 
c      write(*,'(a,3f17.7)') '   new pos ',pos
         

c     Apply antenna offsets and compute the geodetic and cartesian 
c     coordinates and partials 
             
      if( newcrd .or. newant .or.kvflg.gt.0 ) then
         
        call cortyp( pos,l1z,l1n,l1e,l2z,l2n,l2e
     .             , semi,finv,shft,evec0,gdlat,gdhgt,sitepart0 )  
c        write(*,'(a,3f15.7)') '   evec0 ',(evec0(i,1),i=1,3)  
      endif


      return
      end

