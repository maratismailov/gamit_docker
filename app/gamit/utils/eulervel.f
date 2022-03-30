c     Program eulervel

c     Determine the horizontal velocity and uncertainty of a station from Euler vectors

c     Rename of program medvel_new (R. King Jan 1997;  S. McClusky Jun 1998, Jul 2002)

      implicit none

      integer i,j,iel,iplate,plate_num,indx,ierr,jerr,trimlen,iwhere
     .     , pnum

      real*8 lat,lon,rad,pi,rad2deg,fin(3),a,finv,out(3),pos(3),vel(3)
     .     , wvec(3),wsig(3),wcorrel(3,3),svec(3),ssig(3),scorrel(3,3) 
     .     , neu(3),neusig(3),neucor(3,3),epoch
     .     , rate,ratesig,az,azsig,necor,razcor,values(9)
     .     , xrot(100),xrot_sig(100),yrot(100),yrot_sig(100),zrot(100)
     .     , zrot_sig(100),rho_xy(100),rho_xz(100),rho_yz(100)
     .     , geod_pos(3),rot_mat(3,3)

            
      character*1 crd_type
      character*4 plate_name(100),plate
      character*8 site_name 
      character*10 word, pcode
      character*80 apr_file, plate_file, out_file
      character*132 buffer

      pi = 4.d0*atan(1.d0)
      rad2deg = 180.d0/pi

c-----------------------------------------------------------------

       write(*,'(/,a,/,a,/,a,/,a,/,a,//,a/)')
     . 'Utility to calculate station velocities and uncretainties using'
     .,'plate motion euler vectors and their correlation matricies.'
     .,'You have to enter a lat, lon, and plate-motion vector code.'
     .,'Plate motion vectors are defined in the eulervel.dat file.'
     .,'eulervel.dat file entries are in the glorg euler vector format.'
     .,'Available plate-motion vectors codes:' 

*      Open the eulervel.dat file containing the plate-motion vectors
       open(10, file='eulervel.dat', iostat=ierr, status='old')
       if ( ierr .ne. 0 ) then
         print*,'File eulervel.dat does not exist.'
         stop
       endif
       
       plate_num = 1
       iplate = 0 
       jerr = 0  

*      Read the eulervel.dat input data file extract available plate names
*      and their associated euler vectors 
10     do while ( jerr.eq.0 )
         read(10,'(a)', iostat=jerr ) buffer
 
*        See if we can find the plate-motion vector names
         indx = 1
         if( trimlen(buffer).gt.0 .and. buffer(1:1).eq.' ' ) then
           iplate = iplate + 1
           call getword( buffer, plate_name(iplate), indx)  
           write(*,'(i4,a,a)') iplate,'. ', plate_name(iplate)
           iel = 1
         else
           iel = 0
         end if
 
*        Decode the line and add new information if found
         if( iel.gt.0 ) then
             call multiread( buffer, indx, 'R8', ierr, values,
     .                       ' ',9)

           xrot(iplate)     =  values(1)
           xrot_sig(iplate) =  values(2)
           yrot(iplate)     =  values(3)
           yrot_sig(iplate) =  values(4)
           zrot(iplate)     =  values(5)
           zrot_sig(iplate) =  values(6)
           rho_xy(iplate)   =  values(7)
           rho_xz(iplate)   =  values(8)
           rho_yz(iplate)   =  values(9)
cd           print*,'input ',plate_name(iplate),xrot(iplate),yrot(iplate),
cd     .     zrot(iplate),xrot_sig(iplate),yrot_sig(iplate),zrot(iplate),
cd     .     rho_xy(iplate),rho_xz(iplate),rho_yz(iplate)
         endif
       enddo   

cd     uncomment this line to display 'test' plate case
cd     .,'   (0.01 +/- 0.001 rad/My on equator and prime meridian:  test'
              
*     See what type of coordinates the User wants to input     
      write(*,'(/,a,/,a,/,a)') 
     . 'Are station coordinates:', 
     . 'geodetic (g), spherical (s), globk.apr file (a), plate file (p)'
      read(*,'(a)') crd_type 
      if( crd_type.eq.'g' .or. crd_type.eq.'G' ) then
        word = 'geodetic  '
      elseif ( crd_type.eq.'s' .or. crd_type.eq.'S' ) then
        word = 'geocentric'
      elseif ( crd_type.eq.'a' .or. crd_type.eq.'A' ) then
        word = 'globk'
      elseif ( crd_type.eq.'p' .or. crd_type.eq.'P' ) then
        word = 'plate'
      else  
         write(*,'(/,2a)') 'Invalid coordinate type (',crd_type,')'
         write(*,'(a)') 'Start over and enter g, s, a or p'
         stop
      endif

*     Get the coordinates or .apr file name, and plate-motion code from the user
      if ( word .ne. 'globk' .and. word .ne. 'plate' ) then
        write(*,'(/,a,a,a)') 
     .  'Enter ',word(1:trimlen(word)),
     .  ' N-lat E-lon (deg) & plate-motion code'
        read(*,*) lat,lon,plate_num 
      elseif ( word .eq. 'globk' ) then
        write(*,'(/,a,a,a)') 
     .  'Enter ',word(1:trimlen(word)),
     .  ' .apr file name '
        read(*,'(a)',iostat=ierr)  apr_file
        if( ierr.ne.0) then
          print *,'Error reading apr file',ierr
          stop   
        endif
        open(15, file=apr_file, iostat=ierr, status='old')
        if ( ierr .ne. 0 ) then
          print*,'File: ',apr_file,' does not exist.'
          stop   
        else
          print *,'Opened: ',apr_file
        endif                  
        write(*,*) 'Enter plate-motion code number'
        read(*,*,iostat=ierr)  plate_num
        if( ierr.ne.0 ) then
          print *,'Error reading plate-motion code',ierr
          stop
        endif
        iwhere = index(apr_file,'.')
        out_file = apr_file(1:iwhere)//"vel"
        open(20, file=out_file, status='unknown') 
        write(20,15)
 15     format( /,'# VELOCITY ESTIMATES FROM EULERVEL  ',/,
     .     '#  Long.     Lat.',8x,'E & N Rate ',4x,
     .     ' E & N Adj. ',4x,' E & N +-',2x,
     .     ' RHO ',6x,' H Rate   H adj.    +-',2x,'SITE',/,
     .     '#',2x,'(deg)    (deg)',3x,3(7x,'(mm/yr)'),17x,
     .     '(mm/yr)' )
      elseif ( word .eq. 'plate' ) then
        write(*,'(/,a,a,a)') 
     .  'Enter ',word(1:trimlen(word)),
     .  ' .apr file name and plate-motion code'
        read(*,*) plate_file
        open(15, file=plate_file, iostat=ierr, status='old')
        if ( ierr .ne. 0 ) then
          print*,'File: ',plate_file,' does not exist.'
          stop
        endif
        iwhere = index(plate_file,'.')
        out_file = plate_file(1:iwhere)//"vel"
        open(20, file=out_file, status='unknown') 
        write(20,15)
      endif

      if (plate_num .gt. iplate) then
         write(*,'(a,i4)') 'Invalid plate code (',plate_num,')'
         stop
      endif

*     Get the relative plate-motion angular velocity vector 
      wvec(1) =  xrot(plate_num)/rad2deg
      wvec(2) =  yrot(plate_num)/rad2deg
      wvec(3) =  zrot(plate_num)/rad2deg   
      wsig(1) =  xrot_sig(plate_num)/rad2deg
      wsig(2) =  yrot_sig(plate_num)/rad2deg 
      wsig(3) =  zrot_sig(plate_num)/rad2deg  
      wcorrel(1,1) = 1.00000 
      wcorrel(2,1) = rho_xy(plate_num)
      wcorrel(3,1) = rho_xz(plate_num)
      wcorrel(1,2) = wcorrel(2,1)
      wcorrel(2,2) = 1.00000
      wcorrel(3,2) = rho_yz(plate_num)
      wcorrel(1,3) = wcorrel(3,1)
      wcorrel(2,3) = wcorrel(3,2)
      wcorrel(3,3) = 1.00000

      jerr = 0 

*     Loop over the input globk.apr file entries or input station coordinates
20    do while ( jerr.eq.0 ) 
        if ( word .eq. 'globk' .or. word .eq. 'plate' ) then
          read(15,'(a)', iostat=jerr ) buffer
 
*         See if we can find the plate-motion vector names
          indx = 1
          if( trimlen(buffer).gt.0 .and. buffer(1:1).eq.' ' ) then
            call getword( buffer, site_name, indx)  
            iel = 1
          else
            iel = 0
          end if
 
*         Decode the line and add new information if found
          if( iel.gt.0 ) then
            call multiread( buffer, indx, 'R8', ierr, values,
     .                    ' ',7)
            pos(1)     =  values(1)
            pos(2)     =  values(2)
            pos(3)     =  values(3)
            vel(1)     =  values(4)
            vel(2)     =  values(5)
            vel(3)     =  values(6)
            epoch      =  values(7)
cd            print*,'input ',site_name,pos(1),pos(2),pos(3)
cd     .      ,vel(1),vel(2),vel(3),epoch
            if ( word .eq. 'plate' ) then
              call getword( buffer, plate, indx)
              pnum = 1
              pcode = 'not_found'
50            do while ( pcode .eq. 'not_found' )
                print*,'plate_name plate',plate_name(pnum), plate
                if ( plate_name(pnum) .eq. plate ) then
                  write(*,'(a,a4)') 'Valid plate code found (',plate,')'
                  pcode = 'found'
*     Get the relative plate-motion angular velocity vector 
                  wvec(1) =  xrot(pnum)/rad2deg
                  wvec(2) =  yrot(pnum)/rad2deg
                  wvec(3) =  zrot(pnum)/rad2deg   
                  wsig(1) =  xrot_sig(pnum)/rad2deg
                  wsig(2) =  yrot_sig(pnum)/rad2deg 
                  wsig(3) =  zrot_sig(pnum)/rad2deg  
                  wcorrel(1,1) = 1.00000 
                  wcorrel(2,1) = rho_xy(pnum)
                  wcorrel(3,1) = rho_xz(pnum)
                  wcorrel(1,2) = wcorrel(2,1)
                  wcorrel(2,2) = 1.00000
                  wcorrel(3,2) = rho_yz(pnum)
                  wcorrel(1,3) = wcorrel(3,1)
                  wcorrel(2,3) = wcorrel(3,2)
                  wcorrel(3,3) = 1.00000
                endif
                if ( pnum .eq. iplate ) then
                    write(*,'(a,a4)') 'Invalid plate code (',plate,')'
                    stop
                endif
                pnum = pnum + 1
              enddo
            endif  

*           Get the stations geocentric lat and lon 
            rad= sqrt(pos(1)**2+pos(2)**2+pos(3)**2)
            lon= atan2(pos(2),pos(1))
            lat= asin(pos(3)/rad)
            lon= lon*180.d0/pi
            lat= lat*180.d0/pi 
cd            print*,'lat,lon,rad ',lat,lon,rad 

*           Get the stations geodetic lat and lon 
            call XYZ_to_GEOD( rot_mat, pos, geod_pos )
            geod_pos(1) = 90.d0-geod_pos(1)*180.d0/pi
            geod_pos(2) = geod_pos(2)*180.d0/pi
          else
            goto 20
          endif

*       Convert the station coordinates to Cartesian
        elseif ( word .ne. 'globk' ) then
          jerr = 1           
          a = 6378137.d0
          fin(1) = lat/rad2deg
          fin(2) = lon/rad2deg
          fin(3) = 0.d0   
          if ( crd_type.eq.'g' .or. crd_type.eq.'G' ) then 
            finv = 298.257223563d0  
            call gdetic( fin,out,a,finv )
          elseif ( crd_type.eq.'s' .or. crd_type.eq.'S' ) then
            out(1) = fin(1)
            out(2) = fin(2)
          endif
          pos(1) = a*cos(out(1))*cos(out(2))
          pos(2) = a*cos(out(1))*sin(out(2))
          pos(3) = a*sin(out(1))
          lat = out(1)*rad2deg
          lon = out(2)*rad2deg
        endif

*       Get the Cartesian station velocity and uncertainties
                       
        call velpole( pos, wvec, wsig, wcorrel, svec, ssig, scorrel )
        write(*,'(//,a)')
     .  'Geocentric coordinates   lat     lon (deg)      x        y      
     .   z (m)'
        write(*,'(1x,a8,12x,2f9.2,4x,3f10.0,/)') site_name,lat,lon,pos
        write(*,'(a,1x,a4,2a,/,a,6f11.6)')
     .  'from',plate_name(plate_num),' motion'
     .  ,'  wx        sig        wy         sig        wz        sig'
     .  , '   (rad/My)  ', (wvec(i),wsig(i),i=1,3)
        write(*,'(//,a,/,6f8.4)')  
     .  'Cartesian velocities and sigmas (m/yr)',(svec(i),ssig(i),i=1,3)
        write(*,'(/,a,/,3(3f12.8,/))')
     .  'Correlation matrix XYZ',((scorrel(i,j),j=1,3),i=1,3)
                  
*       Convert the Cartesian velocities to local horizontal and output
        call xyzneu (lat, lon, svec, ssig, scorrel, neu, neusig, neucor) 

*       Convert the local north and east to magnitude and direction         
*       Compute the value and sigma of the magntitude
        necor = neucor(1,2)
        call magdir ( neu,neusig,necor,rate,ratesig,az,azsig,razcor )
        write(*,'(/,a,/,1x,a8,16x,11(f8.2))')
     .  'Local velocities (mm/yr)  North     +/-     East    +/-      Up
     .     +/-    Rate     +/-    Az      +/-    Correl.'
     .  ,site_name,(neu(i)*1000.,neusig(i)*1000.,i=1,3)
     .  ,rate*1000.,ratesig*1000.,az*rad2deg,azsig*rad2deg,razcor  
        write(*,'(a,/,3(3f7.3,/))')
     .  'Correlation matrix NEU',((neucor(i,j),j=1,3),i=1,3)  

*       Write output globk.vel velocity file
        if ( word .eq. 'globk' .or. word .eq. 'plate' ) then 
          write(20,220)   geod_pos(2), geod_pos(1),
     .    neu(2)*1.d3, neu(1)*1.d3, 0.d0, 0.d0,
     .    neusig(2)*1.d3, neusig(1)*1.d3, 0.001d0,
     .    neu(3)*1.d3, 0.d0, neusig(3)*1.d3, site_name
 220      format(2(1x,f8.3),1x,6(1x,f7.2),1x,f6.3,2x,
     .                 3(1x,f7.2), 1x,a8)
        endif

      enddo

      stop
      end
       

