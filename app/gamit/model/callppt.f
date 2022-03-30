C Was Program to test calling piercept subroutine & to get lat long of piercept
C Want this program to make the conversion between inertial and earthfixed sat coords and then calculate ppt..

      Subroutine callppt(satcoord,recoord,unvec,slat,slon,debug,ionht)
      implicit none

C-------------------------------------------------------      
C Declare passed variables:
C In:
C satcoord	Satellite coordinates (cartesian)/km
C recoord	receiver coordinates (cartesian)/km
C Out:
C unvec		unit vector satellite to site, geocentric cartesian
C slat, slon    geocentric (spherical) coordinates for ppcoord
      Real*8 slat, slon
      Real*8 satcoord(3), recoord(3) 
      Real*8 unvec(3)
C Other variables
C ionht		height of 'single layer' ionosphere/km
C eradius	Earth's radius /km
C calcrad       Calculated radius of pierce point/m as a check 
      Real*8 ionht, eradius, calcrad
C ppcoord	pierce point coordinates  (cartesian) km, then m.
      Real*8 ppcoord(3)
C Declare variable passed for xyz_to_geod
      Real*8 rot_mat(3,3),geod_pos(3)
      Integer*4 k,IRC,debug
      Real*8 pi

      pi=4.0D0*atan(1.0D0)

C Define variables
C can change ionht to represent different thin layer height. but now passed through from iondel.f
C      ionht=450.0D0
      eradius=6371.0D0
C----------------------------------------------     


C Call subroutine to calculate cartesian coords of pierce point
      call pierce_pt(satcoord, recoord, ionht, eradius, 
     1 ppcoord,unvec, IRC,debug)
C need error message if IRC not 1

C      Print*, 'MODEL\callppt satcoord, recoord, ppcoord '
C      Print*, satcoord, recoord, ppcoord
      Do k=1,3
       ppcoord(k)=ppcoord(k)*1000.0D0
      Enddo
C      Print*, 'MODEL\callppt  Cartesian coords X,Y,Z /m '
C      Print*, ppcoord

Call sub to convert coordinates to spherical coordinates
      Call xyz2sph (ppcoord,slat,slon,calcrad )
        if(debug.ge.3) then
      Print*,'MODEL/callppt slat/rad, slon/rad',slat,slon
        endif
      slat=slat*180.0D0/pi
      slon=slon*180.0D0/pi
        if(debug.ge.3) then
      Print*,'MODEL/callppt slat, slon/deg,calcrad',slat, slon,calcrad
        endif


      End
