      Subroutine write_kmlf( luout,icall,site,x )

c       icall = 1  write header
c              =2   write a site coordinate record

      implicit none
      
      include '../../libraries/includes/const_param.h'    

      integer*4 luout,icall

      real*8 x(3),geod_pos(3),latd,lond,ht,rot(3,3)
 
      character*8 site

             

      if( icall.eq.1 ) then
        write(luout,10)
10      format('<?xml version="1.0" encoding="UTF-8"?>',/
     .         ,'<kml xmlns="http://earth.google.com/kml/2.0">',/
     .         ,'<Document>',/
     .         ,'  <name>My GPS Sites.kml</name>',/
     .         ,'  <Style id="GPSSiteStyle">',/
     .         ,'    <IconStyle id="GPSSiteIconStyle">',/
     .         ,'      <Icon>',/
     .         ,'        <href>root://icons/palette-3.png</href>',/
     .         ,'          <x>192</x>',/
     .         ,'          <y>64</y>',/
     .         ,'          <w>32</w>',/
     .         ,'          <h>32</h>',/
     .         ,'      </Icon>',/
     .         ,'    </IconStyle>',/
     .         ,'  </Style>',/
     .         ,'  <Folder>',/
     .         ,'    <name>My GPS Sites</name>',/
     .         ,'    <open>1</open>')   

      elseif( icall.eq.2 ) then
c       convert to geodetic coordinates 
        call XYZ_to_GEOD(rot,x,geod_pos)
        latd = 90.d0 - geod_pos(1)*180.d0/pi 
        lond = geod_pos(2)*180.d0/pi
        if( lond.gt.180.d0 ) lond = lond - 360.d0
        ht = geod_pos(3)
        write(luout,20) site,site,lond,latd,ht   
20      format('    <Placemark>',/
     .            ,'      <name>',a8,'</name>',/
     .            ,'      <description>',a8,'</description>',/
     .            ,'      <styleUrl>#GPSSiteStyle</styleUrl>',/
     .            ,'      <Point>',/
     .            ,'        <coordinates>',2(f13.8,','),f13.8
     .            ,         '</coordinates>',/
     .            ,'      </Point>',/
     .            ,'    </Placemark>')
      else
        write(luout,30)
30      format('  </Folder>',/
     .        ,'</Document>',/
     .        ,'</kml>')

      endif

      return
      end
