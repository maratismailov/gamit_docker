      Subroutine write_lfile( luout,site,x )

c     Write a GAMIT sperhical l-file   R. King 080228

      implicit none

      integer*4 luout,latd,latm,lond,lonm

      real*8 x(3),dlat,dlon,rad,lats,lons
                        
      character *1 nors,eorw 
      character* 4 siteid  
      character *8 site  
      character*12 sname
      character*60,lfmt1

      data lfmt1/     
     .'(a4,1x,a12,a1,i2,1x,i2,1x,f8.5,1x,a1,i3,1x,i2,1x,f8.5,f13.4)'/  
 
      siteid = site(1:4)
      sname = site//'        '
      call carsph(x,dlat,dlon,rad)
      call degdms(dlat,nors,latd,latm,lats) 
      if( nors.eq.'-' ) then
        nors = 'S'
      else
        nors = 'N'
      endif       
      call degdms(dlon,eorw,lond,lonm,lons)   
      if( eorw.eq.'-' ) then
        eorw = 'W'
      else
        eorw = 'E'
      endif                  
      sname(1:8) = site
      write(luout,lfmt1) site(1:4),sname,nors,latd,latm,lats
     .                   , eorw,lond,lonm,lons,rad

      return
      end





