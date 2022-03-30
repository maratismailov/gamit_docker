      Subroutine write_velf( luout,site,x,xdot )      

c     Write a line of a GLOBK velocity file. R. King 080228


      implicit none
      
      include '../../libraries/includes/const_param.h'    

      integer*4 luout

      real*8 x(3),xdot(3),lond,latd,ev,nv,hv,geod_pos(3)
     .     , neudot(3),rot(3,3)
      real*4 zero,one
                
      character*8 site
      
      data zero/0./,one/1.0/
               
c  Convert Cartesian coordinates to geodetic position, and velocities to ENU    
        
      call XYZ_to_GEOD(rot,x,geod_pos)
      latd = 90.d0 - geod_pos(1)*180.d0/pi 
      lond = geod_pos(2)*180.d0/pi

c  Convert E N U velocities to Cartesian 

      call mmply(rot,xdot,neudot,3,3,1)
      nv = neudot(1)*1.d3
      ev = neudot(2)*1.d3  
      hv = neudot(3)*1.d3

      write(luout,'
     .  (1x,f10.5,1x,f10.5,1x,4f8.2,2f8.2,f7.3,2x,f8.2,2f8.2,1x,a8)')
     .       lond,latd,ev,nv,zero,zero,one,one,zero,hv,zero,one,site
     
      return
      end

               
