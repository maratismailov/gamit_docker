      Subroutine read_glistf( luin,site,x,finished )

c     Read the next line of a GLOBK vel or org file and convert the position and 
c     velocities to Cartesian    R. King 080228

      implicit none
             
      include '../../libraries/includes/const_param.h'    
            
      integer*4 luin,ioerr,nblen
 
      real*8 x(3),lond,latd,ht,geod_pos(3)
                
      character*8 site
      character*128 line
      
      logical finished


      read(luin,'(a)',iostat=ioerr) line  
        if( nblen(line).lt.66 ) then
          finished = .true. 
          return
        elseif( ioerr.eq.0 ) then
          if( line(1:1).eq.' ' ) then
            read(line(61:68),'(a8)',iostat=ioerr) site
            read(line,*,iostat=ioerr) lond,latd,ht 
            ht = ht*1.d3
            if( ioerr.ne.0 ) then
              write(*,'(a,i4)') 'Error decoding glist-file line ',ioerr
              stop   
            endif
          endif
        else
          write(*,'(a,i4)') 'Error reading glist-file line ',ioerr
          stop
        endif      

                      
c  Convert geodetic position and ENU velocities to Cartesian   

      geod_pos(1) = (90.d0 -latd)*pi/180.d0
      geod_pos(2) = lond*pi/180.d0
      geod_pos(3) = ht    
      call GEOD_to_XYZ(geod_pos,x)     
                 
      return
      end

