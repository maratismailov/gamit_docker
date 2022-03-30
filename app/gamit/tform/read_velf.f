      Subroutine read_velf( luin,site,x,xdot,finished )

c     Read the next line of a GLOBK vel or org file and convert the position and 
c     velocities to Cartesian    R. King 080228

      implicit none
             
      include '../../libraries/includes/const_param.h'    
            
      integer*4 luin,nfmt,ioerr
 
      real*8 x(3),xdot(3),neudot(3),lond,latd,ht,ev,nv,hv
     .          , geod_pos(3),loc_coord(3),rot(3,3),rottr(3,3)
                
      character*8 site
      character*128 line
      
      logical finished,comment
                 
c     function
      integer*4 nblen
                
      comment = .true.            
      do while( comment ) 
        read(luin,'(a)',iostat=ioerr) line
cd          print *,'read_velf ioerr line ',ioerr,line
          if( ioerr.eq.-1 ) then
            finished = .true. 
            return
          elseif( ioerr.eq.0 ) then
            if( line(1:1).eq.' ' ) then
              comment = .false.      
              nfmt = nblen(line)
              if( nfmt.eq.109 ) then
                read(line(102:109),'(a8)',iostat=ioerr) site
              elseif( nfmt.eq.113) then
                read(line(106:113),'(a8)',iostat=ioerr) site
              elseif( nfmt.ne.0 ) then
                write(*,'(a,i4)')'Cannot identify vel-file format
     .             , nfmt=',nfmt
                stop
              endif
              read(line,*,iostat=ioerr) lond,latd,ev,nv
              ht = 0.d0  
              if( nfmt.eq.109 ) then
                read(line(76:84),*,iostat=ioerr) hv  
              elseif( nfmt.eq.113) then
                read(line(80:88),*,iostat=ioerr) hv  
              elseif( nfmt.eq.0 ) then
                continue
              endif
              if( ioerr.ne.0.and.ioerr.ne.-1 ) then
                write(*,'(a,i4)')' Error decoding H rate in vel-file '
     .             ,ioerr
                stop 
              endif
            endif
          else
            write(*,'(a,i4)') 'Error reading vel-file line ',ioerr
            stop
          endif  
       enddo    

                      
c  Convert geodetic position and NEU velocities to Cartesian   

c     arguments are colat and lon in radians, ht in meters
      geod_pos(1) = (90.d0 -latd)*pi/180.d0
      geod_pos(2) = lond*pi/180.d0
      geod_pos(3) = ht
      call GEOD_to_XYZ(geod_pos,x)
      
c  Convert E N U velocities to Cartesian 

c     rotation matrix from XYZ to NEU (use transpose)
      call XYZ_to_GEOD(rot, x, loc_coord)   
      call transp(rot,rottr,3,3)    
      neudot(1) = nv*1.d-3
      neudot(2) = ev*1.d-3
      neudot(3) = hv*1.d-3
      call mmply(rottr,neudot,xdot,3,3,1)
     
      return
      end

                 


