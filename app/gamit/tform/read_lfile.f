      Subroutine read_lfile( luin,site,x,date,values_ok,finished )

c     Read the next line of a old-sylte GAMIT l-file and convert the position and 
c     velocities to Cartesian    R. King 080228

      implicit none
             
      include '../../libraries/includes/const_param.h'    
                                 
      integer*4 luin,ioerr,latd,latm,lond,lonm,nblen,i
 
      real*8 x(3),rad,dlat,dlon,rlat,rlon,lats,lons,date
                
      character*1 latflag,lonflag
      character*8 site  
      character*12 sname 
      character*60 lfmt1 
      character*128 line
      
      logical values_ok,finished
          
      data lfmt1/     
     .'(a4,1x,a12,a1,i2,1x,i2,1x,f8.5,1x,a1,i3,1x,i2,1x,f8.5,f13.4)'/  

      date = 0.d0                 
      do i=1,3
        x(i) = 0.d0
      enddo             
                              
      read(luin,'(a)',iostat=ioerr) line   
      if( ioerr.eq.-1 .or. nblen(line).le.0 ) then
        finished = .true.    
      elseif( ioerr.eq.0 ) then
        read( line,lfmt1,iostat=ioerr) site(1:4),sname
     .    , latflag,latd,latm,lats,lonflag,lond,lonm,lons,rad     
        if(ioerr.ne.0 ) then
           write(*,'(a,i4)') 'Error decoding L-file values ',ioerr
           write(*,'(2a)') 'LINE: ',line  
           values_ok = .false.
        else
          read(line(75:83),'(f9.4)',iostat=ioerr) date
c         take no action on missing date, just leave at zero
        endif    
        site(5:8) = '_GPS'
        call dmsdeg(latflag,latd,latm,lats,dlat)   
        call dmsdeg(lonflag,lond,lonm,lons,dlon)   
        rlat = dlat*pi/180.d0
        rlon = dlon*pi/180.d0
        x(1) = rad*dcos(rlat)*dcos(rlon) 
        x(2) = rad*dcos(rlat)*dsin(rlon)
        x(3) = rad*dsin(rlat)  
c       check for reasonableness
        rad = dsqrt(x(1)**2 + x(2)**2 + x(3)**2)
        if (rad.gt.6.d6 .and. rad.lt.7.d6 ) then 
          values_ok = .true.
        else
           values_ok = .false.
        endif
      else
        write(*,'(a,i4)') 'Error reading l-file line ',ioerr
        stop
      endif      

      return
      end

                 


