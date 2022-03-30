      program interp_log_tab
c
c     Program to interpolate a postseismic log values obtained from tsview.
c
c Input Variables
c ---------------
c file_name -- name of input file containing log values.
c interp_int -- interpolation interval (days).
c interp_start_date -- interpolation start date (YYYY MM DD).
c interp_stop_date -- interpolation stop date (YYYY MM DD).
c
c Working Variables
c------------------
c
c site -- 8 char site name
c lat lon -- latitude and longitude of the site (decimal degrees)
c v_n, v_e -- north and east linear velocity of site (mm/yr)
c a_n, a_e, a_n_sig a_e_sig -- north and east log function amplitude and sigma (mm)
c tau -- log function time constant
c wrms_e, wrms_n -- wrms of the log function fit to the time series data (mm)
c data_start, data_stop -- start and stop times of time series data for site (yyyy mm dd hh mm) 
c log_date -- log function strating epoch (yyyy mm dd hh mm)  
 
c Output Variables
c-----------------
c
c n_pos e_pos
c n_vel e_vel
c n_pos_sig e_pos_sig
c n_vel_sig e_vel_sig
c
      implicit none

c     include '../../kf/includes/const_param.h'
 
      real*8 lat, lon, v_e, v_n, tau(3),    
     .    a_n(3), a_n_sig(3), a_e(3), a_e_sig(3), 
     .    wrms_e, wrms_n, values(30),dpy

      real*8 interp_start_jd, interp_stop_jd, data_start_jd, 
     .       data_stop_jd, log_date_jd, jd, sec, delta_t, n_pos(5),
     .       n_vel(5), e_pos(5), e_vel(5), n_pos_sig, e_pos_sig,
     .       n_vel_sig, e_vel_sig,data_start(5),data_stop(5)
c
      integer*4 i,j,indx,jerr,ierr,iel,trimlen,ioerr,date(5),ihr,imin,
     .          interp_start(5), interp_stop(5), interp_int,
     .          log_date(5),iarg,iclarg
c     .          data_start(5),data_stop(5)
c
      character*8  site_name, arg
      character*80 file_name
      character*200 buffer
c
c initialize some defaults
      interp_start(1) = 1999 
      interp_start(2) = 08      
      interp_start(3) = 18       
      interp_start(4) = 0
      interp_start(5) = 0
      interp_stop(1)  = 2006
      interp_stop(2)  = 04
      interp_stop(3)  = 01
      interp_stop(4)  = 0
      interp_stop(5)  = 0
      do i = 1,5
        data_start(i) = 0.d0
        data_stop(i)  = 0.d0
        log_date(i)   = 0
      enddo
      jerr = 0 
      dpy = 365.25d0

c read command line input

c Get the run-string arguments

      iarg = iclarg(1,file_name)
      open(unit=1,file=file_name,status='old',iostat=ioerr) 
      if( ioerr.ne.0 ) then
        print *,'Error opening input file ',file_name
        if ( file_name .eq." " ) then
          call help_file 
          stop
        endif
        print *,'Error opening input file ',file_name,ioerr
      else
        print *,'Opened inputfile ',file_name
      endif

      iarg = iclarg(2,arg)
      read(arg,'(i4)') interp_int 

      iarg = iclarg(3,arg)
      read(arg,'(i4)') interp_start(1)

      iarg = iclarg(4,arg)
      read(arg,'(i4)') interp_start(2) 

      iarg = iclarg(5,arg)
      read(arg,'(i4)') interp_start(3) 

      iarg = iclarg(6,arg)
      read(arg,'(i4)') interp_stop(1)

      iarg = iclarg(7,arg)
      read(arg,'(i4)') interp_stop(2) 

      iarg = iclarg(8,arg)
      read(arg,'(i4)') interp_stop(3) 

c Print the input to the screen

      print *,'Interp Int. Start YYYY MO DD HR MN Stop YYYY MO DD HR MN'
      print *,interp_int,interp_start,interp_stop

c      call rcpar(1,file_name)
c      call rcpar(2,interp_int)
c      call rcpar(3,arg) 
c     interp_start(1) = arg
c      call rcpar(4,arg) 
c      interp_start(2) = arg
c      call rcpar(5,arg) 
c      interp_start(3) = arg
c      call rcpar(6,arg) 
c      interp_stop(1) = arg
c      call rcpar(7,arg) 
c      interp_stop(2) = arg
c      call rcpar(8,arg) 
c      interp_stop(3) = arg

c     convert the interpolation start and stop year, month,day, hr, min, sec into JD  
      sec=0.d0
      call ymdhms_to_jd(interp_start, sec, interp_start_jd)
      call ymdhms_to_jd(interp_stop, sec, interp_stop_jd)
c check if user knows how to use program.
      if ( file_name .eq." " ) then
        call help_file 
        stop
      endif

c      open(unit=1,file=file_name,status='old',iostat=ioerr) 
c      if( ioerr.ne.0 ) then
c        print *,'Error opening input file: ',file_name,ioerr
c        stop
c      else
c        print *,'Opened input file: ',file_name
c      endif

c loop over input data
10    do while ( jerr.eq.0 )
         read(1,'(a)', iostat=jerr ) buffer
 
*        See if we can find the site name
         indx = 1
         if( trimlen(buffer).gt.0 .and. buffer(1:1).eq.' ' ) then
            call getword( buffer, site_name, indx)
            iel = 1
         else
            iel = 0
         end if
 
*       Decode the line and add new information if found
         if( iel.gt.0 ) then
            call multiread( buffer, indx, 'R8', ierr, values,
     .                       ' ',26)
            lon =           values(1)
            lat =           values(2)
            v_e =           values(3)
            v_n =           values(4)
            log_date(1)  =  values(5)
            log_date(2)  =  values(6)
            log_date(3)  =  values(7)
            tau(1) =        values(8)
            a_e(1) =        values(9)
            a_e_sig(1) =    values(10)
            a_n(1) =        values(11)
            a_n_sig(1) =    values(12)
            tau(2) =        values(13)
            a_e(2) =        values(14)
            a_e_sig(2) =    values(15)
            a_n(2) =        values(16)
            a_n_sig(2) =    values(17)
            tau(3) =        values(18)
            a_e(3) =        values(19)
            a_e_sig(3) =    values(20)
            a_n(3) =        values(21)
            a_n_sig(3) =    values(22)
            wrms_e =        values(23)
            wrms_n =        values(24)           
            data_start(1) = values(25)
            data_stop(1) =  values(26)
c            data_start(2) = values(26)
c            data_start(3) = values(27)
c            data_stop(2) =  values(29)
c            data_stop(3) =  values(30)
c
            write(*,*) " "
            write(6,15)site_name,lat,lon,v_e,v_n,log_date,tau(1),a_e(1),
     .      a_e_sig(1),a_n(1),a_n_sig(1),tau(2),a_e(2),a_e_sig(2),
     .      a_n(2),a_n_sig(2),tau(3),a_e(3),a_e_sig(3),a_n(3),
     .      a_n_sig(3),wrms_e,wrms_n,data_start(1),data_stop(1)

15          format('Site: ',a8,1x,2(f8.3,1x),2(f7.2,1x),i4,1x,4(i2,1x),/
     .      ,f5.0,1x,f6.1,1x,f6.1,1x,f6.1,1x,f6.1,/
     .      ,f5.0,1x,f6.1,1x,f6.1,1x,f6.1,1x,f6.1,/
     .      ,f5.0,1x,f6.1,1x,f6.1,1x,f6.1,1x,f6.1,/
     .      ,f4.1,1x,f4.1,/
     .      ,f8.3,1x,f8.3)
c     .      ,i4,1x,4(i2,1x),i4,1x,4(i2,1x))
            write(*,*) " "
 
c     convert the data start and stop year, month, day, hr, min, sec into JD  
c            data_start(4) = 0; data_start(5) = 0
c            call ymdhms_to_jd(data_start, sec, data_start_jd)
            call decyrs_to_jd(data_start(1),data_start_jd)
c            data_stop(4) = 0; data_stop(5) = 0
c            call ymdhms_to_jd(data_stop, sec, data_stop_jd)
            call decyrs_to_jd(data_stop(1),data_stop_jd)
  
c     convert the log function start  year, month, day, hr, min, sec into JD
c            log_date(4) = 0; log_date(5) = 0
            call ymdhms_to_jd(log_date, sec, log_date_jd)
c
c Loop thru interpolation time range stepping by the interpolation interval
c Compute site position and velocity @ specified time points.

            jd = interp_start_jd
            write(6,20)
20          format(' Site      Lon       Lat     Date          Delta T',
     .      '      North T1       North T2    North T3   North All ',
     .      'North Total North Sig     East T1     East T2     ',
     .      'East T3    East All   East Total  East Sig')
            write(6,25)
25          format('           deg      deg    YYYY MM DD HR MN   days',
     .      '           mm      ',
     .      '    mm         mm           mm         mm         mm    ',
     .      '        mm         mm        mm            mm      ',
     .      '    mm        mm')
30          do while ( jd .le. interp_stop_jd ) 
               delta_t = jd - log_date_jd 
c               print*,'Delta T: ',delta_t, ' days'    
          
c Compute site position and velocity if interpolation time > data start time and > log start date
               if ( jd .ge. log_date_jd .and. jd .ge. data_start_jd 
     .             .and. jd .le. data_stop_jd ) then
                  n_pos(1) =  a_n(1)*log(1+(delta_t/tau(1)))
                  n_pos(2) =  a_n(2)*log(1+(delta_t/tau(2)))
                  n_pos(3) =  a_n(3)*log(1+(delta_t/tau(3)))
                  n_pos(4) = n_pos(1) + n_pos(2) + n_pos(3)
                  n_pos_sig =  a_n_sig(1)*log(1+(delta_t/tau(1)))
     .                      +  a_n_sig(2)*log(1+(delta_t/tau(2)))  
     .                      +  a_n_sig(3)*log(1+(delta_t/tau(3)))
                  n_pos(5) = delta_t/dpy * v_n + n_pos(4)
                  n_vel(1) = a_n(1)/tau(1)*dpy/(1+(delta_t/tau(1)))
                  n_vel(2) = a_n(2)/tau(2)*dpy/(1+(delta_t/tau(2))) 
                  n_vel(3) = a_n(3)/tau(3)*dpy/(1+(delta_t/tau(3)))
                  n_vel(4) = n_vel(1) + n_vel(2) + n_vel(3)
                  n_vel_sig = a_n_sig(1)/tau(1)*dpy/(1+(delta_t/tau(1)))
     .                     +  a_n_sig(2)/tau(2)*dpy/(1+(delta_t/tau(2))) 
     .                     +  a_n_sig(3)/tau(3)*dpy/(1+(delta_t/tau(3)))
                  n_vel(5) = v_n + n_vel(4)
                  e_pos(1) =  a_e(1)*log(1+(delta_t/tau(1)))
                  e_pos(2) =  a_e(2)*log(1+(delta_t/tau(2)))
                  e_pos(3) =  a_e(3)*log(1+(delta_t/tau(3)))
                  e_pos(4) = e_pos(1) + e_pos(2) + e_pos(3)
                  e_pos_sig =  a_e_sig(1)*log(1+(delta_t/tau(1)))
     .                      +  a_e_sig(2)*log(1+(delta_t/tau(2)))
     .                      +  a_e_sig(3)*log(1+(delta_t/tau(3)))
                  e_pos(5) = delta_t/dpy * v_e + e_pos(4)
                  e_vel(1) = a_e(1)/tau(1)*dpy/(1+(delta_t/tau(1)))
                  e_vel(2) = a_e(2)/tau(2)*dpy/(1+(delta_t/tau(2))) 
                  e_vel(3) = a_e(3)/tau(3)*dpy/(1+(delta_t/tau(3)))
                  e_vel(4) = e_vel(1) + e_vel(2) + e_vel(3)
                  e_vel_sig = a_e_sig(1)/tau(1)*dpy/(1+(delta_t/tau(1)))
     .                     +  a_e_sig(2)/tau(2)*dpy/(1+(delta_t/tau(2))) 
     .                     +  a_e_sig(3)/tau(3)*dpy/(1+(delta_t/tau(3)))
                  e_vel(5) = v_e + e_vel(4)

                  call jd_to_ymdhms(jd,date,sec)

                  write(6,85)site_name,lon,lat,date,delta_t,n_pos(1),
     .            n_pos(2),n_pos(3),n_pos(4),n_pos(5),n_pos_sig,
     .            e_pos(1),e_pos(2),e_pos(3),e_pos(4),e_pos(5),e_pos_sig
85                format (1x,a8,1x,f7.3,1x,f8.3,1x,i4,1x,4(i2,1x),1x,
     .            f5.0,1x,'|',1x,2(f11.3,1x,f11.3,1x,f11.3,1x,f11.3,
     .            1x,f11.3,1x,f9.3,1x,'|'),1x,'POS')               
                  write(6,95)site_name,lon,lat,date,delta_t,n_vel(1),
     .            n_vel(2),n_vel(3),n_vel(4),n_vel(5),n_vel_sig,
     .            e_vel(1),e_vel(2),e_vel(3),e_vel(4),e_vel(5),e_vel_sig
95                format (1x,a8,1x,f7.3,1x,f8.3,1x,i4,1x,4(i2,1x),1x,
     .            f5.0,1x,'|',1x,2(f11.3,1x,f11.3,1x,f11.3,1x,f11.3,
     .            1x,f11.3,1x,f9.3,1x,'|'),1x,'VEL')
cd                print *,site_name, date, n_pos, n_pos_sig, e_pos, 
cd     .          e_pos_sig, n_vel, n_vel_sig, e_vel, e_vel_sig

               endif

               jd = jd + interp_int

            enddo
         else
            goto 10
         endif
      enddo
      stop
      end  

      subroutine help_file
      write(*,*) " "
      write(*,*) "HELP for interp_log_tab"
      write(*,*) " "
      write(*,*) "usage: interp_log_tab <File> <Inter> <Start> <Stop>"
      write(*,*) " "
      write(*,*) "where:"
      write(*,*) "  File: Name of the file containing log parameters"
      write(*,*) "  Inter: Interpolation interval in days"              
      write(*,*) "  Start: YYYY MM DD to start interpolation"
      write(*,*) "  Stop: YYYY MM DD to stop interpolation"                                                 
      write(*,*) " "                                   
      write(*,*) "Note: Interpolation will not be performed outside the"                                                 
      write(*,*) "the valid start/stop data range given in the input"
      write(*,*) "file."
      write(*,*) " "
      write(*,*) "Example: "
      write(*,*) "interp_log_tab log_file.inp 7 1999 08 17 2000 01 01"
      write(*,*) " "
      return
      end
