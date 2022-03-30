       subroutine svbody_coords( iuy,satcrd,sun,xhat,yhat,zhat
     .                       ,prn,svantbody,yatt,ievent,iepoch )

c
c** This subroutine computes the GPS - fixed coordinate system and takes
c   into account the yaw attitude of the satellite as read from the y table.
c   All vector inputs must be in a geocentric frame such that the Z axis
c   of the frame is toward the northern hemisphere (e.g., J2000.).
c  
c  IN: 
c      iuy       :  unit for yaw table (0 if no table)              I*4      
c      satcrd    :  satellite coordinates                           R*8(6)
c      sun       :  sun coordinates                                 R*8(6)
c      prn       :  prn of satellite                                I*4
c      svantbody :  antenna/body descriptor                         C*20
c      yatt      :  y-file yaw attitude of satellite                R*8      
c      iepoch    :  for debug                                       I*4 
c 
c  OUT:
c       xhat,yhat,zhat  :  unit vectors of sv axes (in sv system)
c
c  P. Tregoning
c  3 December 1997  (renamed from gps_coords to svbody_coords by R King 28 October 2014)

      implicit none

      include '../includes/dimpar.h'

      real*8 satcrd(6),sun(6),xhat(3),yhat(3),zhat(3)
     .      ,rad2deg,cosa,sina,sc_x_t(3),sc_y_t(3),amag3,dot 
     .      ,ideal_yaw,beta,yatt,sc_W(3),model_yaw,yaw_error,yaw
      integer*4 iuy,prn,ievent,i
      logical debug            
c      for debug only
      integer*4 iepoch  
      character*20 svantbody
          
c      function

      data rad2deg/57.29577951308232d0/,debug/.false./
       external amag3,dot
c
      save rad2deg,sc_X_t,sc_Y_t
c                                       
cd      if( iepoch.eq.1027 ) then
cd        debug = .true.
cd      else
cd         debug = .false.
cd      endif
                 
      cosa = 0.d0
      sina = 0.d0
      yaw_error=0.0d0                      

c If orbit-normal mode, compute the axes
      if( ievent.eq.3 ) then
        do i=1,3
           xhat(i) = satcrd(i+3) /
     .         dsqrt(satcrd(4)**2+satcrd(5)**2+satcrd(6)**2)
        enddo 
        if(debug) print *,'orbit normal yatt,xhat ',yatt,(xhat(I),I=1,3)
        return

c Otherwise, compute tne nominal axes and yaw angle
      else
        call yaw_attit( satcrd,sun,ideal_yaw,beta,sc_x_t,sc_y_t
     .                , zhat,sc_W,svantbody)  
        if( debug ) print *,'svantbody sun satcrd ideal_yaw '
     .     ,svantbody,sun,satcrd,ideal_yaw

c RWK140125: Reverse the sign for the Block IIR satellites, which have their
c            X-axis defined in the opposite sense from the others
c RWK 170505: New yaw routine orbits/kouba_gps.f removes the sign reversal
c            to be consistent with current IGS standards
c        if( svantbody(1:9).eq.'BLOCK IIR') then
c          ideal_yaw = -ideal_yaw
c       this reassignment to avoid changing yatt for second pass 
c          model_yaw = yatt  + 180.d0      
c          ideal_yaw = ideal_yaw + 180.d0     
c          do i=1,3
c            sc_x_t(i) = -sc_x_t(i)
c          enddo           
c        else
c           model_yaw = yatt
c        endif
        model_yaw = yatt 
c       get both values positive between 0 and 360 for difference test below
        if( model_yaw.gt.360.d0 ) model_yaw = model_yaw - 360.d0  
        if( model_yaw.lt.0.d0 ) model_yaw = model_yaw + 360.d0
        if( ideal_yaw.gt.360.d0) ideal_yaw = ideal_yaw - 360.d0 
        if( ideal_yaw.lt.0.d0 ) ideal_yaw = ideal_yaw + 360.d0
       
c       if no input yaw file, set modeled  yaw = ideal yaw
        if( iuy.eq.0 ) model_yaw = ideal_yaw
c       now compute the yaw error (being yatt - ideal_yaw)
        yaw_error = model_yaw - ideal_yaw
        yaw_error = yaw_error-nint(yaw_error/360.0d0)*360.0d0  
c** debug to look for floating point exception
cd      if( dabs(yatt).lt.0.1d0 .or. dabs(yatt).gt.360.d0 .or.
cd     .    dabs(ideal_yaw).lt.0.1d0 .or. dabs(ideal_yaw).gt.360.d0 ) then
cd          print *,yatt,ideal_yaw,yaw_error
cd      endif
cd      if( dabs(yaw_error).lt.1.d-12) print *,'yaw error',yaw_error
       
        if(debug) print *,'ideal_yaw,yatt,yaw_error'
     .      ,ideal_yaw,yatt,yaw_error

c** Now adjust the X and Y coordinate by the yaw error. Here we can assume that the
c   yaw angle from the yaw table (yatt) already contains the correct adjustments 
c   for any yaw manoeuvres etc etc. Therefore, irrespective of what the satellite
c   has done, if there is a difference between yatt and ideal_yaw the sc axes xhat
c   and yhat must be adjusted. 
        if(dabs(yaw_error).gt.1.0d-3)then      
          cosa=dcos(yaw_error/rad2deg)
          sina=dsin(yaw_error/rad2deg)
          do i=1,3
            xhat(i)=sc_X_t(i)*cosa+sc_Y_t(i)*sina
            yhat(i)=sc_Y_t(i)*cosa-sc_X_t(i)*sina 
          enddo                  
        if(debug) print *,'model yaw prn antbody yatt scX scY xhat yhat'
     .          ,       prn,svantbody,model_yaw,sc_X_t,sc_Y_t,xhat,yhat
         if(debug) print *,'ideal_yaw model_yaw yaw_error '
     .                   ,ideal_yaw,model_yaw,yaw_error
        else
          do i=1,3
            xhat(i) = sc_X_t(i)
            yhat(i) = sc_Y_t(i)   
          enddo 
        if(debug) print *,'ideal yaw prn antbody yatt scX scY xhat yhat'
     .        ,         prn,svantbody,ideal_yaw,sc_X_t,sc_Y_t,xhat,yhat
        endif
cd       if( iepoch.gt.2 ) stop 

c      endif if on orbit normal or nominal yaw
       endif

       return
       end

c  PT971202: the code below is the original code from Bar-Sever and computes
c            the partial derivatives - used to estimate the max yaw rate. 
c            Since we have no plans to estimate this parameter in GAMIT I have
c            just commented it all out but left it here in case someone wants
c            to do this at a later stage. Good luck to you if you do! There is
c            some more code in yaw_ctrl (part of YAWTAB now).


c** Now compute the coordinate partials if necessary.
c      if (.not. part) return
c
cc
cc** There is no dependence on velocity.
c      do i=1,3
c         do j=1,3
c            Dsc_XDv(i,j)=0.0d0
c            Dsc_YDv(i,j)=0.0d0
c            Dsc_ZDv(i,j)=0.0d0
c          enddo
c      enddo
cc
cc** Compute the partial of Z direction wrt s/c position.
cc   (Z=R/|R| ==> Zij=-Dij/|R| + Ri*Rj/|R|**3 = (-Dij+Zi*Zj)/|R| ,
cc   where R is the s/c position vector and D is Kronecker's Delta).
cc       do i=1,3
cc          Dsc_ZDr(i,i)=-(1.0d0-sc_Z(i)*sc_Z(i))/earth_sc_dist
cc          do 100 j=1,i-1
cc             Dsc_ZDr(i,j)=sc_Z(i)*sc_Z(j)/earth_sc_dist
cc             Dsc_ZDr(j,i)=Dsc_ZDr(i,j)
cc          enddo
c
c
cc** Compute the partials of the X and Y directions.
c
cc If the Earth-probe-Sun angle is very small, zero the partials and return.
c     if (EPS_angle .lt. 1.0d-6 .or.
*                                                ! singularity
c    &     EPS_angle .gt. 180.0d0-1.0d-6) then
c         do i=1,3
c            do j=1,3
c               Dsc_XDr(i,j)=0.0d0
c               Dsc_YDr(i,j)=0.0d0
c               dSC_Xdyaw(i)=0.0d0
c               dSC_Ydyaw(i)=0.0d0
c             enddo
c          enddo
c         return
c      endif
c
cc* Start with compute the partials of the X and Y directions before yaw
cc  adjustment, i.e., of sc_X_t and sc_Y_t.
cc
cc Compute the sc_Y_t partial.
c      a=1.0d0/(sin(EPS_angle/rad2deg)*earth_sc_dist*sc_sun_dist)
c      do i=1,3
c         ymat(i,i)=(1.0d0-sc_Y_t(i)*sc_Y_t(i))*a
c         do j=1,i-1
c            ymat(i,j)=-sc_Y_t(i)*sc_Y_t(j)*a
c            ymat(j,i)=ymat(i,j)
c         enddo
c       enddo
c      do i=1,3
c         Dsc_Y_tDr(i,1)=ymat(i,2)*sun(3)-ymat(i,3)*sun(2)
c         Dsc_Y_tDr(i,2)=ymat(i,3)*sun(1)-ymat(i,1)*sun(3)
c         Dsc_Y_tDr(i,3)=ymat(i,1)*sun(2)-ymat(i,2)*sun(1)
c      enddo
cc
cc Compute the sc_X_t partial.
c      do j=1,3
c         call cross(sc_Y_t,Dsc_ZDr(1,j),Dsc_X_tDr(1,j))
c         call cross(Dsc_Y_tDr(1,j),sc_Z,vec)
c         do i=1,3
c            Dsc_X_tDr(i,j)=Dsc_X_tDr(i,j)+vec(i)
c         enddo
c      enddo
cc
cc* Now account for the yaw adjustmet.
c      do i=1,3
c         do j=1,3
c            Dsc_XDr(i,j)=Dsc_X_tDr(i,j)*cosa+Dsc_Y_tDr(i,j)*sina
c            Dsc_YDr(i,j)=Dsc_Y_tDr(i,j)*cosa-Dsc_X_tDr(i,j)*sina
c          enddo
c       enddo
cc
cc** Now compute the partials of sc_X and sc_Y wrt the max_yaw_rate parameter.
c      dsina=cosa*Dyaw_angleDp/rad2deg
c      dcosa=-sina*Dyaw_angleDp/rad2deg
c      do i=1,3
c         Dsc_XDp(i)=sc_X_t(i)*dcosa+sc_Y_t(i)*dsina
c         Dsc_YDp(i)=sc_Y_t(i)*dcosa-sc_X_t(i)*dsina
c      enddo

cc PT 950508: now compute d( inert coord)/d(yaw rate)
c      call matmpy(Dsc_xDr,Dsc_xDp,dSC_Xdyaw,3,3,1)
c      call matmpy(Dsc_yDr,Dsc_yDp,dSC_Ydyaw,3,3,1)
c
cc
cc** Bye.
