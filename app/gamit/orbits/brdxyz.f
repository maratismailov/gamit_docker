Copyright (c) Massachusetts Institute of Technology, 1987. All rights reserved.

       subroutine brdxyz( iwkn,sow,iwknbc,bephem,satprm,ibcerr,nprn )
c
c     Compute the satellite earth-fixed position (X,Y,Z) from
c     the broadcast ephemeris
c
c     Written by P. Morgan, February, 1987 based on a routine of
c     same name by R. King
c
c     Modifications by K. Feigl (Jun 87), Y. Bock (Nov 87), and R. King (Oct 95)
c     to pass times and arrays better.

c     Input:
c        iwkn   : week number at requested epoch
c        sow    : seconds-of-week at requested epoch
c        iwknbc : week number of input BC message
c        bephem : elements of input BC message (bephem(1) = seconds-of-week ref epoch)

c     Output:
c        satprm : position and velocity at requested epoch
c        ibcerr : error flag (0 if ok, -1 if invalid ephemeris)

c     The meanings of the transfered parameters can be found in
c     the print statments.

      implicit none
                     
      integer*4 iwkn,sow,iwknbc,ibcerr,iter,nprn,len,rcpar,i
                          
      character*80 prog_name
      character*256 message

      real*8
     1       bephem(16),satprm(6)
     2     , xetoe,xem0,xedn,xeart,xeecc
     3     , xei0,xeidt,xeom0,xew,xecuc,xecus,xecrc,xecrs
     4     , xecic,xecis,xeomd
     5     , gm,erate,t,xn,a
     6     , anom,e1,e2,enom,cosv,sinv,v,phi,xinc
     7     , di,du,dr,u,r,xp,yp,asc,x,y,z,xdot,ydot,zdot
     8     , term,xpdot,ypdot,asctrm,secdif
      
      logical debug /.false./

C GM in units of km**3/sec**2
cc     rwk 160804: changed value of gm to EGM08 and change km to meters
cc    data gm/398600.50d+00/,erate/7.292115147d-05/                   
      data gm/398600.4415d+09/,erate/7.292115147d-05/                   
              

c     get calling program name and m-file name for report_stat
      len = rcpar(0,prog_name)


      ibcerr = 0
     
      if( debug)  print *,'BRDXYZ bephem ',bephem
                 

c     unpack the arrays
c     ephemeris parameters (distance units all in meters)
      xetoe  =bephem(1)
      xem0   =bephem(2)
      xedn   =bephem(3)
      xeart = bephem(4) 
      xeecc  =bephem(5)
      xei0   =bephem(6)
      xeidt  =bephem(7)
      xeom0  =bephem(8)
      xeomd  =bephem(9)
      xew    =bephem(10)
      xecuc  =bephem(11)
      xecus  =bephem(12)
      xecrc  =bephem(13)
      xecrs  =bephem(14)
      xecic  =bephem(15)
      xecis  =bephem(16)

c     compute the time from reference epoch

      t = secdif(iwkn,sow,iwknbc,xetoe)

c     compute the semi-major axis and check for reasonableness

      a = xeart*xeart        
      if (dabs(a) .lt. 1.d7 .or. dabs(a) .gt. 1.d9) then
          write (message,'(a,i3,a,f7.0,a,d12.4)')
     .       'Bad BC ephemeris record   prn=',nprn,' sow='
     .       ,bephem(1),' a=',a
           call report_stat('FATAL',prog_name,'orbits/brdxyz'
     .                      ,' ',message,0)
         ibcerr = 1
         return
      endif

c     compute the mean motion

      xn = dsqrt(gm/(a*a*a)) +xedn

c     debug
c

c     write(6,401) t,xetoe,xem0,xn,a,
c    1xeecc,xei0,xeidt,xeom0,xeomd,xew
c     write(10,401) t,xetoe,xem0,xn,a,
c    1xeecc,xei0,xeidt,xeom0,xeomd,xew
c401   format(/,1x,'rewrite of basic quantities',//,1x,
c     2't = ttag - xetoe , in seconds, is ',d22.15,/,1x,
c     3'xetoe the epoch time for the ephemeris is ',d22.15,/,1x,
c     4'xem0, the mean anomaly is ',d22.15,/,1x,
c     5'xn is ',d22.15,/,1x, 'a the semi-major axis is ',d22.15,/,1x,
c     6'xeec the orbit eccentricity is ',d22.15,/,1x,
c     7'xeio the orbit inclination is ',d22.15,/,1x,
c     8'xeidt the rate of change of orbit inclination is ',d22.15,/,1x,
c     9'xeom0 the longitude of the ascending node is ',d22.15,/,1x,
c     1'xeomd the rate of change of the ascending node is ',d22.15,/,1x,
c     2'xew the argument of perigee is ',d22.15)
c     write(6,402) xecuc,xecus,xecrc,xecrs,xecic,xecis
c     write(10,402) xecuc,xecus,xecrc,xecrs,xecic,xecis
c 402  format(/,1x,'xecuc the cuc polynomial parameter is ',d22.15,/,1x,
c     1'xecus the cus polynomial parameter is ',d22.15,/,1x,
c     2'xecrc the crc polynomial parameter is ',d22.15,/,1x,
c     3'xecrs the crs polynomial parameter is ',d22.15,/,1x,
c     4'xecic the cic polynomial parameter is ',d22.15,/,1x,
c     5'xecis the cis polynomial parameter is ',d22.15,/,1x)

c     compute the mean anomaly and solve Kepler's eqn
c     (successive approx. ok for low eccentricity)

      anom= xem0 + xn*t
      e1= anom
      iter= 0
  410 e2= anom + xeecc*dsin(e1)
      if( dabs(e2-e1).lt.1.d-10 ) goto 415
      e1= e2
      iter = iter + 1
      if( iter.gt.7 ) goto 997
      goto 410
  415 enom= e2

c     compute the angular elements

      cosv= ( dcos(enom)- xeecc)/(1.d0 - xeecc*dcos(enom) )
      sinv= dsqrt(1.d0 - xeecc**2)*dsin(enom)
     1                        / (1.d0 - xeecc*dcos(enom))
      v= datan2(sinv,cosv)
      phi= v + xew
      du = xecus*dsin(2.d0*phi) + xecuc*dcos(2.d0*phi)
      dr = xecrc*dcos(2.d0*phi) + xecrs*dsin(2.d0*phi)
      di = xecic*dcos(2.d0*phi) + xecis*dsin(2.d0*phi)
      u  = phi + du
      r  = a*(1.d0 - xeecc*dcos(enom)) + dr
      xinc= xei0 + di + xeidt*t
      xp= r*dcos(u)
      yp= r*dsin(u)
      asctrm=xeomd-erate
      asc= xeom0 + asctrm*t - erate*xetoe
c
c    write(6,888) t,xn,anom,e1,enom,cosv,sinv,v
c    1          ,  phi,du,dr,di,u,xinc,asc
c     write(10,888) t,xn,anom,e1,enom,cosv,sinv,v
c    1          ,  phi,du,dr,di,u,xinc,asc
c 888  format (/,1x,'dt the time difference is ',d22.15,/,1x,
c     1'xn is ',d22.15,/,1x,
c     2'anom is ',d22.15,/,1x,
c     3'  e1 is ',d22.15,/,1x,
c     4'enom is ',d22.15,/,1x,
c     5'cosv is ',d22.15,/,1x,
c     6'sinv is ',d22.15,/,1x,
c     7'   v is ',d22.15,/,1x,
c     8' phi is ',d22.15,/,1x,
c     9'  du is ',d22.15,/,1x,
c     1'  dr is ',d22.15,/,1x,
c     2'  di is ',d22.15,/,1x,
c     3'   u is ',d22.15,/,1x,
c     4'xinc is ',d22.15,/,1x,
c     5' asc is ',d22.15)

c   compute the coordinates of the satellite in space
c   rectangular coordinates and convert to km for the t-file

      x= xp*dcos(asc) - yp*dcos(xinc)*dsin(asc)
      y= xp*dsin(asc) + yp*dcos(xinc)*dcos(asc)
      z= yp*dsin(xinc)
      term=(xn*a)/dsqrt(1.d0-xeecc*xeecc)
      xpdot=-dsin(u)*term
      ypdot=(xeecc+dcos(u))*term
      xdot= xpdot*dcos(asc) - ypdot*dcos(xinc)*dsin(asc)
     .     - xp*dsin(asc)*asctrm - yp*dcos(xinc)*dcos(asc)*asctrm
      ydot= xpdot*dsin(asc) + ypdot*dcos(xinc)*dcos(asc)
     .     + xp*dcos(asc)*asctrm - yp*dcos(xinc)*dsin(asc)*asctrm
      zdot= ypdot*dsin(xinc)
cd      write(6,900) nprn,iwkn,sow,iwknbc,xetoe,t,x,y,z,xdot,ydot,zdot
cd      write(10,900) nprn,iwkn,sow,iwknbc,xetoe,t,x,y,z,xdot,ydot,zdot
cd 900  format (
cd     . /,1x,'PRN#, Week number, sow, Week number (be), xetoe, delt ',
cd     .        i3,1x,2(i4,f8.0),f8.0,
cd     1   1x,'Satellite position (m) & velocity (m/sec):',/,
cd     2   5x,'x    = ',d22.15,
cd     3   5x,'y    = ',d22.15,
cd     4   5x,'z    = ',d22.15,/,
cd     5   5x,'xdot = ',d22.15,
cd     6   5x,'ydot = ',d22.15,
cd     7   5x,'zdot = ',d22.15)

c     put the the orbital parameters into array for easy passing
      satprm(1)=x
      satprm(2)=y
      satprm(3)=z
      satprm(4)=xdot
      satprm(5)=ydot
      satprm(6)=zdot           
      do i=1,6
        satprm(i) = satprm(i)/1.d3
      enddo

      return

 997  write(message,998)iter
 998  format('Error, the iterative solution of Keplers equations does',
     .' not occur in (',2i2,') iterations. Review input nav file')
      call report_stat('FATAL',prog_name,'orbits/brdxyz',' ',message,0)

      end
