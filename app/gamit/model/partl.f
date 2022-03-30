      Subroutine PARTL( ichan,rhat,drate,svec,freqtr,norbpart
     .                , npart,sitepart,atpart,polepart,ut1part,svantpart
     .                , tmpart,iepoch )


C     Calculate the partial derivatives of phase with respect to
C     site coordinates, and clock, atmospheric, and orbital parameters.
C     Written by R.King 14 March 1987 from code in old COMPAR and OUTPUT.
C
C     Modified by Simon McClusky to calculate earth orientation parameter partials.
C     April 1994.  
c
c     Modified by R. King to add partials wrt SV antenna offsets.  September 1998

      implicit none

      include '../includes/dimpar.h'

      integer*4 norbpart,npart,ichan,i,j,k 
      integer*4 iepoch

      real*8 parts(3),sitepart(3,3),rhat(3),tmpart(maxlab,maxsat)
     .     , svec(maxyt2),freqtr,freqr1,atpart,vlight,drate
     .     , polepart(3,4),ut1part(3,2),svantpart(3,3)


      logical debug 

      data vlight/299792.458d0/,debug/.false./
                            
c Input:  
c
c   ichan          channel on X-file
c   rhat(3)        unit vector from station in direction of satellite
c   drate          delay rate (s/s)
c   svec           vector of position(3), velocity(3), and partials of position 
c                   (3*(6+norbpart)) for satellite from T-file
c   freqtr         frequency of L1 transmission   
c   norbpart           number of non-gravitational parameters for satellites (3, 6, or 9)
c   npart          total number of partials for C-file 
c   sitepart(3,3)  partials of inertial cartesian station position wrt lat,long,rad
c   atpart         partial of delay (s) wrt atmospheric zenith delay (s)
c   polepart(3,4)  partials of inertial station position wrt pole position (x,y) and rates
c   ut1part(3,2)   partials of inertial station position wrt UT1 and rate
c   svantpart(3,3) partials of inertial satellite position wrt SV antenna offsets (S/C x,y,z)

c Output:
c
c   tmpart(maxlab,maxsat) matrix of partials of delay wrt model parameters 
c                  i= parameter, j= channel (satellite) 
c       i= 1       station latitude      (radians)
c          2       station longitude     (radians)
c          3       station radius        (km)
c          4       station zenith delay  (m)
c          5       station clock         (s)
c          6-11    satellite initial position and velocity (km, km/s)
c          next norbpart  satellite non-gravitational force parameters (dimensionless)
c          next 3     satellite antenna offsets (x,y,z)
c          next 6     Earth orientation parameters (xpole, ypole, UT1, and their rates)
c
c          Case of BERNE (9-parameter) radiation pressure model
c             12-20  non-gravitational force parameters
c             21-23  SV antenna offsets
c             24-29  EOPs
                         
c Actual cases (necessary to delineate until norbpart added to C-file header): 
c   Station partials only:                        npart =  5
c   Station + orbit + EOP with norbpart=3 :           npart = 20
c                              norbpart=6 :           npart = 23
c                              norbpart=9 :           npart = 26 
c   Station + orbit + SVant + EOP with norbpart = 6 : npart = 26
c                                      norbpart = 9 : npart = 29
 
c This routine uses the partials of Cartesian coordinates (station or satellite) wrt 
c the parameters and computes the partial of phase wrt parameters using the chain rule

c
c     d(tau)           d(tau)     d(X)
c     ------------  =  ------  *  ------------
c     d(parameter)     d(X)       d(parameter)                                   
c
c  to compute the partial of delay wrt parameters and then converting to L1
c  phase by multiplying by the nominal L1 frequency, and where d(tau)/d(X) is 
c  the unit vector from the station to the satellite (Xsat - Xsite)  and 
c  d(X)/d(parameter) is d(Xsat)/d(parameter) or -d(Xsite)/d(parameter) 
c  depending on whether the parameter is a satellite or station parameter.  
c  The partials wrt to clock epoch or zenith delay can be computer directly 
c  from the delay rate and zenith delay without need of the chain rule.

c       Station Coordinates

      do i=1,3
        parts(i)= 0.D0
      enddo
      do  i=1,3
        do j=1,3
          parts(i) = parts(i) - sitepart(j,i) * rhat(j)
        enddo
      enddo
      do i=1,3
        tmpart(i,ichan)= ( parts(i) / vlight ) * freqtr
      enddo

c       Atmospheric Delay

c     units are meters (since 930914--cm prior to that)
      tmpart(4,ichan) = freqtr * atpart / vlight / 1.0D+03


c       Receiver Clock Epoch (rate and acceleration no longer supported)   

      freqr1=freqtr*(1.d0-drate)
      tmpart(5,ichan) = freqr1

c       Satellite Orbital Parameters

      if (npart.eq.5) goto 999
c      If npart = 5, then no satellite or Earth orientation partials partials


c       Integrated orbital partials
          
      k=0
      do i = 6, 6+norbpart-1
c     do i = 6, 15
        tmpart(i,ichan)= 0.d0
        do j=1,3
          k=k+1
          tmpart(i,ichan) = 
     .       tmpart(i,ichan)+(svec(6+k)/vlight)*rhat(j)*freqtr
        enddo 
        if( debug ) then 
           if( iepoch.eq.2460 )  
     .       print *,'Chn ',ichan,' Orb part ',i,tmpart(i,ichan)
        endif
      enddo


c       Satellite antenna offsets
       
      k = 0
      do i = 6+norbpart, 6+norbpart+2
        tmpart(i,ichan) = 0.d0  
        k=k+1       
        do j=1,3  
c          if( iepoch.eq.51 .and. ichan.eq. 3 ) then        
c            print *,'ichan i j k svantpart rhat '
c     .             , ichan,i,j,k,svantpart(j,k),rhat(j)
c          endif
c         change km to m
          tmpart(i,ichan) =  tmpart(i,ichan) + 
     .          svantpart(j,k)*1.d-03/vlight*rhat(j)*freqtr
        enddo 
c          print *,'i tmpart ',tmpart(i,ichan)
      enddo     
c      if( iepoch.eq.51 .and. ichan.eq. 3 ) stop
     
 

c       Earth orientation parameter partials with respect to phase.

c**NOTE rwk980905:  Pole partials are added but ut1 is subtracted.  I cannot
c                   see a reason for the sign change.  I would have expected both
c                   to be subtracted (consistent with site vector being subtracted).
c                   These check out in a solution, however, so there must be a sign
c                   change in the computations in eopart.
               
c     pole partials
      k = 0               
      do i = 6+norbpart+3, 6+norbpart+6
c*      do i = 6+norbpart, 6+norbpart+5
        k= k+1                    
        tmpart(i,ichan) = 0.d0 
        do j=1,3
          tmpart(i,ichan) = tmpart(i,ichan) + 
     .                  polepart(j,k)/vlight * rhat(j) * freqtr
        enddo 
      enddo
c     ut1-tai partials  
c     pole partials
      k = 0               
      do i = 6+norbpart+7, 6+norbpart+8
c*      do i = 6+norbpart+4, 6+norbpart+5
        k= k+1 
        tmpart(i,ichan) = 0.d0 
        do j=1,3
          tmpart(i,ichan) = tmpart(i,ichan) -
     .                   ut1part(j,k)/vlight * rhat(j) * freqtr
        enddo 
      enddo

999   return
      end

