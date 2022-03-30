      Subroutine GETNAV( debug,icall,iwkn,sow,gnss,iprn,nsats,sats
     .                 , xsat,svclock,satok,iwkne,xetoe )
            
c     Get satellite position and clock from a RINEX navigation file
c     (velocity computed but not passed)

c     Written by R. King   6 Sept 89
c     Modified by K. Feigl Oct 90 to read whole table into memory.
c     Modified by R. King Dec 2015 to use the last value prior values 
c           rather than the closest ones. 
     
      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'

c  Input:
c     debug   : If .true., print stuff to screen
c     If icall = 0 read whole file into memory.
c     If icall > 0 return the ephemeris requested iprn for t=(iwkn,sow)
c     iwkn,sow: GPS week number and seconds-of-week requested
c     gnss    : GNSS requested
c     iprn    : PRN requested
c     nsats   : number SVs in sats array
c     sats(maxsat): PRNs for channels 1-nsats 

      integer*4 icall,iwkn,iprn,nsats,sats(maxsat)
      real*8 sow      
      logical debug

c   Output
c     xsat(3)  : coordinates of SV at t=(iwkn,sow)
c     svclock: SV clock offset at t=(iwkn,sow) (sec)
c     satok: False if bad or missing SV; otherwise true             
c     iwkne:  GPS week number of nav-file entry used
c     sowe :  Seconds-of-week of nav-file entry used
                
      real*8 xsat(3),svclock
      logical satok
      
c   Navigation ephemeris values (see lib/reade.f for description)
            
      integer maxeph
      parameter (maxeph = 10000)  
      integer*4 iwkne
      integer*4 mprns(maxeph),iwknes(maxeph)   
c     number of navigation epochs for any SV (maximum is maxeph)
      integer*4 in
      real*8 trans_sow,sowe
      real*8 bephem(16),        bclock(6)
      real*8 bephems(16,maxeph),bclocks(6,maxeph)
      real*8 bepheml(16),       bclockl(6)
      real*8 subfr1 (8),xetoe
      real*8 subfr1s(8,maxeph)

      real*8  xem0,xedn,xeart,xeecc
     3       ,xei0,xeidt,xeom0,xew,xecuc,xecus,xecrc,xecrs
     3       ,xecic,xecis,xeomd
     4       ,xetoc,xeaf0,xeaf1,xeaf2,xeadc,xetdg
     5       ,erate,t,xn,a
     6       ,anom,e1,e2,enom,cosv,sinv,v,phi,xinc
     7       ,di,du,dr,u,r,xp,yp,asc,dt
     8       ,term,xpdot,ypdot,acheck
     9       ,frame_time,gm
        
c   Other local variables
      integer*4 iflag,iter,irec,i,j,ii
      real*8 secdif,tdiff,tdiff_last,xsatdot(3)
      character*1 gnss,sys 
      character*80 prog_name
      character*256 message  
c     this used for debug: set false after printing the first-epoch values
      logical first_epoch,converged,matched

      data  gm     /398600.80d+09/
     .,     erate  /7.292115147d-05/
           
c  Variables Saved from initial read

      save bephems,bclocks,subfr1s,mprns,iwknes,in,first_epoch

       call rcpar(0,prog_name)      

cd      print *,'GETNAV icall iprn ',icall,iprn
      if (icall .eq. 0) then
                         
        first_epoch = .true.
        in = 0
        iflag = 0
c       rewind to be sure
        rewind (unav)

        do while  (iflag.ne.-1 .and. in.lt.maxeph )
          in = in + 1                            
          if (in .lt. maxeph) then
cd            print *,'calling READE icall in ',icall,in
            call reade( unav,icall,sys 
     .                , iflag,trans_sow,mprns(in),iwknes(in)
     .                , bephems(1,in),bclocks(1,in),subfr1s(1,in) )       
            if( sys.ne.gnss.or.iflag.gt.0 ) in = in - 1
c           the header is read, now set icall to read the rest (RINEX only)
            icall = 1
          else
            write(message,'(a,i6)') 
     .             'Too many ephemeris records in efile: ',in
            call report_stat('WARNING',prog_name,'getnav',' '
     .                      , message,0)
            if( uinfor.gt.0) 
     .         write (uinfor,*) 'GETNAV: too many ephemeris records'
          endif
        enddo 
        write(message,'(a,i6,a)') ' Read in ',in,' ephemeris records.'
        if( uinfor.gt.0 ) write(uinfor,'(a)') message
        if( debug ) write (uscren,'(a)') message
        if( debug.and.mprns(i).eq.7  ) then        
          print *,'printing all PRN 7 records'
          do i=1,in
            if( mprns(i).eq.7) 
     .        print *,'i mprns bephems(1) bephems(5) '
     .         ,i, mprns(i),bephems(1,i),bephems(5,i)
          enddo
        endif   
        return
      endif 
       
 

c  icall > 1: Loop over records looking for the closest time before the 
c             time requested for GPS, Beidou, and Galileo (cannot use 
c             nav-file for Glonass)
c    
      matched = .false.
      tdiff_last = 99999.
      i= 0  
      irec = 0 
      satok = .true.
      do while ( .not.matched .and. i.lt.(maxeph-1)  ) 
        i= i + 1 
c       check only if the satellites match
        if( mprns(i).eq.iprn )  then     
          irec = i
          tdiff = dabs(secdif( iwkn,sow,iwknes(i),bephems(1,i)))
          if( tdiff.lt.tdiff_last ) then               
            tdiff_last = tdiff
          else 
c           set the match when tdiff starts increasing
            matched = .true.  
          endif
        endif
      enddo
      if( irec.eq.0 ) then
        write(message,'(a,a1,i2)') 'No nav-file values for ',gnss,iprn
        call report_stat('WARNING',prog_name,'getnav',' ',message,0)
      else                    
        xetoe = bephems(1,irec)
        iwkne = iwknes(irec)
        do j = 1,8
          subfr1(j) = subfr1s(j,irec)
        enddo
        do j = 1,6
          bclock(j) = bclocks(j,irec)
        enddo
        do j = 1,16
          bephem(j) = bephems(j,irec)
        enddo     
        if (debug) then 
          write(*,'(a,a1,i2,a,i5,f11.3,a,i5,f11.3,a,i4)') 
     .       'Nav-file match for ',gnss,iprn,' at '
     .      ,iwkn,sow,' is ',iwkne,xetoe,' irec ',irec                       
        endif        
      endif
    
c Assign the clock values

      xetoc=bclock(1)
      xeaf0=bclock(2)
      xeaf1=bclock(3)
      if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .    gnss.eq.'J' ) then
        xeaf2=bclock(4)
        xeadc=bclock(5)
        xetdg=bclock(6)
      elseif( gnss.eq.'R' ) then
        frame_time = bclock(4)
      endif

c  Assign the ephemeris values
                        
      if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .  gnss.eq.'J'.or.gnss.eq.'I' ) then
        xetoe  =bephem(1)
        xem0   =bephem(2)
        xedn   =bephem(3)
        xeart  =bephem(4)          
        xeecc  =bephem(5) 
cd        if(iprn.eq.7) print *,'assign xetoe xeecc ',xetoe, xeecc 
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
      elseif( gnss.eq.'R') then
        call report_stat('FATAL',prog_name,'getnav',' '
     .   ,'Glonass nav orbit not accurate enough, must use sp3 file',0)
      else                                         
        write(message,'(a,a1,a)') 'GNSS (',gnss,') not recognized'
        call report_stat('FATAL',prog_name,'getnav',' ',message,0)
      endif 
                      
c  Get the clock offset

      dt= sow - xetoc
      svclock= xeaf0 + xeaf1*dt + xeaf2*dt**2
      if( debug ) 
     .  print *,'SV clock sow xetoc dt xeaf0 xeaf1 xeaf2 svclock '
     .                  , sow,xetoc,dt,xeaf0,xeaf1,xeaf2,svclock  
c     if( gnss.eq.'R') then
c         for glonass need to subtract the leap second
c          call timcon(1,iwkn,sow,iyr,mon,iday,imin,sec,utcoff)
c          svclock = svclock - utcoff
c         endif
c       endif

c  Get the SV position  

      t= sow -xetoe
      if( t.gt.302400.d0 ) t = t - 604800.d0
      if( t.lt.-302400.d0 ) t = t + 604800.d0
      if(debug) print *,'  At t-ref ',t,':',xsat
      acheck = xeart*xeart
      if (dabs(acheck) .lt. 1.d6 .or. dabs(acheck) .gt. 1.d8) then 
        write(message,'(a,f12.6)') 'GETNAV: semimajor A (m) = ',acheck 
        if(debug) print *,message
        if( uinfor.gt.0 ) 
     .     write (uinfor,*) 'GETNAV: semimajor A (m) = ',acheck
        call report_stat('WARNING',prog_name,'getnav',' ',message,0)
        satok = .false.
        return
      else
        xn = dsqrt(gm/xeart**6) +xedn
        a= (gm/xn**2)**0.3333333333333333d0
* MOD TAH 000714: Changed to use real semimajor axis.  Also "retarded"
*         t by 67 ms to account at least partially for light travel 
*         time.  (Strictly calc should be iterated but this gets us
*         closer than ignoring the effect.
        a = acheck
        t = t - 0.067d0
      endif
c
ccc      if(debug.and.first_epoch) then 
      if(debug.and.iprn.eq.18 ) then
        print 401,iprn,t,xetoe,xem0,xn,a
     .     , xeecc,xei0,xeidt,xeom0,xeomd,xew
401    format(/,1x,'Read from nav-file:',//,1x,
     1 'iprn satellite number is ',i2,/,1x,
     2 't = sowk - xetoe , in seconds, is               ',1pe22.15,/,1x,
     3 'xetoe the epoch time for the ephemeris is       ',1pe22.15,/,1x,
     4 'xem0, the mean anamoly is                       ',1pe22.15,/,1x,
     5 'xn is                                           ',1pe22.15,/,1x,
     . 'a the semi-major axis is                        ',1pe22.15,/,1x,
     6 'xeec the orbit eccentricity is                  ',1pe22.15,/,1x,
     7 'xeio the orbit inclination is                   ',1pe22.15,/,1x,
     8 'xeidt the rate of change of orbit inclination   ',1pe22.15,/,1x,
     9 'xeom0 the longitude of the ascending node is    ',1pe22.15,/,1x,
     1 'xeomd the rate of change of the ascending node  ',1pe22.15,/,1x,
     2 'xew the arguement of perigee is                 ',1pe22.15)
        print 402, xecuc,xecus,xecrc,xecrs,xecic,xecis
 402   format(/,1x,
     . 'xecuc the cuc polynomial parameter is           ',1pe22.15,/,1x,
     1 'xecus the cus polynomial parameter is           ',1pe22.15,/,1x,
     2 'xecrc the crc polynomial parameter is           ',1pe22.15,/,1x,
     3 'xecrs the crs polynomial parameter is           ',1pe22.15,/,1x,
     4 'xecic the cic polynomial parameter is           ',1pe22.15,/,1x,
     5 'xecis the cis polynomial parameter is           ',1pe22.15,/,1x)
      endif

      anom= xem0 + xn*t
c     solve kepler's eqn (successive approx. ok for low eccentricity)
      e1= anom
      iter= 0        
      converged = .false.
      do while (.not.converged )
        e2= anom + xeecc*dsin(e1)
        if( dabs(e2-e1).lt.1.d-10 ) then
          converged = .true.
        else
          e1= e2
          iter = iter + 1
          if( iter.gt.10 ) then
           satok = .false.              
           write(message,'(a,i2,a)') 'Keplers eqn for PRN ',iprn
     .                  , ' did not converge in < 10 iterations'
           call report_stat('WARNING',prog_name,'getnav',' ', message,0)
           return
          endif
        endif
      enddo
      enom= e2
      cosv= ( dcos(enom)- xeecc)/(1.d0 - xeecc*dcos(enom) )
      sinv= dsqrt(1.d0 - xeecc**2)*dsin(enom)
     .                        / (1.d0 - xeecc*dcos(enom))
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
      asc= xeom0 + (xeomd - erate)*t - erate*xetoe

      if(debug.and.first_epoch) then
        print 888,t,xn,anom,e1,enom,cosv,sinv,v
     .             ,  phi,du,dr,di,u,xinc,asc
 888    format (/,1x,'dt the time difference is ',1pe22.15,/,1x,
     1   'xn is ',1pe22.15,/,1x,
     2   'anom is ',1pe22.15,/,1x,
     3   '  e1 is ',1pe22.15,/,1x,
     4   'enom is ',1pe22.15,/,1x,
     5   'cosv is ',1pe22.15,/,1x,
     6   'sinv is ',1pe22.15,/,1x,
     7   '   v is ',1pe22.15,/,1x,
     8   ' phi is ',1pe22.15,/,1x,
     9   '  du is ',1pe22.15,/,1x,
     1   '  dr is ',1pe22.15,/,1x,
     2   '  di is ',1pe22.15,/,1x,
     3   '   u is ',1pe22.15,/,1x,
     4   'xinc is ',1pe22.15,/,1x,
     5   ' asc is ',1pe22.15)
        first_epoch = .false.
      endif
     
      xsat(1) = xp*dcos(asc) - yp*dcos(xinc)*dsin(asc)
      xsat(2) = xp*dsin(asc) + yp*dcos(xinc)*dcos(asc)
      xsat(3) = yp*dsin(xinc)
      term=(xn*a)/dsqrt(1.d0-xeecc*xeecc)
      xpdot=-dsin(u)*term
      ypdot=(xeecc+dcos(u))*term
      xsatdot(1)= xpdot*dcos(asc) - ypdot*dcos(xinc)*dsin(asc)
      xsatdot(2)= xpdot*dsin(asc) + ypdot*dcos(xinc)*dcos(asc)
      xsatdot(3)= ypdot*dsin(xinc)
                           
cd      print *,'test iprn sow irec satok ',iprn,sow,irec,satok 
      return
      end


