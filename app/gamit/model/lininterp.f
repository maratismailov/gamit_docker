      subroutine lininterp( a,b,val_time,valtab,obstime,nrow,ncol
     .                    , valinterp,exttim )

c  performs a linear interpolation to calculate the value
c  at a particular time. Written specifically to interpolate
c  the tabulated 6-hourly atmospheric pressure loading 
c  displacements but could work equally well for any other
c  circumstance.
c
c  INPUT:
c    a,b         : dimensions of the matrices in the calling program
c    val_time    : array of time tags of tabulated values
c    valtab      : array of corresponding tabulated values
c    obstime     : time to which the tabulated values should be
c                  interpolated
c    nrow        : number of rows of tabulated values     
c    ncol        : number of independent variables to be
c                  interpolated 

c
c  OUTPUT:
c     valinterp  : interpolated values 
c     exttim     : if requested interpolation is outside the range
c                  of the table, this variable will have the time
c                  (neg or pos) by which the range is exceeded;
c                  if within range, extflg = 0.
c                  
c
c  P. Tregoning
c  8 January 2004

      implicit none

      integer*4 a,b,ncol,nrow,t1,t2,i

      real*4 val_time(a),valtab(a,b),exttim,valinterp(ncol)
      real*8 obstime,dval,dt,tstep,delta,eps
      character*256 message

      logical found
c      logical debug
c      data debug/.false./

c debug  
c      integer*4 j
c      print*,'a,b,obstime',a,b,obstime
c      print*,'val_time',val_time
c      print*,'nrow and ncol',nrow,ncol
c      print*,'valtab:'
c      do i = 1,nrow
c        print*,(valtab(i,j),j=1,ncol)
c      enddo
                                    
c     eps is a small number tolerance to allow use of the end values
c     currently set to be 90s if times are in days
      data eps/.001d0/
                  
c      if( obstime.gt.220.082 ) debug = .true.
c      if( obstime.gt.220.0826 ) debug = .false.

c   see if the requested time is outside the range of the array   
      exttim = 0.        
c      if( debug ) then
c        print *,'obstime nrow val_time(n) eps '
c     .     ,obstime,nrow,val_time(nrow),eps
c        endif
      if( obstime.lt.(val_time(1)-eps) )  then
        exttim = obstime - val_time(1)  
        do i=1,ncol
         valinterp(i) = valtab(1,i)
        enddo
        return
      elseif( obstime.gt.(val_time(nrow)+eps) ) then
        exttim =  obstime - val_time(nrow) 
c        print *,'exttim ',exttim 
        do i=1,ncol
          valinterp(i) = valtab(nrow,i)
        enddo
        return
      endif

c   run through the time array to find which two times the requested
c   time lies between
      found = .false.
      t1 = 0
      do while (.not.found)
        t1 = t1+1
        if(t1+1.gt.nrow)then
          write(message,'(a,f15.7,a,2f15.7)')
     .     'Requested time',obstime,' outside tabulated values. Min/max'
     .    ,val_time(1),val_time(nrow)
          call report_stat('FATAL','MODEL','lininterp',' '
     .                  ,message,0)
        endif
c PT050415: there can be a roundoff problem if obstime is infinitessimally
c           less that the first val_time(t1) .... ! 
c RWK070111: Modified to allow interpolation within eps of last value
        if((dabs(obstime-val_time(t1)).lt.0.01d0.or.
c     .   obstime.ge.val_time(t1)).and.obstime.lt.val_time(t1+1))then
     .    obstime.ge.val_time(t1)).and.obstime.lt.(val_time(t1+1)+eps)) 
     .          then
          found = .true.
        endif
c PT061201: trap the case that the obs time is exactly the last time in the array
c RWK070111: Now allow it to be eps later than last time
c        if(dabs(obstime-val_time(t1+1)).lt.0.00001d0)then 
        if(dabs(obstime-(val_time(t1+1)+eps)).lt.0.00001d0)then 
          found = .true.
        endif
      enddo    

      t2 = t1 + 1

c      print*,'a,b,ncol,nrow,obstime',a,b,ncol,nrow,obstime
c      print*,'lininterp: input values'
c      do i=1,nrow
c        print*,(valtab(j,i),j=1,ncol)
c      enddo
c      print*,'lininterp: times:',(val_time(i),i=1,nrow)
c      print*,'lininterp: t1 and t2:',t1,t2

c  now linearly interpolate all the parameters between these two times
      do i = 1,ncol
        dt = obstime-val_time(t1)
        tstep = val_time(t2)-val_time(t1)
        dval = valtab(t2,i)-valtab(t1,i)

        delta = dval * dt / tstep

c        print*,'dval,dt,tstep,delta',dval,dt,tstep,delta
c        print*,'val at t1 and t2',valtab(t1,i),valtab(t2,i)

        valinterp(i) = valtab(t1,i)+delta
c        print*,'val:',valinterp(i)
      enddo

      return
      end

