      subroutine recod(isit1,isit2,it,dt,omc,record)
c
c     record all effective observation.
c     1. set up time resolution threshold as 0.5 year
c     2. sorting record table by site index and time
cc    This sorting doesn't seem to happen here since we now 
cc        minic=1 hardwired in readrv.f
c
      include 'solvem.fti'
c     integer isit1,isit2,i1,i2,it,i,ia,ib,itime,it2
      integer isit1,isit2,i1,i2,it,i,ia,ib,it2
c     integer id,ic,j
      integer id,j
c     integer record(4000,4)
      real*8 record(4000,4),itime,ic
      real*8 dt,omc,tlim

c     catch it if the site number is larger than ok
      if (isit1.gt.nsit.or.isit2.gt.nsit) goto 50

c     geodetic velocity
      if (it.eq.24) then
         i1 = 0
         i2 = isit1
c        itime = nint(dt*2.0d0)
         itime = dt*2.0d0
         goto 30
      endif

c     baseline vetors or vector rates 
      if (it.ge.25.and.it.le.30) then
c        itime = nint(dt*2.0d0)
         itime = dt*2.0d0
         i1 = isit1
         i2 = isit2
         omc = 0.0d0
         goto 30
      endif
c
c     I don't understand the need for this
c     tlim = 2.0d0
      tlim = 1.0d0

c     terrestrial rate observations need 1 year
      if (it.gt.10.and.it.le.20) tlim = 365.24d0
     
c     itime = nint(dt*tlim)
      itime = dt*tlim  
     
c     time resolution threshold is defined as 0.5 years 
c     (default) or 0.5 days if the m/yr option is specified
c     ie, iomode(18)=1

      if (iomode(18).gt.0) then
c        round to nearest 0.5 days ~= 1.0d-3 years
c        itime = anint(itime*1.0d3)/1.0d3
c        itime = anint(itime*365.25d0)/365.25d0
         itime = anint(itime)*1.0d0
c        print*, 'recod: itime ',itime
      else
c        default to nearest integer (nearest 0.5 years)
         itime = anint(itime)*1.0d0
      endif

      i1 = min(isit1,isit2)
      i2 = max(isit1,isit2)
c  
c      at the present time miniv is hardwired =1
      if (jobs.eq.0.or.miniv.eq.1) goto 30
c
c     arrange records by site index
      id = 1
      do 10 i = 1,jobs
         ia = nint(record(i,1))
         ib = nint(record(i,2))
         ic = record(i,3)
         it2 = nint(record(i,4))
         
c        test

         if (i1.lt.ia) then
            id = i
            goto 40
         endif
         if (i1.eq.ia) then
            if (i2.lt.ib) then
               id = i
               goto 40
            endif
            if (i2.eq.ib.and.itime.lt.ic) then
               id = i
               goto 40
            endif
            if (it.eq.it2.and.i2.eq.ib.and.itime.eq.ic) goto 100
         endif
 10   continue
      
 30   jobs = jobs+1
      record(jobs,1) = real(i1)
      record(jobs,2) = real(i2)
      record(jobs,3) = itime
      record(jobs,4) = real(it)
      anorm(jobs) = omc
      goto 100
c
 40   jobs = jobs+1
      do i = jobs,id+1,-1
         do j = 1,4
            record(i,j) = record(i-1,j)
         enddo
         anorm(i) = anorm(i-1)
      enddo
      record(id,1) = real(i1)
      record(id,2) = real(i2)
      record(id,3) = itime
      record(id,4) = real(it)
      anorm(id) = omc
      goto 100
c
 50   print*,' Site index overflow:'
      print*,isit1,' ',isit2,' > ',nsit
c
 100  continue
      return
      end
c     

