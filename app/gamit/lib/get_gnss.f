      Subroutine get_gnss( gnss,nprn,asvid,isvid,nobtyp,obs,illi,issi )

c     Remove satellites for GNSS other than the one requested from the RINEX arrays
c     Renamed and modified from lib/gps_only.f   R. King 27 August 2014

      implicit none
              
      include '../includes/dimpar.h'
      include '../includes/makex.h'

      integer*4 i,j,k
                                
c        GNSS system requested 
      character*1 gnss

c        number of SVs --reset on output
      integer*4 nprn

c        character id for SV (G/blank=GPS, R=GLONASS)
      character*1 asvid(maxchn)

c        numberial id for SV (PRN for GPS)
      integer*4 isvid(maxchn)
                                      
c        number of observation types for requested GNSS SVs
      integer*4 nobtyp
               
c        observations
      real*8 obs(maxchn,maxobt)

c        loss-of-lock and signal strength indicators
      integer*4 illi(maxchn,maxobt),issi(maxchn,maxobt)

c        index to requested GNSS SVs 
      logical lgnss(maxchn)

c        temporary variables and arrays for sorting
      character*1 asvtmp(maxchn)
      integer*4 ntmp,isvtmp(maxchn),illitmp(maxchn,nobtyp)
     .        , issitmp(maxchn,maxobt)
      real*8 obstmp(maxchn,maxobt)
                         
              
cd      print *,'GET_GNSS obs(1-5) '
cd      do i=1,nprn
cd        write(*,*) i,(obs(i,j),j=1,5)
cd      enddo
      do i=1,nprn
         if( asvid(i).eq.' ' ) asvid(i) = 'G'
      enddo

      ntmp = 0                       
      do i=1,nprn
        if( asvid(i).eq.gnss ) then
          lgnss(i) =  .true.
           ntmp = ntmp + 1
        else
          lgnss(i) = .false.
        endif
      enddo              

c       if data found for only the system requested, nothing to do
      if( ntmp .eq. nprn ) then
         return                     

      else
                    
         j = 0
         do i=1,nprn
            if( lgnss(i) ) then
               j = j +1
               isvtmp(j) = isvid(i)
               asvtmp(j) = asvid(i)
               do k=1,nobtyp
                 obstmp(j,k) = obs(i,k)
                 illitmp(j,k) = illi(i,k)
                 issitmp(j,k) = issi(i,k)
               enddo
            endif
         enddo
 
         do i=1,ntmp    
           asvid(i) = asvtmp(i)
           isvid(i) = isvtmp(i)
           do k=1,nobtyp
             obs(i,k) = obstmp(i,k)
             illi(i,k) = illitmp(i,k)
             issi(i,k) = issitmp(i,k)
           enddo
         enddo     
         do i=ntmp+1,nprn
           asvid(i) = ' '
           isvid(i) = 0
         enddo 
         nprn = ntmp

      endif

cd      print *,'End GET_GNSS obs '
cd      do i=1,nprn
cd        write(*,*) i,(obs(i,j),j=1,5)
cd      enddo
      return
      end
