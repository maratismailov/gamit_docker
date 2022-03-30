c
      Subroutine FILOMC( istat,phi,phi2,iphi,omc
     .                 , isnr,ierrfl,iertag,el,elev)
c
c     fill phase and wl arrays with omc data
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'models.h' 

      logical okdata,lgood,iskip_dec,iskip_cut     
      integer*4 istat,jsat
      integer*4 iee,iga
      real*8 facwl
      integer isnr(maxdat,maxsat)
      real*8  omc(maxdat,maxsat),phi(maxsit,maxsat),phi2(maxsit,maxsat)
      integer iphi(maxsit,maxsat),ierrfl(maxsat)
      real*8  elev(maxsit,maxsat),el(maxsat)
      integer iertag(maxsit,maxsat) ,i
     
      logical debug/.false./

      if (iepoch.lt.istart) goto 100

      facwl=17.0d0/137.0d0

      do 30 jsat=1,nsat
c
c Bob's pushing approach

         iee = iertag(istat,jsat)
         if (iee.ne.10) iertag(istat,jsat) = ierrfl(jsat)

c** rwk 190813: iphi was never initialized so obsevations that are
c**             not present (e.g. L1-only) but not flagged, get used.
c**             With earlier code, the array was zero but not now.
c**             No, though true, this doesn't seem to be the explanation
         iphi(istat,jsat) = 0                  
         if(debug) print *,'FILOMC 0 iphi(4,,i) '
     .     ,iepoch,(iphi(4,i),i=1,32)



c Simon's decimation approach
         iskip_dec = .false.
         if (mod((iepoch-istart),idecim).ne.0) iskip_dec = .true.
         if(debug)   print *,'FILOMC iepoch istat jsat elevcut el '
     .              ,iepoch,istat,jsat,elevcut(istat)*convd,el(jsat) 
  
c Check to ensure data point is above specified cutoff angle if requested
         iskip_cut = .false.
         if ( el(jsat) .lt. elevcut(istat)*convd ) iskip_cut = .true.

c Check to see if data is flagged good for use in the solution
         okdata = .false.
               if (     lgood(ierrfl(jsat))
     .            .and. isnr(1,jsat).ge.minsnr
     .            .and. iusesa(jsat).ne.0
     .            .and. iusest(istat).ne.0
     .            .and. .not. iskip_dec
     .            .and. .not. iskip_cut  ) okdata = .true.          
cd          if(debug.and.istat.eq.4.and.jsat.eq.25) then
          if(debug) then 
            print *,
     .       'FILOMC istat jsat ierrfl isnr iusesa iusest okdata '
     .             ,  istat,jsat,ierrfl(jsat),isnr(1,jsat),iusesa(jsat)
     .             , iusest(istat),okdata
            print *,'FILOMC 1 iphi(4,,i) ',iepoch,(iphi(4,i),i=1,32)
           endif 

c If data is good store in arrays for later use
         if (okdata) then
c
            elev(istat,jsat) = el(jsat)
            iphi(istat,jsat) = 1
            if(l2flag.ge.-1) then
              phi(istat,jsat) = omc(1,jsat)
            elseif(l2flag.eq.-2) then
              phi(istat,jsat) = omc(2,jsat)
            endif
c
            if(l2flag.gt.0) then
               phi2(istat,jsat) = omc(2,jsat)
            else
               phi2(istat,jsat) = 0.d0
            endif
            wl1(istat,jsat) = omc(1,jsat)-omc(2,jsat)
            wl1(istat,jsat) = wl1(istat,jsat)
     .            -(omc(3,jsat)+omc(4,jsat))*facwl
            if(lquick.ge.1) then
               iga = intsam(istat)/inter
               if( iga.lt.1 ) iga = 1
c  **  add bias at gap only if usual quick algorithm (lquick = 2 means
c      full solution but with implicit biases)
               if( lquick.eq.1 ) then
                  if (igaps(istat,jsat).ge.iga) iertag(istat,jsat) = 10
               endif
            endif
            igaps(istat,jsat) = 0

c If data not good zero out data arrays, and set no data flag
         else                
            if(debug) 
     .          print *,'FILOMC 2 iphi(4,,i) ',iepoch,(iphi(4,i),i=1,32)

            elev(istat,jsat) = 0.d0
            iphi(istat,jsat) = 0
            phi(istat,jsat) = 0.d0
            phi2(istat,jsat) = 0.d0
            wl1(istat,jsat) = 0.d0
            igaps(istat,jsat) = igaps(istat,jsat)+1
         endif    
         if(debug) then                     
           print *,'FILOMC 3 iphi(4,,i) ',iepoch,(iphi(4,i),i=1,32)
           print *,'FILOMC iepoch istat jsat phi '
     .      ,iepoch,istat,jsat,phi(istat,jsat)
         endif
 30   continue
    

 100  continue   
      if(debug)
     .    print *,'FILOMC 3 iphi(4,,i) ',iepoch,(iphi(4,i),i=1,32)
      return
      end
