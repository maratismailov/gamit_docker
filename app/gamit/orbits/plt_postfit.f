      Subroutine plt_postfit( iop,pltnam )                             
     
c   Create files from which to plot the postfit residuals
c 
c   Written by S. McClusky  Sept 1994; modified by R. King Dec 2000

        implicit none
c
      include '../includes/dimpar.h'
      include '../includes/orbits.h'
      include 'orbfit.h'
                            
      logical first_epoch

      character*3 in_sys, out_sys
      character*28 pltf
      character*80 pltnam

      integer*4 iop,iepoch,jsat,kepoch,lu,indx,indx_end,i,j
      integer*4 num_epochs(maxsat)
   
      real*8 rotmat(3,3)
      real*8 az_diff(maxsat,maxepc)
      real*8 len_diff(maxsat,maxepc)

      real*8 dt,km_m
      real*8 xyzsat_pos(3), llhsat_pos(3), dxyz(3), out_comp(3)
      real*8 orbllh(3,maxsat,maxepc),dneu(3,maxsat,maxepc),pi          

      real*8 maxdneu(3,maxsat),maxdxyz(3,maxsat),maxdrac(3,maxsat)

* MOD TAH 200328: Added actual time to residal output.
      real*8 MJD  ! Modified JD + fractional part of day for time
                  ! of each residual.

c*    real*8 maxdvel(3,maxsat)
      data km_m/1000.d0/
      parameter ( pi            = 3.1415926535897932D0 )

c   If ground-track values needed, compute them now

      if( iop.eq. 4 ) then
                              
        in_sys = 'XYZ'
        out_sys = 'NEU'
                                        
        do iepoch = iepstart,iepstop

          do jsat = 1,nsat 
                     
            if( lobs(jsat,iepoch) ) then
              do i = 1,3
                xyzsat_pos(i) = ref_pos_vel(i,jsat,iepoch)*1000.d0
                dxyz(i) = omc(i,jsat,iepoch) 
              enddo
c             convert X,Y,Z to Lat, Long, Height, and delta X,Y,Z, to delta N,E,U...
              call rotate_geod( dxyz,out_comp,in_sys,out_sys,xyzsat_pos
     .                     , llhsat_pos,rotmat ) 
c             convert Lat, Long, Height from radians/meters to degrees/kilometres.
              do i = 1,3
                if (i.eq.1) then
                 orbllh(i,jsat,iepoch)=90.d0 - (llhsat_pos(i)*180.d0/pi)
                elseif(i.eq.2) then
                  orbllh(i,jsat,iepoch) = llhsat_pos(i)*180.d0/pi
                elseif(i.eq.3) then
                  orbllh(i,jsat,iepoch) = llhsat_pos(i)/1000.d0
                endif
                dneu(i,jsat,iepoch) = out_comp(i) 
              enddo

c             compute Azimuth and Length of horizontal orbit difference components......
              az_diff(jsat,iepoch) = datan2(out_comp(2),out_comp(1)) 
              az_diff(jsat,iepoch) = az_diff(jsat,iepoch)*180.d0/pi
              if ( az_diff(jsat,iepoch).lt.0.d0)
     .           az_diff(jsat,iepoch) =  360.d0+az_diff(jsat,iepoch)
              len_diff(jsat,iepoch)=dsqrt(out_comp(2)**2+out_comp(1)**2)
            endif
c           endif on valid obs
          enddo
c         end loop on SVs
        enddo
c       end loop on epochs
      endif
c    --end if on ground-track values


c   Count the epochs and find the maximum value for each component

      do j=1,nsat
        num_epochs(j) = 0
        do i=1,3
          maxdxyz(i,j)=0.d0
          maxdrac(i,j)=0.d0
          maxdneu(i,j)=0.d0
        enddo
      enddo
      do iepoch = iepstart,iepstop 
         do jsat=1,nsat
           if( lobs(jsat,iepoch) ) then 
             num_epochs(jsat) = num_epochs(jsat) + 1
             do i=1,3
               if (dabs((omc(i,jsat,iepoch))).gt.maxdxyz(i,jsat))
     .               maxdxyz(i,jsat) = dabs((omc(i,jsat,iepoch)))
               if (dabs((drac(i,jsat,iepoch))).gt.maxdrac(i,jsat)) 
     .               maxdrac(i,jsat) = dabs((drac(i,jsat,iepoch)))
               if (iop.eq.4 ) then
                 if (dabs((dneu(i,jsat,iepoch))).gt.maxdneu(i,jsat))
     .               maxdneu(i,jsat) = dabs((dneu(i,jsat,iepoch)))
               endif
             enddo
           endif
         enddo
      enddo


c   Convert time interval from day into min

      dt=delt(1)/60
      
c   Prepare plot output files
             
      indx = 1
c     find next nonblank character of pltfile descriptor
      do while (pltnam(indx:indx) .eq. ' ' .and. indx .lt. 25)
         indx = indx+1
      enddo     
c     find the next blank character of pltfile descriptor
      indx_end = indx
      do while (pltnam(indx_end:indx_end) .ne. ' '
     .            .and. indx_end .lt. 25)
c     increment end indx
        indx_end = indx_end + 1
      enddo   

c   Loop over all SVs and write the plot files

      do jsat =1,nsat

c       create the SV-specific filename and open the file
        write(pltf,'("plt_",a,".",i2)') pltnam(indx:(indx_end-2))
     .                               ,isat(jsat)
        if ( pltf(indx_end+4:indx_end+4) .eq. " " ) then
          write(pltf(indx_end+4:indx_end+4),'("0")')
        endif                                       
        lu = luplt(jsat)
        open(unit=lu,file=pltf,status='unknown')

c       write the label text

        if (iop.eq.1) then
          write(lu,'("Delta X (m)")')
          write(lu,'("Delta Y (m)")')
          write(lu,'("Delta Z (m)")')
        endif
        if (iop.eq.2) then
          write(lu,'("D-radial (m)")')
          write(lu,'("D-along (m)")')
          write(lu,'("D-cross (m)")')
        endif
        if (iop.eq.3) then
          write(lu,'("V-X (km/s)")')
          write(lu,'("V-Y (km/s)")')
          write(lu,'("V-Z (km/s)")')
        endif
        if (iop.eq.4) then
          write(lu,'("Latitude  (deg)")')
          write(lu,'("Longitude (deg)")')
          write(lu,'("Height (km)")')
          write(lu,'("Delta N (m)")')
          write(lu,'("Delta E (m)")')
          write(lu,'("Delta U (m)")')
          write(lu,'("VECT AZ (deg)")')
          write(lu,'("VECT LEN (m)")')
        endif
        write(lu,'("Number epochs (Time Interval=",f6.2," min)")') dt
        write(lu,'(i4)') num_epochs(jsat)

c       write time series series
c MOD TAH Output more digits (f6.2 to f9.5)
        first_epoch = .true.
        do kepoch = 1,nepoch
        
          if( lobs(jsat,kepoch) ) then
* MOD TAH 200328: Get the MJD+fractional day of this epoch and 
*           add to end of line so that users know what time 
*           the residual refers.
*           delt(1) -- seconds, jdb(1) PEP JD (2400001+MJD)
*           tb(1) -- Seconds of day for start.
            mjd = (jdb(1)-2400001) + 
     .            (tb(1)+(kepoch-1)*delt(1))/86400.d0
            if (iop.eq.1) then
              if (first_epoch) then
                write( lu, '(f9.5," Max DX")') maxdxyz(1,jsat) * km_m
                write( lu, '(f9.5," Max DY")') maxdxyz(2,jsat) * km_m
                write( lu, '(f9.5," Max DY")') maxdxyz(3,jsat) * km_m
                first_epoch = .false.
              endif
              write(lu,'(i4,3f11.5,1x,F13.5)') kepoch
     .                        ,(omc(i,jsat,kepoch)*km_m,i=1,3), mjd
            endif
            if (iop.eq.2) then  
              if (first_epoch) then
                write( lu, '(f9.5," Max R")') maxdrac(1,jsat) * km_m
                write( lu, '(f9.5," Max A")') maxdrac(2,jsat) * km_m
                write( lu, '(f9.5," Max C")') maxdrac(3,jsat) * km_m 
                first_epoch = .false.
              endif
              write(lu,'(i4,3f11.5,1x,F13.5)') 
     .               kepoch,(drac(i,jsat,kepoch)*km_m,i=1,3),mjd
            endif
            if (iop.eq.3) then             
              call report_stat('FATAL','ORBFIT'
     .       ,'orbits/plt_postfit',' ','Velocity plots not yet coded',0)
c              if (first_epoch) then
c                write( lu, '(f9.5)') maxdvel(1,jsat) * km_m
c                write( lu, '(f9.5)') maxdvel(2,jsat) * km_m
c                write( lu, '(f9.5)') maxdvel(3,jsat) * km_m 
c                first_epoch = .false.
c              endif
c              write(lu,'(i4,3f11.5)') kepoch,(dvel(i,jsat,kepoch),i=1,3)
            endif
            if (iop.eq.4) then   
              if (first_epoch) then
                write( lu, '(f9.5)') maxdneu(1,jsat) * km_m
                write( lu, '(f9.5)') maxdneu(2,jsat) * km_m
                write( lu, '(f9.5)') maxdneu(3,jsat) * km_m   
                first_epoch = .false.
              endif
              write(lu,'(i4,8f14.7,1x,F13.5)') kepoch
     .           ,(orbllh(i,jsat,kepoch)*km_m,i=1,3)
     .           ,(dneu(i,jsat,kepoch)*km_m,i=1,3)
     .           ,az_diff(jsat,kepoch),len_diff(jsat,kepoch)*km_m, mjd
            endif           

          endif
c         endif on valid obs                

        enddo
c       end loop on epochs
      
      enddo
c---- end loop on SVs

      return
      end
