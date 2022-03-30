      Subroutine write_summary (ntfiles,tfnam )
     
c     Calculate the post-fit rms and write the orbital summaries to the rms and fit files

* MOD TAH 180106: Added full 3-D RMS to output summary
 
      implicit none

      include '../includes/dimpar.h'
      include '../includes/orbits.h'
      include 'orbfit.h'
                

      integer*4 ntfiles,iepoch,jsat,nobstot,islot2,i,j,k,jbad_sat

      real*8 sumxyz(3,maxsat),sumrac(3,maxsat)
     .     , sum2xyz(3,maxsat),sum2rac(3,maxsat),rmstot(maxsat)
     .     , rmsxyz(3,maxsat),rmsrac(3,maxsat)
     .     , mean_rmsxyz(3),mean_rmsrac(3),postfit_sum2,mean_rmstot
     .     , prefit_rms,chidf,big,sigma,fract,dot,spanhrs
     .     , adjtmp,mts_mas,pre_prm,post_prm,km_m

      real*8 sum23D(maxsat)   ! Sum of 3D difference squured
      real*8 rms3d(maxsat)    ! Computed 3D RMS
      real*8 mean_rms3D       ! Mean value of 3D RMS

      character*19 dattim 
      character*80 tfnam(maxtfil)  
      character*256 message

      data big/999.99/,mts_mas/15.d0/,km_m/1000.d0/
                  
c     Prepare and write out the rms and fit file headers
                                                                   
      call runtim( dattim )
      write(lurms,'(a,a19)') 'Summary of rms differences: ORBFIT run '
     .                       ,dattim 
      write(lurms,'(a,a14)') 'Reference T-file: ',tfnam(1)(1:14)
      write(lurms,'(a,20a14)') 'External T-files: '
     .                         ,(tfnam(i)(1:14),i=2,ntfiles)  
      spanhrs = nepoch*delt(1)/3600.d0           
      write(lurms,'(a,f7.2,a)') 'Reference T-file span: '
     .                          ,spanhrs,' hrs'

* MOD TAH 200121: List any bad satellites
      if( nbad_sat.gt.0 ) then
* MOD TAH 200706: Increased for,at length to allow for more than 10 
*         bad satellites.
          write(lurms,'(a,I3,1x,a,1x,50I3)') 'Deleted Satellites', 
     .          nbad_sat,'List', ibad_sat(1:nbad_sat) 
      endif
      write(lurms,'(/,2a,/2a)')
     . 'PRN   Total     delta-X   delta-Y   delta-Z '
     .,'  d-Radial  d-Along   d-Cross      d-3D'
     .,'=============================================================='
     .,'====================='
      write(lufit,'(a,a19)') 'Parameters adjustments: ORBFIT run '
     .                       ,dattim 
      write(lufit,'(a,a14)') 'Reference T-file: ',tfnam(1)(1:14)
      write(lufit,'(a,20a14)') 'External T-files: '
     .                         ,(tfnam(i)(1:14),i=2,ntfiles)  
      spanhrs = nepoch*delt(1)/3600.d0
      write(lufit,'(a,f7.2,a)') 'Total overlapped time range: '
     .                          ,spanhrs,' hrs'


c     Since we want the postfit residuals for every SV for every arc, 
c     do the calculations explicitly
       
      nobstot = 0
      do i=1,nsat
        nobstot = nobstot + nobs(i)
      enddo
      call report_stat('STATUS','ORBFIT','orbits/write_summary',' '
     .                , 'Calculating residuals',0)
      prefit_sum2 = 0.d0
      do iepoch = iepstart,iepstop
        do jsat = 1, nsat
           if( lobs(jsat,iepoch) ) then
             do i=1,3
                prefit_sum2 = prefit_sum2 + omc(i,jsat,iepoch)**2
             enddo
           endif
        enddo
      enddo   
      prefit_rms = dsqrt(prefit_sum2/dble(nobstot)) * km_m
      write(lufit,'(a,f8.5)') 'Prefit  RMS (m): ',prefit_rms
c     postfit residuals 
      if( nparam.gt.0 ) then 
        do iepoch = iepstart,iepstop
          do jsat = 1, nsat
            if( lobs(jsat,iepoch) ) then
              do i=1,3
                do j=1,nparam  
                  omc(i,jsat,iepoch) = omc(i,jsat,iepoch) - 
     .                            adjust(j)*part(i,j,jsat,iepoch)
                enddo
              enddo 
            endif
          enddo
        enddo
        postfit_sum2 = 0.d0
        do iepoch = iepstart,iepstop
          do jsat = 1, nsat
             if( lobs(jsat,iepoch) ) then
               do i=1,3
                  postfit_sum2 = postfit_sum2 + omc(i,jsat,iepoch)**2
               enddo
             endif
          enddo
        enddo 
c       since weight is 1.0 (m), sqrt(chi2/df) same as rms
        chidf = dsqrt(postfit_sum2/dble(nobstot))

        write(lufit,'(a,f8.5)') 'Postfit RMS (m): ',chidf * km_m
* MOD TAH 090227: Set a minimum here in case we are fitting exactly
        if( chidf.lt.1.d-6 ) then  ! Less than 1 mm in km
            chidf = 1.d-6
        endif

c       Write parameters adjustments
             
        write(lufit,'(/,3a,/)') '    PARAMETER                  '
     .  , '       A PRIORI   ADJUST          POSTFIT    SIGMA  '
     .  , 'FRACT'
        do i=1,nparam
          islot2 = mod(islot(i),100)
          if( islot(i).le.3 ) then
c           convert translations from km to m
            pre_prm = apr_prm(i) * km_m
            adjtmp = adjust(i) * km_m    
            sigma = dsqrt(amat(i,i)) * chidf * km_m  
          elseif( islot(i).eq.4. or. islot(i).eq.5 .or. islot(i).eq.6
     .             .or. islot(i).eq.8 .or. islot(i).eq.9  ) then
c           convert X- or Y-rotations (and Z-inertial) from arc-seconds to mas
            pre_prm = apr_prm(i) * 1000.d0 
            adjtmp = adjust(i) * 1000.d0        
            sigma = dsqrt(amat(i,i)) * chidf * 1000.d0
          elseif( islot(i).eq.10 ) then
c           convert terrestrial Z-rotation from time seconds to mas
            pre_prm = apr_prm(i) * 1000.d0 * mts_mas
            adjtmp = adjust(i) * 1000.d0 * mts_mas 
            sigma = dsqrt(amat(i,i)) * chidf * 1000.d0 * mts_mas
          elseif( islot(i).eq.11 .or. islot(i).eq.12 ) then
c            convert X- or Y- rotation rates from arc-sec/day to mas/day
            pre_prm = apr_prm(i) * km_m 
            adjtmp = adjust(i) * 1000.d0 
            sigma = dsqrt(amat(i,i)) * chidf * 1000.d0
          elseif( islot(i).eq.13 ) then
c            convert Z rotation rate from time-sec/day to mas/day       
            pre_prm = apr_prm(i) * 1000.d0 * mts_mas
            adjtmp = adjust(i) * 1000.d0 * mts_mas  
            sigma = dsqrt(amat(i,i)) * chidf * 1000.d0 * mts_mas
          elseif( islot(i).eq.7)  then    
c           convert scale to ppb
            pre_prm = apr_prm(i) * 1.d9
            adjtmp = adjust(i) * 1.d9    
            sigma = dsqrt(amat(i,i)) * chidf * 1.d9
          elseif( islot(i).ge.100 ) then
c           orbital parameters XYZ
            if( islot2.ge.1.and.islot2.le.3 ) then
c             convert ICs from km and km/s to m and m/s
              pre_prm = apr_prm(i) * km_m 
              adjtmp = adjust(i) * km_m   
              sigma = dsqrt(amat(i,i)) * chidf * km_m
c           orbital parameters Velocities
            elseif( islot2.ge.4.and.islot2.le.6 ) then
c             convert ICs from km/s to mm/s
              pre_prm = apr_prm(i) * km_m *1000
              adjtmp = adjust(i) * km_m  * 1000 
              sigma = dsqrt(amat(i,i)) * chidf * km_m *1000
              call sub_char(prmnam(i),'(m/s)','(mm/s)')
            else
c             no conversion for radiation pressure parameters
              pre_prm = apr_prm(i)
              adjtmp = adjust(i) 
              sigma = dsqrt(amat(i,i)) * chidf
            endif
          endif   
          post_prm = pre_prm + adjtmp  
          fract = adjtmp/sigma  
c         add dividers between sections
          if( islot(i).gt.100 .and. islot2.eq.1 ) write(lufit,'(1x)')  
          write(lufit,'(a30,f16.5,f9.5,f17.5,f9.5,f7.1)')  
     .            prmnam(i),pre_prm,adjtmp,post_prm,sigma,fract
        enddo 

      else
c        no estimation, use the prefit residuals
         chidf = dsqrt(prefit_sum2/dble(nobstot))
         write(lufit,'(a,f8.5,1x,a)') 'Postfit RMS (m): ',chidf * km_m,
     .        'No Parameters'

      endif


c     Get residuals and statistics for each SV

      bad_rmstot = 0.0
      jbad_sat = 0
                      
      do jsat = 1,nsat
        sum23D(jsat) = 0.d0
        do i=1,3           
          sumxyz(i,jsat) = 0.d0
          sumrac(i,jsat) = 0.d0 
          sum2xyz(i,jsat) = 0.d0
          sum2rac(i,jsat) = 0.d0
        enddo                               
        do k = iepstart,iepstop  
c         see if an observation at this epoch for this SV
          if( lobs(jsat,k) ) then    
c           get radial, along-track,cross-track    
cd            print *,'WRITE_SUMMARY jsat epoch ref_pos_vel ',
cd     .        jsat,k,(ref_pos_vel(i,jsat,k),i=1,6)
            call xyz2rac( ref_pos_vel(1,jsat,k),omc(1,jsat,k)
     .                , drac(1,jsat,k))
            do i=1,3                           
              sumxyz(i,jsat)  = sumxyz(i,jsat)  + omc(i,jsat,k)
              sum2xyz(i,jsat) = sum2xyz(i,jsat) + omc(i,jsat,k)**2
              sumrac(i,jsat)  = sumrac(i,jsat)  + drac(i,jsat,k)
              sum2rac(i,jsat) = sum2rac(i,jsat) + drac(i,jsat,k)**2
              sum23D(jsat)    = sum23D(jsat)    + omc(i,jsat,k)**2
            enddo
          endif
        enddo
        do i=1,3
          rmsxyz(i,jsat) = dsqrt(sum2xyz(i,jsat)/dble(nepoch))
          rmsrac(i,jsat) = dsqrt(sum2rac(i,jsat)/dble(nepoch))
        enddo
        rms3D(jsat) = sqrt(sum23D(jsat)/nepoch)
c       total rms 
        rmstot(jsat) = min( big
     .                  ,dsqrt(dot(rmsxyz(1,jsat),rmsxyz(1,jsat))/3))
c       write rms values for each SV to file
        write(lurms,'(i3,8(1x,f9.5))') isat(jsat),rmstot(jsat)*km_m,
     .       (rmsxyz(i,jsat)*km_m,i=1,3),(rmsrac(i,jsat)*km_m,i=1,3),
     .        rms3D(jsat)*km_m

cd        print *,'WRITE_SUMMARY nsat jsat isat ',nsat,jsat,isat(jsat)
c       Check if any SV have bad misfits. If so identify the worst.
        if ( max_fit_tol .gt. 0.0 .and. 
     .     rmstot(jsat)*km_m .gt. max_fit_tol .and.
     .     rmstot(jsat)*km_m .gt. bad_rmstot ) then
          jbad_sat = isat(jsat) 
          bad_rmstot = rmstot(jsat)*km_m
cd          print *,'jbad_sat ',jbad_sat
        endif
      enddo
       
c     Write the mean values at the bottom of the rms file
                         
      do i=1,3
        mean_rmsxyz(i) = 0.d0
        mean_rmsrac(i) = 0.d0
      enddo
      mean_rms3D  = 0.0d0
      mean_rmstot = 0.d0

      do jsat = 1,nsat
        do i = 1,3
           mean_rmsxyz(i) = mean_rmsxyz(i) + rmsxyz(i,jsat)**2
           mean_rmsrac(i) = mean_rmsrac(i) + rmsrac(i,jsat)**2
        enddo
        mean_rmstot = mean_rmstot + rmstot(jsat)**2
        mean_rms3D  = mean_rms3D  + rms3D(jsat)**2
      enddo
      do i = 1,3
        mean_rmsxyz(i) = dsqrt(mean_rmsxyz(i)/dble(nsat))*km_m
        mean_rmsrac(i) = dsqrt(mean_rmsrac(i)/dble(nsat))*km_m
      enddo
      mean_rmstot = dsqrt(mean_rmstot/dble(nsat))*km_m  
      mean_rms3D  = sqrt(mean_rms3D/nsat)*km_m

      write(lurms,'(2a)') '============================================'
     .           ,'======================================='
      write(lurms,'(a,8(f9.5,1x),/)') 
     .      'MEAN',mean_rmstot,mean_rmsxyz,mean_rmsrac, mean_rms3D
c     write the total in the STATUS file for sh_gamit grep'ing
      write(message,'(a,f9.5)') 'Overall fit (rms) to external orbit ='
     .                          ,mean_rmstot
      call report_stat('STATUS','ORBFIT','orbits/write_summary',' '
     .                , message,0)
c    Does the solution needs iteration.
      if ( jbad_sat .ne. 0 ) then
        nbad_sat = nbad_sat + 1
        ibad_sat(nbad_sat) = jbad_sat
        write(message,'(a,i3,a,f9.4,a,f6.4)') 'Satellite: '
     .  ,ibad_sat(nbad_sat),' Misfit: ',bad_rmstot
     .  ,' Larger than allowable max_fit_tol: ',max_fit_tol
        call report_stat('STATUS','ORBFIT','orbits/write_summary',' '
     .                , message,0) 
* MOD TAH 200121: Write out only the good solution (rewind lurms) 
!       write(lurms,'(a)') message
        rewind(lurms)
        return
      endif

c     Write the new G-file
 
      if( nparam.gt.0 ) then 
cd       print *,'WRITE_SUMMARY calling write_g nsat isat satnam 29-32'
cd     .      ,nsat,(isat(i),i=29,32),(satnam(i),i=29,32)
        call write_g
      endif
     
      return
      end
