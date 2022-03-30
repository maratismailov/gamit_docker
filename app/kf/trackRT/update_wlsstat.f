      subroutine update_wlsstat(ep, iNL1, iNL2, na)

      implicit none

*     Routine to update the widelane sums when the the
*     number of integer ambiquities is changes


      include '../includes/const_param.h'
      include '../../libraries/includes/const_freq.h'
      include 'trackRT.h'         ! Common block

* PASSED 
      integer*4 ep            ! Counter for number of epochs of data processed.
      integer*4 iNL1, iNL2    ! Change in L1 and L2 cycles
      integer*4 na            ! Ambiquity being changed

* LOCAL
      real*8 dmw, dex  ! Change on MW and EX widelanes

****  Update the sums
      dmw = iNL1-iNL2
      dex = iNL1-iNL2*fL1/fL2
      wls_sqr(1,na) = wls_sqr(1,na) 
     .              - 2*dmw*wls_sum(1,na)
     .              + dmw**2*wls_num(na)
      wls_sqr(2,na) = wls_sqr(2,na) 
     .              - 2*dex*wls_sum(2,na)
     .              + dex**2*wls_num(na)

      wls_sum(1,na) = wls_sum(1,na) - dmw*wls_num(na)
      wls_sum(2,na) = wls_sum(2,na) - dex*wls_num(na)

****  Update the current expected values for the curr_sd ex and mw 
*     wide lanes
      curr_sdmw(na) = curr_sdmw(na) - dmw
      curr_sdex(na) = curr_sdex(na) - dex


      if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
     .write(*,220) ep, na, iNL1, iNL2, curr_sdmw(na), curr_sdex(na),
     .             ambiq_all(1,na),  ambiq_all(2,na)
 220  format('WLFIXUPD ',i6,' NA ', i3,' iNL1L2 ',2I5,
     .       ' WideLanes ',2F13.3,' L1L2 Cycles ',2F13.1)

      return
      end


      subroutine acc_wlsstat(ep, dmw, dex, na, opt)

      implicit none

*     Routine to accumulate the Wide statistic 

      include 'trackRT.h'         ! Common block

* PASSED 
      integer*4 ep            ! Counter for number of epochs of data processed.
      integer*4 na            ! Ambiquity being changed

      real*8 dmw, dex  ! Change on MW and EX widelanes

      character*(*) opt

****  Sum the statistics
      if( opt(1:1).eq.'R' ) then   ! Re-set, this is new estimate starting
         wls_sum(1,na) = dmw
         wls_sum(2,na) = dex
         wls_sqr(1,na) = dmw**2
         wls_sqr(2,na) = dex**2
         wls_num(na)   = 1
      else                         ! Normal sum
         wls_sum(1,na) = wls_sum(1,na) + dmw
         wls_sum(2,na) = wls_sum(2,na) + dex
         wls_sqr(1,na) = wls_sqr(1,na) + dmw**2
         wls_sqr(2,na) = wls_sqr(2,na) + dex**2
         wls_num(na)   = wls_num(na)   + 1
      end if

      return
      end


      subroutine remapdd( ep, old_sv, new_sv ) 

      implicit none

*     Routine to map the sums from one satellite to another
*     in accumulating the statistics of the double difference
*     widelanes.

      include 'trackRT.h'         ! Common block

* PASSED 
      integer*4 ep              ! Counter for number of epochs of data processed.
      integer*4 old_sv(max_site), new_sv(max_site)  ! Old and new reference satellites

* LOCAL
      integer*4 i,j   ! Counters
      integer*4 new_na  ! Ambiquity number for new SV

      real*8 mean(2)    ! Mean MW and EX widelanes at new SV
      real*8 omean(2)   ! Old mean value before update



****  The process is not quite rigorous here:  We compute the mean
*     value of the new-sv relative to the old-sv and then apply
*     mean value to all the effected differnces

****  Loop over stations to do this
      do i = 2, num_site
*        Find the ambiquity number do the new-sv minus old-sv
         new_na = 0
         do j = 1, num_ambs
            if( bf_ents(1,j).eq.i .and. bf_ents(2,j).eq.new_sv(i) .and.
     .          new_sv(i).ne.old_sv(i) ) then
                new_na = j
                if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
     .          write(*,140) site_names(i), old_sv(i), new_sv(i), ep
 140            format('REF PRN Switch at ',a,' from PRN ',I2.2,
     .                 ' to PRN ',I2.2,' at Ep ',I6)
               exit
            endif
        end do

****    Now compute mean value and remove from other ambiquities
   
        if( new_na.gt.0 ) then
            if( wls_num(new_na).gt.0 ) then
                mean(1) = wls_sum(1,new_na)/wls_num(new_na)
                mean(2) = wls_sum(2,new_na)/wls_num(new_na)
            else
                mean(1) = 0
                mean(2) = 0
            endif
*           Now map everyone over to the new ambiquity
            do j = 1, num_ambs
               if( bf_ents(1,j).eq.i ) then
*                 Correct station, so remove this mean
*                 from the data.  If this is the original
*                 ref_svs, replace it with new svs values
                  if( bf_ents(2,j).ne.old_sv(i) .and. 
     .                bf_ents(2,j).ne.new_sv(i) ) then
C                     print *,'SWITCHB ',ep, j, new_na, mean,  
C    .                   wls_sum(:,j), wls_sqr(:,j),  wls_num(j) 
                      omean(1) = wls_sum(1,j)/wls_num(j)
                      omean(2) = wls_sum(2,j)/wls_num(j)

                      wls_sqr(1,j) = wls_sqr(1,j) 
     .                             - 2*wls_sum(1,j)*mean(1)
     .                             +  mean(1)**2*wls_num(j)
                      wls_sqr(2,j) = wls_sqr(2,j) 
     .                             - 2*wls_sum(2,j)*mean(2)
     .                             +  mean(2)**2*wls_num(j)

                      wls_sum(1,j) = wls_sum(1,j) - mean(1)*wls_num(j)
                      wls_sum(2,j) = wls_sum(2,j) - mean(2)*wls_num(j)
C                     print *,'SWITCHA ',ep, j, new_na, mean, 
C     .                   wls_sum(:,j), wls_sqr(:,j),  wls_num(j) 
                  elseif( bf_ents(2,j).eq.old_sv(i) ) then
*                     Flip sign and replace
                      wls_sqr(1,j) = wls_sqr(1,new_na) 
                      wls_sqr(2,j) = wls_sqr(2,new_na) 
                      wls_sum(1,j) = -wls_sum(1,new_na) 
                      wls_sum(2,j) = -wls_sum(2,new_na) 
                      wls_num(j)   = wls_num(new_na)
                  end if
               end if
            end do
*           Now reset the new_na values
            wls_sqr(1,new_na) = 0
            wls_sqr(2,new_na) = 0
            wls_sum(1,new_na) = 0 
            wls_sum(2,new_na) = 0 
            wls_num(new_na)   = 0
         elseif( new_sv(i).ne.old_sv(i) ) then
            write(*,420) ep, new_sv(i), old_sv(i), site_names(i)
 420        format('**WARNING** No ambiquity at epoch ',i6,
     .             ' New SV ',i2,' Old SV ',i2,' Site ',a)
         endif
      end do

***** Thats all
      return
      end

      subroutine res_wlsRT( ddmw, ddex, eNL1, eNL2, iNl1, iNL2)

      implicit none

*     Subroutine to get estimates of L1 and L2 numbers of cycles from
*     MW and EX widelane (assimption is EX should be zero).


* PASSED 
      real*8 ddmw, ddex  ! MW and EX widelanes
      real*8 eNL1, eNL2  ! Estimates of NL1 and NL2

      integer*4 iNL1, iNL2  ! Integer values

* LOCAL
      integer*4 L1, L2      ! Loop test values
      real*8  rmsmin        ! Min RMS
      real*8  wgh           ! Wgh to MW-MW
      real*8  errmw, errex  ! Error in MW and EX widelines

****  Get estimate of number of L1 and L2 cycles
      wgh = 0.25
      eNL1 = 77*ddmw/17 - 60*ddex/17
      eNL2 = 60*ddmw/17 - 60*ddex/17

*     Now try to match to to match better with small
*     search
      rmsmin = 1.e20
      do L1 = nint(eNL1)-3, nint(eNL1)+3
         do L2 = nint(eNL2)-2, nint(eNL2)+2  
            errmw = (ddmw-(L1-L2))**2
            errex = (ddex-(L1-L2*77.d0/60.d0))**2
            if( errmw*wgh+errex.lt.rmsmin ) then
                rmsmin = errmw*wgh+errex
                iNL1 = L1
                iNL2 = L2
            end if
         end do
      end do

****  That all
      return
      end

 
                
      subroutine rep_wlsstat(ep)

      implicit none

*     Subroutine to update WLS ambiquity and L1/L2 ambqiquity estimates

      include '../includes/const_param.h'
      include 'trackRTObs.h'
      include 'trackRT.h'         ! Common block

* PASSED
      integer*4 ep   ! Epoch counter

* LOCAL
      integer*4 i,j   ! Loop counters

      real*8 mean(2), rms(2)  ! Mean and RMS of WLs
      real*8 eNL1, eNL2  ! Estimates of NL1 and NL2
      integer*4 iNL1, iNL2  ! Integer values


***** Loop over all the ambiquities
      do i = 1, num_ambs
         if( wls_num(i).gt.0 ) then
            do j = 1,2
               mean(j) = WLS_sum(j,i)/WLS_num(i)
               rms(j) = sqrt(abs(WLS_sqr(j,i)-mean(j)**2*WLS_num(i))/
     .              max(1,WLS_num(i)-1))
            end do
            call res_wlsRT( mean(1), mean(2), eNL1, eNL2, iNl1, iNL2)
            if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
     .      write(*,220) i, ep, site_names(bf_ents(1,i)), 
     .           bf_ents(2:5,i), mean, rms, WLS_num(i),
     .           eNL1, eNL2, iNl1, iNL2,  eNl1-iNL1, eNl2-iNL2,
     .           ambiq_all(:,i)
 220        format('BF_AMB ',i4,' Ep ',i6,1x,a4,1x,' PRN ',I2.2,3x,
     .        ' Range ',2I6,' FX ',o3,'  Means ',2F8.2,
     .        ' RMS ',2F6.2,' # ',I6,' eN12 ',2F8.2,
     .        ' iN12 ',2I4,' Err ',2F6.2,1x,2F14.1)
C        else
C           if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
C    .      write(*,240) i, ep, site_names(bf_ents(1,i)), 
C    .           bf_ents(2:6,i), WLS_num(i), ambiq_all(:,i)
C240        format('BF_AMB ',i4,' Ep ',i6,1x,a4,1x,' PRN ',I2.2,3x,
C    .        ' Range ',2I6,' FX ',o3,' Last ',I6,' WLS_num ',i2,
C    .        ' Ambiq ',2F14.1)

         end if
      end do

***** Thats all
      return
      end

         
CTITLE DEC_WLSTAT 

      subroutine dec_wlstat(ep, ndd, na)

      implicit none

*     Routine to remove the contribution of the any epoch
*     from the MW and EX statistics when data are edited 

      include '../includes/const_param.h'
      include '../../libraries/includes/const_freq.h'
      include 'trackRTObs.h'
      include 'trackRT.h'         ! Common block

* PASSED 
      integer*4 ep            ! Counter for number of epochs of data processed.
      integer*4 na            ! Ambiquity being changed
      integer*4 ndd           ! Double difference number

* LOCAL
      integer*4 j, js      ! Counters

      real*8 ddmw, ddex    ! Change on MW and EX widelanes
      real*8 mw(4), ex(4)  ! One-way values
      real*8 mean(2), rms(2)  ! For debug test

****  Form up the double difference wide lanes
      do j = 1,4
         js = wls_obn(j,na)
         if( js.gt.0 ) then
            mw(j) = (RT_omcraw(1,js)- RT_omcraw(2,js)-
     ,             dfsf*(RT_omcraw(3,js)*fL1+RT_omcraw(4,js)*fL2)/
     .                                              vel_light) 
            ex(j) = (RT_omcraw(1,js)- RT_omcraw(2,js)*(fL1/fL2))
         else
            write(*,120) na, j, wls_obn(:,na)
 120        format('ERROR: Amb ',i4,' Ent ',i2,' WLS_OBN ',4i4,
     .             ' No DD but being edited?')
            mw(j) = 0.0
            ex(j) = 0.0
         endif
      end do

****  Now form DD version
      ddmw = (mw(1)-mw(2))-(mw(3)-mw(4))
      ddex = (ex(1)-ex(2))-(ex(3)-ex(4))


      wls_sum(1,na) = wls_sum(1,na) - ddmw
      wls_sum(2,na) = wls_sum(2,na) - ddex
      wls_sqr(1,na) = wls_sqr(1,na) - ddmw**2
      wls_sqr(2,na) = wls_sqr(2,na) - ddex**2
      wls_num(na)   = wls_num(na)   - 1
      bf_ents(4,na) = bf_ents(6,na) ! Copy last valid epoch back
      if( ep.ne.0 ) then

         do j = 1,2
            mean(j) = WLS_sum(j,na)/WLS_num(na)
            rms(j) = sqrt(abs(WLS_sqr(j,na)-mean(j)**2*WLS_num(na))/
     .           max(1,WLS_num(na)-1))
         end do
         if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
     .   write(*,220) na, ep, site_names(bf_ents(1,na)), 
     .        bf_ents(2:5,na), mean, rms, WLS_num(na)
 220     format('BF_AMB ',i4,' Ep ',i6,1x,a4,1x,' G ',I2.2,
     .        ' Range ',2I7,' Flg ',o4,' Means ',2F8.2,
     .        ' RMS ',2F6.2,' # ',I6)
      end if

***** Thats all
      return
      end

                


            

