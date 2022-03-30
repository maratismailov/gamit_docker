
      program simsnr

      implicit none 
 
*     This program will fit SNR values and estimate corrections
*     to phase data based on the oscillations of the SNR.
 
      include 'simsnr.h'
      include '../includes/const_param.h'
      
*   ierr    - IOSTAT error
*   rcpar   - Reads the runstring
*   len_run - Length of the runstring read
*   win     - Lenght of Box-car smoothing window (in epochs)  (not used)
*   i,j     - Loop counters
 
 
      integer*4 ierr, rcpar, len_run, i,j, len_if
      
      character*256 line
 
****  Get the name of the input snr file (generated by svsnr)
      len_if = rcpar(1, infile)
      if( len_if.le.0 ) then
          call proper_runstring('simsnrr','simsnr.hlp',1)
      end if
 
      open(100,file=infile, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',infile,1,'simsnr')
      
      len_run = rcpar(2, line)
      if( len_run.gt.0 ) then
           read(line,*, iostat=ierr) num_rf,(h(i),r(i),an(i),
     .                               i=1,num_rf)
           do i = 1, num_rf
              an(i) = an(i)*pi/180
           end do
      else
           num_rf = 2
           h(1) = 10.7
           h(2) =  0.7
           r(1) =  0.25
           r(2) =  0.25
           an(1) = 0.d0
           an(2) = 0.d0
      end if
 
****  Tell user what is happen
      write(*,120) infile(1:len_if), num_rf
 120  format('* simsnr: Fitting SNR data in ',a,/,
     .       '* Reflectors: ',i4,/,
     .       '*   Height (L1 lambda)    R    Angle (deg)')
 
      do i = 1, num_rf
          write(*,140) i,h(i), r(i), an(i)*180/pi
 140      format('*',i3,F10.2,1x,F10.3,1x,f10.2)
      end do
      
****  Now read all of the input file
 
      call read_snr(100)
 
      write(*,220) num_epochs, num_sat
 220  format('* ',i5,' epochs of data with ',i3,' satellites found')
  
*     Clear the phase adjusments
      do i = 1, num_epochs
          do j = 1, num_sat
              phs_l1a(i,j) = 0.0
              phs_L2a(i,j) = 0.0
            
          end do
      end do
 
      call sim_elev
      
      end
      
 
CTITLE READ_SNR
 
      subroutine read_snr( unit )

      implicit none 
 
*     Routine to read the SNR file generated by svsnr.
 
      include 'simsnr.h'
      include '../includes/const_param.h'
 
* PASSED Variables
*         unit      - Unit number
 
      integer*4 unit
 
* LOCAL Variables
 
*   ierr, jerr      - IOSTAT errors.
*   i,j             - Loop counter
*   trimlen         - Length of string
*   chan            - Channel number for the prn read
*   prn             - PRN number
*   date(5)         - Date of measurement
*   code_type       - Indicates if X-correlation (1) or code
*                     correlation (0)

      integer*4 ierr, jerr, i, trimlen, chan, prn, date(5),
     .          code_type
      
*   sectag          - Seconds tage
*   sod             - Seconds of day (not used)
*   jdread          - Julian date computed from read date
*   lastjd          - Last jd read 
*   azread, elread  - Azimith and elevation angles read (degrees)
*   snr1, snr2      - Signal-to-noise ratios read from file

      real*8 sectag, sod, jdread, lastjd, azread, elread, 
     .       snr1, snr2
 
*   line            - Line read from file
 
      character*256 line
 
****  Start reading the file
      num_epochs = 0
      num_sat = 0
      lastjd = 0
      ierr = 0
 
      do while ( ierr.eq.0 )
 
          read(unit,'(a)', iostat=ierr ) line
          if( ierr.eq.0 .and. line(1:1).eq.' ' .and.
*                                             ! Process the line
     .        trimlen(line).gt.0 ) then
 
              read(line,120,iostat=jerr) date, sectag, sod, prn,
     .                azread, elread, snr1, snr2, code_type
 120          format(i5,4i3,F6.1,f10.1,1x,'PRN ',i2.2,1x, 
     .               2f10.3,1x,2f6.1,3x,i2)
 
              call ymdhms_to_mjd( date, sectag, jdread)
 
****          See if the epoch has changed
              if( abs(jdread-lastjd).gt.1.d0/86400.d0 ) then
                  num_epochs = num_epochs + 1
                  lastjd = jdread
                  epoch(num_epochs) = jdread
                  do i = 1, max_sat
                      flags(num_epochs,i) = 0
                  end do
              end if
 
****          Now see what channel the Prn is in
              call get_chan( prn, prn_list, num_sat, chan)
 
****          OK save the values in the appropriate locations
              az(num_epochs,chan) = azread*pi/180.d0
              el(num_epochs,chan) = elread*pi/180.d0
 
              snr_L1o(num_epochs,chan) = snr1
              snr_L2o(num_epochs,chan) = snr2
          
 
*             Set flags bit to show that we have an SNR measurement
*             at this time
              call sbit(flags(num_epochs,chan),1,1)
*             Set L1 bit to show that it is code measurement
              call sbit(flags(num_epochs,chan),2,1)
*             Check to see of code measurement at L2, if it is
*             set bit to show this.  (If not a code measurement then
*             L1 SNR gain will be removed before fiting the L2 gain curve).
              if( code_type.eq.0 ) then
                  call sbit(flags(num_epochs,chan),3,1)
              end if
          end if
      end do
 
****  Thats all
      return
      end
 
CTITLE GET_CHAN
 
      subroutine get_chan( prn, prn_list, num_sat, chan )
 
      implicit none 

*     This routine will return the chan number for a given PRN
*     or add the PRN to the list
 
*   prn     - PRN number to be matched
*   prn_list(*) - List of PRN's already found
*   num_sat     - Number of PRN's already found
*   chan        - The returned channel number
 
      integer*4 prn, prn_list(*), num_sat, chan
 
* LOCAL Variables
 
*         i     - Loop counter
 
      integer*4 i
 
****  See if we ready have
      i = 1
      do while ( prn_list(i).ne.prn .and. i.le.num_sat )
          i = i + 1
      end do
 
      if ( prn_list(i).ne. prn ) then
          num_sat = num_sat + 1
          prn_list(num_sat) = prn
          chan = num_sat
      else
          chan = i
      end if
 
****  Thats all
      return
      end
 
CTITLE sim_elev
 
      subroutine sim_elev

      implicit none 
 
*     This routine will fit a polynomial in sin(elev) to the L1 and
*     L2 SNR values
 
*     The basic L1 model is
*     SNR_L1 = A sin(el) + B sin(el)**2 + ..
*     For L2 (X-correlation) we first remove the L1 SNR model (scale to 1
*     at zenith) and then fit a second set of coefficients to sin(el).  For
*     the L2 code tracking channels an over-all scale factor is estimated.
 
      include 'simsnr.h'
      include '../includes/const_param.h'
 
* LOCAL Variables
 
*   i,j,k       - Loop counters
*   ipivot(max_poly)    - Pivot elements for inversion  (not used)
*   date(5)     - Date of measurement
*   code_type   - Code type (set to 0 for X-correlation, 1 for code
*               - correlation)
 
      integer*4 i,j,k, date(5), code_type, jerr
 
*       kbit        - Checks bit status
 
      logical kbit
 
*   se          - Sin(elevation)
*   apart(max_poly)     - Paritial derivatives for estimates (not used)
*   scale(max_poly)     - Scale factors for inversion (not used)
*   est_snr     - Estimated SNR from polynomial fit
*   est_gain        - est_snr divided by SNR at zenith (to get gain
*               - curves)  (not used)
*   sectag      - Seconds tag of date
*   sod         - Seconds of day.
 
 
 
      real*8 se,est_snr,
     .    sectag, sod, ce, l1_main, l2_main, 
     .    dL1(2), dL2(2), l1r, l1i, l2r, l2i, d1, d2
 
      
*****  Now compute the observed L1 gains, and set up L2 by removing L1 gain
      do i = 1, num_epochs
         call mjd_to_ymdhms( epoch(i), date, sectag)
         sod = date(4)*3600.d0 + date(5)*60.d0 + sectag
         do j = 1, num_sat
*                                                 ! SNR good
              if( kbit(flags(i,j),1) ) then
                  se = sin(el(i,j))
                  ce = cos(el(i,j))
                  est_snr = 0
                  L1_main = 100*se
                  L2_main = 100*se**2
*
*                 Now get the geometry for two reflectors
                  L1r = L1_main
                  L1i = 0.d0
                  L2r = L2_main
                  L2i = 0.d0
                  
                  do k = 1, num_rf
                     d1 = (h(k)/sin(el(i,j)+an(k)))*
     .                    (1.d0 - cos(2*(el(i,j)+an(k))))
                     d2 = d1*60/77
                     
                     dL1(1) = r(k)*L1_main*ce**4*cos(2*pi*d1+pi)
                     dL1(2) = r(k)*L1_main*ce**4*sin(2*pi*d1+pi)
                     dL2(1) = r(k)*L2_main*ce**4*cos(2*pi*d2+pi)
                     dL2(2) = r(k)*L2_main*ce**4*sin(2*pi*d2+pi)
                     
                     L1r = L1r + dL1(1)
                     L1i = L1i + dL1(2)
                     L2r = L2r + dL2(1)
                     L2i = L2i + dL2(2)
                     
                  end do

****              Now get total amp and phs
                  snr_L1o(i,j) = sqrt(L1r**2 + L1i**2)
                  phs_L1a(i,j) = atan2(l1i,L1r)/(2*pi)
                  
                  snr_L2o(i,j) = sqrt(L2r**2 + L2i**2)
                  phs_L2a(i,j) = atan2(l2i,L2r)/(2*pi)
                  
                  code_type = 1
                  if( kbit(flags(i,j),3) ) code_type = 0
                  write(*,120,iostat=jerr) date, sectag, sod, 
     .                prn_list(j), az(i,j)*180/pi, el(i,j)*180/pi, 
     .                snr_L1o(i,j), snr_L2o(i,j),  code_type,
     .                phs_L1a(i,j), phs_L2a(i,j)       
 120              format(i5,4i3,F6.1,f10.1,1x,'PRN ',i2.2,1x, 
     .               2f10.3,1x,2f6.1,3x,i2,2x,2F10.5)
                                                   
              end if
          end do
      end do
 
 
*     Thats all
      return
      end
 

 