CTITLE REMAP_RX

      subroutine remap_rx

      implicit none

*     Routine to remap phase data to common reference frequencies.
*     Remapping is needed for Glonass and may also work when GPS
*     and other CDMA GNSS systems are combined.

      include '../includes/xfile_def.h' 
      include 'track_com.h'
      include '../../libraries/includes/freq_def.h'   ! GAMIT Frequency table

*


* LOCAL VARIABLES
      integer*4 ns   ! Station loop
     .,         ep   ! Epoch loop
     .,          j   ! Channel loop

      integer*4 pn   ! PRN number
     .,         lv   ! Satellite number from list of PRNs

      integer*4 PtoL   ! Function to return list number given the PRN


*     Set the reference frequencies based on the GNSS systems to 
*     be processed.
      if ( index(tr_gnss,'G').gt.0 ) then
         fR1 = gps_f1
         fR2 = gps_f2
      elseif ( index(tr_gnss,'R').gt.0 ) then
         fR1 = glonass_f1
         fR2 = glonass_f2
      elseif ( index(tr_gnss,'E').gt.0 ) then
         fR1 = galileo_f1
         fR2 = galileo_f5
      elseif ( index(tr_gnss,'C').gt.0 ) then
         fR1 = beidou_f2
         fR2 = beidou_f7
      else
         write(*,120) trim(tr_gnss) 
 120     format('** DISASTER ** Unknown GNSS system requesed. ',a,
     .          ' input with TR_GNSS command')
         stop 'TRACK: Unknown GNSS'
      endif

***** Set the new frequencies ratios for MW-WL, EX-WL, LC etc
      call set_combs

****  Now cycle through all the data, remapping the phase to new reference
*     frequencies.  Non-integer ambiguities will be accounted for when 
*     ambiguities are resolved.
      do ns = 1, num_site
         do ep = 1, num_epochs
            do j = 1, num_chan_se(ns,ep)
*              Get the satellite 
               pn = ctop_cse(j,ns,ep)
               lv = PtoL(pn)
               if ( lv.gt.0 ) then
               L1o_all_cse(j,ns,ep) = L1o_all_cse(j,ns,ep)*fR1/fL1(lv)
               L2o_all_cse(j,ns,ep) = L2o_all_cse(j,ns,ep)*fR2/fL2(lv)
               endif
            end do
         end do
      end do

***** Thats all
      return
      end




 

