CTITLE  TIME_SPACE
 
      subroutine time_space(tic_space, indx, intv)
c
c     routine to convert the tic spacing into increments of
c     either years, months, days, hours or minutes
c
c Include files
c -------------
*                         ! the parameter file
      include 'plot_param.h'
c
*                         ! the plot common block
      include 'plot_com.h'
c
c Variables
c ---------
c tic_space -- the tic spacing desired
c indx -- the index in the tax_data (time_axis data) array for
c     the quanity
c intv -- the number of intervals to be spaced bewteen tic marks
c mult -- the number of intervals in the next highest unit spacing
c
      real*4 tic_space
 
c
      integer*4 indx, intv, mult
 
c
c.... See into which interval tic_space falls
      intv = int(tic_space/365.0)
      indx = 1
*                                   ! try month spacing
      if( intv.lt.1 ) then
         intv = int(tic_space/30.)
         indx = 2
*                                   ! try day spacing
         if( intv.lt.1 ) then
            intv = int(tic_space)
            indx = 3
*                                   ! try hour spacing
            if( intv.lt.1 ) then
               intv = int(tic_space*24.0)
               indx = 4
*                                    ! use minutes
               if( intv.lt.1 ) then
                  intv = int(tic_space*1440.0)
                  indx = 5
                                     ! use seconds
                  if(  intv.lt.1 ) then
                      intv = int(tic_space*86400.)                   
                      indx = 6
                  end if
               end if
            end if
         end if
      end if
c
c.... Now trim the spacing so that we will hit multiples of the
c     the spacing of the next higest unit
      if( indx.gt.1 .and. intv.gt.0 ) then
         mult = tax_data(1,indx)/intv
         intv = tax_data(1,indx)/mult
      end if
c
c MOD JLD 850111 To use half the number of tic marks do decrease crowding
c     for the hour plots
      if (indx .eq. 4) intv = intv * 2
c
      return
      end
 
