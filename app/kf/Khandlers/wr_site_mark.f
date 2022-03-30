CTITLE WR_SITE_MARKV_APR
 
      subroutine wr_site_markv_apr(iout, site)
 

      implicit none
 
*     Routine to output the apriori markov statistics for each site
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
 
* Variables
* ---------
* iout -- the output device lu
* site -- the number of the site to be output
 
      integer*4 iout, site
 
*   i,j     - Loop counters
      integer*4 j
 
*
* local variables
* ---------------
* are_markov -- logical indicating if there markov values
* ifmar  -- logical function which is true is markov process on
 
 
      logical are_markov, ifmar
 
***** Firstly check to see if any markov elements for this site
      are_markov = .false.
      are_markov = are_markov .or. ifmar(clk_mar_apr(1,site),3)
      are_markov = are_markov .or. ifmar(atm_mar_apr(1,site),3)
      are_markov = are_markov .or. ifmar(atm_az_mar_apr(1,site),2)
 
***** If something Markov output
*                             ! output the values
      if( are_markov ) then
 
         write(iout,100) site_names(site)
 100     format(1x,a,' Markov process values')
 
*                                                 ! output values
         if( ifmar(clk_mar_apr(1,site),3) ) then
            write(iout,110) (clk_mar_apr(j,site),j=1,3)
 110        format(9x,' Clock process:         ',3(f8.5,1x))
         end if
 
         if( ifmar(atm_mar_apr(1,site),3) ) then
            write(iout,120) (atm_mar_apr(j,site),j=1,3)
 120        format(9x,' Atmosphere process:    ',3(f8.5,1x))
         end if
 
         if( ifmar(atm_az_mar_apr(1,site),2) ) then
            write(iout,130) (atm_az_mar_apr(j,site),j=1,2)
 130        format(9x,' Azimuthal atm process: ',2f8.5)
         end if
 
      end if
 
      end
 
