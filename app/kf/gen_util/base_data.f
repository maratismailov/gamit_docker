ctitle
 
      integer*4 function base_data(data_type, num_sites)

      implicit none 
 
c
c     Routine to compute the maximum number of observables possible
c     at each epoch of an experiment.
c
c     RESTRICTION: ONLY ONE DELAY TYPE CAN BE TURNED ON AT ONE TIME:
c
c Variables
c ---------
c base_data -- the maximum number of observables per epoch for
c     this SOLVK run
c data_type -- a bit masked integer with specification of data to be
c     used. Bit 1 -- group delays;
c               2 -- phase delays;
c               3 -- SB delays;
c               4 -- delay rates
c     CURRENTLY; ONLY ONE DELAY TYPE CAN BE TURNED ON AT ONE TIME:
c num_sites   -- number of sites in this experiment.
c
*   nstep       - Number of data types per epoch in this
*               - solution
      integer*4 data_type, num_sites, nstep
 
c
c.... compute maximum number observables for each epoch
c     see if delay being used.
*     Use the NUM_STEP subroutine to get the data per epoch.
 
      call num_step( data_type, nstep )
      base_data = (num_sites*(num_sites-1)/2)* nstep
 
c.... Thats all
      return
      end
 
c........................................................................
