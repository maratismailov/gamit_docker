ctitle
 
      subroutine bak_file_map(ibepoch, imar_val, ipost_res,
     .   ibelev, ibak_recs, ibak_sou, ibak_unw, num_mar, data_type,
     .   num_site, bak_recl )

      implicit none 
 
c
c     This subroutine sets up the dynamic mapping of the contents
c     of the bak file.  This file is used for plotting the results
c     of the SOLVK solution.
c     The record length of bak_file is even so that we can use VIS
c     in the plotting routine to loop through the array
c
c
c Variables
c ---------
c ibepoch   -- the start index of bepoch
c imar_val  -- the start index of mar_val
c ipost_res -- the start index of post_res
c ibelev    -- the start index of belev
c ibak_recs -- the start index of bak_recs
c ibak_sou   -- the start index of bak_sou
c ibak_unw   -- the start index of bak_unw
c
c num_mar   -- the number of markov parameters
c data_type -- the bit masked word with data types used
c num_site  -- the number of sites in the experiment
c
c bak_recl -- he record length of the bak_file (depends on the
c     parameters estimated during the solution.)
c
      integer*4 ibepoch, imar_val, ipost_res, ibelev,
     .    ibak_recs, ibak_sou, ibak_unw, num_mar, data_type, num_site,
     .    bak_recl
 
c
c The layout of a record of BAK_FILE is:
c     bac_recl =
c   [ 2     + num_mar*2      + 2*nbl       +   num_site  + nbl      +
c     ^----^ ^--------------^ ^----------^  ^----------^  ^--------^
c     epoch   markov           residuals     elevation     unweight
c             (r*4 value and   (r*4 value    angles        flags
c              sigma)          and sigma)
c
c    + 1 + nbl   ]
c     ^---------^
c      source
c      and record
c      numbers
c
c     The residuals, elevation angles and record numbers are arranged
c     in station order.  If a station is not used at a particular epoch
c     its entries will be zero.
c
c Local variables
c ---------------
c nbl -- number of baseline data per epoch possible (one data per
c     baseline per data type)
c base_data -- function to compute nbl
c
      integer*4 nbl, base_data
 
c
c.... Get number of baseline data possible per epoch
      nbl = base_data(data_type,num_site)
c
c.... Set up the pointers
      ibepoch    = 1
      imar_val   = ibepoch    + 2
      ipost_res  = imar_val   + 2*num_mar
      ibelev     = ipost_res  + 2*nbl
      ibak_unw   = ibelev     + num_site
      ibak_sou   = ibak_unw   + nbl
      ibak_recs  = ibak_sou   + 1
c
      bak_recl   = ibak_recs  + nbl - 1
c
c.... Make the record length an even number
      bak_recl   = 2*int((bak_recl+1)/2)
c
c.... Thats all
      return
      end
 
c .....................................................................
 
