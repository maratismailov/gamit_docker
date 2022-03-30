CTITLE    .................................................................
 
      subroutine first_tic0(vmin, ref_val, first_tic, idate,
     .   intv,indx, tax_data)
c
c     Routine to compute the position of the first tic mark for
c     time axes.
c
c Variables
c ---------
c vmin -- the value of the min value on the axis
c ref_val -- the reference value for the minimum value
c first_tic -- position of first tick mark
c idate -- the date of first tic mark
c intv  -- the interval to be spaced on the axis
c indx  -- the indx of the lable type in tax_data
c tax_data -- an array of values which give the turn over values
c     for each type of label spacing, and the reset values when
c     a turn over ocurrs (see &plobd)
c
      real*4 vmin, first_tic
 
c
      real*8 ref_val
 
c
      integer*4 intv, indx, tax_data(2,6), idate(6)
 
c
c Local variables
c ---------------
c first -- the true of vmin
c
      real*8 first
 
c
c comjd -- computed julian date
c
      real*8 comjd
 
c
c Functions
c ---------
c fjldy_8 -- real*8 function to compute julian date
c
      real*8 sectag
      integer*4 i, iterm
 
c
c.... Compute time a mimimum
      first = vmin + ref_val 
      sectag = 0.0d0
c
C     call epoc_8(idate(2),idate(3),idate(1),idate(4),idate(5),
C    .   first)
      call mjd_to_ymdhms(first,idate,sectag)
      idate(6) = nint(sectag)
c
      iterm = 6
c

c
c.... Now loop by the spacing until we get inside the plot
!      comjd = 2378496.500d0
!      comjd = first - 3*365d0
      comjd = first
      if( indx.eq.1) then
         comjd = comjd - 365*intv
      else if( indx.eq.2) then 
         comjd = comjd - 30*intv
      else if( indx.eq.3) then 
         comjd = comjd - 1*intv
      else if( indx.eq.4) then
         comjd = comjd - 1./24*intv
      else if( indx.eq.5) then
         comjd = comjd - 1./1440.0*intv
      else if( indx.eq.6) then
         comjd = comjd - 1./86400.0*intv
      endif
      
            
      call mjd_to_ymdhms(comjd, idate, sectag)
      idate(6) = nint(sectag)
      
c.... Now reset the date so that we get a value before the first tic
      do i = indx+1,6
         idate(i) = tax_data(2,i)
      end do
c
      do while (comjd.lt.first)
         call inc_cax(idate, intv,indx, tax_data, comjd)
      end do
c      
      first_tic = comjd - ref_val

c
      return
      end
 
