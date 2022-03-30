CTITLE LAGRANGE_INTP
 
      subroutine Lagrange_intp( INF, PTS, epoch, value, dim)

      implicit none 
 
 
*--------------------------------------------------------------------
*     4 point Lagrange (equally spaced) interpolation routine.
*     Routine also computes the derivative of the function.
*     The routine will simualtaneosly interpolate several
*     difference type of data.  The number of types of data
*     is given by variable DIM.
*
*     The derivative is returned in units of units of the
*     tabluar points per unit of the spacing entries (usally days)
*
*     REF: Abramowitz, M., and I.A. Stegun, Handbook of mathematical
*         functions, Dover Publ. Inc., New York, 1972. Pages 879 and 883.
*         (pp1046)
*
*     T.Herring                   11:02 AM WED., 14 Jan., 1987
*--------------------------------------------------------------------
 
*   dim         - Number of types of data to interpolate.
*   i,j         - Loop counters
*   st_index    - the start index in the table to be used.  This
*               - is the position in PTS of the first point
*               - imediately before EPOCH.
 
 
      integer*4 dim, i,j, st_index
 
*   coeff(-1:2) - The 4 coefficients used for the interpolation
*   dt          - Difference in time bewteen the first data point
*               - to be used in the interpolation and the epoch
*               - In units of the spacing of the table.
*   epoch       - Epoch at which we want intepolated value
 
*   INF(3)      - Gives information about the spacing of
*               - the data points.  The values are
*               - 1 -- Julian date of first point.
*               - 2 -- Spacing of the data in days
*               - 3 -- Number of points in table
*   PTS(dim,1)  - The tabluar values of the quantity to be
*               - interolated.
 
*   value(dim,2)- Interpolated value and its derivative at
*               - 'epoch'.
 
      real*8 coeff(-1:2), dt, epoch, INF(3), PTS(dim,1), value(dim,2)
 
***** START, Determine the index of the first data point imediately
*     before the epoch, and get the time difference
 
      st_index = (epoch-INF(1))/INF(2) + 1
*                                               ! Make sure we use valid entry
      if( st_index-1.lt.1      ) st_index = 2
      if( st_index+2.gt.inf(3) ) st_index = inf(3) - 2
 
      dt       = (epoch - (INF(1)+(st_index-1)*INF(2)) )/INF(2)
 
***** Now compute the coefficients for the intepolation
 
      coeff(-1) = - dt*(dt-1)*(dt-2)/6
      coeff( 0) =  (dt**2-1)*(dt-2)/2
      coeff( 1) = - dt*(dt+1)*(dt-2)/2
      coeff( 2) =  (dt**2-1)*dt/6
 
*     Get interpolated value
 
*                         ! Loop over different types
      do i = 1,dim
 
          value(i,1) = 0.d0
          do j = -1,2
              value(i,1) = value(i,1) + coeff(j)*PTS(i,st_index+j)
          end do
      end do
 
***** Now do the derivative
 
      coeff(-1) = -(3*dt**2-6*dt+2)/6
      coeff( 0) =  (3*dt**2-4*dt-1)/2
      coeff( 1) = -(3*dt**2-2*dt-2)/2
      coeff( 2) =  (3*dt**2-1     )/6
 
*     Compute the derivative
 
      do i = 1, dim
 
          value(i,2) = 0.d0
          do j = -1,2
              value(i,2) = value(i,2) + coeff(j)*PTS(i,st_index+j)
          end do
 
*         Get into the correct units
*                                          ! divide by data spacing
          value(i,2) = value(i,2)/INF(2)
 
      end do
 
***** Thats all
      return
      end
 
