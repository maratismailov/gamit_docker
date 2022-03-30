CTITLE CHECK_TRANS
      subroutine check_trans( p1, p2, nc, nt, coeff )
 
      implicit none 

 
*     This routine forms up the coefficients and pointers for the
*     transission elements in the state transission matrix.  Only
*     the none zero, non-unit values are saved.  All off-diagonal
*     terms are assumed to be multiplied by time in years.  This
*     routine is called once for each pair of parameters to be checked.
*
*     Routine assumes that sucessive values of p2 are called for
*     the onw value of p1.
 
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   p1,p2   - Parameter numbers of the paramters being tested
*           - (if both are non zero then it means we need an
*           - element in the transission matrix.)
*   nc      - Pointer to current array index in trans_col and
*           - trans_coeff.  Trans_col gives the columns of the
*           - state transission matrix.
*   nt      - pointer to current array index in trans_row.
*           - Trans_row points to the rows of the state
*           - transission matrix, the array index to trans_col
*           - with the first entry for the column numbers for
*           - for this row, and the number of column numbers for
*           - this row.
 
      integer*4 p1,p2, nc, nt
 
*   coeff   - Coefficient for the state transission element.
*           - (assumed to be multipled by delta t when used, and
*           - one added if diagonal element.)
 
 
      real*4 coeff
 
***** If both of the parameters are being estimated
 
*                                         ! Yes, both are being
      if( p1.ne.0 .and. p2.ne.0 ) then
*                                         ! estimated, so there is a
*                                         ! transission element.
*                         ! Next column in the transsion matrix
          nc = nc + 1
          call check_glb_max('transission columns', nc, max_trans_col)
          trans_col(nc)      = p2
          trans_coeff(nc)    = coeff
 
*         See if we have changed row
          if( p1.ne.trans_row(1,nt) .or.
*                                             ! Yes we have changed row
     .        nt.eq.0               ) then
*                                             ! or first transission
*                                             ! so save current row,
*                                             ! pointer to first element
*                                             ! in trans_col.
              nt = nt+1
              call check_glb_max('transission rows', nt, max_trans_row)
 
              trans_row(1,nt) = p1
              trans_row(2,nt) = nc
          end if
 
*         Now increment number of columns in transission matrix for this
*         parameter
          trans_row(3,nt) = trans_row(3,nt) + 1
 
      end if
 
***** Thats all
      return
      end
 
