CTITLE VAR_COMP
 
      subroutine var_comp(tran,covar, vars, temp_var, nres,ndata,
     .   option)
 

      implicit none 
 
c
c     routine to calculate the variances of a results of transformation
c     tran, using covariance matrix covar.  Results is returned in
c     vars.
c
c Variables
c ---------
c tran -- matrix transformation from data to results
c covar -- covariance matrix of the data
c vars -- the computed variances of the results
c temp_var -- temporary storage during calculations
c nres -- dimension of the results
c ndata -- dimension of the data.
c option -- determines if the complete covariance matrix will
c     will be computed (if option.ne.0 ) or if just variances
c     will be computed (option=0)
c
 
      integer*4 nres, ndata, option
 
*   i,j,k   - Loop counters
      integer*4 i,j
 
c
      real*8 tran(nres,1), covar(ndata,1), vars(nres,1),
     .    temp_var(ndata,1)
 
c
c     The routine operates by performing:
c
c                           T
c     vars = tran*covar*tran
c
c.... First do covar*tran(T) and save in temp_var
*                         ! row
      do i = 1, ndata
*                         ! column
         do j = 1, nres
c
c....       Dot the row and column of covar and tran
            call dvdot(temp_var(i,j),covar(1,i),1,tran(j,1),nres,ndata)
         end do
      end do
c
c.... Now complete the multiplication; see option to see if
c     we will compute matrix or just variances
      do i = 1,nres
c
*                                  ! just do variances
         if( option.eq.0 ) then
            call dvdot(vars(i,1),tran(i,1),nres,temp_var(1,i),1,ndata)
*                                  ! we will do whole covariance matrix
         else
            do j = 1,nres
               call dvdot(vars(i,j),tran(i,1),nres, temp_var(1,j),1,
     .            ndata)
            end do
         end if
*                 ! loop over over rows
      end do
c
      return
      end
 
CTITLE CHK_CORXYZ
      subroutine chk_corxyz(cov_xyz)

      implicit none 

*     Routine to set small correlations to zero exactly
      real*8 cov_xyz(3,3)
      integer*4 i,j
      real*8 rho

      do i = 1,3
         do j = 1,3
            if( i.ne.j ) then
                rho = cov_xyz(i,j)/sqrt(cov_xyz(i,i)*
     .                                  cov_xyz(j,j))
                if( abs(rho).lt.1.d-4) then
                    cov_xyz(i,j) = 0.d0
                end if
            end if
         end do
      end do

      return
      end
