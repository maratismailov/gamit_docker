CTITLE CHECK_EST
 
      subroutine check_est( glb_parn, ndimg, nc, num,
     .           loc_parn, ndiml, estimated, ndime, num_est )

      implicit none  
 
*     Routine to scan the parmeter number in glb_parn, and to copy
*     the value to loc_parn if the parameter number appears in
*     the array estimated.
 
 
*   estimated(ndime,1)  - list of parameter numbers which
*                       - we want to copy
*   glb_parn(ndimg,num) - Parameter number from the global
*                       - solution
*   i,j,k               - Loop counters
*   loc_parn(ndiml,1)   - Parmeter number array to be created
*                       - if parameter is in the estimated list
*   ndimg       - First dimension of glb_parn
*   nc          - Number of elements in glb_parn to check
*   ndiml       - First dimension of loc_parn
*   ndime       - First dimension of estimated
*   num         - Number of items in glb_parn
*   num_est     - Number of items in estimated list
 
      integer*4 ndime, ndimg, num, ndiml   
      integer*4 estimated(ndime,*), glb_parn(ndimg,num), i,j,k,
     .    loc_parn(ndiml,*), nc, num_est
 
****  First clear the loc_parn array
 
      call clear_loc_parn( loc_parn, ndiml*num )
 
***** Now loop over values in glb_parn, seeing if we can find them in
*     estimated
 
      do i = 1, num
          do j = 1, nc    
 
*             Scan estimated
              do k = 1, num_est
                  if( estimated(1,k).eq.glb_parn(j,i) ) then
                      loc_parn(j,i) = estimated(1,k)
                  end if
              end do
          end do
      end do
 
***** Thats all
      return
      end
 
