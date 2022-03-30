      subroutine keep_loc( parn, ndim, nc, ltog, cnum, gnum)

      implicit none 
 
*     This routine to scan the list of parameters in parn
*     which are arranged by site eg. and remove the elements
*     which are not for local quantities eg sites.  Nominally
*     meant for sites, sources and satellites.
*
 
*   ndim            - Dimension of first element of parn
*   nc              - nunber of entries to check
*   cnum            - Number of local stations/sources/svs
*   gnum            - number of global  ""      ""     ""
*   parn(ndim,gnum) - Parmaeter numbers arranged by site
*                   - source or sv
*   ltog(cnum)      - Local to global site relation
 
      integer*4 ndim, nc, cnum, gnum, parn(ndim,gnum), ltog(cnum)
 
*   i,j             - Counters
*   iel             - Simple pointer
 
 
      integer*4 i,j, iel
 
*     Start by looping over local sites and setting parn element
*     negative for these.  When this is finished all non-neg elements
*     are set tp zero.
 
      do i = 1, cnum
          iel = ltog(i)
          do j = 1,nc
              parn(j,iel) = -parn(j,iel)
          end do
      end do
 
*     Now that all local parameters are negative, set all non-neg
*     values to zero.
 
      do i = 1, gnum
          do j = 1, nc
              if( parn(j,i).ge.0 ) then 
                  parn(j,i) = 0
              else
                  parn(j,i) = -parn(j,i)
              end if
          end do
      end do
 
****  Thats all
      return
      end
 
 
