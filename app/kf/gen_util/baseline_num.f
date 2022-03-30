CTITLE 'BASELINE_NUMBER'
 
 
      integer*4 function baseline_number(i,j)

      implicit none 
 
*     J.L. Davis                   7:51 PM  TUE., 14  APR., 1987
 
*     Gives the baseline number for the site pair i, j
 
*          i,j          - Site numbers
*   ,   max_site        - Max of i,j
*   ,   min_site        - Min of i,j
 
      integer*4 i,j, max_site, min_site
 
      max_site = max(i,j)
      min_site = min(i,j)
 
      baseline_number = max_site * (max_site - 3) / 2 + min_site + 1
 
      end
 
