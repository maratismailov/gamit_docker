CTITLE set_freqs 

      subroutine set_freqs(nsat)

      implicit none

*     GAMIT routine called by GETFIL and based on kf/ctogobs routine to assign 
*     the frequency quantities needed by CVIEW for each SV.  These are the same 
*     for all SVs except for GLONASS  

* It's not clear to me why the /clean routines pass 'nsat' (and also nobs, ncfls)
* rather than putting in common in cview.h where the PRN #s are stored, but I've
* left this scheme to avoid changing too many routines. --rwk 151229

* INCLUDES
                                                   
      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
           
*  Frequencies in cview.hread from C-files or set from the GNSS indicated on the x- or RINEX file
*     
      integer*4 nsat,i

      do i=1,nsat  
  
        gear(i) =  fL2(i)/fL1(i) 
        faclr(i) = (fL2(i)/fL1(i))/(1.d0 - (fL2(i)/fL1(i))**2) 
        facwl(i) = (fL1(i)-fL2(i))/(fL1(i)+fL2(i))  
c       this one in cview.h but apparently not used
c       facl2 (0.220779220779 for GPS)

      enddo      
 
cd      print *,'SET_FREQS fL1 fL2 faclr facwl '
cd     .       ,fL1(1),fL2(1),faclr(1),facwl(1)
cd      stop
                                                             
      return
      end   

                          



