      Subroutine QUATERNION

*     Compute the quaternion from the yaw angle
*     R. King 190927

*    (Not yet written)
                    
      implicit none

      include '../includes/dimpar.h'          
      include 'orbex.h'

      integer*4 i,j,iepcy 

* Function
      real*8 timdif 

      logical debug/.true./
                   
      print *,'QUATERINON noepoch nosat ',noepoch,nosat
      print *,' Set temporarily set to 0.5 '

*  Loop over the y-file values, selecting the span requested
        
      noepoch = 0 
      do iepcy = 1,nyepoch
       
        if( 
     .   (timdif(jdy(iepcy),ty(iepcy),obx_jdstart,obx_tstart).ge.-0.1d0)
     .       .and.
     .   (timdif(jdy(iepcy),ty(iepcy),obx_jdstop,obx_tstop).le.0.1d0) )
     .       then        
          noepoch = noepoch+ 1 
          jdobx(noepoch) = jdy(iepcy)
          tobx(noepoch) = ty(iepcy)   
          if(debug.and.iepcy.gt.240.and.iepcy.lt.243) then 
            print *,'QUATERNION iepcy noepoch obx_jdstart obx_tstart '
     .                ,iepcy,noepoch,obx_jdstart,obx_tstart
            print *,'   obx_jdstop obx_tstop ',obx_jdstop,obx_tstop
            print *,'   jdy ty ',jdy(iepcy),ty(iepcy)
          endif 
          do j=1,nosat
            do i=1,4
              quatern(i,j,noepoch) = 0.5d0
            enddo
          enddo
        endif                     

      enddo 
      if(debug) print *,'QUATERNION noepoch ',noepoch                      
      if(debug) print *,'QUATERNION 1,1,1 ',quatern(1,1,1)       

      return
      end 






