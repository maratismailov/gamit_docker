      program test_svcode
      implicit none

*     Test conversion code
      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'

      
      integer*4 qnum_svs   ! number of satellites to test
      integer*4 rcpar, lenrun, i,j , na
      integer*4 qparn(max_svs_elem)
      integer*4 old_code, new_code 
      integer*4 svcode_to_code

      character*128 runstring

      data qparn / 1,1,1, 1,1,1, 1,1,0,1,0, 1,1,1,1,1,1, 0,0,0, 1,1,1/

      cglb_vers = 106

      lenrun = rcpar(1,runstring) 
      read(runstring,*) qnum_svs
      write(*,110) qnum_svs
 110  format('Generating codes for ',i4,' satellites',/,
     .       'SVS  PARN  Old Code   New code')
      do i = 1, qnum_svs
         do j = 1,max_svs_elem
            if( qparn(j).ne.0 ) then
                na = na + 1
                old_code = 51 + j*256 + i*65536
                new_code = svcode_to_code( old_code) 
                write(*,120) i,j,old_code, new_code
 120            format(i4,1x,i3,1x,O10,1x,O10)
            end if
         enddo
         print *,' '
      enddo 

      end


      
