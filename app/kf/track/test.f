      program test

      implicit none


      include '../includes/xfile_def.h' 

      character*2 datatypes(xf_maxdat)

      datatypes(1) = 'L1'    
      datatypes(2) = 'L2'    
      datatypes(3) = 'C1'    
      datatypes(4) = 'P2'    
      datatypes(5) = 'P1'    
      datatypes(6) = 'S1'    
      datatypes(7) = 'S2'    
      datatypes(8) = 'C2'    
      datatypes(9) = 'L5'

      xf_ndat = 9 

      call assign_dattyp (datatypes )
      print *,'XF_DATTYP ',xf_dattyp(1:6,1)

      end
