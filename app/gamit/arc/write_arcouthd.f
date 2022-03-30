      Subroutine write_arcouthd

c     Write the satellite-independent headers of the print file
c     R. King 21 May 2012 from code in arc.

      implicit none

      include '../includes/dimpar.h'  
      include '../includes/global.h'
      include '../includes/arc.h'
  
      write(iarh,'(/,a)')
     .  ' Time Type   Integration Interval  Tabular Interval'
      write(iarh,'(a)')
     .  ' ---------   --------------------  ----------------'
      write(iarh,'(1x,a5,19x,f8.4,12x,f6.1)')
     .   time_type,    diint,                delt

 
      write(iarh,'(/,a)')
     .  ' Reference Frame       Precession  Nutation  Gravity Model'
      write(iarh,'(a)') 
     .  ' --------------------  ----------  --------  -------------'
      write(iarh,'(1x,a20,7x,a5,5x,a5,10x,a5)')
     .   frame_name,            precmod,    nutmod,   gravmod
                                               

      write(iarh,'(/,a)') 
     .   ' Solar Rad Model  Earth Rad Model  Antenna Rad Model'
      write(iarh,'(a)')
     .   ' ---------------  ---------------  -----------------'
      write(iarh,'(11x,a5,13x,a5,14x,a5)')
     .     srpmod,          eradmod,         antradmod
      
      return
      end
