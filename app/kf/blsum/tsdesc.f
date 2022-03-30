CTITLE TSDESC

      subroutine tsdesc( unit )

      implicit none

*     Routine to write out the extended header to the tsfiles.
*

      integer*4 unit   ! Unit number to be used.

      integer*4 i      ! Loop counter

      character*128 SFD(28)   ! 28 field descriptor lines
      character*256 Header    ! Header line immediately after End Field Descriptor.


      SFD( 1) = 'Start Field Description'
      SFD( 2) = 'YYYYMMDD      Year, month, day for the ' //
     .          'given position epoch'
      SFD( 3) = 'HHMMSS        Hour, minute, second for ' //
     .          'the given position epoch'
      SFD( 4) = 'JJJJJ.JJJJJ   Modified Julian day for the' // 
     .          ' given position epoch'
      SFD( 5) = 'X             X coordinate, Specified Reference' //
     .          ' Frame, meters'
      SFD( 6) = 'Y             Y coordinate, Specified Reference' //
     .          ' Frame, meters'
      SFD( 7) = 'Z             Z coordinate, Specified Reference' //
     .          ' Frame, meters'
      SFD( 8) = 'Sx            Standard deviation of the X' //
     .          ' position, meters'
      SFD( 9) = 'Sy            Standard deviation of the Y' // 
     .          ' position, meters'
      SFD(10) = 'Sz            Standard deviation of the Z' //
     .          ' position, meters'
      SFD(11) = 'Rxy           Correlation of the X and Y position'
      SFD(12) = 'Rxz           Correlation of the X and Z position'
      SFD(13) = 'Ryz           Correlation of the Y and Z position'
      SFD(14) = 'Nlat          North latitude, WGS-84 ellipsoid,' //
     .          ' decimal degrees'
      SFD(15) = 'Elong         East longitude, WGS-84 ellipsoid,' //
     .          ' decimal degrees'
      SFD(16) = 'Height (Up)   Height relative to WGS-84 ellipsoid, m'
      SFD(17) = 'dN            Difference in North component from ' //
     .          'NEU reference position, meters'
      SFD(18) = 'dE            Difference in East component from' //
     .          ' NEU reference position, meters'
      SFD(19) = 'du            Difference in vertical component ' //
     .          'from NEU reference position, meters'
      SFD(20) = 'Sn            Standard deviation of dN, meters'
      SFD(21) = 'Se            Standard deviation of dE, meters'
      SFD(22) = 'Su            Standard deviation of dU, meters'
      SFD(23) = 'Rne           Correlation of dN and dE'
      SFD(24) = 'Rnu           Correlation of dN and dU'
      SFD(25) = 'Reu           Correlation of dEand dU'
      SFD(26) = 'Soln          "rapid", "final", "suppl/suppf", ' //
     .      '"campd", or "repro" corresponding to products  generated'
      SFD(27) = ' with rapid or final orbit products, ' //
     .          'in supplemental processing, campaign data ' // 
     .          'processing or reprocessing'
      SFD(28) = 'End Field Description' 
      Header  = '*YYYYMMDD HHMMSS JJJJJ.JJJJ         X             Y '//
     .          '            Z            Sx        Sy       Sz    ' //
     .          ' Rxy   Rxz    Ryz            NLat         Elong   ' //
     .          '      Height         dN        dE        dU       ' //
     .          '  Sn       Se       Su      Rne    Rnu    Reu  Soln'
       
****  Write out the lines
      do i = 1, 25
         write(unit,'(a)') trim(SFD(i))
      end do
      write(unit,'(a,a)') trim(SFD(26)), trim(SFD(27))
      write(unit,'(a)')   trim(SFD(28))

      write(unit,'(a)') trim(Header)

***** Thats all
      return
      end
