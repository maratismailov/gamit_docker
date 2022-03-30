CTITLE ASSIGN_DATTYP

      subroutine assign_dattyp( datatypes ) 

      implicit none

*     Routine to assign the mapping from the data read in the rinex file
*     to the xf_dattyp array.

      include '../includes/xfile_def.h' 

* Assignments:
* GPS       L1  L2   C1/P1   C2/P2   (L5/C5 posssible)
* Glonass   L1  L2   C1      C2
* Galileo   L1  L5   C1      C5      (L6/C6  L7/C7  L8/C8 possible)
* MOD TAH 200223: Replaced L7 with L6 and visa versa
* Beidou    L2  L6   C2      C6      (L7/C7 possible).
*

* PASSED
      character*(*) datatypes(xf_maxdat)

* LOCAL
      integer*4 i    ! Loop counter
      character*256 datastr

*     Convert the array into a string so we can use index
*     to get element location
      do i = 1,xf_ndat
         write(datastr((i-1)*2+1:i*2),'(a2)') datatypes(i)
      end do

      xf_dattyp = 0

***** Get the satellite system type
*     GPS 
         if( index(datastr,'L1').gt.0 ) 
     .       xf_dattyp(1,1) = (index(datastr,'L1')-1)/2+1
         if( index(datastr,'L2').gt.0 ) 
     .       xf_dattyp(2,1) = (index(datastr,'L2')-1)/2+1
         if( index(datastr,'P1').gt.0 ) 
     .       xf_dattyp(3,1) = (index(datastr,'P1')-1)/2+1
         if( index(datastr,'P2').gt.0 ) 
     .       xf_dattyp(4,1) = (index(datastr,'P2')-1)/2+1
         if( index(datastr,'C1').gt.0 ) 
     .       xf_dattyp(5,1) = (index(datastr,'C1')-1)/2+1
         if( index(datastr,'C2').gt.0 ) 
     .       xf_dattyp(6,1) = (index(datastr,'C2')-1)/2+1
*     Glonass
         if( index(datastr,'L1').gt.0 ) 
     .       xf_dattyp(1,2) = (index(datastr,'L1')-1)/2+1
         if( index(datastr,'L2').gt.0 ) 
     .       xf_dattyp(2,2) = (index(datastr,'L2')-1)/2+1
         if( index(datastr,'P1').gt.0 ) 
     .       xf_dattyp(3,2) = (index(datastr,'P1')-1)/2+1
         if( index(datastr,'P2').gt.0 ) 
     .       xf_dattyp(4,2) = (index(datastr,'P2')-1)/2+1
         if( index(datastr,'C1').gt.0 ) 
     .       xf_dattyp(5,2) = (index(datastr,'C1')-1)/2+1
         if( index(datastr,'C2').gt.0 ) 
     .       xf_dattyp(6,2) = (index(datastr,'C2')-1)/2+1
*     Galileo
         if( index(datastr,'L1').gt.0 ) 
     .       xf_dattyp(1,3) = (index(datastr,'L1')-1)/2+1
         if( index(datastr,'L5').gt.0 ) 
     .       xf_dattyp(2,3) = (index(datastr,'L5')-1)/2+1
         if( index(datastr,'P1').gt.0 ) 
     .       xf_dattyp(3,3) = (index(datastr,'P1')-1)/2+1
         if( index(datastr,'P5').gt.0 ) 
     .       xf_dattyp(4,3) = (index(datastr,'P5')-1)/2+1
         if( index(datastr,'C1').gt.0 ) 
     .       xf_dattyp(5,3) = (index(datastr,'C1')-1)/2+1
         if( index(datastr,'C5').gt.0 ) 
     .       xf_dattyp(6,3) = (index(datastr,'C5')-1)/2+1

*      Biedou
* MOD TAH 200223: Replaced L7/C7 with L6/C6 and assigned
*        slot 3/4 to C2/C6 (no P range values for Beidou(.
         if( index(datastr,'L2').gt.0 ) 
     .       xf_dattyp(1,4) = (index(datastr,'L2')-1)/2+1
         if( index(datastr,'L6').gt.0 ) 
     .       xf_dattyp(2,4) = (index(datastr,'L6')-1)/2+1
         if( index(datastr,'C2').gt.0 ) 
     .       xf_dattyp(3,4) = (index(datastr,'C2')-1)/2+1
         if( index(datastr,'C6').gt.0 ) 
     .       xf_dattyp(4,4) = (index(datastr,'C7')-1)/2+1


      end

