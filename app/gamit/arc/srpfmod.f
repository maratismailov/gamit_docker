      Subroutine srpfmod (antbody,sbmass,angb,srpxacc,srpyacc,srpzacc)

C Calculate accelerations due to solar radiation pressure along the spacecraft 
c X and Z axes using fourier coefficients supplied by Marek Ziebart (UCL).
C E Petrie Dec 2010
c R King  Oct 2014: Temporarily limit this to GPS and keep 'iblock' identifier; 
c               later generalize to use c*20 antbody for all GNSS
c R King  October 2015: Use c*20 antbody to allow all GNSS

C  In:    antbody: satellite body type 
C         sbmass : satellite mass
C         angb   : sun-satellite-earth angle (radians)
C
C  Out:  srpxacc: accelerations along the spacecraft x axis /km s-2
C	 srpyacc: accelerations along the spacecraft y axis /km s-2
C	 srpzacc: accelerations along the spacecraft z axis /km s-2

c
      implicit none
                     
      character*20 antbody
      integer*4 i,err
      real*8 sbmass,angb,srpxacc,srpyacc,srpzacc
C     Dimension to fifty - will be ok if only 17 coeffs.
      integer*4 nx(50),ny(50),nz(50),j,n(50)
      real*8 a0x,a0y,a0z,bx(50),ax(50),bz(50),ay(50),by(50),az(50)
      character*80 line
      logical first, first2
      data first /.true./ first2 /.true./
      save first,a0x,a0y,a0z,ax,ay,az,bx,by,bz,first2
      
      srpxacc=0.d0
      srpyacc=0.d0
      srpzacc=0.d0
cc      if (antiblock.eq.3) then      
      if( antbody(1:9).eq.'BLOCK IIA' ) then
C       calculate accels in km s-2 using UCL IIA coefficients and default satellite mass of 972.9kg
        call UCLSRPeditIIA (angb,srpxacc,srpyacc,srpzacc)
C            convert to satellite mass from svnav.dat (is 972.9kg, but proofs against changes)
             srpxacc=srpxacc*972.9d0/sbmass
             srpyacc=srpyacc*972.9d0/sbmass
             srpzacc=srpzacc*972.9d0/sbmass
cc      else if (iblock.eq.4.or.iblock.eq.5.or.iblock.eq.6) then 
      elseif( antbody(1:11).eq.'BLOCK IIR-A' .or.
     .        antbody(1:11).eq.'BLOCK IIR-B' .or. 
     .        antbody(1:11).eq.'BLOCK IIR-M' ) then 
C       Use the UCL coefficients up to 50, X, Y. Z body fixed refernce system.
C       These are officially for SVN 41, 43, 44, 46, 51 (IIR-A 2-6 as 1 failed on launch)
C        read in the coefficients only for the first time the routine is called for one of the block 4,5,6 satellites
C        Coefficients are based on a mass of 1100kg and mean solar irradiance of 1368W m-2
         if (first) then  
            first=.false.
C set up 
          a0x=0.d0
          a0y=0.d0
          a0z=0.d0
              do j=1,50
                 nx(j)=0.d0
                 ny(j)=0.d0
                 nz(j)=0.d0
                 ax(j)=0.d0
                 ay(j)=0.d0
                 az(j)=0.d0
                 bx(j)=0.d0
                 by(j)=0.d0
                 bz(j)=0.d0
              enddo
C         this read in code works but is not very failproof if the files are changed.
C          X coefficients
            open (unit = 5678, file = 
C     .      "/chandler/home/petrie/orbitrefs/data/MarekSRPfourierX.txt",
C     .      "/home/users/nejp5/gg/tables/MarekSRPfourierX.txt",
C           Files should be linked to day directory by sh_gamit
     .       "MarekSRPfourierX.txt",
     .      status='old',iostat=err)
C          Read six text header lines
            do i=1,6
              read (5678,*) line
C             print*,i, line
            enddo

C           Read a0
            read(5678,*) a0x
C           print*, a0x

            read (5678,*) line
C           print*, line
C           Read coefficients
            do j=1,50
              read(5678,*) nx(j),ax(j),bx(j)
C           print*, nx(j),ax(j),bx(j)
            enddo
            close (5678)
C           Print*,err,'iostat'
C           Print'(1x,i10,/)',nx

Cd           Print*,'ax'
Cd           Print*,ax
Cd           Print*,'bx'
Cd           Print*,bx
C       Y coefficients
           open (unit = 5678, file = 
C     .     "/chandler/home/petrie/orbitrefs/data/MarekSRPfourierY.txt",
C     .      "/home/users/nejp5/gg/tables/MarekSRPfourierY.txt",
C           Files should be linked to day directory by sh_gamit
     .      "MarekSRPfourierY.txt",
     .      status='old',iostat=err)
C          Read six text header lines
           do i=1,6
           read (5678,*) line
C           print*,i, line
           enddo

C          Read a0
           read(5678,*) a0y
C           print*, a0y

           read (5678,*) line
C           print*, line
C          Read coefficients
           do j=1,50
            read(5678,*) ny(j),ay(j),by(j)
C            print*, ny(j),ay(j),by(j)
           enddo
           close (5678)
C           Print*,err,'iostat'
C           Print'(1x,i10,/)',ny

C           Print*,'ay'
C           Print*,ay
C           Print*,'by'
C           Print*,by

C        Z coefficients
            open (unit = 5678, file = 
C     .      "/chandler/home/petrie/orbitrefs/data/MarekSRPfourierZ.txt",
C     .      "/home/users/nejp5/gg/tables/MarekSRPfourierZ.txt",
C           Files should be linked to day directory by sh_gamit
     .      "MarekSRPfourierZ.txt",
C     .      "testinput.txt",
     .       status='old',iostat=err)
C           Read six text header lines
            do i=1,6
            read (5678,*) line
C            print*,i, line
            enddo

C           Read a0
            read(5678,*) a0z
C            print*, a0z

            read (5678,*) line
C            print*, line
C           Read coefficients
            do j=1,50
             read(5678,*) nz(j),az(j),bz(j)
C             print*, nz(j),az(j),bz(j)
            enddo
            close (5678)
C            Print*,err,'iostat'
C            Print'(1x,i10,/)',nz
C           Print*,'ay'
C           Print*,ay
C           Print*,'by'
C           Print*,by
C            Print*,'az'
C            Print*,az
C            Print*,'bz'
C            Print*,bz
C           End of reading in the X,Y,Z fourier coefficients (first time only).
         endif
C Debug check save worked
C           Print*,'ax'
C           Print*,ax
C           Print*,'bx'
C           Print*,bx
C           Print*,'ay'
C           Print*,ay
C           Print*,'by'
C           Print*,by
C            Print*,'az'
C            Print*,az
C            Print*,'bz'
C            Print*,bz


C   Calculate the accelerations
C       X axis
        srpxacc=a0x/2.d0
        do i=1,50
         srpxacc=srpxacc+(ax(i)*dcos(float(i)*angb)+
     .                 bx(i)*dsin(float(i)*angb))
C        print*,'float(i)',float(i)
        enddo
C Y axis
        srpyacc=a0y/2.d0
        do i=1,50
         srpyacc=srpyacc+(ay(i)*dcos(float(i)*angb)+
     .                 by(i)*dsin(float(i)*angb))
C        print*,'float(i)',float(i)
        enddo
C Z axis
        srpzacc=a0z/2.d0
        do i=1,50
         srpzacc=srpzacc+(az(i)*dcos(float(i)*angb)+
     .                 bz(i)*dsin(float(i)*angb))
        enddo
C Convert to km s-2
C        Print*,'srpxacc',srpxacc
C        Print*,'srpzacc',srpzacc
        srpxacc=srpxacc*1.0d-3
        srpyacc=srpyacc*1.0d-3
        srpzacc=srpzacc*1.0d-3
C Allow for different input mass from svnav.dat
        srpxacc=srpxacc*1100.d0/sbmass
        srpyacc=srpyacc*1100.d0/sbmass
        srpzacc=srpzacc*1100.d0/sbmass
C        Print*,'srpxacc',srpxacc
C        Print*,'srpzacc',srpzacc
C        End of if for IIR-A, IIR-B, IIR-M)
c***RWK 151019: I can't find where iblock is ever 99 for setup,
c               so skip this code until we sort this out
cc      elseif (iblock.eq.99) then
CC Set up the fourier coefficients (nb, if use these in future need to check what 
C the default mass used for calculation was - have put 1100kg in the code)
CC Coefficients for IIR_B_X
C a0		
       a0x = 2.43338094794787d-10
C        b(n)			a(n)
        bx(1) = -1.07050014256139d-07
        ax(1) = 2.73060613405594d-10
        bx(2) = -3.13819577323511d-10
        ax(2) = -2.57874457104631d-10
        bx(3) = -2.25322507615340d-09
        ax(3) = -2.11475735358375d-10
        bx(4) = -1.94820171428050d-10
        ax(4) = 5.73115877541095d-11
        bx(5) = 1.38078042445690d-09
        ax(5) = 1.57085323611020d-11
        bx(6) = -1.25454403672251d-11
        ax(6) = 1.00292678688194d-10
        bx(7) = -9.23949007654900d-11
        ax(7) = -2.72124739682274d-11
        bx(8) = 4.38227256709194d-11
        ax(8) = -3.35736158987152d-11
        bx(9) = -7.21528755064295d-11
        ax(9) = 2.62913750621447d-11
        bx(10) = -2.48972568123080d-11
        ax(10) = -1.23004756680251d-11
        bx(11) = 4.97289063344050d-11
        ax(11) = -3.17027478788134d-11
        bx(12) = 3.44254300015349d-11
        ax(12) = 3.77010081520677d-11
        bx(13) = -1.32732584044218d-11
        ax(13) = 1.36705821493913d-11
        bx(14) = -6.13451559908228d-12
        ax(14) = 6.13224213274386d-12
        bx(15) = 5.39857538503087d-11
        ax(15) = 7.41190831968585d-12
        bx(16) = -1.59963590787475d-11
        ax(16) = -1.92499323577743d-11
        bx(17) = -4.92939194089531d-11
        ax(17) = -5.84155930264133d-12
C  n	 b(n)	                a(n)
C  1	-1.07050014256139E-07	2.73060613405594E-10
C  2	-3.13819577323511E-10	-2.57874457104631E-10
C  3	-2.25322507615340E-09	-2.11475735358375E-10
C  4	-1.94820171428050E-10	5.73115877541095E-11
C  5	1.38078042445690E-09	1.57085323611020E-11
C  6	-1.25454403672251E-11	1.00292678688194E-10
C  7	-9.23949007654900E-11	-2.72124739682274E-11
C  8	4.38227256709194E-11	-3.35736158987152E-11
C  9	-7.21528755064295E-11	2.62913750621447E-11
C  10	-2.48972568123080E-11	-1.23004756680251E-11
C  11	4.97289063344050E-11	-3.17027478788134E-11
C  12	3.44254300015349E-11	3.77010081520677E-11
C  13	-1.32732584044218E-11	1.36705821493913E-11
C  14	-6.13451559908228E-12	6.13224213274386E-12
C  15	5.39857538503087E-11	7.41190831968585E-12
C  16	-1.59963590787475E-11	-1.92499323577743E-11
C  17	-4.92939194089531E-11	-5.84155930264133E-12
C  18	0	0
C  19	0	0
C  20	0	0
C  21	0	0
C  22	0	0
C  23	0	0
C  24	0	0

CC Coefficients for IIR_B_Z
C a0		
       a0z =4.75367049250363d-10
C n	b(n)	a(n)
        bz(1) = 2.47844542710212d-10
        az(1) = -1.08372331619382d-07
        bz(2) = 3.29701497325601d-11
        az(2) = 4.13334132815567d-10
        bz(3) = 1.65833544007341d-10
        az(3) = -3.61512366093177d-10
        bz(4) = 1.73228126301479d-11
        az(4) = -4.26988412194273d-10
        bz(5) = 3.10435026926945d-11
        az(5) = 1.30177780901669d-09
        bz(6) = 5.34888708648183d-12
        az(6) = -9.81924654714910d-11
        bz(7) = -9.97896093090592d-12
        az(7) = 4.24338572656987d-10
        bz(8) = -4.55052766891820d-12
        az(8) = -1.35691362845386d-10
        bz(9) = 8.06492692901475d-12
        az(9) = 4.40729318633969d-10
        bz(10) = 8.66280346177117d-12
        az(10) = 3.50036828005455d-11
        bz(11) = 6.24370255801238d-13
        az(11) = 4.90513519009512d-10
        bz(12) = -2.19308807022689d-12
        az(12) = -9.10002332636802d-11
        bz(13) = -1.58384699717405d-11
        az(13) = 4.25632570186904d-10
        bz(14) = -4.51222576045757d-13
        az(14) = 1.88111028053916d-11
        bz(15) = 1.06395957461653d-11
        az(15) = 4.05941885091675d-10
        bz(16) = 1.10950943473286d-12
        az(16) = -7.45730521956724d-11
        bz(17) = 8.11004747777902d-13
        az(17) = 4.00354024787598d-10
C  n	b(n)	a(n)
C  1	2.47844542710212E-10	-1.08372331619382E-07
C  2	3.29701497325601E-11	4.13334132815567E-10
C  3	1.65833544007341E-10	-3.61512366093177E-10
C  4	1.73228126301479E-11	-4.26988412194273E-10
C  5	3.10435026926945E-11	1.30177780901669E-09
C  6	5.34888708648183E-12	-9.81924654714910E-11
C  7	-9.97896093090592E-12	4.24338572656987E-10
C  8	-4.55052766891820E-12	-1.35691362845386E-10
C  9	8.06492692901475E-12	4.40729318633969E-10
C  10	8.66280346177117E-12	3.50036828005455E-11
C  11	6.24370255801238E-13	4.90513519009512E-10
C  12	-2.19308807022689E-12	-9.10002332636802E-11
C  13	-1.58384699717405E-11	4.25632570186904E-10
C  14	-4.51222576045757E-13	1.88111028053916E-11
C  15	1.06395957461653E-11	4.05941885091675E-10
C  16	1.10950943473286E-12	-7.45730521956724E-11
C  17	8.11004747777902E-13	4.00354024787598E-10
C  18	0	0
C  19	0	0
C  20	0	0
C  21	0	0
C  22	0	0
C  23	0	0
C  24	0	0
C Calculate the accelerations:
C X axis
        srpxacc=a0x/2.d0
        do i=1,17
         srpxacc=srpxacc+(ax(i)*dcos(float(i)*angb)+
     .                 bx(i)*dsin(float(i)*angb))
C        print*,'float(i)',float(i)
        enddo
C Z axis
        srpzacc=a0z/2.d0
        do i=1,17
         srpzacc=srpzacc+(az(i)*dcos(float(i)*angb)+
     .                 bz(i)*dsin(float(i)*angb))
        enddo
C Convert to km s-2
C        Print*,'srpxacc',srpxacc
C        Print*,'srpzacc',srpzacc
        srpxacc=srpxacc*1.0d-3
        srpzacc=srpzacc*1.0d-3
        srpxacc=srpxacc*1100.d0/sbmass
        srpzacc=srpzacc*1100.d0/sbmass
C        Print*,'srpxacc',srpxacc
C        Print*,'srpzacc',srpzacc
      else
        if (first2) then
          first2=.false.
cd          Print*,'ARC/srpfmod iblock',iblock,
cd     .     ' no fourier coeffs for this block'
        endif
      endif
C        Print*,'iblock',iblock
C        Print*,'srpxacc',srpxacc
C        Print*,'srpyacc',srpyacc
C        Print*,'srpzacc',srpzacc
      return
      end
