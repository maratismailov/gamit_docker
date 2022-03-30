      Program coeffs
C ejp Jan 2011 
c testing how to read in coefficients. Doesn't test for changes in format.
C test only, not part of main code.
      implicit none

      integer*4 iblock,i,err
      real*8 angb,srpxacc,srpzacc
C     Dimension to fifty - will be ok if only 17 coeffs.
      integer*4 nx(50),ny(50),nz(50),j,n(50)
      real*8 a0x,a0z,bx(50),ax(50),bz(50),ay(50),by(50),az(50)
      character*80 line


C set up 
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

C       Use the UCL coefficients up to 50, X, Y. Z body fixed rerefernce system.
C       These are officially for SVN 41, 43, 44, 46, 51 (IIR-A 2-6 as 1 failed on launch)
        open (unit = 5678, file = 
     .  "/chandler/home/petrie/orbitrefs/data/MarekSRPfourierX.txt",
C     .  "testinput.txt",
     .   status='old',iostat=err)
C       Read six text header lines
        do i=1,6
        read (5678,*) line
cd        print*,i, line
        enddo

C       Read a0
        read(5678,*) a0x
cd        print*, a0x

        read (5678,*) line
cd        print*, line
C       Read coefficients
        do j=1,50
         read(5678,*) n(j),ax(j),bx(j)
cd         print*, n(j),ax(j),bx(j)
        enddo
        close (5678)
cd        Print*,err,'iostat'
cd        Print'(1x,i10,/)',n

cd        Print*,'ax'
cd        Print*,ax
cd        Print*,'bx'
cd        Print*,bx
        end
