      program trsat

c     convert satellite ephemeris units from
c     GLOBAL convention to GAMIT G-file convention.

c
c     Kurt Feigl Jan 31, 1991
c
      integer*4 ioerr
      integer maxel
      parameter (maxel = 9)
      real*8 satels(maxel)
      character*6 satcode

c extract prn.ext glba_r.prt 6 | grep -v '*' | ~/gu/util/sccs/trsat >! gfile

c     then you have to add the time tag to the output

c     here is an extract command file to get the right stuff

c field 1 "Inert.  X     (m)" 1 CH readline 1 2
c field 2 "Inert.  X     (m)" 1 R8 readline 0 1
c field 3 "Inert.  Y     (m)" 1 R8 readline 0 1
c field 4 "Inert.  Z     (m)" 1 R8 readline 0 1
c field 5 "Inert.  dX/dT (mm/s)" 1 R8 readline 0 1
c field 6 "Inert.  dY/dT (mm/s)" 1 R8 readline 0 1
c field 7 "Inert.  dZ/dT (mm/s)" 1 R8 readline 0 1
c field 8 "Solar Rad. X  (unit)" 1 R8 readline 0 1
c field 9  "Solar Rad. Y  (unit)" 1 R8 readline 0 1
c field 10 "Solar Rad. Z  (unit)" 1 R8 readline 0 1


c outform 1 "(1x, a6)"
c outform 2 "(1x, D20.14)"
c outform 3 "(1x, D20.14)"
c outform 4 "(1x, D20.14)"
c outform 5 "(1x, D20.14)"
c outform 6 "(1x, D20.14)"
c outform 7 "(1x, D20.14)"
c outform 8 "(1x, d20.14)"
c outform 9  "(1x, d20.14)"
c outform 10 "(1x, d20.14)"

c run


      icount = 1

      ioerr = 0
      do 100 while (ioerr .eq. 0)
         read(unit   = 5,
     +      fmt    = '(1x,a6,9(2x,d20.14))',
     +      err    = 500,
     +      end    = 1000,
     +      iostat = ioerr)
     +      satcode,(satels(i),i=1,maxel)

c           change underscores to spaces
            if (satcode(4:4) .eq. '_') satcode(4:4) = ' '
c           change leading zero to space
            if (satcode(5:5) .eq. '0') satcode(5:5) = ' '

c           change PRN numbers greater than 30 (eclipsing)
c           back to their rightful values
            read  (satcode(5:6),'(i2)') iprn
            if (iprn .gt. 30) iprn = iprn - 30
            write (satcode(5:6),'(i2)') iprn

c           convert meters to kilometers
            satels(1) = satels(1) * 1.e-3
            satels(2) = satels(2) * 1.e-3
            satels(3) = satels(3) * 1.e-3
c           convert mm/s to km/ss
            satels(4) = satels(4) * 1.e-6
            satels(5) = satels(5) * 1.e-6
            satels(6) = satels(6) * 1.e-6
c           remaining 3 parameters (7-9) are
c           dimensionless in both conventions.

            write (6,'(a6)') satcode
            write (6,'(8(d20.14,/),d20.14)') (satels(i),i=1,maxel)
 100  continue
 500  continue
      call ferror (ioerr,6)
1000  continue

      stop
      END

