************************************************************************
******************** PARAMETERS FOR DISPLACE.F *************************
************************************************************************
      
      include '../includes/const_param.h'

      integer MAXSTA, MAXFLT, one, zero

      parameter ( one = 1 )
      parameter ( zero = 0 )

      parameter (MAXSTA=20000, MAXFLT=100)

      integer makegrid, numsta, numflts, xnum, ynum

      character*8 sname(MAXSTA), scale(MAXFLT)
      character*256 infile, outfile

      real slong(MAXSTA), slat(MAXSTA), 
     +      flong(MAXFLT), flat(MAXFLT), 
     +      fstrike(MAXFLT), depth(MAXFLT), 
     +      dip(MAXFLT), length(MAXFLT), width(MAXFLT), 
     +      u1(MAXFLT), u2(MAXFLT), u3(MAXFLT), 
     +      fx, fy,
     +      xmin, ymin, xinc, yinc, 
     +      origlong, origlat, 
     +      x(MAXSTA), y(MAXSTA), 
     +      tx(MAXSTA), ty(MAXSTA), 
     +      rx(MAXSTA), ry(MAXSTA), 
     +      maxlon, minlon, maxlat, minlat,
     +      parallel1, parallel2, 
     +      disple(MAXSTA), displn(MAXSTA), displu(MAXSTA),
     +      deltae(MAXSTA), deltan(MAXSTA), deltau(MAXSTA)

C     common / dis_real / slong, slat, 
C    +      flong, flat, 
C    +      fstrike, depth, 
C    +      dip, length, width, 
C    +      u1, u2, u3, 
C    +      fx, fy,
C    +      xmin, ymin, xinc, yinc, 
C    +      origlong, origlat, 
C    +      x, y, 
C    +      tx, ty, 
C    +      rx, ry, 
C    +      maxlon, minlon, maxlat, minlat,
C    +      parallel1, parallel2, 
C    +      disple, displn, displu,
C    +      deltae, deltan, deltau

C     common / dis_int/ makegrid, numsta, numflts, xnum, ynum
C     
C     common / dis_ch/  sname, scale, infile, outfile

* MOD TAH 031026: Add information to fit globk model
* MOD TAH 170415: Added more eq_def information
      real*8 glat, glong, gdepth, grad  ! Depth in km, grad (m) 
      real*8 eq_pos(3)    ! XYZ coordinates of earthquake
      real*8 neq          ! Sum of 1/d^2 partials
      real*8 bh, bv       ! B-vector for hamp and vamp

      integer*4 nume,     ! Number of values in estimate
     .     eq_date(5)     ! YMDMH of earthquake (eq_def line)
   
      character*8 eq_code    ! Earthqake code
      

      common / gdisp / glat, glong, gdepth, grad,
     .                 eq_pos, neq, bh, bv, 
     .                 nume, eq_date
      common / gdisp_ch / eq_code





