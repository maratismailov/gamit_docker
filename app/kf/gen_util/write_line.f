CTITLE WRITE_LINE
 
      subroutine write_line(iout, parn, first,second,units, est,
     .   adjust, sigma, type)

      implicit none 
 
 
c
c     This routine will output the parameter estimate line to the
c     the output device.
c
c Variables
c ---------
c iout -- the output LU
c parn -- the parameter number
c first -- the first label on the line
c second -- the second label
c units -- the units label
c est   -- the estimated parameter
c adjust -- the adjustment to the parameter
c sigma  -- the sigma on the adjustment
c type -- indicates if we need to change the estimated value in someway
c     options  action
c     type  0  do nothing
c     type  1  convert est to deg or hrs, min, secs
c
      integer*4 iout, parn, type
 
c
      character*(*) first, second, units
 
c
      real*8 est, adjust, sigma
 
c
c Local variables
c ---------------
c parnum_lab -- label for parameter number.  May be a number,
c      'loc.' to indicate that this results is computed, or blank.
c      If parn is greater than zero, field will contain parameter number.
c      If parn if less than zero, field will contain loc.
c      If parn is zero, field will be left blank.
c deg -- deg part of est
c min -- minutes part
c sec  -- the seconds part
c
      character*8 parnum_lab
 
c
      real*8 deg, min, sec
      integer*4 ideg, imin
c
c
c.... Firstly check to see if parameter number, 'loc. ' or blank label
c     should be given
*                               ! set label to loc.
      if( parn.lt.0 ) then
         parnum_lab = 'Loc. '
      else
*                               ! set label to parameter number
         if( parn.gt.0 ) then
            write(parnum_lab,'(i5,". ")') parn
*                               ! parn is zero so leave field blank
         else
            parnum_lab = '     '
         end if
      end if
c
c.... Check for type 0
*                              ! just normal output of est
      if( type.eq.0 ) then
c
c....    Here we have two options.  If est is zero then donot output
c        any value
*                              ! do not output estimate itself
         if( est.eq.0 ) then
c
            write(iout,100) parnum_lab(1:7),first,second,units,
     .         adjust, sigma
  100       format(a,a,1x,a,a,t57,2(1x,f12.5))
c
*                              ! OK output estimate
         else
c
            write(iout,120) parnum_lab(1:7),first,second, units,
     .         est, adjust, sigma
  120       format(a,a,1x,a,a,t39,f18.5,2(1x,f12.5))
c
         end if
c
*                ! type 0
      end if
c
c.... check for type 1
*                             ! convert estimate to deg or hrs, min, sec
      if( type.eq.1 ) then
c
         call mas_to_dms(est, deg,min, sec)
         ideg = nint(deg)
         imin = nint(min)
c
         if( ideg.ne.0 ) then
            write(iout,200) parnum_lab(1:7), first,second, units,
     .         ideg, imin, sec, adjust, sigma
 200        format(a,a,1x,a,a,t39,i3,I3,f11.7,2(1x,f12.4))
         else
            write(iout,210) parnum_lab(1:7), first,second, units,
     .         ideg, imin, sec, adjust, sigma
 210        format(a,a,1x,a,a,t39,i2,1x,I3,f11.7,2(1x,f12.4))
         end if
c
      end if
c
c.... Thats all
      return
      end
 
