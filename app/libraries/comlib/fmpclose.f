
CTITLE FMPCLOSE

      integer*4 function fmpclose( dcb )

      implicit none

*     This routine will close the the fortran unit connected to
*     dcb buffer.  The unit number is in dcb(1)

      integer*4 dcb(16)

* LOCAL VARIABLES
* ierr - IOSTAT error on close

      integer*4 ierr

****  close the unit number
      close ( dcb(1), iostat=ierr )
      fmpclose = -ierr

****  Thats all
      return
      end 
