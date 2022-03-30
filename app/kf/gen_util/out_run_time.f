ctitle
 
      subroutine out_run_time(iout, label, run_time)

      implicit none 
 
c
c     routine to output a label with runtime attached
c
c Variables
c ---------
c iout -- the output LU number
c label -- a label to written with runtime
c run_time -- the SOLVK runtime
c
      integer*4 iout, run_time(1), j
 
c
 
      character*(*) label
 
c
c.... output the line
      write(iout,100) label, (run_time(j),j=1,6)
  100 format(/,a,' Runtime ',i4,'/',i2,'/',i2,1x,i2,':',i2.2,'.',i2.2)
c
      return
      end
 
c......................................................................
