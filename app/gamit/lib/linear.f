      Subroutine LINEAR( xin,m,dint,x1,x2,ya,out)

c     Purpose:  To linearly interpolate a 1-d vector 'ya' of
c               size m and dimension maxel.  x1 and x2 are
c               miminimum and maximum range, dint is the
c               interval, and xin is the desired point of
c               interpolation.
                                                         
c     R. King 20 July 2018, replacing the earlier subroutine that
c     omitted the array limits and used integer tabular points and 
c     interval.

      implicit none
                
      include '../includes/dimpar.h'
             
      integer*4 m,indx,len,rcpar,i

      real*8 xin,x1,x2,dint,ya(maxel),fract,out
                           
      character*80 prog_name
      character*256 message          

c        Get the program name for report_stat
      len = rcpar(0,prog_name)
       
c        Find the index in the array and the fraction beyond it
      indx = int( (xin-x1)/dint )  + 1
       if( indx.gt.m.or.indx.gt.maxel ) then 
          write(message,'(2(a,i3),a)')  'Interpolation index (',indx
     .    ,') > array index (',m,') or dimension (',maxel,')' 
        call report_stat('WARNING',prog_name,'lib/linear',' ' 
     .                    ,message,0)
      endif
      fract =  (xin-x1)/dint - float(indx-1)
cd       print *,'LINEAR ya ',(ya(i),i=1,m)      
cd       print *,'LINEAR m x1 xin dint indx fract '
cd     .               , m,x1,xin,dint,indx,fract 

c        Interpolate the value
      out = ya(indx) + fract*(ya(indx+1)-ya(indx))    

      return
      end


