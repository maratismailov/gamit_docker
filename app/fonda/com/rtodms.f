      subroutine rtodms(arad,ideg,iminu,asec,job)
c	
c     transform a value from radian to degree, minute and second
c
c     job  :  1 = from radian to degree, minute, second
c             2 = from degree to degree, minute, second
c
cmk in the southern hemisphere this code results in 
cmk negative degrees, minutes and seconds.
      implicit real*8 (a-h,o-z)
      integer ideg,iminu,job 

      a1 = arad
      if (job.eq.2) goto 10

      rtod = 45.0d0/datan(1.0d0)
      a1 = a1*rtod

 10   ideg = int(a1)
      a2 = (a1-dble(ideg))*60.0d0
      iminu = int(a2)
      asec = (a2-dble(iminu))*60.0d0

      return
      end
