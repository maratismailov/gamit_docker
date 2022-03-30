      subroutine get_refractn(idatf)
c
c     get refraction corrections.
c     unit:  second
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      integer idatf,isit
      integer ib,ic,i0,j
      integer lift_arg,match_name,count_arg
      character*8 sitnam,r_apr
      character*128 line
c
c     initialize
      do isit = 1,maxsit
         defe(isit) = 0.0d0
         defn(isit) = 0.0d0
      enddo
      fac = dtor/3.6d3
c
c     read refraction correction file
      do 20 isit = 1,maxsit
         read (idatf,'(a)',end=200) line
         ic = count_arg(line)
         if (ic.le.0) goto 20
         ib = lift_arg(line,sitnam,1)
         if (ib.le.0) goto 20
         i0 = match_name(nsit,ib,sname,sitnam)
         if (i0.le.0) then
            print*,' GET_REFRACTN name mismatch:',sitnam
            goto 20
         endif
         do 10 j = 2,3
            ib = lift_arg(line,r_apr,j)
            if (ib.le.0) goto 10
            read (r_apr,*) ref
            if (j.eq.2) defe(i0) = ref*fac
            if (j.eq.3) defn(i0) = ref*fac
 10      continue
 20   continue

 200  continue

      return
      end
