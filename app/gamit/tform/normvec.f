      program normvec

c     3-d vector norm and uncertainty

      real*8 x ,y, z, c(3,3),jac(3),work(3),val(9)
      real*8 length,siglength,ddot

      integer i

      write(6,'(a)') 'Input vector of form'
      write(6,'(a)') '  x y z sigx sigy  sigz corxy corxz coryz'
      write(6,'(a)') '  <CNTL> Z to stop'

 10   continue
      read(5,*,end=999) (val(i),i=1,9)

      x = val(1)
      y = val(2)
      z = val(3)
      c(1,1)  = val(4)**2
      c(2,2)  = val(5)**2
      c(3,3)  = val(6)**2
      c(1,2)  = val(4)*val(5)*val(7)
      c(1,3)  = val(4)*val(6)*val(8)
      c(2,3)  = val(5)*val(6)*val(9)
      c(2,1)  = c(1,2)
      c(3,1)  = c(1,3)
      c(3,2)  = c(2,3)

      length = dsqrt(x*x + y*y + z*z)

      jac(1) = x/length
      jac(2) = y/length
      jac(3) = z/length

      call dgemv ( 'N', 3, 3
     .           , 1.d0, c, 3
     .           , jac, 1
     .           , 0.d0, work, 1)

      siglength = ddot(3,jac,1,work,1)
      siglength = dsqrt(siglength)

      write(6,'(/,2f15.4)') length,siglength

      goto 10
 999  continue

      end
