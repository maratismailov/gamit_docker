      subroutine propcv(b,cb,g,x,cx,mp,np,m,n)
c     propagate a covariance x = G * b 
c                            Cx = G * Cb * G**T


      real*8 b(np),x(mp),cb(np,np),cx(mp,mp),g(mp,np)
      integer max
      parameter (max = 1000)
      
      real*8 w(max,max)

      integer m,n,mp,np

      if (mp .gt. max .or. np .gt. max) then
         write (*,*) 'PROPCV: needs redimensioning'
         stop
      endif

c     calculate x = G * b
      call dotd(g,b,x, mp,np, np,1, mp,1, m,n,1, 1)

c     calculate w = Cb * G**T
      call dotd(Cb,g,W, np,np, mp,np, max,max, n,n,m, 3)

c     calculate Cx = G * W
      call dotd(g,w,cx, mp,np, max,max, mp,mp, m,n,m, 1)

      return
      end


