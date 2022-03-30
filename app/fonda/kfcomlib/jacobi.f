      subroutine jacobi(a,n,np,d,v,nrot,b,z)

      implicit none

*     Eigenvalue and eigen vector determination routine
*     from numerical recipes.

*     Computes all the eigen vectors and values for a real symmetric
*     matrix A of size NxN and dimensioned to NPxNP.  The upper
*     triangular portion of the matrix is destroyed in the operation.
*     D(N) returns the eigen values; V(nxn) returns the eigenvectors.
*     NROT is the number of rotations needed.  B and Z are work arrays.

*     MODIFIED from original by passing in the work arrays and explicit
*     declaration of some but not all variables

* n, np  - size and dimension of matrix
* nrot   - number of rotations needed.

      integer*4 n, np, nrot, ip, iq, i, j

* a(np,np) - Matrix whose eigenvalues are to be found
* d(np)    - Eigenvalues (not sorted)
* v(np,np) - Eignenvectors
* b(np), z(np) - Work arrays.

      real*8 a(np,np),d(np),v(np,np),b(np),z(np),sm,tresh,g,h,t,tau
     .     , theta,c,s

*     Start the computation.
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.
11      continue
        v(ip,ip)=1.d0
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.d0
13    continue
      nrot=0

* MOD TAH 971022: Changed 50 to np/2
      do 24 i=1,np/2+1
        sm=0.
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.)return
        if(i.lt.4)then
          tresh=0.2d0*sm/n**2
        else
          tresh=0.
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip)))
     *         .and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5d0*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.d0+theta**2))
                if(theta.lt.0.d0)t=-t
              endif
              c=1./sqrt(1.d0+t**2)
              s=t*c
              tau=s/(1.d0+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
      write(*,*) '**WARNING**', np/2, ' iterations should never happen'
      return
      end
