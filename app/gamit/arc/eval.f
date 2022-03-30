Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      subroutine eval (s,l2)

c     W.B.Smith, March 1966
c     Rick Abbot  October 1984, modification to do GPS satellites only

      implicit none

      include '../includes/dimpar.h'   
      include '../includes/arc.h'


      integer*4 idim,l2,i,k

      real*8 s,whc1,hh,sbfn,g,h,w,whc,amx,bmx,sta,emx,e1

      logical debug/.false./

      dimension w(6),h(6),hh(6),whc(6),e1(maxyt2)

      data hh(2),hh(3),hh(4),hh(5),hh(6)/2.d0,3.d0,4.d0,5.d0,6.d0/
      data w(1),w(2),w(3),w(4),w(5),w(6)/
     1 3.1559193121691931d-01, 2.2833333333333333d0,
     2 3.7500000000d0,4.2500000000d0,3.0000000000d0,1.0000000000d0/
c                                                             

      if(debug) print *,'EVAL l1 l2 n1 n2 hc ',l1,l2,n1,n2,hc 

c        y and f predictions from data computed in last successful
c        step
100   h(1) = hc
      if(debug) print *,'EVAL 100 s hc ',s,hc  
      whc1=w(1)*hc
      do i=2,6
        h(i)=hc/hh(i)
      enddo
      do k = n1,n2
        i=k          
        fnord(i,1) = fnord(i,3) + h(1)*(anord(i) + h(2)*(bnord(i) 
     .     + h(3)*(cnord(i)+h(4)*(dnord(i)+h(5)*enord(i)))))   
        if(debug.and.i.eq.1) print *,'h a b c d e f nord '
     .    , h(1),anord(i),bnord(i),cnord(i),dnord(i),enord(i),fnord(i,1)
        y(i,1) = y(i,4) + h(1)*(fnord(i,3)+h(2)*(anord(i)+h(3)*(bnord(i)
     .        + h(4)*(cnord(i) +h(5)*(dnord(i)+h(6)*enord(i))))))
      enddo
c
c        y is corrected        iteration 1
      idim = 1
      do k = n1,n2
        i = k
        fnord(i,2) = sbfn(i,idim,s+hc)
        e1(i)=fnord(i,2)-fnord(i,1)
        y(i,2)=y(i,1)+whc1*e1(i)
      enddo

c        y is corrected        iteration 2
      idim = 2
      do k = n1,n2
        i = k
        fnord(i,2) = sbfn(i,idim,s+hc)
        e1(i) = fnord(i,2) - fnord(i,1)
        y(i,3) = y(i,1) + w(1) * e1(i)*hc
      enddo         
      if(debug) print *,'EVAL idim hc s fnord1 fnord2 e1 y2 '
     .           ,idim,hc,s,fnord(1,1),fnord(1,2),e1(1),y(1,2)



      if (k2 - 19)1010,710 ,700
700   if (n2.gt.meq) go to 1010
c
c        interval control logic
c
710   g = epsi/dabs(hc)
      amx = dabs(y(1,3) - y(1,2))
      bmx = dabs(y(1,2) - y(1,1))
      emx = dabs(e1(1))
      if (1-neq)720,740,740
720   do i=2,meq
        amx = dmax1(amx,dabs(y(i,3) - y(i,2)))
        bmx = dmax1(bmx,dabs(y(i,2) - y(i,1)))
        emx = dmax1(emx,dabs(e1(i)))
      enddo
      amx=amx+1.d-35
      bmx=bmx+1.d-30

c        stability test is performed if bmx is g.t. epsi/16
740   sta=epsi/16.d0   
      if(debug) then
        write (*,'(1x,8i5)') neq,meq,n1,n2,l1,k2,l4,l5
        write (*,'(1x,6(d12.5,1x))') amx,bmx,emx,g,epsi,hc
      endif
      if(bmx-sta) 760,760,750
750   continue
      if (amx - (1.25d-01*bmx))760,760,900

c        accuracy test
760   if (emx - g)800,1000,900

c        tests satisfied. attempt made to double interval
800   if(l5)1000,810,810
810   if(l4)1000,820,820
820   if( 260.d0*emx-g)830,1000,1000
830   if(25.6d0*amx-bmx)840,840,1000
840   hc = hc*2.0d0                           
      if(debug) print *,'hc doubled ',hc
      l4=-5
      l2 = -1
      go to 1900

c        test unsatisfied. interval halved.
  900 continue
      hc = hc/2.0d0                      
      if(debug) print *,'hc halved ',hc 
      l4=-5
      if (dabs(hc) - hmn)920,920,910
910   l2 = -1
      go to 1900    
920   call report_stat('FATAL','ARC','eval',' '
     .,'Interval to minimum, something wrong with ICs or model',0)
c*920   write (6,930)
c*930   format (//10x,60hinterval to minimum    stop, go to next integrati
c*     1on job               )
      l2=1
      go to 1900
c
1000  if (k2.eq.19) go to 1010
      if ((n2.gt.meq).or.(meq.eq.neq)) go to 1010
      n1 = meq + 1
      n2 = neq
      go to 100
c
c        a,b,c,d,e are corrected. y,f,s are updated. k2 and l4 are
c        advanced
1010  s=s+hc
      l4 = l4+ 1
      k2 = k2+ 1
      if(l5)1020,1030,1020
1020  l5=-1
1030  continue
      whc(1)=1.0d0
      do i=2,6
        whc(1)=whc(1)*hc
        whc(i)=w(i)/whc(1)
      enddo
      do i=1,neq
        fnord(i,3) = fnord(i,2)
        y(i,4) = y(i,iptr3)
        anord(i) = anord(i) 
     .   + h(1)*(bnord(i)+h(2)*(cnord(i)+h(3)*(dnord(i)+h(4)*enord(i))))
     .   + whc(2)*e1(i)
         bnord(i) = bnord(i) 
     .    + h(1)*(cnord(i)+h(2)*(dnord(i)+h(3)*enord(i))) + whc(3)*e1(i)
         cnord(i) = cnord(i) + h(1)*(dnord(i) + h(2)*enord(i)) 
     .    +  whc(4)*e1(i)
         dnord(i) = dnord(i) + h(1)*enord(i) + whc(5)*e1(i)
         enord(i) = enord(i) + whc(6)*e1(i) 
      enddo
1900  continue         
      if(debug) print *,'EVAL nord (1) a b c d e f '
     .   , anord(1),bnord(1),cnord(1),dnord(1),enord(1),fnord(1,3)
cd      stop

      return
      end
