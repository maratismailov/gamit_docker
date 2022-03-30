Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
c
      subroutine oline(id,istat,jstat,free_fix,sigs,baseln,dbase,mode)
c
c     write a single - line record to o-file
c     mode = 1: (x,y,z)
c     mode = 2: (n,e,u)
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      character*1 bfixcd,unders
      character*4 buf4a,buf4b,upperc,free_fix
      character*3 buf3
      character*16 buf16

      real*8 sigs(3,3),baseln,dbase

      integer id(maxsit),istat,jstat,mode,iascii,is1,iyear,ijk

c     underscore for printing
      iascii=95
      unders=char(iascii)

c       get 4 char code for first station from c-file name
        is1=id(istat)
        if (is1.lt.1) is1=1
        buf16 = cfiln(is1)
        buf4a = buf16(2:5)
        call uppers(buf4a)
c       get 4 char code for 2nd  station from c-file name
        is1=id(jstat)
        if (is1.lt.1) is1=1
        buf16 = cfiln(is1)
        buf4b = buf16(2:5) 
        call uppers(buf4b)
c       year number  
        iyear = itor(3,is1)
        call check_y2k(iyear)
c       bias fixing code letter
        if ( free_fix.eq.'fixd' ) then
            write(bfixcd,'(a1)') upperc('x')
        else
            write(bfixcd,'(a1)') upperc('r')
        endif
c       get 3 char day number string from c-file name
c       warning, this does not yet work for multi-session solutions
        buf16 = cfiln(istat)
        buf3 = buf16(8:10)
        if (mode.eq.1)
     .   write (15,20) buf4a,unders,buf4b,iyear,buf3,bfixcd,
     .   (sigs(ijk,1),sigs(ijk,2),ijk=1,3),baseln,dbase,
     .   (sigs(ijk,3),ijk=1,3)
 20    format (a4,a1,a4,1x,i4,'.',a3,a1,1x,
     .         'X',f14.4,1x,'+-',1x,f8.4,1x,
     .         'Y',f14.4,1x,'+-',1x,f8.4,1x,
     .         'Z',f14.4,1x,'+-',1x,f8.4,1x,
     .         'L',f14.4,1x,'+-',1x,f8.4,1x,
     .         ' Correlations (X-Y,X-Z,Y-Z) = ',3f10.5)

       if (mode.eq.2)
     .   write (15,30) buf4a,unders,buf4b,iyear,buf3,bfixcd,
     .   (sigs(ijk,1),sigs(ijk,2),ijk=1,3),baseln,dbase,
     .   (sigs(ijk,3),ijk=1,3)

 30    format (a4,a1,a4,1x,i4,'.',a3,a1,1x,
     .         'N',f14.4,1x,'+-',1x,f8.4,1x,
     .         'E',f14.4,1x,'+-',1x,f8.4,1x,
ckf  .         'U',f15.4,1x,'+-',1x,f8.4,1x,
     .         'U',f14.4,1x,'+-',1x,f8.4,1x,
     .         'L',f14.4,1x,'+-',1x,f8.4,1x,
     .         ' Correlations (N-E,N-U,E-U) = ',3f10.5)

       return
       end
