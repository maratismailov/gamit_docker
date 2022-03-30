      subroutine put_fix(wcmd, lcmd, chi)
c
c     decomposite fix_para command
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      character*(*) wcmd
      character*10 code
      integer i, j, i0, ic, ib, match_name
      integer lcmd,count_arg,lift_arg,idx(9)
c  
c     decomposite command line
      ic = count_arg(wcmd)
c     not enough arguments
      if (ic.le.1) goto 200
c     pointer to coordinate and velocity
      ib = lift_arg(wcmd,code,1)
      if (ib.le.0) goto 200
      do i = 1,6
         idx(i) = 0
      enddo
      do i = 1,ib
         if (code(i:i).eq.'x') idx(1) = 1
         if (code(i:i).eq.'y') idx(2) = 1
         if (code(i:i).eq.'z') idx(3) = 1
         if (code(i:i).eq.'u') idx(4) = 1
         if (code(i:i).eq.'v') idx(5) = 1
         if (code(i:i).eq.'w') idx(6) = 1
      enddo

c     print header for chi2 output
      print*,''
      print*,' Site        chi2      delta chi2'
      do 20 i = 2,ic
         ib = lift_arg(wcmd,code,i)
         if (ib.le.0) goto 20
         if (code(1:3).eq.'all') then
            chio = chi
            do 10 j = 1,nsit
               call reduce_para(j,idx,chi)
 10         continue
            print*,' chi2 after fix components at ',
     .         'ALL sites : ',chi,chi-chio
            goto 200
         endif
         i0 = match_name(nsit,ib,sname,code)
         if (i0.le.0.or.i0.gt.nsit) then
            print*,' PUT_FIX name mismatch:',code
            goto 20
         endif
         chio = chi
         call reduce_para(i0,idx,chi)
c        print*,' chi2 after fix site : ',i0,chi,chi-chio
         print'(1x,a8,1x,f12.2,1x,f12.2)',
     .      sname(i0),chi,chi-chio
         fixcnt = fixcnt+1
         fixsit(fixcnt) = i0
 20   continue

 200  continue

      return
      end
