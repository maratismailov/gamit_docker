      subroutine get_quake_list(idatf)
c
c     get earthquake time and influenced site list
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
c
      integer idatf,i1,idx,iq,isit,ia
      integer ib,ic,i0,lsit,j,lbad,iread
      integer lift_arg,match_name,count_arg
      dimension apr_val(6),cjac1(9),tmp1(9),tmp2(9),tmp3(9)
      character*8 sitnam,e_apr
      character*128 line
c
c     read earthquake number
      read (idatf,*) nquake
      if (nquake.le.0) goto 200
c
c     read site name and constraint
      idx = 1
      isit = 0
      call zero1d(1,maxprm,gvm)
      do 50 iq = 1,nquake
         iq_ind(iq) = idx
         read (idatf,*) q_time,lsit
         quake_time(iq) = q_time
         i1 = 0
         lbad = 0
         do 20 iread = 1,lsit
            read (idatf,'(a)',end=200) line
            ic = count_arg(line)
            if (ic.le.0) goto 20
            ib = lift_arg(line,sitnam,1)
            if (ib.le.0) goto 20
            i0 = match_name(nsit,ib,sname,sitnam)
            if (i0.le.0) then
               print*,' GET_QUAKE_LIST name mismatch:',sitnam
               lbad = lbad+1
               goto 20
            endif
            isit = isit+1
            quake_sit(isit) = i0
            do 10 j = 2,4
               ib = lift_arg(line,e_apr,j)
               if (ib.le.0) then
                  apr_val(j-1) = 1.0d1
                  goto 10
               endif
               read (e_apr,*) apr_val(j-1) 
 10         continue
            call zero1d(1,6,tmp1)
            tmp1(1) = apr_val(1)**2
            tmp1(3) = apr_val(2)**2
            tmp1(6) = apr_val(3)**2
            if (tmp1(1).lt.1.0d-6) tmp1(1) = 1.0d-6
            if (tmp1(3).lt.1.0d-6) tmp1(3) = 1.0d-6
            if (tmp1(6).lt.1.0d-6) tmp1(6) = 1.0d-6
            a1 = slat(i0)
            a2 = slon(i0)
            a3 = srad(i0)
            call getjac(a2,a1,a3,cjac1,6)
            call atwa(3,3,cjac1,tmp1,tmp2,tmp3,1)
            call cholsk(tmp2,tmp3,1,3,ia)
            ia = (isit-1)*6
            do j = 1,6
               gvm(ia+j) = tmp2(j)
            enddo
            i1 = i1+1
            if (i1.ge.lsit-lbad) goto 30
 20      continue
 30      idx = idx+lsit-lbad
 50   continue

 200  iq_sit = isit
      print*,' Identified sites with coseismic displacements =',
     .   iq_sit
      print*,' quake number:',nquake

      return
      end
