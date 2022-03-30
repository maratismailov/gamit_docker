      subroutine cmodel(ierr)
c
c     get coodinate change model
c
c     coordinate mode : (cmode)
c         1. directly read network file (realized in getnet.f)
c         2. add earthquake induced coordinate jump
c         3. ......
c
c     unit:
c         x, y, z : meter
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'
      
      integer idx,iq,isit,ierr
      integer ib,i0,lsit,lbad,iread
      integer lift_arg,match_name
      character*8 sitnam
      character*128 line
c
      ierr = 0
      if (cmode.eq.2) then
         open (18,file=quakfil,status='old',err=1000)
c
c        read earthquake number
         read (18,*) nquake
         if (nquake.le.0) goto 200
c
c        read site name
         idx = 1
         isit = 0
         do 50 iq = 1,nquake
            iq_ind(iq) = idx
            read (18,*) q_time,lsit
            quake_time(iq) = q_time
            lbad = 0
            do 20 iread = 1,lsit
               read (18,'(a)',end=200) line
               read (line,*) ce,cn,cu
               ib = lift_arg(line,sitnam,4)
               if (ib.le.0) goto 20
               i0 = match_name(nnet,ib,sname,sitnam)
               if (i0.le.0) then
                  print*,' CMODEL name mismatch:',sitnam
                  lbad = lbad+1
                  goto 20
               endif
               isit = isit+1
               quake_sit(isit) = i0
               quake_ce(isit) = ce
               quake_cn(isit) = cn
               quake_cu(isit) = cu
 20         continue
 30         idx = idx+lsit-lbad
 50      continue

 200     iq_sit = isit
         close (18)
      endif

 1000 continue
      return
      end
