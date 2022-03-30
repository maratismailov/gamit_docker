      subroutine sensitive(key,job,chi)
c
c*aut FERHAT Gilbert
c*aut GRGS
c
c*ver version 1.00 July, 11 1995
c*
c*rol compute the "sensibilite" pour un relevement
c
c************************************************************************
c
      character*64  sub_file
c
      character*8 sname2(maxsit)
      character*8 sitnam
c
c    
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
      integer idatf,iobs,i,j,ifile,isit1,ie1,it
      integer illsit,len,nblen,n1,nlen,match_name
      integer iuse(maxsit),vuse(maxsit)
      character*8 name1
      character*64 subfil
      equivalence (iuse,map)
      equivalence (vuse,map(maxsit+1))
c
c
      rewind(14)
c     read site number
      read (14,*) nsit
c
c     read site name(id)
      do i = 1,nsit
         if (fcode.le.2) then
            read (idatf,'(a8)') sitnam
            sname(i) = sitnam
         endif
      enddo
c
      read(14,'(a)',err=140,end=200) sub_file
      rewind(14)
      close(14)
c
      print *,' subfile= ',sub_file
      len = nblen(sub_file)
      open (31,file=sub_file(1:len),status='old',err=130)
c
         read (31,*) n1
         do 20 j = 1,n1
            read (31,'(a8)') name1
            nlen = nblen(name1)
            isit1 = match_name(nsit,nlen,sname,name1)
            if (isit1 .le. 0) goto 20
            if (iuse(isit1).lt.minic) iuse(isit1) = iuse(isit1)+1
            if (iuse(isit1).ge.minic) vuse(isit1) = vuse(isit1)+1
 20      continue
         read (31,*) ie1
         read (31,*) ie1,it,iobs
ccc         if (it.eq.31.or.it.eq.32) jaux = jaux+6
         if (it.eq.2.and.iobs.ge.3) then
             do i = 1, nb_exp
cc              read(31,120,end=200) an,iobs1,iobs2,obs3
                read(31,'(a)',end=200) buffer
                read(buffer,*) an,iobs1,iobs2,obs3,
     &                 sigma_angle
                nlen1 = lift_arg(buffer,sitea,6)
c               print*,'1 nlen1 ',nlen1
                nlen1 = lift_arg(buffer,point_viseb,7)
c               print*,'2 nlen1 ',nlen1
                nlen1 = nblen(sitea)
c               isite1 = match_name(nsit,nlen1,sname,site)
                isite1 = match_name(nb_pts,nlen1,sname,sitea)
c               print*,' nsit  =',nsit
c               print*,' nlen1 =',nlen1
c               print*,' sname =',sname
c               print*,' isite1=',isite1
c               print*,' sitea =',sitea
                if (isite1 .le. 0) goto 300
c               print*, ' ff1.5'
                nlen1 = nblen(point_viseb)
c               isite2=match_name(nsit,nlen1,sname,point_vise)
                isite2=match_name(nb_pts,nlen1,sname,point_viseb)
                if (isite2 .le. 0) goto 300
c               print*, ' ff2'
c 120           format(2x,f9.4,4x,f9.3,5x,f5.1,9x,a8,2x,a8)
 120            format(2x,f9.4,3x,i3,3x,i2,3x,f11.8,3x,
     &           f6.4,2x,a8,2x,a8)
c
c               write(*,*) an,iobs1,iobs2,obs3
c    &          ,sigma_angle,sitea,point_viseb
c               write(*,120) an,iobs1,iobs2,obs3
c    &          ,sigma_angle,sitea,point_viseb
               do kk = 1, nb_pts
                  if (sitea.eq.station(kk)) then                     
                      write(28,*) phi(kk),xlambda(kk)
c                      write(*,*)  phi(kk),xlambda(kk)
                  endif
               enddo
               do kk = 1, nb_pts
                  if (point_viseb.eq.station(kk)) then
                      write(28,*) phi(kk),xlambda(kk)
                      write(28,*) '> '
c                     write(*,*) phi(kk),xlambda(kk)
c                     write(*,*) '> '
                  endif   
               enddo
             enddo
         endif
c
    



         close(31)

      return
      end
