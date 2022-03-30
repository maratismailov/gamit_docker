c
      subroutine gmt_sketch(ifile)
c
c**************************************************************************
c
c*aut Gilbert FERHAT
c*aut GRGS
c
c*ver Oct  9, 1995 v0_00
c
c*par               input parameters
c
c*par ifile       : #28 GMT sketch file
c*par sub_file    : sequential input sub_file (#31) is read in the file list
c*par iomode(1)   : = 1 if priori value file (#13) is specified
c
c*par               output parameters
c
c*par sketch_file : sketch file for GMT (#28)
c
c*rol 1st : read the a priori coordinates file
c*rol 2nd : read the FIRST observation data file specified in the file list
c*rol 3rd : create a GMT sketch file
c
c***************************************************************************
c 
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
      integer       ifile,kk,iphi_int,iphi_2nd,lambda_int,lambda_2nd
      integer       nb_files,len,nb_pts,i,no_exp,itype_observ,nb_exp
      integer       isite1,isite2
      integer       nlen1,nblen  
      integer       match_name, lift_arg  
      character*1   nord,est,star
c     character*8   station(900),site,point_vise
      character*8   station(900),sitea,point_viseb
      character*80  skipp
c     character*65  buffer
      character*128 buffer
      character*22  full_name
      character*64  sub_file
c     sub_file    : sequential input file(#31) specified in the file list  
      real*8        alti, annee, xlambda_3rd
      real*8        phi(900), xlambda(900), sigma_angle, obs3
      integer       iobs1, iobs2,nombre_pts
c
c*************************************************************************
c     open apriori file
c*************************************************************************
c
c     print *, 'starting with gmt_sketch subroutine'
c     file #13 is the apriori value file
      rewind (13)
      kk = 0
 10   continue
      read(13,'(a1)',end=30) star
      if (star.ne.' ') goto 10        
      if (star.eq.' ') then
          backspace 13
          kk = kk + 1
          read(13,20,end=30) station(kk),full_name(1:12),nord,
     &       iphi_int,iphi_2nd,phi_3rd,est,lambda_int,
     &       lambda_2nd,xlambda_3rd,alti,zer1,zer2,zer3,
     &       annee,zer4,zer5,zer6
 20       format (1x,a8,1x,a12,1x,a1,2(i2,1x),f8.5,1x,
     .            a1,i3,1x,i2,1x,f8.5,f13.4,3f8.4,f9.3,3f8.4)
c
c
c         write(*,20) station(kk),full_name(1:12),nord,
c    &       iphi_int,iphi_2nd,phi_3rd,est,lambda_int,
c    &       lambda_2nd,xlambda_3rd,alti,zer1,zer2,zer3,
c    &       annee,zer4,zer5,zer6
c    &       lambda_2nd,xlambda_3rd,alti,zero,zero,zero,
c    &       annee,zero,zero,zero
c
          phi(kk) = dble(iphi_int) + dble(iphi_2nd)/60.0d0 +
     &              phi_3rd/3600.0d0
c
          if (nord.eq.'N') phi(kk) =   phi(kk)
          if (nord.eq.'S') phi(kk) = - phi(kk)
c
          xlambda(kk) = dble(lambda_int) + dble(lambda_2nd)/60.0d0 
     &                 + xlambda_3rd/3600.0d0
          if (est.eq.'E') xlambda(kk) =   xlambda(kk)
          if (est.eq.'W') xlambda(kk) = - xlambda(kk)
c
      endif
c
      goto 10
 30   continue
c     print*, 'nb de points lus ',kk     
      nombre_pts=kk
c
c*********************************************************************
c     open observation data file
c*********************************************************************
c
c     the file list #14 was opened before in the subroutine GETSIT.F
c     lines read before in file #14 are 
c           -  1st : the number of stations
c           -  and then all the stations names
      read(14,*) nb_files       ! nb_files : number of sequential input files
c                                            specified in the file list
c     We just read the first file name
c     and want to draw the sketch just for this one
c
c     get sub-file name (#31)
      read(14,'(a)',err=140,end=200) sub_file
c     print *,' subfile= ',sub_file
      len = nblen(sub_file)
      open (31,file=sub_file(1:len),status='old',err=130)
c
c     read(31,40) nb_pts
      read(31,'(i)') nb_pts
c     print*, '  nb_pts= ',nb_pts
c40   format(2x,i3)
      do i = 1, nb_pts
         read(31,'(a10)') skipp
      enddo
      read(31,'(a10)') skipp
c     print*, '  nb_total_exp= ',skipp
c   
      if (iomode(5).gt.0)  then 
         write(28,*) '> '
 100     continue
         read(31,'(a)',end=200) buffer
         read(buffer,*) no_exp,itype_observ,nb_exp
c
c        write(*,*)  no_exp,itype_observ,nb_exp
c        write(*,110)  no_exp,itype_observ,nb_exp
 110     format(4x,i3,4x,i1,4x,i4)
c
c        itype_observ = 2,3  : SOLVEM code for horizontal angle 
         if (itype_observ.eq.2.or.itype_observ.eq.3) then
             do i = 1, nb_exp
cc              read(31,120,end=200) an,iobs1,iobs2,obs3
                read(31,'(a)',end=200) buffer
                read(buffer,*) an,iobs1,iobs2,obs3,
     &                 sigma_angle
                nlen1 = lift_arg(buffer,sitea,6)
c               print*,'1 nlen1 ',nlen1
c               print*,'1 nlen1 ',nlen1
                nlen1 = lift_arg(buffer,point_viseb,7)
c               print*,'2 nlen1 ',nlen1
c               print*,'2 nlen1 ',nlen1
                nlen1 = nblen(sitea)
c               isite1 = match_name(nsit,nlen1,sname,site)
                isite1 = match_name(nb_pts,nlen1,sname,sitea)
cccc            isite1 = match_name(nsit,nlen1,sname,sitea)
cc              print*,' nsit  =',nsit
cc              print*,' nlen1 =',nlen1
c               print*,' sname =',sname
cc              print*,' isite1=',isite1
cc              print*,' sitea =',sitea
                if (isite1 .le. 0) goto 300
cc              print*, ' ff1.5'
c               nlen1 = nblen(point_viseb)
c               isite2=match_name(nsit,nlen1,sname,point_vise)
                isite2=match_name(nb_pts,nlen1,sname,point_viseb)
cccc            isite2=match_name(nsit,nlen1,sname,point_viseb)
cc              print*,' nsit  =',nsit
c               print*,' nlen1 =',nlen1
c               print*,' sname =',sname
c               print*,' isite2=',isite2
cc              print*,' point_viseb =',point_viseb
                if (isite2 .le. 0) goto 300
c               print*, ' ff2'
c 120           format(2x,f9.4,4x,f9.3,5x,f5.1,9x,a8,2x,a8)
 120            format(2x,f9.4,3x,i3,3x,i2,3x,f11.8,3x,
     &           f6.4,2x,a8,2x,a8)
c
c               write(*,*) an,iobs1,iobs2,obs3
c    &          ,sigma_angle,sitea,point_viseb
c
c               print*, ' avec format'
c               write(*,120) an,iobs1,iobs2,obs3
c    &          ,sigma_angle,sitea,point_viseb
c              do kk = 1, nb_pts
               do kk = 1, nombre_pts
                  if (sitea.eq.station(kk)) then                     
                      write(28,*) phi(kk),' ',xlambda(kk)
     &               ,' ',station(kk)
c                      write(*,*)  phi(kk),' ',xlambda(kk)
                  endif
               enddo
c              do kk = 1, nb_pts
               do kk = 1, nombre_pts
                  if (point_viseb.eq.station(kk)) then
                      write(28,*) phi(kk),' ',xlambda(kk),
     &                ' ',station(kk)
                      write(28,*) '> '
cc                    write(*,*) phi(kk),xlambda(kk)
cc                    write(*,*) '> '
c                     write(*,*) phi(kk),xlambda(kk)
c                     write(*,*) '> '
                  endif   
               enddo
             enddo
         endif
c
c        itype_observ = 4  : SOLVEM code for baseline length
         if (itype_observ.eq.4) then
             do i = 1, nb_exp
                read(31,'(a80)',end=200) skipp
c               write(*,*) skipp
             enddo
         endif
         goto 100
      else
        stop
      endif

c
c     troubled stop
 130  continue
      print *, 'Cannot open sequential input file '
      print *, 'specified in your file list'
      stop
 140  continue
      print *, 'Cannot read the name of the sequential'
      print *, 'specified in the file list'
      stop
c
cc -----------------------------------------------------------
c
 200  continue 
      rewind (31)
      close(31)
      close(28)
c     rewind (14)
      rewind (14)
c      close(14)
c     print *, 'close 31 and 28, will rewind 14'
c     ios = 10
c      open(14,file='jr58a.new.fls',status='old',err=210,iostat=ios)
c     read the file list (#14)
      print *,'opening file list'
      read(14,'(a20)') skipp
c     print *,'   ',skipp
      backspace 14
c     read(14,40) nb_pts
      read(14,'(i)') nb_pts
c     print *,'nb pts= ',nb_pts
c     print *,'GMT_sketch.f finished'
      do i = 1, nb_pts
        read(14,'(a10)') skipp
c        print *, 'station ',skipp
c         print *,'i= ',i
c         backspace 14
      enddo
c      close(14)
c     
c     opening file list
c     read(14,'(a20)') skipp
c     print *,'   ',skipp
c     backspace 14
c     read(14,40) nb_pts
      return
c  --------------------------------------------------------
c
 210  continue
      print *,' suicide due to file list'
      stop
 300  print*,' mismatch site name at GMT_SKETCH: ',
     .   sitea,' ',point_viseb,nlen1,isite1,isite2
      stop
      end
