      subroutine remedy(idatf,record,illsit)
c
c     check every site, identify "ill" site.
c     all records have been sorted in subroutine recod.
c     at each site, at every epoch need two observations to 
c        determine 2-d position
c
c     unit:  omc (second or mm)
c
c

      include 'solvem.fti'
      integer isit,i1,i2,iuse(maxsit),vuse(maxsit),idatf
      integer loop,inew,illsit,tuse(maxsit),jnew,i,mc
cmk   integer irepe(maxsit),record(4000,4),isfix(maxsit,2)
      integer irepe(maxsit),isfix(maxsit,2)
      real*8 record(4000,4)
      logical sil(maxsit)
      equivalence (iuse,map)
      equivalence (vuse,map(maxsit+1))
c
c     temporary assuming all data are fine in list mode
      if (idatf.eq.14) then
         do isit = 1,nsit
            iuse(isit) = minic+1
            vuse(isit) = miniv+1
         enddo
         goto 100
      endif     
c
c     find out fixed site index
      do 70 isit = 1,nsit
         isfix(isit,1) = 0
         isfix(isit,2) = 0
         i1 = (isit-1)*6
         i2 = fix(i1+1)*fix(i1+2)*fix(i1+3)
         if (i2.gt.0) then
            isfix(isit,1) = 1
            iuse(isit) = minic
         endif
         i2 = fix(i1+4)*fix(i1+5)*fix(i1+6)
         if (i2.gt.0) then
            isfix(isit,2) = 1
            vuse(isit) = miniv
         endif
 70   continue
c
c     step 1: sort out observations by sites
      call sorto(record,sil,isfix,1)

      print*,'  Checking for sites with less ',
     .       'than 2 connections ...'
c
c     step 2: identify the site with less than 2 connections
      loop = 0
      illsit = 0
 150  inew = 0
      loop = loop+1
      do 20 isit = 1,nsit
         sil(isit) = .false.
         if (iuse(isit).lt.minic.or.vuse(isit).lt.miniv) then
            sil(isit)=.true.
            inew = inew+1
         endif
         if (isfix(isit,1).eq.0) iuse(isit) = 0
         if (isfix(isit,2).eq.0) vuse(isit) = 0
         irepe(isit) = 0
         tuse(isit) = 0
 20   continue        
c
      call sorto(record,sil,isfix,2)
c 
      jnew = 0
      do 30 isit = 1,nsit
         sil(isit) = .false.
         if (iuse(isit).lt.minic.or.vuse(isit).lt.miniv) then
            sil(isit)=.true.
            print*,sname(isit),' has less than the required obs!'
            jnew = jnew+1
         endif
 30   continue
c
      print*,' Loop (',loop,') ..... ',jnew,' ill site(s) found.'
      if (inew.lt.jnew) goto 150
      illsit = inew
      if (nsit-illsit.lt.2) stop ' No enough sites!'
      if (inew.eq.0) go to 100
c
c     step 3: add sites with weak observations to strengthen network
      loop = 0
      print*,' Add sites with weak observations to strengthen ',
     .   'the network. Adding fixed sites ...'
      mc = max(minic,1)
 110  inew = 0
      loop = loop+1
      do 40 isit = 1,nsit
         if(.not.sil(isit)) goto 40
         if (iuse(isit).ge.mc.or.vuse(isit).gt.0) then
            sil(isit) = .false.
            inew = inew+1
         endif
         if (isfix(isit,1).eq.0) iuse(isit) = 0
         if (isfix(isit,2).eq.0) vuse(isit) = 0
         irepe(isit) = 0
         tuse(isit) = 0
 40   continue        
c
      call sorto(record,sil,isfix,3)
c 
      jnew = 0
      do 50 isit = 1,nsit
         if(.not.sil(isit)) goto 50
         if (iuse(isit).ge.mc.or.vuse(isit).gt.0) then
            sil(isit) = .false.
            jnew = jnew+1
         endif
 50   continue   
c
      illsit = illsit-jnew
      print*,' Loop (',loop,') ..... ',illsit,' ill site(s) found.'
      if (inew.lt.jnew) goto 110
      if (nsit-illsit.lt.2) stop ' No enough sites!'
c
c     list the ill sites
      jnew = 0
      do 80 isit = 1,nsit
c        I don't understand this? kurt
ckf      if (.not.sil(isit)) goto 80
ckf      jnew = jnew+1
ckf      tuse(jnew) = isit
ckf      iuse(isit) = 0
ckf      vuse(isit) = 0
         if (sil(isit)) then
            jnew = jnew+1
            tuse(jnew) = isit
            iuse(isit) = 0
            vuse(isit) = 0
            write (6,*) tuse(jnew),'  ',sname(isit),isit
         endif
 80   continue
ckf   write (6,'(20i4)') (tuse(i),i=1,jnew)
      if (iomode(10).gt.0) write (24,'(10(1x,i3,1x,a4))') 
     .   (tuse(i),sname(tuse(i))(1:4),i=1,jnew)
c        
 100  continue


      return
      end
c     

      subroutine sorto(record,sil,isfix,mode)
c
c     mode = 1:  sort out observations by sites
c     mode = 2:  identify ill sites
c     mode = 3:  identify sites with weak data
c
cmk   !This is particularly difficult code to understand!
c
      include 'solvem.fti'
cmk   integer irec,i1,i2,it2,i0,it,iuse(maxsit),vuse(maxsit)
      integer irec,i1,i2,i0,it,iuse(maxsit),vuse(maxsit)
cmk   integer tuse(maxsit),mode,iold1,iold2
      integer mode,iold1,iold2
cmk   integer irepe(maxsit),record(4000,4),isfix(maxsit,2)
      integer irepe(maxsit),isfix(maxsit,2)
      real*8 record(4000,4),it2,tuse2(maxsit)
      real*8 omc,vest
      logical sil(maxsit)
      equivalence (iuse,map)
      equivalence (vuse,map(maxsit+1))
c
      i0 = 0
      iold1 = 0
      iold2 = 0
      do 40 irec = 1,jobs
         i1 = int(record(irec,1))
         i2 = int(record(irec,2))
         it2 = record(irec,3)
         it = int(record(irec,4))
c        it > 30 modes have no contributions
         if (it .gt. 30) goto 40
c        it > 20 modes 
         if (it .gt. 20 .and. it. le. 30) then
            if (mode.eq.2.and.sil(i2)) goto 40
            if (mode.eq.3.and..not.sil(i2)) goto 40
            if (isfix(i2,1).eq.1.and.isfix(i2,2).eq.1) goto 40
            if (it.ge.25.and.it.le.30) then
               irepe(i1) = irepe(i1)+3
               irepe(i2) = irepe(i2)+3
c               if (it-it/2*2.gt.0) then
                  iuse(i1) = iuse(i1)+3
                  iuse(i2) = iuse(i2)+3
c               endif
               vuse(i1) = vuse(i1)+3
               vuse(i2) = vuse(i2)+3
               tuse2(i1) = it2
               tuse2(i2) = it2
               go to 40
            endif
            if (i2.ne.i0.and.it2.ne.tuse2(i2)) then
               irepe(i2) = irepe(i2)+3
               if (it.eq.21.or.it.eq.23) iuse(i2) = iuse(i2)+3
               if (it.eq.25.or.it.eq.27) iuse(i2) = iuse(i2)+3
               if (it.eq.25.or.it.eq.27) iuse(i1) = iuse(i1)+3
               if (it.eq.22.or.it.eq.24) vuse(i2) = vuse(i2)+3
               if (it.eq.26.or.it.eq.28) vuse(i2) = vuse(i2)+3
               if (it.eq.26.or.it.eq.28) vuse(i1) = vuse(i1)+3
               if (it.gt.10.and.it.le.20) vuse(i2) = vuse(i2)+3
               if(irepe(i2).eq.1) then
                  i0 = i2
                  tuse2(i2) = it2
                  irepe(i2) = 0
               endif
            endif
            go to 40
         endif

c        azimuth observation
         if (it.eq.1) then
            if (mode.eq.2.and.sil(i2)) goto 40
            if (mode.eq.3.and..not.sil(i2)) goto 40
            if (isfix(i2,1).eq.1.and.isfix(i2,2).eq.1) goto 40
            if (it2.ne.tuse2(i2)) irepe(i2) = 0
            if (it2.eq.tuse2(i2).and.irepe(i2).ge.2) goto 40
            irepe(i2) = irepe(i2)+1
            tuse2(i2) = it2
            if(iuse(i2).lt.minic) then
               iuse(i2) = iuse(i2)+1
            else
               vuse(i2) = vuse(i2)+1
            endif
            go to 40
         endif

c        baseline length, zenith height and leveling
cmk      Also angle and direction, but is this right?
         if (it. ge. 2 . and. it .le. 6) then
            if (mode.eq.2.and.(sil(i1).or.sil(i2))) goto 40
            if (mode.eq.3.and..not.sil(i1).and..not.sil(i2)) goto 40
            if (mode.eq.3.and.sil(i1).and.sil(i2)) goto 40
            if (i1.eq.iold1.and.i2.eq.iold2) then

               if (it2.ne.tuse2(i1)) then
                  vest = (anorm(irec)-omc)/dble(nint(it2-tuse2(i1)))
cmk               !magic number! what is 1.0d2?
                  if (mode.eq.1.and.dabs(vest).gt.1.0d2) 
     .               print*,' V = ',i1,i2,vest
                  if (dabs(vest).gt.1.0d2) goto 40
               endif

               if (it2.eq.tuse2(i1).and.irepe(i1).ge.1) goto 40
               if (it2.ne.tuse2(i1).and.irepe(i1).ge.2) goto 60

            else
               if (it2.ne.tuse2(i1).or.irepe(i1).ge.2) irepe(i1) = 0
               if (it2.ne.tuse2(i2).or.irepe(i2).ge.2) irepe(i2) = 0
               iold1 = i1
               iold2 = i2
               omc = anorm(irec)
            endif

            if (isfix(i1,1).eq.1.and.isfix(i1,2).eq.1) goto 60
            if (it2.eq.tuse2(i1).and.irepe(i1).ge.2) goto 60
            irepe(i1) = irepe(i1)+1
            tuse2(i1) = it2

            if(iuse(i1).lt.minic) then
               iuse(i1) = iuse(i1)+1
            else
               vuse(i1) = vuse(i1)+1
            endif

 60         if (isfix(i2,1).eq.1.and.isfix(i2,2).eq.1) goto 40
            if (it2.eq.tuse2(i2).and.irepe(i2).ge.2) goto 40
            if (i1.eq.iold1.and.i2.eq.iold2.and.irepe(i2).ge.2) goto 40
            irepe(i2) = irepe(i2)+1
            tuse2(i2) = it2

            if(iuse(i2).lt.minic) then
               iuse(i2) = iuse(i2)+1
            else
               vuse(i2) = vuse(i2)+1
            endif

            go to 40
         endif
c
c        baseline length rate, zenith height rate and leveling rate
         if (it. ge. 12 . and. it .le. 16) then
            if (mode.eq.2.and.(sil(i1).or.sil(i2))) goto 40
            if (mode.eq.3.and..not.sil(i1).and..not.sil(i2)) goto 40
            if (mode.eq.3.and.sil(i1).and.sil(i2)) goto 40
            if (i1.eq.iold1.and.i2.eq.iold2) then
               if (it2.ne.tuse2(i1)) then
                  vest = (anorm(irec)-omc)
                  if (mode.eq.1.and.dabs(vest).gt.1.0d2) 
     .               print*,' V = ',i1,i2,vest
                  if (dabs(vest).gt.1.0d2) goto 40
               endif
               if (it2.eq.tuse2(i1).and.irepe(i1).ge.1) goto 80
            else
               if (it2.ne.tuse2(i1).or.irepe(i1).ge.1) irepe(i1) = 0
               if (it2.ne.tuse2(i2).or.irepe(i2).ge.1) irepe(i2) = 0
               iold1 = i1
               iold2 = i2
               omc = anorm(irec)
            endif
            if (isfix(i1,1).eq.1.and.isfix(i1,2).eq.1) goto 80
            if (it2.eq.tuse2(i1).and.irepe(i1).ge.1) goto 80
            irepe(i1) = irepe(i1)+1
            tuse2(i1) = it2
            if(iuse(i1).lt.minic) then
               iuse(i1) = iuse(i1)+1
            else
               vuse(i1) = vuse(i1)+1
            endif
 80         if (isfix(i2,1).eq.1.and.isfix(i2,2).eq.1) goto 40
            if (it2.eq.tuse2(i2).and.irepe(i2).ge.1) goto 40
            if (i1.eq.iold1.and.i2.eq.iold2.and.irepe(i2).ge.1) goto 40
            irepe(i2) = irepe(i2)+1
            tuse2(i2) = it2
            if(iuse(i2).lt.minic) then
               iuse(i2) = iuse(i2)+1
            else
               vuse(i2) = vuse(i2)+1
            endif
            go to 40
         endif
c
 40   continue
      return
      end

