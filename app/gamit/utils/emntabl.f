c     Create Nutation, Sun, and Moon Tables for GAMIT processing
c     from a PEP N-body Ephemeris.   R. King November 1984

      integer*4 kkmax,kksmax
      parameter(kkmax=1200)
      parameter(kksmax=150)

      integer*4 jd1,jd2,span,name(6,10),jdbd0(10),ibody,inut,isun
     .        , imoon,iout,intrd,intnut,intmon,intsun
     .        , jds1,jds2,nmoon,jm,jdbd,idir,nrec,ivl,mvl,nvel,itrt
     .        , jdbd1,jdbd2,nbdy2,npg,jvlbd,kkm,kkn,kks,kk,i,j,k
      data      jd1/0/,jd2/99999999/,span/0/
     .        , ibody/10/inut/11/,isun/12/,imoon/13/,iout/8/
     .        , intrd/0/,intnut/-1/,intmon/-1/,intsun/4/

      integer*2 nplnt,npl(10),ncp(10),intb(10),nbdy,nbdy1,intbd
     .        , kbd(40),nskip,kskip
      data nplnt/0/,nskip/0/,kskip/9999/
c
      real*8 mass1(10),relftb(10),beta(6,10),mass(10)
     .     , frct,merc(6,10),body(6,5,8),mon(6,40)
     .     , tti,convds,ltvel,tts,aultsc,fract,aukm,aukmvl
     .     , izero(100),moon(6,kkmax),earth(6,kksmax),tmp3(3)
      data   convds/4.8481368110953599d-6/,ltvel/2.99792458d5/
     .     , aultsc/499.00478370d0/,fract/0.d0/,izero/100*0.d0/

      real*4 t1(kkmax),dpsi(kkmax),deps(kkmax),units,psid(40),epsd(40)
     .     , librt(40,3),ts(kksmax),epsbd
      data units/1./

      character*1 ans,anssun,ansmon,ansnut,ansrot,dump,lowerc
      character*5 frame
      character*40 nbody,soltab,luntab,nutabl
      data nbody/'nbody740'/,soltab/'soltab'/,luntab/'luntab'/
     .     ,nutabl/'nutabl'/
      character*80 titsun,titmon,titnut
      character*128 title

      logical first,skphdr
      data    first/.true./,skphdr/.false./

      real*8 a2000(3,3)
      DATA A2000/.9999257079523629D0, .0111789381264276D0,
     1           .0048590038414544D0, -.0111789381377700D0,
     2           .9999375133499888D0, -.0000271579262585D0,
     3           -.0048590038153592D0, -.0000271625947142D0,
     4           .9999881946023742D0/

      write(6,'(a)')'Program to write GAMIT tables from PEP N-body file'
      write(6,'(a)') 'Do you want to write a Sun table (y/n)?'
      read(5,'(a)') anssun
      if( lowerc(anssun).eq.'y') then
        open( unit=isun,file=soltab,form='formatted',status='unknown'
     .    , err=9020 )
        write(6,'(a)') 'Enter title of Sun table'
        read(5,'(a80)') titsun
      endif
      write(6,'(a)') 'Do you want to write a Moon table (y/n)?'
      read(5,'(a)') ansmon
      if( lowerc(ansmon).eq.'y') then
        open( unit=imoon,file=luntab,form='formatted',status='unknown'
     .    , err=9030 )
        write(6,'(a)') 'Enter title of Moon table'
        read(5,'(a80)') titmon
      endif
      write(6,'(a)') 'Do you want to write a nutation table (y/n)?'
      read(5,'(a)') ansnut
      if( lowerc(ansnut).eq.'y') then
        open( unit=inut,file=nutabl,form='formatted',status='unknown'
     .    , err=9040 )
        write(6,'(a)') 'Enter title of nutation table'
        read(5,'(a80)') titnut
      endif
      write(*,'(a)') 'Enter jdstart jdstop (integers - free format)'
      read(5,*) jd1,jd2
      write(*,'(a)') 'Do you want to print out values? (y/n)'
      read(5,'(a)') dump
      if ( lowerc(dump).eq.'y') then
         write(*,'(a)') 'Print every Nth record - Enter N'
         read(5,*) nskip
      endif

c         Initialization

      kk= 0


c         Open the N-body file and read the header records

      open( unit=ibody,file=nbody,form='unformatted',status='old'
     .    , err=9010 )
  100 read(ibody) title
      write(*,'(/,1x,a,/,1x,a128)') 'Title on N-body file:',title
      read(ibody) nbdy,(npl(i),i=1,nbdy),(ncp(i),i=1,nbdy)
     .          , (intb(i),i=1,nbdy),jdbd1,jdbd2,(jdbd0(i),i=1,nbdy)
     .          , ((beta(i,j),i=1,6),j=1,nbdy)
     .          , ((name(i,j),i=1,6),j=1,nbdy)
     .          , nmoon,nbdy1,intbd,jvlbd,epsbd,(kbd(i),i=1,40),itrt,npg
     .          , (mass1(i),i=1,nbdy),(relftb(i),i=1,nbdy)

c        Set reference frame of N-body

  110  write(*,'(a)') 'Input frame of N-body ephemeris (B1950 or J2000)'
       read(5,'(a5)') frame
       call uppers(frame)
       if( frame.ne.'B1950' .and. frame.ne.'J2000') then
          write(*,'(a)') 'Frame invalid'
          goto 110
       endif
       write(*,'(a)') 'Input frame is ',frame
       if( frame.ne.'J2000') then
          write(*,'(a)') 'Do you want to rotate to J2000'
          read(5,'(a)') ansrot
          if( lowerc(ansrot).eq.'y') frame = 'J2000'
        endif


c        Compute masses and constants

      if( mass1(10).eq.0.d0) mass1(10) = 82.3005883d0
      do i=1,nbdy
        mass(i) = 0.d0
        if( mass1(i).gt.1.d0 ) mass(i) = 1.d0/mass1(i)
      enddo
      if( lowerc(anssun).eq.'y' .or. lowerc(ansmon).eq.'y' ) then
         write(*,'(/,a,f13.8,a)') 'AULTSC= ',aultsc
     .        ,'  Do you want to change it? (y/n)'
         read(5,'(a)') ans
         if( lowerc(ans).eq.'y') then
            write(*,'(a)') 'Enter new value of AULTSC:'
            read(5,*) aultsc
         endif
         aukm = ltvel*aultsc
         aukmvl = aukm/86400.d0
         write(*,'(a)') 'Do you want velocities for earth & moon? (y/n)'
         nvel = 0
         read(5,'(a)') ans
         if( lowerc(ans).eq.'y') nvel = 1
      endif


c       Write out the N-body headers

      write(*,'(/,a,6x,2i10,i5)') 'JDBD1,JDBD2,NBDY=',jdbd1,jdbd2,nbdy
      idir = 1
      if(jdbd1.gt.jdbd2) idir = -1
      if( .not.skphdr) then
        write(*,'(/,a)')'   NPL     NPC           NAME            MASS1'
        write(*,'(/,(2x,2i5,2x,6a4,4x,1pd22.15))')  ( npl(k),ncp(k)
     .       , (name(j,k),j=1,6),mass1(k),k=1,nbdy )
        write(*,'(/,a,10i5,/,3(5x,10i5,/))') 'kbd= ',kbd
        write(*,'(a,3x,6i3,1x,1pd10.3)')
     .          'nmoon,nbdy1,intbd,jvlbd,intrd,npg,epsbd'
     .       ,   nmoon,nbdy1,intbd,jvlbd,intrd,npg,epsbd
      endif

      write(*,'(a,2i10)') 'JD1 JD2 requested ',jd1,jd2
      jds1 = max0 ( min0(jdbd1,jdbd2), jd1-20 )
      jds2 = min0 ( max0(jdbd1,jdbd2), jd2+20 )
      if( jd1.lt.0 ) jds1 = -jd1-20
      if( jd2.lt.0 ) jds2 = -jd2+20
      jdbd1 = jds1
      jdbd2 = jds2
      if( idir.le.0 ) then
         jdbd1 = jds2
         jdbd2 = jds1
      endif
      write(*,'(a,2i10)') 'JDBD1 JDBD2 selected: ',jdbd1,jdbd2
      nrec = 0
      if( jdbd1.ge.jdbd2 ) goto 9001
      nbdy2 = nbdy1 -1

  400 read(ibody,end=9000)jdbd,frct,ivl,mvl
     .    , ((merc(i,j),i=1,ivl),j=1,10)
     .    , (((body(i,j,k),i=1,ivl),j=1,5),k=1,nbdy2)
     .    , ((mon(i,j),i=1,mvl),j=1,40)
     .    , (psid(j),epsd(j),j=1,40), ((librt(j,i),i=1,3),j=1,40)

      if( (jdbd-jdbd1)*idir.lt.0 ) goto 400
      if( first) write(*,'(a,i8,/)') 'Start of copy at ',jdbd
      first = .false.

c       Fill arrays for output

      if( lowerc(anssun).eq.'y' ) then
        tts = 0.d0
        do j=1,5
          jm= 8*j - 7
          kks = kks + 1
          if( kks.gt.kksmax ) goto 9060
          ts(kks)= dfloat(jdbd-2400000) + tts
          do i=1,3
            earth(i,kks) = (body(i,j,2)-mass(10)*mon(i,jm))*aukm
            if(nvel.gt.0) earth(i+3,kks) =
     .                    (body(i+3,j,2)-mass(10)*mon(i+3,jm))*aukmvl
          enddo
          if( lowerc(ansrot).eq.'y' ) then
            call matmpy( a2000,earth(1,kks),tmp3,3,3,1)
            do i=1,3
            earth(i,kks) = tmp3(i)
            enddo
          endif
        tts = tts + 4.d0
        enddo

      endif

      if( lowerc(ansmon).eq.'y' .or. lowerc(ansnut).eq.'y' ) then
        tti = 0.d0
        do j=1,40
          kk = kk + 1
          if( kk.gt.kkmax ) goto 9050
          t1(kk) = dfloat(jdbd-2400000) + tti
          if( lowerc(ansnut).eq.'y' ) then
            dpsi(kk) = psid(j)/convds
            deps(kk) = epsd(j)/convds
          endif
          if( lowerc(ansmon).eq.'y' ) then
            do i=1,3
              moon(i,kk) = mon(i,j)*aukm
              if( nvel.gt.0 ) moon(i+3,kk)= mon(i+3,j)*aukmvl
            enddo
            if( lowerc(ansrot).eq.'y' ) then
               call matmpy( a2000,moon(1,kk),tmp3,3,3,1)
               do i=1,3
                 moon(i,kk) = tmp3(i)
               enddo
             endif
          endif
        tti = tti + 0.5
        enddo
      endif

      nrec = nrec + 1
      if( lowerc(dump).eq.'y' ) then
        kskip = kskip + 1
        if( kskip.ge.nskip ) then
           kskip = 0
           write(*,'(/,a,i10,f12.7,a)') 'JDBD, FRACT ',jdbd,fract
           write(iout,'(/,a,i10,f12.7,a,2i3)') 'JDBD, FRACT ',jdbd,fract
     .           ,'  IVL, MVL= ',ivl,mvl
        endif
        if( lowerc(anssun).eq.'y' ) then
           write(iout,'(a,/,(1x,6d20.11))')
     .        'Earth',((body(i,j,2),i=1,6),j=1,5)
        endif
        if( lowerc(ansmon).eq.'y' ) then
           write(iout,'(a,/,(1x,6d20.11))')
     .        'Moon',((mon(i,j),i=1,6),j=1,40)
        endif
        if( lowerc(ansnut).eq.'y' ) then
            write( iout,'(a,(1x,3e15.8,2x,2e15.8))')
     .             'Psi Eps',(psid(i),epsd(i),i=1,40)
        endif
      endif

      if( (jdbd-jdbd2)*idir.ge.0 ) goto 9001
      goto 400

 9000 print *,'End-of-file reached on N-BODY'
 9001 write(6,'(/,a,i8,a,i10,/)') 'End of run at ',jdbd
     .     , ' Records read from N-body file',nrec

      if( lowerc(ansnut).eq.'y' ) kkn = kk
      if( lowerc(ansmon).eq.'y' ) kkm = kk

      write(6,'(a,3i5)')'Number of (Sun,Moon,Nutation) values selected='
     .                 , kks,kkm,kkn


c       Write the output tables

      if( lowerc(anssun).eq.'y' )
     . call tabsun(isun,titsun,intsun,1.d0,1.d0,frame,nvel,kks,ts,earth)
      if( lowerc(ansmon).eq.'y' )
     . call tabmon(imoon,titmon,intmon,1.d0,1.d-3,frame,nvel,kk,t1,moon)
      if( lowerc(ansnut).eq.'y' )
     .   call tabnut( inut,titnut,intnut,units,kk,t1,dpsi,deps )
      stop

 9010 write(*,'(a,a40)') 'Error opening nbody file = ',nbody
      stop
 9020 write(*,'(a,a40)') 'Error opening Sun file = ',soltab
      stop
 9030 write(*,'(a,a40)') 'Error opening Moon file = ',luntab
      stop
 9040 write(*,'(a,a40)') 'Error opening Nutation file = ',nutabl
      stop
 9050 write(6,'(a,i5,a,i5)') 'Number of values selected=',kk
     .                     ,  '  gt kmax=',kkmax
      stop
 9060 write(6,'(a,i5,a,i5)') 'Number of values selected=',kks
     .                     ,  '  gt kksmax=',kksmax
      stop
      end

c
      subroutine tabsun( iout,title,int,units,outunt,frame,nvel,n,t,x )

c     Write a GAMIT-format output Sun table from an array of values.

c     iout      output unit number
c     title     title for output file
c     int       interval (days) of values in output table
c     units     input units (km)   (output is km)
c     outunt
c     nvel      flag for coordinates-only (=0) or velocities (=1)
c     n         number of values in input arrays
c     t         array of input dates (NB: PEP JD, not MJD)
c     x         array of input coordinates

      integer*4 kksmax
      parameter( kksmax=150 )

      integer*4 ix(6),npjd,ifix,jd1,jd2,nvel,iout,int,npr,i,j,n

      real*8 x(6,kksmax),units,outunt,xx

      real*4 t(kksmax)

      character*80 title
      character*24 varfmt
      character*5  frame
      data         varfmt/'(1x,i5,6i11)            '/

      jd1 = ifix(t(1)) + 2400000
      jd2 = ifix(t(n)) + 2400000
      npr = 3
      if( nvel.gt.0 ) npr = 6

      write(iout,10) title,varfmt,nvel,jd1,jd2,npr,int,outunt,frame
   10 format(a80,/,a24,4x,i2,1x,i7,1x,i7,1x,i2,1x,i2,1x,1pe15.0,1x,a5)

      do  j = 1,n
        do i=1,npr
          xx = x(i,j)*units/outunt
          ix(i) = idint(dsign(dabs(xx)+0.5d0,xx))
        enddo
        npjd = ifix(t(j))
        write(iout,varfmt) npjd,(ix(i),i=1,npr)
      enddo

      return
      end


      subroutine tabmon( iout,title,int,units,outunt,frame,nvel,n,t,x )

c     Write a GAMIT-format output Moon table from an array of values.

c     iout      output unit number
c     title     title for output file
c     int       interval (days) of values in output table
c     units     input units (km)   (output is km)
c     outunt
c     nvel      flag for coordinates-only (=0) or velocities (=1)
c     n         number of values in input arrays
c     t         array of input dates (NB: PEP JD, not MJD)
c     x         array of input coordinates

      integer*4 kkmax
      parameter( kkmax=1200 )

      integer*4 ix(6),npjd,ifix,jd1,jd2,nvel,iout,int,npr,i,j,n

      real*8 x(6,kkmax),units,outunt,xx

      real*4 t(kkmax)

      character*80 title
      character*24 varfmt
      character*5  frame
      data         varfmt/'(1x,i5,6i11)            '/

      jd1 = ifix(t(1)) + 2400000
      jd2 = ifix(t(n)) + 2400000
      npr = 3
      if( nvel.gt.0 ) npr = 6

      write(iout,10) title,varfmt,nvel,jd1,jd2,npr,int,outunt,frame
   10 format(a80,/,a24,4x,i2,1x,i7,1x,i7,1x,i2,1x,i2,1x,1pe15.0,1x,a5)

      do  j = 1,n
        do i=1,npr
          xx = x(i,j)*units/outunt
          ix(i) = idint(dsign(dabs(xx)+0.5d0,xx))
        enddo
        npjd = ifix(t(j))
        write(iout,varfmt) npjd,(ix(i),i=1,npr)
      enddo

      return
      end


      subroutine tabnut( iout,title,int,units,n,t,x,y )

      integer*4 kkmax
      parameter( kkmax=1200 )

c
c     iout      output unit number
c     title     title for output file
c     int       interval (days) of values in output table
c     units     input units (km)   (output is 1.e-4 sec)
c     n         number of values in input arrays
c     t         array of input dates (NB: PEP JD, not MJD)
c     x         array of input values of psi
c     y         array of input values of eps

      integer*4 ix(4),iy(4),npjd(4),ifix,npr,jd1,jd2,int,iout
     .        , i,j,k,n
      data npr/4/

      real*4 t(kkmax),x(kkmax),y(kkmax),units,outunt
      data outunt/1.e-4/

      real*8 xx,yy

      character*80 title
      character*24  varfmt
      data          varfmt/'(1x,i5,8i8,8x,i2)       '/

      jd1 = ifix(t(1)) + 2400000
      jd2 = ifix(t(n)) + 2400000

      write(iout,10) title,varfmt,jd1,jd2,npr,int,outunt
   10 format(a80,/,a24,11x,i7,1x,i7,1x,i2,1x,i2,1x,1pe15.0)

      i = 0
   20 k = 0
   30 i = i + 1
      if( i.gt.n ) goto 50
      k = k + 1
      xx = x(i)*units/outunt
      yy = y(i)*units/outunt
      ix(k) = idint(dsign(dabs(xx)+0.5,xx))
      iy(k) = idint(dsign(dabs(yy)+0.5,yy))
      npjd(k) = ifix(t(i))
      if( k.ne.4 ) goto 30
      write(iout,varfmt) npjd(1),(ix(j),iy(j),j=1,4)
      goto 20
   50 continue
      if( k.eq.0 ) goto 70
      write(iout,varfmt) npjd(1),(ix(j),iy(j),j=1,k)
      write(6,60) k
   60 format(/,1x,'Need to edit nutation table to add ',i1, 'in col 80')

   70 return
      end















