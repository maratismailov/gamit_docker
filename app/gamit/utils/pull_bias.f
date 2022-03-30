      Program pull_bias

c    K. Feigl, B. Oral, R. King  1990-1992

c**   compile a list of unfixed widelanes

      integer   i,numfid,iqfile,nqfile,inum,j,num
      parameter (num=100)
      character qfile(100)*80,text*80,key*80
      character ion_con*5,nor_rms*12
      character fids(num)*4,save(100)*15
      character infile*80
      logical fcheck

c*    added following lines to read the input file  containing  sites names to be ignored
c*    also modified parameter (num=9000) and character fids(num)*4
c      open(unit=15,file="/p1data/gps/util/sccs/pull_bias.in "
c     .     ,status='unknown')
c      open(unit=15,file='~/gu/util/sccs/pull_bias.in',status='unknown')

c       Input the list of sites to exclude
          
      infile = ' ' 
      infile(1:14) = './pull_bias.in'
      if( .not.fcheck(infile)) then
         infile(1:22) = '../tables/pull_bias.in'
         if( .not.fcheck(infile) ) then
           write(6,'(/,a)')
     .       ' No pull_bias.in file available, use all sites'
           goto 5
         endif
      endif
      print *,' Opening ',infile
      open(unit=15,file=infile,status='unknown')
      numfid = 0
      do 1 i=1,num
        read(15,'(a4)',end=2) fids(i)
        if( fids(i).eq.'    ' ) goto 2
        numfid = numfid + 1
    1   continue
    2 write(6,'(/,i3,a)') numfid,' sites excluded:'
      write(6,'(i5,5x,a4)') (i,fids(i),i=1,numfid)
      close(15)

c      Input Q-file names

    5 write(6,'(/,a)')
     .    'Enter Q-files, one line at a time, ending with EOF or blank'
      i = 1
    6 read(5,'(a)',end=7) qfile(i)
      if( qfile(i).eq.' ' ) goto 7
      i = i+1
      goto 6
    7 nqfile = i-1

      do iqfile = 1, nqfile
      open(1,file=qfile(iqfile),status='old')

      key(1:14) = '              '
      do while (key(1:14).ne.'L1,L2 estimate')
         read(1,'(2x,2a)') key(1:14),text
      enddo

*     decode ionosphere constraint
      write(ion_con,'(a)') text(16:20)

*     decode normalized RMS
      read(1,'(a)') (text,i=1,3)
      write(nor_rms,'(a)') text(37:48)
      print*,qfile(iqfile)(1:14),' ppm: ',ion_con,' RMS: ',nor_rms

*     skip to widelane fixed section
      key(1:13) = '             '
      do while (key(1:13).ne.'L2-L1 biases '
     .    .and. key(1:13).ne.'   == summary')
         read(1,'(2x,a)') key(1:13)
      enddo

      do while (key(1:5).ne.'Label')
         read(1,'(8x,a)') key(1:5)
      enddo

*     find unfixed wide lanes
      inum = 0
      do while (key(5:5).ne.' ')
         read(1,'(2x,2a)') key(1:9),text
         if (key(4:4).eq.'*' .and. key(7:9).eq.'L21') then

*           write only regionals
            do j = 1, numfid

               if    (text(1:4).eq.fids(j)
     .           .or. text(6:9).eq.fids(j)) goto 10
            enddo
            inum = inum+1
            print '(1x,a9,1x,a)',key(1:9),text
            save(inum) = text(1:15)
10          continue
         endif
      enddo

*     find unfixed regional narrowlanes not included in above set
      do i = 1, 2
         do while (key(1:4).ne.'KEYS')
            read(1,'(1x,a4)') key(1:4)
         enddo
         key(1:9) = '         '
      enddo

      do while (key(5:9).ne.'B1L1 ')
         read(1,'(2x,2a)') key(1:9),text
      enddo

      do while (key(5:9).eq.'B1L1 ')
         if (key(4:4).eq.'*') then

*           check for fiducials and widelane frees
            do j = 1, numfid
               if    (text(1:4).eq.fids(j)
     .           .or. text(6:9).eq.fids(j)) goto 20
            enddo
            do j = 1, inum
               if (text(1:15).eq.save(j)) goto 20
            enddo
            print '(1x,a9,1x,a)',key(1:9),text
20          continue
         endif
         read(1,'(2x,2a)') key(1:9),text
      enddo

      close(1)
      enddo

      end

