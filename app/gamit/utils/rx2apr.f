	Program RX2APR

c	Program to covert grep'ed XYZ coordinates from RINEX files
c	appropriate format for input into gapr_to_l
c
c	usage: RX2APR flnm
c
c     files:
c	flnm = file(s) of xyz coordinates from rinex files
c     flnmapr  = file of formatted xyz coordinates and velocities
c
c     variables:
c     xvel,yvel,zvel unknown station velocities; set to zero
c     xpos,ypos,zpos station positions from rinex file

c  written by Feng Shen at MIT  Mar 1994
c       
      implicit none

	character*128 flnm,flnmapr
	character*3 day
	character*4 yr
	integer*4 len_run,rs,rcpar,ierr,trimlen,iyr,iday
	logical file_stat
	real*4 epoch 
c
	len_run=1
	file_stat=.true.
	rs = 0
c	check runstring
c	get filename

       rs=rs+1
      len_run=rcpar(rs,flnm)
      if(rs .eq. 1 .and. len_run .eq. 0) then
        print *,' '
        print *,'To run RX2APR enter: rx2apr xyz.rnx <year> <day>'
        print *,' '
        print *,'  where xyz.rnx is produced by the UNIX command '
        print *,' '
        print *,'    grep POSITION *.93o > xyz.rnx '
        print *, ' '
        print *,'  and <year> <day> specify the reference epoch for'
        print *,'  the output file xyz.rnx.apr'
        stop
      endif

c	get year

	    rs=rs+1
	    len_run=rcpar(rs,yr)
	    if(rs .eq. 2 .and. len_run .eq. 0) then
	        print*,'RX2APR: incomplete runstring'
	        print*,'To run program: rx2apr <filename> <year> <day>'
	        stop
	    endif

c	get day

	    rs=rs+1
	    len_run=rcpar(rs,day)
	    if(rs .eq. 3 .and. len_run .eq. 0) then
	        print*,'RX2APR: incomplete runstring'
	        print*,'To run program: rx2apr <filenames> <year> <day>'
	        stop
	    endif

	    read(yr,'(i4)') iyr
	    read(day,'(i4)') iday
          epoch=float(iyr) + (float(iday)/365.25)

c	    open file to be edited

	    open(unit=10,file=flnm,status='old',iostat=ierr)
	    if(ierr .ne. 0) then
	       print*,'Error opening input file',flnm,ierr
	       file_stat=.false.
	    endif
c
	    if(file_stat .eqv. .true.) then

c	    open new file

	    flnmapr=flnm(1:trimlen(flnm))//'.apr'
	    if(ierr .eq. 0) then
	        open(unit=20,file=flnmapr,status='unknown',iostat=ierr)
	    endif
	    if(ierr .eq. 0) then
              call xyzedit (10,20,ierr,flnm,epoch)
	    endif
      close (20)
	   endif
      write(*,'(/,2a)') 'Created file ',flnmapr
	end
	subroutine xyzedit (nunit,munit,ierr,flnm,epoch)
c	variable definition:
c		nunit=input file number
c		munit=output file number
c		ierr=error indicator on file reads
c		inline=128 character character string
c		flnm=name of file to read from
c	   	xpos,ypos,zpos=xyz coordinates of local stations from rinex files
c		xvel,yvel,zvel=velocities of stations; for new stations, these are
c				 unknown so they are set to be zero
c		xave,yave,zave=average position of sations if occupied more than once
c			        during the campaign
c		ncount=number of occupations at a station in a campaign
c
	integer*4 nunit,munit,ierr,ncount,indref,indxb,indxe,line_count
	real*4 xvel,yvel,zvel,epoch,a,c
	real*8 xpos,ypos,zpos,xave,yave,zave
	character*(*) flnm
	character*128 inline
	character*4 sta1,sta2,post
      logical eof /.false./

c	define ellipsoid axes; used to delineate really bad outliers

	a=6.38816e6
	c=6.34618e6

c	define velocities

	xvel=0.
	yvel=0.
	zvel=0.

c	define name tail necessary for translation using program glbtol

	post='_GPS'

c	read first line of file

      read(nunit,'(a)',iostat=ierr) inline
      line_count = 1
      if(ierr .gt. 0) then
	      print*,'error reading ',flnm
	      print*,'last line read ',inline
      endif
      if(ierr .lt. 0) then
         print*,'err EOF reading file ',flnm
      endif

c	determine index of station name and xyz coordinates in character string inline

c	indref=index(inline,'O:')
      indref=index(inline,'o:')
      indxb=indref-11
      indxe=indref-6

c	find station name

	sta1=inline(indxb:indxe)


c	find station coordinates
	read (inline (indref+3:),*) xpos,ypos,zpos
c	print*,'sta1 x y z ',sta1,xpos,ypos,zpos


c	initialize averages

      xave=xpos
      yave=ypos
      zave=zpos
      ncount=1

      do while (ierr .eq. 0)
	     read(nunit,'(a)',iostat=ierr) inline
            line_count = line_count + 1

	     if(ierr .gt. 0) then
	        print*,'error reading ',flnm
	        print*,'last line read ',inline
c	        reading error, get out of loop
	    endif
	    if(ierr .lt. 0) then
           if( line_count.le.2 ) then
	          write(*,'(a,i5,a)') 'Unexpected EOF after '
     .                        ,  line_count,' records'
             stop
           else
             eof = .true.
           endif
	    endif

	    sta2=inline(indxb:indxe)
	    read (inline (indref+3:),*) xpos,ypos,zpos
c            print*,'xpos= ,ypos=, zpos=', xpos,ypos,zpos

	    if(sta2 .ne. sta1 .or. ierr .lt. 0) then

c	      rewrite lines to format appropriate for a *.apr file
		  xave=xave/float(ncount)
         yave=yave/float(ncount)
         zave=zave/float(ncount)
c         print*,' xavefinal=,yave=, zave=',xave,yave,zave

         write(20,100) sta1,post,xave,yave,zave,xvel,yvel,zvel,epoch
100      format(1x,2a4,3(1x,f13.4),3(1x,f7.4),1x,f9.4)
          if (eof) goto 99

c	        reset averages

	        xave=xpos
	        yave=ypos
	        zave=zpos
	        ncount=1
c                print*,'xaveset=,yaveset=,zaveset=',xave,yave,zave
	    else
	        xave=xave+xpos
	        yave=yave+ypos
	        zave=zave+zpos
	        ncount=ncount+1

c                print*,'ncount xaveadd=, yaveadd=,zaveadd='
c     .                , nount,xave,yave,zave

c	        check if positions are non-zero; if any are zero, remove them from average

	        if(xpos .eq. 0 .or. ypos .eq. 0 .or. zpos .eq. 0) then
	            xave=xave-xpos
	            yave=yave-ypos
	            zave=zave-zpos
	            ncount=ncount-1
	        endif
	    endif
	    sta1=sta2
      end do
   99 write(*,'(/,a,i3,a)') 'Normal end after reading ',line_count
     .      ,' records'
	close (munit)
	return
	end

