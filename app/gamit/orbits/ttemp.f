	subroutine ttemp(luo,tfname,spfmt,nprn,nrec,isw,lunit,nintrs)
c
c purpose:	read NGS data records assuming the header has been
c		skipped, then store those records in t-file format
c               for interpolation by gsatel.f
c
c by Simon McClusky MIT Jan 1994
c
        implicit none
c
      include '../includes/dimpar.h'
	character*3 spfmt(2)
        character*7 stat
        character*24 tfname
        integer iyr,mon,iday,ihr,min,lunit,i,j,k,maxepc
	integer luo,isw,nprn(2),nrec(2),ioerr,ipksat,nintrs(2)
        parameter (maxepc=673)
        real*8 torb(3,maxepc),sec
c	data sc/1.e3/
c
c set up temporary tfile information
c
        lunit=luo+10
        stat='unknown'
c
c open temporary tfile
c
        call topens(tfname,stat,lunit,ioerr)
c
c read ngs file data records (SP1 format)
c
	do i=1,nrec(isw)
	if (spfmt(isw).eq.'sp1') then
		read(luo,900) iyr,mon,iday,ihr,min,sec
900	format(4x,i4,4i3,1x,f10.7)
			do j=1,nprn(isw)
			read(luo,910) (torb(k,j),k=1,3)
910	format(5x,3f13.5)
		        enddo
c
c write out temporary tfile data
c
            WRITE(lunit,ERR=993,IOSTAT=IOERR)
     .      ((torb(k,IPKSAT),k=1,3),IPKSAT=1,nprn(isw))

        nintrs(isw) = 3
c
c read ngs file data records (SP3 format)
c
	else
		read(luo,930) iyr,mon,iday,ihr,min,sec
930	format(3x,i4,4i3,f12.8)
			do j=1,nprn(isw)
			read(luo,940) (torb(k,j),k=1,3)
940	format(4x,3f14.6)
			enddo
            WRITE(lunit,ERR=993,IOSTAT=IOERR)
     .      ((torb(k,IPKSAT),k=1,3),IPKSAT=1,nprn(isw))

        nintrs(isw)=3
	endif
	enddo
	return
993     write(*,995)tfname
995     format(1x,'Error writing out temporary tfile ',a8)
	end
