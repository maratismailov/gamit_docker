	subroutine rdspdt(luo,statim,endtim,spfmt,nprn,isw,tag,orb,nrec)
c
c purpose:	read NGS data records assuming the header has been
c		skipped, then store those records within a specified
c		time range (between statim and endtim)
c
c by Peng Fang, pfang@pgga.ucsd.edu, April 93
c
      integer*4 maxsat,maxepc
        parameter (maxsat=24,maxepc=673)
	character*3 spfmt(2)
	integer*4 luo,isw,nprn(2),nrec(2)
      integer*4 julday,iyr,mon,iday,ihr,min,ii,i,j,k
	real*8 sec,statim,endtim,rectim,sc
	real*8 tag(2,maxepc),orb(4,maxsat,maxepc,3)
	data sc/1.e3/
	ii=1
	do i=1,nrec(isw)
	if (spfmt(isw).eq.'sp1') then
		read(luo,900) iyr,mon,iday,ihr,min,sec
900	format(4x,i4,4i3,1x,f10.7)
	rectim=julday(mon,iday,iyr)-2400001+(ihr+(min+sec/60)/60.)/24.
		if (rectim.ge.statim.and.rectim.le.endtim) then
			tag(isw,ii)=rectim
			do j=1,nprn(isw)
			read(luo,910) (orb(isw,j,ii,k),k=1,3)
910	format(5x,3f13.5)
			do k=1,3
			orb(isw,j,ii,k)=orb(isw,j,ii,k)*sc
			enddo
			enddo
			ii=ii+1
		else
			do j=1,nprn(isw)
			read(luo,'(1x)')
			enddo
		endif
	else
		read(luo,930) iyr,mon,iday,ihr,min,sec
930	format(3x,i4,4i3,f12.8)
	rectim=julday(mon,iday,iyr)-2400001+(ihr+(min+sec/60)/60.)/24.
		if (rectim.ge.statim.and.rectim.le.endtim) then
			tag(isw,ii)=rectim
			do j=1,nprn(isw)
			read(luo,940) (orb(isw,j,ii,k),k=1,3)
940	format(4x,3f14.6)
			do k=1,3
			orb(isw,j,ii,k)=orb(isw,j,ii,k)*sc
			enddo
			enddo
			ii=ii+1
		else
			do j=1,nprn(isw)
			read(luo,'(1x)')
			enddo
		endif
	endif
	enddo
	nrec(isw)=ii-1
	return
	end
