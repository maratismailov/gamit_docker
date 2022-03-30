	subroutine rdsphd(lusp,mjd,fmjd,iprn,isw,
     *		delt,nrec,nprn,orbtyp,org,crdsys,spfmt)
c
c purpose:	read in sp1 or sp3 header information
c
	character*3 orbtyp(2),spfmt(2)
	character*4 org(2)
	character*5 crdsys(2)
	integer lusp,mjd(2),nrec(2),nprn(2),iprn(2,85),isw,iyr,mon
     .  , iday,ihr,min,igwk,i
	real*8 sec,delt(2),fmjd(2),dumy

c check the format

	read(lusp,'(a)') spfmt(isw)
	read(lusp,'(a)') spfmt(isw)
	rewind (lusp)

	if (spfmt(isw)(1:2).ne.'##') then
c
c read in sp1 format
	read(lusp,900) iyr,mon,iday,ihr,min,sec,delt(isw),mjd(isw),
     *		fmjd(isw),nrec(isw),orbtyp(isw),org(isw)
900	format(4x,i4,4i3,f11.7,f15.7,i6,f16.14,1x,i6,a1,1x,a3)
	read(lusp,910) nprn(isw),(iprn(isw,i),i=1,34),
     *		crdsys(isw)(4:5),igwk
910	format(3x,i2,1x,34i2,a2,i4)
	spfmt(isw) = 'sp1'

	else

c read in sp3 format
	read(lusp,920) iyr,mon,iday,ihr,min,sec,nrec(isw),
     *		crdsys(isw),orbtyp(isw),org(isw)
920	format(3x,i4,4i3,1x,f11.8,2x,i6,7x,a5,1x,a3,1x,a4)
	read(lusp,930) igwk,dumy,delt(isw),mjd(isw),fmjd(isw)
930	format(3x,i4,1x,f15.8,1x,f14.8,1x,i5,1x,f15.13)
	read(lusp,940) nprn(isw),(iprn(isw,i),i=1,85)
* MOD TAH 060201: Changed 17i3 to 17(1x,I2) to allow for G01 etc.
940	format(4x,i2,3x,17(1x,i2),4(/,9x,17(1x,i2)),15(/))
	spfmt(isw) = 'sp3'

	endif

c convert delt from second to day
c	delt(isw)=delt(isw)/3600/24

        return
        end



