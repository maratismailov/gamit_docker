	subroutine utc2gps(idyoyr,iy,ihr,min,sec)

c	purpose:	convert dayofyear, hr, min, sec
c              from UTC to GPST
c                        after adding hour,min.sec
c	input:		idyoyr, iy, hours, minute, sec
c	output:		idyoyr, iy
c
c       R. King 13 July 94 from Bock/Fang FIXDRV version


	integer*4 idyoyr,iy,ihr,min,imo,iday,jd,julday
	real*8 sec,utcoff,taiutc,tsec

c         print *,'Debug utc2gps',idyoyr,imo,iday,iy
         call monday(idyoyr,imo,iday,iy)
c         print *,idyoyr,imo,iday,iy
         jd=julday(imo,iday,iy)
         utcoff = taiutc(jd) - 19.d0
         tsec = dble(ihr)*3600.d0 + dble(min)*60.d0 + sec
C        add to jd
c         print *,'jd,tsec,utcoff ',jd,tsec,utcoff
         call timinc(jd,tsec,utcoff)
C        convert back to year and day of year
         call dayjul(jd,iy,idyoyr)
         call ds2hms(iy,idyoyr,tsec,ihr,min,sec)
c         print *,idyoyr,iy,ihr,min,sec


	return
	end
