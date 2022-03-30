       subroutine jul2hms(time,yr,mo,day,hr,min,sec)

c
c  function to convert time in PEP JD to year, month, day,
c  hours, minutes, seconds. Year will be full (ie 1995 not 95)
c
c  IN:      time:  time in PEP Julian days   R*8
c                  (=true Julian date + 0.5)
c
c OUT:        as the names imply!            I*4
c           sec:                             R*8
c
c LOCAL:    doy: day of year                 I*4
c           sod: seconds of day              R*8
c            jd: julian day number           I*4
c
c Paul Tregoning 950324
c
       implicit none

       real*8 time,sec,sod
       integer yr,mo,day,hr,min,doy,jd

       jd = idint(time)
       call dayjul(jd,yr,doy)
       call monday(doy,mo,day,yr)
       sod = (time - dble(jd))*86400.d0
       call ds2hms(yr,doy,sod,hr,min,sec)

       return
       end
