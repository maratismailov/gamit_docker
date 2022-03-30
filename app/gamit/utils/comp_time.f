      subroutine comp_time(epo_s,epo_f,obs_int,times)

c  Purpose: Compute the start,end and mean times of the observations and
c           return them in times(2),times(3),times(6)
c
c IN:
c     epo_s, epo_f          - epochs used in solve solution                   I*4
c     obs_int               - observation interval (in seconds)               I*4
c
c OUT:
c     times                 - start,end and mean times of the observations    C*12(6)
c
c P Tregoning
c 11th August, 1995

      implicit none

      integer epo_s,epo_f,obs_int,leap
      real*8 yr_x,doy_x,sod_x,sod_s,doy_s,yr_s,sod_f,doy_f,yr_f
     .      ,sod_m,doy_m,yr_m
      character*12 times(6)

c determine the start time of the xfile
      read(times(2),'(f2.0,1x,f3.0,1x,f5.0)')yr_x,doy_x,sod_x
      yr_s = yr_x
      doy_s = doy_x
      sod_s = sod_x

c  check if the observations start in a leap year
      if(mod(yr_s+1900.d0,4.d0).eq.0.d0)then
        leap = 1.d0
      else
        leap = 0.d0
      endif

c  determine the start time of the obs
      if(epo_s.ne.1)then
        sod_s = sod_x + (epo_s - 1)*obs_int
        if(sod_s.ge.86400.d0)then
          sod_s = sod_s - 86400.d0
          doy_s = doy_s + 1.d0
        endif
        if(doy_s.gt.366.d0+leap)then
          doy_s = doy_s - 365.d0 - leap
          yr_s = yr_s + 1.d0
        endif
      write(times(2),100)int(yr_s),int(doy_s),int(sod_s)
100   format(i2,':',i3.3,':',i5.5)
      endif

c  determine the end time of the obs
      sod_f = sod_x + (epo_f - 1)*obs_int
      doy_f = doy_x
      yr_f = yr_x
       do while (sod_f.ge.86400.d0)
        sod_f = sod_f - 86400.d0
        doy_f = doy_f + 1.d0
      enddo
      if(doy_f.gt.366.d0+leap)then
        doy_f = doy_f - 365.d0 - leap
        yr_f = yr_f + 1.d0
      endif
      write(times(3),100)int(yr_f),int(doy_f),int(sod_f)

c  now compute the mean time
      yr_m = (yr_s + yr_f)/2.d0
      doy_m = (doy_s + doy_f)/2.d0 + (365.d0+leap)*mod(yr_m,1.d0)
      sod_m = (sod_s + sod_f)/2.d0 + 86400.d0*mod(doy_m,1.d0)

      if(sod_m.ge.86400.d0)then
        sod_m = sod_m - 86400.d0
        doy_m = doy_m + 1.d0
      endif

      if(doy_m.gt.366.d0+leap)then
        doy_m = doy_m - 365.d0 - leap
        yr_m = yr_m + 1.d0
      endif

      write(times(6),100)int(yr_m),int(doy_m),int(sod_m)

      return
      end
