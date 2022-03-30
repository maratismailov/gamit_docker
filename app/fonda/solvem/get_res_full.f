      subroutine get_res_full(ifile,idatf,iobs1,it,prerms,pstrms,
     .                        unit,wrms,iefob)
c
c     calculate residuals with full covariance data
c
c     data type: (it)
c         1. astrometric azimuth
c         2. horizontal angle
c         3. horizontal direction
c         4. baseline length
c         5. zenith height
c         6. leveling
c
c        11. astrometric azimuth rate
c        12. horizontal angle rate
c        13. horizontal direction rate
c        14. baseline length rate
c        15. zenith height rate
c        16. leveling rate
c
c        21. 3-D geocentric coordinate
c        22. 3-D geocentric velocity
c        23. 3-D geodetic coordinate
c        24. 3-D geodetic velocity
c        25. 3-D geocentric baseline vector
c        26. 3-D geocentric baseline rate vector
c        27. 3-D geodetic baseline vector
c        28. 3-D geodetic baseline rate vector
c
c        31. 3-D spherical coordinate
c        32. 3-D spherical frame velocity
c        33. 3-D Cartesian coordinate
c        34. 3-D Cartesian frame velocity
c
c     unit:
c         x, y, z : m
c        vx,vy,vz : mm/year
c        err      : mm
c        time     : year
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'

      character*6 unit(3)
      character*128 lines(600),line
      character*8  name1
      integer ifile,idatf,i,lift_arg,nblen,match_name
      integer iobs1,it,iefob,ib,ibd,ie1,k,isit1,is1
      integer il,j,j1,nlen
      real*8  err,obs,omc,oc,cal
      dimension  err(3),obs(3),omc(3),oc(3),cal(3)
      dimension covo(6),temp(6)
      character*1 symb
      logical ill
c
c        full observations are always it > 20
         if (it.gt.20.and.it.le.40) k = iobs1*3
         if (it.eq.32.or.it.eq.34) then
            unit(1) = '(mm/y)'
            unit(2) = '(mm/y)'
            unit(3) = '(mm/y)'
         else
            unit(1) = '(d,m )'
            unit(2) = '(d,m )'
            unit(3) = '(s,mm)'
         endif
         write (ifile,70) (unit(k),k=1,3), unit(3),unit(3)
c
         iefob = 0
c        get data part first
         do 20 ib = 1,iobs1
            read (idatf,'(a128)') lines(ib)
 20      continue
c
c        get covariance submatrix
         small = 1.0d-20
         do 30 ib = 1,iobs1
            do 40 j = 1,3
               j1 = j*(j+1)/2
               if (it.eq.31) read(idatf,110) ibd,(scale(i),i=1,ibd)
               if (it.eq.33) read(idatf,110) ibd,(scale(i),i=1,ibd)
               if (it.eq.34) read(idatf,110) ibd,(scale(i),i=1,ibd)
c              if (it.eq.34) read(idatf,130) ibd,(scale(i),i=1,ibd)
               covo(j1) = scale(ibd)
               if (scale(ibd).gt.small) then
                  err(j) = dsqrt(scale(ibd))
               else
                  err(j) = dsqrt(small)
               endif
               if (j.gt.1) covo(j1-1) = scale(ibd-1)
               if (j.eq.3) covo(j1-2) = scale(ibd-2)
 40         continue
            line = lines(ib)
            read (line,*) time1,obs(1),obs(2),obs(3)
            k = lift_arg(line,name1,5)
            if (k.le.0) goto 30
            nlen = nblen(name1)
            isit1 = match_name(nsit,nlen,sname,name1)
            if (isit1 .le. 0) goto 30
            dt = time1-rtime
            call calc_omc_full(it,isit1,dt,obs,cal,omc,oc,ill)
            is1 = jtoi(isit1)
            iefob = iefob+3
c
c           get weighting matrix
            call nrmscl(covo,obs,scale,3,1,0)
            call cholsk(covo,temp,1,3,ie1)
            call nrmscl(covo,obs,scale,3,2,0)
            call atwa(1,3,omc,covo,wr,temp,1)
            wrms = wrms+wr
            do il = 1,3
               if (ill) then
                  symb = '*'
                  iefob = iefob-1
               else 
                  prerms = prerms+oc(il)**2
                  pstrms = pstrms+omc(il)**2
                  symb = ' '
               endif
c
c              output residuals
               fac = 1.0d0
               if (it.eq.31.and.il.le.2) fac = rtod*3.6d3
               if (it.eq.31.and.il.eq.3) fac = 1.0d3
               if (it.eq.34.or.it.eq.33) fac = 1.0d3
               write (ifile,80) symb,sname(isit1),'        ',
     .            time1,obs(il),cal(il),err(il)*fac,omc(il),oc(il),
     .            (omc(il)/(err(il)*fac))
            enddo
 30      continue     
c
 70   format ('*   Site1     Site2      Time   ',
     .     4x,' Obs.',a6,6x,'Calc.',a6,3x,'Sigma',
     .     a6,4x,'Adj.',a6,6x,'O-C ',a6,2x,'v/s') 
 80   format (a1,1x,a8,2x,a8,1x,f10.3,2(1x,f15.5,1x),f12.5,2f15.5,f9.1)

 50   continue
 90   continue
c
 110  format (1x,i3,'. ',50(5(1x,d23.16),:,/,6x))
 130  format(i5,'. ', 10(1x,d11.4),:/,200(7x,10(1x,d11.4),:/) )
c
      return
      end
