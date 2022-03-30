      subroutine get_estimates(ih,qnum_sites,qnum_sats,apr,adj,num_orb
     .                       ,qnum_par,qparn_sites,qparn_svs,qparn_eop
     .                       ,qnum_eop,gps_utc,fix_eop)

c  Purpose: read from the hfile the apriori values and the adjustments to the
c           site/sat parameters. Also return the number of orbital parameters
c           estimated for each satellite.
c
c IN:   ih          -  hfile unit number                         I*4
c       qnum_sites  - number of sites in hfile                   I*4
c       qnum_sats   - number of satellites in hfile              I*4
c       gps_utc     - gpst - utc                                 R*8
c OUT:
c       apr         - apriori values for all parameters          R*8(lots!)
c       adj         - adjustments for all parameters             R*8(lots!)
c       num_orb     - number of orbital elements estimates/sat   I*4
c       qnum_par    - total number of parameters in hfile        I*4
c       fix_eop     - apriori values of any fixed eop parameters R*8(6)
c
c NOTE: Maximum number of parameters = number of sites*3 + num sats * num elements
c                                      + 6 eop parameters

      implicit none

      include '../includes/dimpar.h'

      integer ih,qnum_sites,qnum_sats,num_orb,i,j,par_num,qnum_par
     .       ,count
     .      ,qparn_sites(3,qnum_sites),qparn_svs(maxorb,qnum_sats)
     .      ,qparn_eop(2,3),qnum_eop
      real*8 apr(qnum_sites*3 + qnum_sats*maxorb+6),gps_utc
      real*8 adj(qnum_sites*3 + qnum_sats*maxorb+6),fix_eop(6)
      character*4 prn
      character*80 line

      count = 0

c  skip 3 lines
      read(ih,'(a)')line
      read(ih,'(a)')line
      read(ih,'(a)')line

c  read in all the site coords
      i=0
      do while (line(7:11).ne.'ORBIT')
        read(ih,'(a)')line
        if(line(6:6).eq.'*'.and.line(7:11).ne.'ORBIT')then
          i=i+1
          count = count + 1
          read(line,'(27x,f24.16,d27.20)')apr(count),adj(count)
            qparn_sites(1,i) = count
          do j=2,3
            count = count + 1
            read(ih,'(27x,f24.16,d27.20)')apr(count),adj(count)
            qparn_sites(j,i) = count
          enddo
        endif
      enddo

c  rewind to the first line of the orbital elements
      backspace(ih)

c  now read in the first line of satellite orbital parameters
      i=1
      read(ih,100)prn,apr(count+i),adj(count+i)
100   format(22x,a4,f24.16,d27.20)
      line(23:26) = prn
      count = count + 1
      qparn_svs(i,1) = count

c  keep reading until the satellite number changes
      do while (line(23:26).eq.prn)
        read(ih,'(a)')line

        if(line(23:26).eq.prn)then
          i=i+1
          count = count + 1
          read(line(27:50),'(f24.16)')apr(count)
          read(line(51:77),'(f27.20)')adj(count)
          qparn_svs(i,1) = count
        endif
      enddo

c rewind to the first line of the 2nd satellite
      backspace(ih)
      num_orb = i

c now read the remaining satellites
      do i =2,qnum_sats
        do j=1,num_orb
          count = count + 1
          read(ih,110)par_num,apr(count),adj(count)
110       format(i5,22x,f24.16,d27.20)
          qparn_svs(j,i) = count
        enddo
      enddo

c now read in the eop parameters
c PT950815: There will be 6 parameters - but not all may have been estimated.

      qnum_eop = 0
      do i=1,3
        do j=1,2
          read(ih,'(a)')line
          read(line(1:5),'(i5)')par_num
          if(par_num.ne.0)then

c we have an eop parameter - check if it is estimated or not
            if(line(6:6).eq.'*')then
              count = count + 1
              qnum_eop = qnum_eop + 1
              read(line,110)par_num,apr(count),adj(count)
              qparn_eop(j,i) = count

c if it is UT1-TAI, convert it to UT1-UTC
              if(line(7:13).eq.'UT1-TAI')then
                apr(count) = apr(count) - gps_utc - 19.d0
              endif
            else

c store the value as a fixed value, to be written out later in the soln/apriori block
              read(line,'(i5,22x,f24.16)')par_num,fix_eop(2*i-2+j)
              qparn_eop(j,i) = 0

c if it is UT1-TAI, convert it to UT1-UTC
              if(line(7:13).eq.'UT1-TAI')then
                fix_eop(2*i-2+j) = fix_eop(2*i-2+j) - gps_utc - 19.d0
              endif
            endif
          endif
        enddo
      enddo

c  now determine the total number of parameters in hfile
      qnum_par = count

c read one more blank line
      read(ih,'(a)')line

      return
      end
