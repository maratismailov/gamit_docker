      subroutine site_id(isnx,idat,qnum_sites,qsite_code
     .                   ,qsite_names,apr)

c Purpose: write out the list of site codes, names, domes # and coords
c          to a sinex file
c
c IN:
c      isnx         - unit number of sinex file               I*4
c      idat         - unit number of sinex.dat file           I*4
c      qnum_sites   - number of sites in hfile                I*4
c      qsite_code   - 4 char site codes                       C*4(qsite_num)
c      qsite_names  - long names of sites                     C*20
c      apr          - apriori site coords                     R*8(qnum_sites*3)
c
c P Tregoning
c 9th August, 1995

       implicit none

       integer i,j,isnx,idat,qnum_sites
       character*4 qsite_code(qnum_sites)
       character*9 domes
       character*20 qsite_names(qnum_sites)
       character*80 line
       real*8 apr(qnum_sites*3),glat(3),glong(3),ght,xyz(3),llr(3)
     .       ,height,semi,finv,rot_mat(3,3)

      parameter (semi = 6378.137d0,finv=298.257222101)

c  write out the block header
      write(isnx,100)
100   format('+SITE/ID',/
     .      ,'*CODE PT __DOMES__ T _STATION DESCRIPTION APPROX_LON_ '
     .      ,'APPROX_LAT_ _APP_H_')

c now loop through all sites, get the domes number, and output it all

      do i=1,qnum_sites

c get the domes number from the file
          call get_domes(idat,qsite_code(i),domes)

c convert the goecentric coords (in radians,m) to geodetic coords (radians and km's)
          llr(1) = apr(i*3-2)
          llr(2) = apr(i*3-1)
          llr(3) = apr(i*3)

          call llr_to_xyz(llr,xyz,rot_mat)

c convert xyz to km's
          do j=1,3
            xyz(j) = xyz(j)/1000.d0
          enddo

          call geoxyz( semi,finv,glat(1),glong(1),ght,geodrad
     .               , xyz(1),xyz(2),xyz(3),2 )
          call raddeg(glat(1),glat)
          call raddeg(glong(1),glong)

          write(line(1:31),109)int(glong(1)),int(glong(2)),glong(3)
     .       ,int(glat(1)),int(glat(2)),glat(3),ght*1000.d0
109       format(i3.3,i3.2,f5.1,1x,i3,i3.2,f5.1,f8.1)

c write out the record
          write(isnx,110)qsite_code(i),domes,qsite_names(i),line(1:31)
110       format(' ',a4,'  A ',a9,' P ',a20,1x,a31)

      enddo

c  close the block
      write(isnx,'(a)')'-SITE/ID'
      call write_dash(isnx)

      return
      end

