      subroutine wrtant_info(isnx,qnum_sites,antenna,ant_off,phs_off
     .                      ,aphs_mod,qsite_code,start_time,end_time)

c Purpose: Write out the blocks giving antenna offsets from ARP to
c          phase centre and ground mark to ARP
c
c IN:  isnx       - sinex file unit #                                 I*4
c      qnum_sites - # sites in hfile                                  I*4
c      antenna    - type and s/#                                      C*16(maxsit,2)
c      ant_off    - antenna offset - mark to ARP (UNE)                R*8(maxsit,3)
c      phs_off    - antenna offset - mark to phase centres (UNE,UNE)  R*8(maxsit,6)
c      aphs_mod   - antenna phase centre model                        C*4(qnum_sites)
c      qsite_code - 4 char site name                                  C*4(qnum_sites)
c      start_time - start time of observations                        C*12
c      end_time   - stop time of observations                         C*12
c
c  P Tregoning
c  10th August, 1995

      implicit none

      include '../includes/dimpar.h'

      integer isnx,qnum_sites,i,j,k,ant_def_num
      real*8 ant_off(maxsit,3),phs_off(maxsit,6),dh(6)
      character*4 aphs_mod(qnum_sites),qsite_code(qnum_sites)
      character*12 start_time,end_time
      character*16 antenna(maxsit,2),ant_out(maxsit)
      logical match

      ant_def_num = 0

c  oepn the SITE/GPS_PHASE_CENTER block
      write(isnx,100)
100   format('+SITE/GPS_PHASE_CENTER',/
     .     ,'*                       _up___ _north _east_ _up___ _north'
     .     ,' _east_',/
     .     ,'*__antenna type__ _S/N_ __L1-ARP(m)_________ __L2-ARP(m)__'
     .     ,'_______ Az_El_Model')


c  loop over all the sites
      do i=1,qnum_sites
          do j=1,3
            dh(j) = phs_off(i,j) - ant_off(i,j)
            dh(j+3) = phs_off(i,j+3) - ant_off(i,j)
          enddo

c PT950909: only write out this information if either there is a serial
c           number OR it is the first occurance of this type of antenna.
c           This will then conform to sinex standards

          if(antenna(i,2)(1:5).ne.'-----')then
            write(isnx,110)(antenna(i,j),j=1,2),(dh(j),j=1,6)
     .                     ,aphs_mod(i)
110         format(1x,a16,' ',a5,6(1x,f6.4),' ',a4)
          else

c now check whether this type of antenna has already been written out
            match = .false.
            if(ant_def_num.eq.0)then
              ant_def_num = ant_def_num + 1
              ant_out(ant_def_num) = antenna(i,1)
              write(isnx,110)(antenna(i,j),j=1,2),(dh(j),j=1,6)
     .                       ,aphs_mod(i)
            endif

            do k=1,ant_def_num
              if(antenna(i,1).eq.ant_out(k))match = .true.
            enddo

            if(.not.match)then
              ant_def_num = ant_def_num + 1
              ant_out(ant_def_num) = antenna(i,1)
              write(isnx,110)(antenna(i,j),j=1,2),(dh(j),j=1,6)
     .                       ,aphs_mod(i)
            endif
          endif

      enddo

c  close the block
      write(isnx,'(a)')'-SITE/GPS_PHASE_CENTER'
      call write_dash(isnx)


c  oepn the SITE/ECCENTRICITY block
      write(isnx,120)
120   format('+SITE/ECCENTRICITY',/
     .     ,'*                                             _up_____'
     .     ,' _north__ _east___',/
     .     ,'*Code PT _Occ T _Data Start_ __Data End__ typ __ARP-benchm'
     .     ,'ark (m)_______')


c  loop over all the sites
      do i=1,qnum_sites

          write(isnx,130)qsite_code(i),start_time,end_time
     .                   ,(ant_off(i,j),j=1,3)
130       format(' ',a4,'  A    1 P ',a12,' ',a12,' UNE'
     .           ,3(1x,f8.4))
      enddo

c  close the block
      write(isnx,'(a)')'-SITE/ECCENTRICITY'
      call write_dash(isnx)


      return
      end


