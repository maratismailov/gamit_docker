      subroutine rdsit_info(ih,idat,qnum_sites,qsite_code,qsite_names
     .                     ,cfiles,rcvr,antenna,aphs_mod,elev_cutoff
     .                     ,nzen,ant_off,phs_off )

c  Purpose: get all necessary site information for the sites in the hfile
c
c  This includes: -  receiver type
c                 -  recvr serial number
c                 -  antenna type
c                 -  antenna serial number
c                 -  height from ground mark to ARP (L1 and L2)
c                 -  height from ARP to L1 and L2 phase centres
c
c
c  NOTE: height ARP to phase centre is going to depend on what phase centre model
c        was used - (ref. B Schupler poster at IUGG 1995 - cutoff angle in determining
c        phase map changed both the map and the phase centre offset!)
c
c  IN:
c        ih             -  unit number of hfile                          I*4
c        idat           -  unit number of sinex.dat                      I*4
c        qnum_sites     -  number of sites in hfile                      I*4
c
c  OUT:
c        qsite_code     -  site in hfile                                 C*4(qnum_sites)
c        qsite_names    -  full name of site                             C*20(qnum_sites)
c        cfiles         -  names of input data files                     C*10(qnum_sites)
c        rcvr           -  type, s/n, firmware version                   C*16(qnum_sites,3)
c        antenna        -  type, s/n                                     C*16(qnum_sites,2)
c        aphs_mod       -  phase centre model used                       C*4(qnum_sites)
c        elev_cutoff    -  cutoff angle of solution                      R*8(qnum_sites)
c        nzen           -  number of zenith delay parameters             I*4(qnum_sites)
c        ant_off        -  antenna offset mark to ARP
c                             (L1_U, L1_N, L1_E, L2_U, L2_N, L2_E)       R*8(qnum_sites,6)
c        phs_off        -  antenna offset mark to phase centres
c                             (L1_U, L1_N, L1_E, L2_U, L2_N, L2_E)       R*8(qnum_sites,6)
c
c Paul Tregoning
c 14th July, 1995

      implicit none

      include '../includes/dimpar.h'
      integer i,j,k,ih,idat,qnum_sites,nzen(qnum_sites),ierr,num_chg
     .        ,count,unused

      character*4 qsite_code(qnum_sites),aphs_mod(qnum_sites)
     .            ,codes(maxsit,2)
      character*10 cfiles(qnum_sites)
      character*20 qsite_names(qnum_sites)
      character*16 rcvr(maxsit,3),antenna(maxsit,2)
      character*256 line

      real*8 elev_cutoff(qnum_sites),ant_off(maxsit,3)
     .      ,phs_off(maxsit,6)

      line = ' '

c blank out the site names
      do i=1,qnum_sites
         qsite_names(i) = ' '
      enddo

c read sinex.dat to get site name changes
      num_chg = 0
      do while (line(1:10).ne.'-SITE_CODE')
        read(idat,'(a)')line
        if(line(1:10).eq.'+SITE_CODE')then
          i=0
          do while (line(1:1).ne.'-')
            read(idat,'(a)')line
            if(line(1:1).ne.'-')then
              i = i+1
              read(line(1:10),'(1x,a4,1x,a4)')codes(i,1),codes(i,2)
            endif
            num_chg = i
          enddo
        endif
      enddo

c  read in all the station names and receiver, antenna type
      count = 0
      unused = 0
      do i=1,qnum_sites
         read(ih,'(a)')line
         if(line(32:33).ne.' 0')then
           count = count + 1
           read(line,100,iostat = ierr)qsite_code(count)
     .                      ,qsite_names(count),(rcvr(count,j),j=1,3)
     .                      ,(antenna(count,j),j=1,2)
100        format(6x,a4,2x,a15,9x,a16,a16,13x,a6,2x,a16,a6)
           write(*,101,iostat = ierr)qsite_names(count),rcvr(count,1)
     .                      ,antenna(count,1)
101        format(2x,a15,2x,a16,2x,a16)
         else
           unused = unused + 1
         endif
       enddo

c  reset the number of sites
      qnum_sites = count
       write(*,'(/,a,i2,a)')' Solution used ',qnum_sites,' sites.'

c  skip 2 lines
       read(ih,'(a)')line
       read(ih,'(a)')line

c  read in all the antenna offset information
       count = 1
       k=0
       do while (count+k.le.qnum_sites+unused)
           read(ih,'(a)')line
          if(line(7:10).eq.qsite_code(count))then
            read(line,110,iostat=ierr)(ant_off(count,j),j=1,3)
     .         ,(phs_off(count,j),j=1,6),aphs_mod(count)
     .         ,elev_cutoff(count),nzen(count)
110         format(30x,3(f7.5,2x),6(f7.5,2x),1x,a4,f7.2,5x,i2)
            count = count + 1
          else
            k = k + 1
          endif
       enddo

c  check each site code, and change the name if necessary
      do i=1,qnum_sites
         do j=1,num_chg
            if(qsite_code(i).eq.codes(j,1))then
              qsite_code(i) = codes(j,2)
            endif
         enddo
      enddo

      return
      end
