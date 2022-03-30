      subroutine get_domes(idat,site,domes)

c Purpose: to retrieve the domes number of 'site' from the file domes.dat
c          If there is no entry for 'site' in sinex.dat then it returns
c          '---------' instead.
c
c IN:
c      idat   - unit number of sinex.dat         I*4
c      site   - 4 char site code                 C*4
c OUT:
c      domes  - domes number of site             C*9
c
c P Tregoning
c August, 1995

      implicit none

      integer idat,ierr
      character*4 site
      character*9 domes
      character*80 line

      ierr = 0

c  find the start of the domes block
      do while (line(1:6).ne.'+DOMES')
        read(idat,'(a)')line
      enddo

c  now read file until we get a site match or the end of file
      do while (line(1:1).ne.'-')
        read(idat,'(a)')line
        if(line(1:1).ne.'-')then
          if(line(2:5).eq.site)then
             domes = line(8:17)
             rewind(idat)
             return
          endif
        else
          domes = '---------'
        endif
      enddo

      rewind(idat)
      return
      end
