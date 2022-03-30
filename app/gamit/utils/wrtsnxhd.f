      subroutine wrtsnxhd(isnx,in_h,times,num_par,idat,solv_vers
     .                    ,sol_type,agency)

c Purpose: To write all the ancilliary header records of the sinex file
c
c  P Tregoning
c  9th August, 1995
c
c  IN:
c       isnx         -  unit number of sinex file                          I*4
c       in_h         -  input hfile name                                   C*10
c       times        -  runtime,start,stop,ic,eop,mean times (snx form)    C*12(6)
c                       where mean = mean of obs used
c       num_par      -  number of parameters in solution                   I*4
c       idat         -  unit # of sinex.dat                                I*4
c       solv_vers    -  solve version from hfile header                    C*48
c       sol_type     -  soln type indictors col1: constraints (0,1,2)
c                                           col2: coord (X) or vel(V) soln
c                                           col3: EOP estimated? (E or blank)
c                                           col4: biases (r or i)          C*4
c
c OUT:
c       agency       -  name of agency creating the sinex file             C*50

      implicit none

      integer isnx,num_par,idat,yr,mo,day,hr,min,hrdwr
      character*4 sol_type
      character*10 in_h
      character*12 times(6)
      character*19 dattim
      character*48 solv_vers
      character*50 agency
      character*80 line

      rewind(idat)
c  get the runtime
      call runtim(dattim)
      read(dattim,90)day,mo,yr,hr,min
90    format(i2,1x,i2,1x,i4,2(1x,i2))
      yr = yr - 1900
      call conv_date(yr,mo,day,hr,min,0.d0,times(1))

c read the processing agency from the sinex.dat file
      read(idat,'(20x,a50)')agency


c write the first line of a sinex file
      write(isnx,100)agency(1:3),times(1),agency(1:3),times(2),times(3)
     .              ,num_par,sol_type(1:1),sol_type(2:2),sol_type(3:3)
100   format('%=SNX 0.05 ',a3,' ',a12,' ',a3,' ',a12,' ',a12,' P ',i5.5
     .              ,3(' ',a1))

c write next four lines (a bunch of useless lines!)
      call write_dash(isnx)
      write(isnx,110)
110   format('*        1         2        '
     .       ,' 3         4         5         6         7         8',/
     .       ,'*2345678901234567890123456789012345678901234567890123456'
     .       ,'789012345678901234567890')
      call write_dash(isnx)

c now start the file/reference block of information
c first read lines from sinex.dat
      do while (line(1:1).ne.'-')
        read(idat,'(a80)')line
        if(line(1:1).ne.'-') write(isnx,'(a80)')line
      enddo

c find the hardware type in the hfile solve version line
      hrdwr = index(solv_vers,"(")
c now write out the software version, OUTPUT solution line and close block
      write(isnx,120)agency(1:3),in_h,solv_vers(hrdwr:(hrdwr+6))
     .               ,solv_vers(1:(hrdwr-1))
120   format(' OUTPUT',13x,a3,' hfile solution from ',a10,/
     .       ,' HARDWARE',11x,a,/
     .       ,' SOFTWARE',10x,a,/,'-FILE/REFERENCE')

      call write_dash(isnx)

c  write a file comment
      write(isnx,130)
130   format('+FILE/COMMENT',/
     .      ,' Sinex file created from a single hfile gamit solution')

c PT950901: This is unneccesary in a sinex file
c      if(sol_type(1:1).eq.'2')write(isnx,135)
c135   format(' Loose constraints: Orbits - X,Y,Z             10 m',/
c     .      ,'                             vx,vy,vz          10 m/d',/
c     .      ,'                             rad1,rad2,rad3   100 %',/
c     .      ,'                    Sites  - X,Y,Z            100 m',/
c     .      ,'                    EOP    - Xp,yp           3000 mas',/
c     .      ,'                             dXp,dYp          300 mas/d',/
c     .      ,'                             UT1               20 mts',/
c     .     ,'                             LOD                2 mts/d',/)

c  now write the type of solution (ie bias fixed, constrained etc)

      if(sol_type(1:1).eq.'0')then
        write(line(1:24),'(a)')' Constrained solution - '
      else
        write(line(1:24),'(a)')' Loose solution       - '
      endif

      if(sol_type(4:4).eq.'r')then
        write(line(25:48),'(a)')'Biases as real values   '
      else
        write(line(25:48),'(a)')'Biases as integer values'
      endif

      write(isnx,'(a48,/,a)')line(1:48),'-FILE/COMMENT'
      call write_dash(isnx)

c  write Input Acknowledgement block
      write(isnx,140)agency
140   format('+INPUT/ACKNOWLEDGMENTS',/,' ',a50,/
     .      ,'-INPUT/ACKNOWLEDGMENTS')
      call write_dash(isnx)

      return
      end














