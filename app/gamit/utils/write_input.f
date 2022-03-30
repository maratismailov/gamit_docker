      subroutine write_input(isnx,in_h,times,num_par,sol_type,agency)

c Purpose: To write information of the solution input files on the sinex file
c
c  P Tregoning
c  9th August, 1995
c
c  IN:
c       isnx         -  unit number of sinex file                          I*4
c       in_h         -  input hfile name                                   C*10
c       times        -  runtime,start,stop,ic,eop,mean times (snx form)    C*12(6)
c                       where mean is the mean time of obs used in soln
c       num_par      -  number of parameters in hfile soln                 I*4
c       sol_type     -  flags for type of soln: col1 constr/free
c                                               2: coord(X)/vel(V) soln
c                                               3: eop estimated (E)
c                                               4: biases real(r), integer(i) C*4
c       agency       -  name of agency running sinex program                C*50

      implicit none

      integer isnx,num_par,i
      character*4 sol_type
      character*10 in_h
      character*12 times(6)
      character*50 agency

c  write the input/history block
      write(isnx,100)
100   format('+INPUT/HISTORY',/
     .       ,'*O_FM VER_ AGY TIME_STAMP__ DAT '
     .       ,'DATA_START__ DATA_END____ T PARAM S TYPE_')
      write(isnx,110)agency(1:3),times(1),times(2),times(3),num_par
     .               ,(sol_type(i:i),i=1,3)
110   format(' =SNX 0.05 ',a3,' ',a12,' SIO ',a12,' ',a12,' P ',i5.5
     .       ,3(' ',a1))
      write(isnx,120)
120   format('-INPUT/HISTORY')
      call write_dash(isnx)

c  write the input/files block
      write(isnx,130)
130   format('+INPUT/FILES',/
     .      ,'*_Input_file_ID__ _filename__________ __________Descripti'
     .      ,'on_____________________')
      write(isnx,135)times(1),in_h,times(2)(1:2),times(2)(4:6)
135   format(' SIO ',a12,' ',a10,10x,'Daily Solution ( ',a2,'/',a3,' )')
      write(isnx,'(a)')'-INPUT/FILES'
      call write_dash(isnx)


      return
      end
