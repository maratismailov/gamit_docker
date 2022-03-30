      program svpdiff

      implicit none 
 
*     Program to difference position estimates from svpos program.
 
*   err1, err2      - Errros from two files
*   rcpar           - Read runstring
*   len_run         - Length of string
*   trimlen         - Length of string
*   num             - NUmber of output values
*   date(5)         - Date read from file (ymdhm)
*   i,j             - Looop counters
*   numsv1, numsv2  - Numbers of satellites used in file 1 and 2
 
 
      integer*4 err1, err2, rcpar, len_run, trimlen, num, date(5),
     .    i, numsv1, numsv2
 
*   sectag          - Seconds tag read from file
*   ep1, ep2        - JD with fractional day from the two files
*   pos1(4), pos2(4)    - Three adjustments to the positions plus
*                       - Clock offset
*   sig1(4), sig2(4)    - Four sigmas for the pos1/2 values
 
 
      real*8 sectag, ep1, ep2, pos1(4), pos2(4), sig1(4), sig2(4)
 
*   need1, need2        - Indicates next read dhould vbe of file i
 
 
      logical need1, need2
 
*   in1, in2, out           - File names
 
 
      character*128 in1, in2, out
 
*    line    - Line read from file
 
 
      character*256 line
 
****  Read th runstring
      len_run = rcpar(1,in1)
      if( len_run.le.0 ) then
          write(*,*) 'svpdiff: runstring svpdiff <in1> <in2> <out>'
          stop 'Element 1 missing runstring'
      end if
      len_run = rcpar(2,in2)
      if( len_run.le.0 ) then
          write(*,*) 'svpdiff: runstring svpdiff <in1> <in2> <out>'
          stop 'Element 2 missing runstring'
      end if
      len_run = rcpar(3,out)
      if( len_run.le.0 ) then
          write(*,*) 'svpdiff: runstring svpdiff <in1> <in2> <out>'
          stop 'Element 3 missing runstring'
      end if
 
*     OPen the files
      open(100, file=in1, status='old', iostat=err1)
      call report_error('IOSTAT',err1,'open',in1,1,'svpdiff')
      open(101, file=in2, status='old', iostat=err2)
      call report_error('IOSTAT',err2,'open',in2,1,'svpdiff')
      open(200, file=out,  iostat=err2)
      call report_error('IOSTAT',err2,'open',out,1,'svpdiff')
      write(200,120) in1(1:trimlen(in1)), in2(1:trimlen(in2))
 120  format('DD file from ',a,'-',a,/,
     .    '   Date     dN   +-  dE   +-   dU +-  Num SV')
 
***   Start looping over the files
      need1 = .true.
      need2 = .true.
      num = 0
      do while (err1.eq.0 .and.err2.eq.0 )
          if( need1) then
              read(100,'(a)', iostat=err1) line
              if( trimlen(line).gt.0 .and. line(1:1).eq.' ') then
                  read(line, *,iostat=err1) date, sectag,
     .               (pos1(i), sig1(i),i = 1,4), numsv1
                  call ymdhms_to_jd( date, sectag, ep1)
              else
                  ep1 = -1
              end if
          end if
          if( need2) then
              read(101,'(a)', iostat=err1) line
              if( trimlen(line).gt.0 .and. line(1:1).eq.' ') then
                  read(line, * ,iostat=err2)  date, sectag,
     .               (pos2(i), sig2(i),i = 1,4), numsv2
                  call ymdhms_to_jd( date, sectag, ep2)
              else
                  ep2 = -1
              end if
          end if
          if( ep1.eq.ep2 .and. ep1.ne.-1) then
              need1 = .true.
              need2 = .true.
              do i = 1, 4
                 pos1(i) = pos1(i) - pos2(i)
                 sig1(i) = sqrt(sig1(i)**2 + sig2(i)**2)
              end do
              if( err1.eq.0 .and. err2.eq.0 ) then
                  num = num + 1
                  call jd_to_ymdhms( ep1+1.d-9, date, sectag)
                  if( numsv1.eq.numsv2 )
     .            write(200,300) date, sectag,
     .                   (pos1(i), sig1(i), i = 1,4), numsv1, numsv2
 300              format(I5,4i3,f5.1,3(F10.2,1x,f8.2),1x,F15.2,1x,f8.2,
     .                   2i3)
              end if
          else if( ep1.eq.-1 .and. ep2.eq.-1 ) then
              need1 = .true.
              need2 = .true.
          else if( ep1.gt.ep2 .or. ep2.eq.-1 ) then
              need1 = .false.
              need2 = .true.
          else
              need2 = .false.
              need1 = .true.
          end if
      end do
 
****  Close the files
      write(*,400) num
 400  format(I6,' Position differences written to output')
      close(100)
      close(101)
      close(200)
      end
 
 
 
 
