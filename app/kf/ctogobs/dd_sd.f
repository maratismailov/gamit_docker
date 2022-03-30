      program dd_sd

      implicit none
 
*     Program to double difference SD files
 
*   ep1, ep2        - Current epochs form towo files
*   df11,df12, df21,df22    - Data flags from the SD files
*   dfc             - combined data flag
*   err1, err2      - Errros from two files
*   rcpar           - Read runstring
*   len_run         - Length of string
*   trimlen         - Length of string
*   num             - NUmber of output values
*   cand            - Common name for iand/and function
 
      integer*4 ep1, ep2, df11,df12, df21,df22, dfc, err1, err2, rcpar,
     .    len_run, trimlen, num, cand
 
*   el1, el2        - Elevations
*   az1, az2        - Azimuths
*   L11, L12, dL1       - L1 Phase residuals
*   L21, L22, dL2       - L2 Phase residuals
*   P11, P12, dP1       - P1 Phase residuaPs
*   P21, P22, dP2       - P2 Phase residuaPs
 
      real*8 el1, el2, az1, az2, L11, L12, dL1, L21, L22, dL2,
     .    P11, P12, dP1, P21, P22, dP2
 
*   need1, need2        - Indicates next read dhould vbe of file i
 
      logical need1, need2
 
*   in1, in2, out           - File names
 
      character*128 in1, in2, out, runstring
 
****  Read th runstring
      len_run = rcpar(1,in1)
      if( len_run.le.0 ) then
          write(*,*) 'DD_SD: runstring dd_sd <in1> <in2> <out>'
          stop 'Element 1 missing runstring'
      end if
      len_run = rcpar(2,in2)
      if( len_run.le.0 ) then
          write(*,*) 'DD_SD: runstring dd_sd <in1> <in2> <out>'
          stop 'Element 2 missing runstring'
      end if
      len_run = rcpar(3,out)
      if( len_run.le.0 ) then
          write(*,*) 'DD_SD: runstring dd_sd <in1> <in2> <out>'
          stop 'Element 3 missing runstring'
      end if
 
*     OPen the files
      open(100, file=in1, status='old', iostat=err1)
      call report_error('IOSTAT',err1,'open',in1,1,'DD_SD')
      read(100,'(a)') runstring
      read(100,'(a)') runstring
      open(101, file=in2, status='old', iostat=err2)
      call report_error('IOSTAT',err2,'open',in2,1,'DD_SD')
      read(101,'(a)') runstring
      read(101,'(a)') runstring
      open(200, file=out,  iostat=err2)
      call report_error('IOSTAT',err2,'open',out,1,'DD_SD')
      write(200,120) in1(1:trimlen(in1)), in2(1:trimlen(in2))
 120  format('DD file from ',a,'-',a,/,
     .    '  Epoch  dL1 dL2 dP1 dP2 el1 el2 az1 az2 DF')
 
***   Start looping over the files
      need1 = .true.
      need2 = .true.
      num = 0
      do while (err1.eq.0 .and.err2.eq.0 )
          if( need1) read(100, 200,iostat=err1) el1,az1, ep1,
     .                l11, l21, p11, p21, df11,df21
          if( need2) read(101, 200,iostat=err2) el2,az2, ep2,
     .                l12, l22, p12, p22, df12,df22
200       format(2f9.4,1x,I5,1x,2(F12.4,1x),2(F10.2,1x),2(1x,o12))
          if( ep1.eq.ep2 ) then
              need1 = .true.
              need2 = .true.
              dl1 = l11 - l12
              dl2 = l21 - l22
              dp1 = p11 - p12
              dp2 = p21 - p22
              dfc = cand(cand(df11,df12),cand(df21,df22))
              if( err1.eq.0 .and. err2.eq.0 ) then
                  num = num + 1
                  write(200,300) ep1, dl1, dl2, dp1, dp2,
     .                        el1, el2, az1, az2, dfc
 300          format(1x,I5,1x,2(F12.4,1x),2(F10.2,1x),4f9.4,1x,o12)
              end if
          else if( ep1.gt.ep2 ) then
              need1 = .false.
              need2 = .true.
          else
              need2 = .false.
              need1 = .true.
          end if
      end do
 
****  Close the files
      write(*,400) num
 400  format(I6,' Doubles differences written to output')
      close(100)
      close(101)
      close(200)
      end
 
 
 
 
