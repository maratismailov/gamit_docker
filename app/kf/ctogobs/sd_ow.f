      program sd_ow

      implicit none
 
*     Program to create Single difference files from one-way DPH 
*     files from autcln.

      integer*4 mbf
      parameter ( mbf = 100 )
  
 
*   ep1, ep2        - Current epochs form towo files
*   df11,df12, df21,df22    - Data flags from the SD files
*   dfc             - combined data flag
*   err1, err2      - Errros from two files
*   rcpar           - Read runstring
*   len_run         - Length of string
*   trimlen         - Length of string
*   num             - NUmber of output values
*   cand            - Common name for iand/and function
 
      integer*4 ep1, ep2, df1,df2, dfc, err1, err2, rcpar,
     .    len_run, trimlen, num, cand, je1, je2, jec, sv1, sv2,
     .    bf_mask, ng1, ng2, ngc
 
*   el1, el2        - Elevations
*   az1, az2        - Azimuths
*   L11, L12, dL1       - L1 Phase residuals
*   L21, L22, dL2       - L2 Phase residuals
*   P11, P12, dP1       - P1 Phase residuaPs
*   P21, P22, dP2       - P2 Phase residuaPs
 
      real*8 el1, el2, az1, az2, l11, l12, dl1,
     .       l21, l22, dl2, lc1,lc2, dlc, lg1, lg2, dlg,
     .       rc1, rc2, drc, wl1, wl2, dwl, nl1, nl2, dnl

      real*8 sion, vion, mion, rmsi, mean_ion(mbf)

      integer*4 nbf
 
*   need1, need2        - Indicates next read dhould vbe of file i
 
      logical need1, need2, kbit
 
*   in1, in2, out           - File names
 
      character*128 in1, in2, out, runstring
      character*4 zion
 
****  Read th runstring
      len_run = rcpar(1,in1)
      if( len_run.le.0 ) then
          write(*,*) 'SD_OW: runstring sd_ow <in1> <in2> <out> [ZION]'
          stop 'Element 1 missing runstring'
      end if
      len_run = rcpar(2,in2)
      if( len_run.le.0 ) then
          write(*,*) 'SD_OW: runstring sd_ow <in1> <in2> <out> [ZION]'
          stop 'Element 2 missing runstring'
      end if
      len_run = rcpar(3,out)
      if( len_run.le.0 ) then
          write(*,*) 'SD_OW: runstring sd_ow <in1> <in2> <out> [ZION]'
          stop 'Element 3 missing runstring'
      end if

*     See if zero ion delay option passed
      zion = 'NO'
      len_run = rcpar(4,zion)
      call casefold(zion)
      
 
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
      if( zion(1:1).eq.'Z' ) then
         write(200,120) in1(1:trimlen(in1)), in2(1:trimlen(in2))
 120     format('* SD file from ',a,'-',a,' Ion average removed'/,
     .    '* Epoch  dL1 dL2 dP1 dP2 el1 el2 az1 az2 DF')
      else
         write(200,125) in1(1:trimlen(in1)), in2(1:trimlen(in2))
 125     format('* SD file from ',a,'-',a,/,
     .    '* Epoch  dL1 dL2 dP1 dP2 el1 el2 az1 az2 DF')
      endif

****  See if we need to find mean iondelay
      if( zion(1:1).eq.'Z' ) then           
         need1 = .true.
         need2 = .true.
         num = 0
         ng1 = 0
         ng2 = 0
         ngc = 0
         sion = 0.0d0
         vion = 0.0d0
         nbf = 0
         bf_mask = -1
         call sbit(bf_mask,32,0)
         call sbit(bf_mask,31,0)
         do while (err1.eq.0 .and.err2.eq.0 )
             if( need1) read(100, 200,iostat=err1)  
     .                      ep1, L11, L21, LC1, LG1, 
     .                      RC1, WL1,NL1, sv1, 
     .                      el1, az1,je1, df1
             if( need2) read(101, 200,iostat=err2)  
     .                      ep2, L12, L22, LC2, LG2, 
     .                      RC2, WL2,NL2, sv2, 
     .                      el2, az2,je2, df2
             if( ep1.eq.ep2 ) then
                 need1 = .true.
                 need2 = .true.
                 dlg = lg1 - lg2
                 dfc = ior(df1,df2)
                 if( cand(dfc,bf_mask).eq.0 ) then
                     jec = 0
*                    See if new bias flags
                     if( kbit(dfc,31) .or. kbit(dfc,32) ) then
*                        New bias flag.  Compute mean for previous
                         if( ngc.gt.1 ) then
                             nbf = nbf + 1
                             mion = sion/ngc
                             rmsi = sqrt((vion - mion**2*ngc)/(ngc-1))
                             write(*,150) nbf, mion, rmsi, ngc
 150                         format('* Mean iondelay for segment ',i4,
     .                               F8.2,' RMS ',  F8.3,' Num ',i5)
                             mean_ion(nbf) = mion
                         end if
                         vion = 0.d0
                         sion = 0.d0
                         ngc = 0.d0
                     end if
                 else
                     jec = 1
                 end if

                 if( jec.eq.0 ) then
                     ngc = ngc + 1
                     sion = sion + dlg
                     vion = vion + dlg**2
                 end if
             else if( ep1.gt.ep2 ) then
                 need1 = .false.
                 need2 = .true.
             else
                 need2 = .false.
                 need1 = .true.
             end if
         end do

****     Now compute the mean iondelay and its variance
         if( ngc.gt.1 ) then
             mion = sion/ngc
             rmsi = sqrt((vion - mion**2*ngc)/(ngc-1))
             nbf = nbf + 1
             mean_ion(nbf) = mion
             write(*,150) nbf,mion, rmsi, ngc
         end if
         rewind(100)
         read(100,'(a)') runstring
         read(100,'(a)',iostat=err1) runstring
         rewind(101)
         read(101,'(a)') runstring
         read(101,'(a)',iostat=err2) runstring
         

      else
         mion = 0.d0
      end if

***   Start looping over the files
      need1 = .true.
      need2 = .true.
      num = 0
      ng1 = 0
      ng2 = 0
      ngc = 0
      nbf = 0
      bf_mask = -1
      call sbit(bf_mask,32,0)
      call sbit(bf_mask,31,0)
      do while (err1.eq.0 .and.err2.eq.0 )
          if( need1) read(100, 200,iostat=err1)  
     .                   ep1, L11, L21, LC1, LG1, 
     .                   RC1, WL1,NL1, sv1, 
     .                   el1, az1,je1, df1
c         write(*,200) ep1, L11, L21, LC1, LG1, 
c    .                   RC1, WL1,NL1, sv1, 
c    .                   el1, az1,je1, df1
c         
 200      format(I5,1x,2(1x,F8.2), 1x,F8.3, 4(1x,F8.2),
     .           1x,i3,1x,2f10.4,1x, i1, 1x, o12, 2(1x,F20.4))
          if( need2) read(101, 200,iostat=err2)  
     .                   ep2, L12, L22, LC2, LG2, 
     .                   RC2, WL2,NL2, sv2, 
     .                   el2, az2,je2, df2
c         write(*,200) ep2, L12, L22, LC2, LG2, 
c    .                   RC2, WL2,NL2, sv2, 
c    .                   el2, az2,je2, df2

c          print *,'errs ', err1, err2
          if( ep1.eq.ep2 ) then
              need1 = .true.
              need2 = .true. 

*             Check the quality of the data
              dfc = ior(df1,df2)
              if( cand(dfc,bf_mask).eq.0 ) then
                  jec = 0
                  if( kbit(dfc,31) .or. kbit(dfc,32) .and. 
     .                zion(1:1).eq.'Z') then
*                     New bias flag.  Compute mean for previous
                      nbf = nbf + 1
                      mion = mean_ion(nbf)
                  end if
              else
                  jec = 1
              end if
              dl1 = l11 - l12 - mion
              dl2 = l21 - l22 - mion*(77.d0/60.0d0)
              dlc = lc1 - lc2
              dlg = lg1 - lg2 - mion
              drc = rc1 - rc2
              dwl = wl1 - wl2
              dnl = nl1 - nl2
              if( jec.eq.0 ) ngc = ngc + 1
              if( je1.eq.0 ) ng1 = ng1 + 1
              if( je2.eq.0 ) ng2 = ng2 + 1
            
              
              if( err1.eq.0 .and. err2.eq.0 ) then
                  num = num + 1
                  write(200,200) ep2, dL1, dL2, dLC, dLG, 
     .                   dRC, dWL,dNL, sv2, 
     .                   el2, az2, jec, dfc


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
      write(*,400) num, ng1, ng2, ngc
 400  format(I6,' Single differences File 1 Good ',i5,
     .          ' File 2 good ',i5,
     .          ' Good SD ',i5)
      close(100)
      close(101)
      close(200)
      end
 
 
 
 
