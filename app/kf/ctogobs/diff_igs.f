      program diff_igs

      implicit none

*     Program to difference the clocks between two IGS clock files.


      integer*4 max_ent, max_code

      parameter ( max_ent = (32+50)*10000 )
      parameter ( max_code = 1024 )

      real*8 epoch(2,max_ent), clk(2,max_ent), sig(2,max_ent)
      character*2 type(2,max_ent)
      character*4 code(2,max_ent), uniq_code(max_code)

* Miscellaneous variables
      integer*4 i, j, k, date(5), ierr, num_ent(2), num, trimlen,
     .          rcpar, len_run,  unit, jerr, ref(2),
     .          start(2), end(2), beg, ns, nlen, num_code, nc
      real*8 sectag, dclk, doy, values(2), epc(2), oclk
      real*8 stats(3,max_code), mean, rms, std

      character*4 ref_code

      character*256 infile(2), outfile, line 

      logical done, found

****  Get the two input files to be compared
      len_run = rcpar(1,infile(1))
      if( len_run.eq.0 ) then
          call proper_runstring('diff_igs.hlp','diff_igs',1)
      end if
      open(101, file = infile(1), iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',infile(1),1,'DIFF_IGS') 

      len_run = rcpar(2,infile(2))
      if( len_run.eq.0 ) then
          call proper_runstring('diff_igs.hlp','diff_igs',1)
      end if
      open(102, file = infile(2), iostat=ierr, status='old') 
      call report_error('IOSTAT',ierr,'open',infile(1),1,'DIFF_IGS') 

      len_run = rcpar(3,ref_code)
      if( len_run.eq.0 ) then
          call proper_runstring('diff_igs.hlp','diff_igs',1)
      end if

      len_run = rcpar(4,outfile)
      if( len_run.eq.0 ) then
          outfile = ' '
      end if
 
 
****  Generate output name 
      nlen = trimlen(infile(1)) 
      if( trimlen(outfile).eq.0 ) then
          outfile = infile(1)(1:3) // infile(2)(1:3) // 
     .              infile(1)(4:nlen-4) //
     .              '.' // ref_code
      end if

      write(*,130) outfile(1:trimlen(outfile))
 130  format('DIFF_IGS: Creating ',a)

      open(200,file=outfile,iostat=ierr, status = 'unknown')
      call report_error('IOSTAT',ierr,'open',outfile,1,'DIFF_IGS') 
 

****  Now read the two input files
      num_code = 0
      do j = 1,2
         ierr = 0
         ns = 0
         unit = 100 + j

         do while ( ierr.eq.0 )
            read(unit,'(a)',iostat=ierr) line
            if( ierr.eq.0 .and. 
     .         (line(1:3).eq.'AR ' .or. line(1:3).eq.'AS') ) then

*               Valid clock line found
                ns = ns + 1
                read(line,220,iostat=jerr) type(j,ns), code(j,ns),
     .               date, sectag, num, (values(k),k=1,2)
 220            format(a2,1x,a4,1x,i4,4i3,1x,f9.6,1x,i2,
     .             2x,2(1x,E19.12))
                call casefold(code(j,ns))
                if( code(j,ns)(2:2).eq.' ' ) code(j,ns)(2:2) = '0'
                call ymdhms_to_mjd(date,sectag, epoch(j,ns))
                clk(j,ns) = values(1)
                if (jerr.eq.-1 .or. 
     .              index(infile(1),'gfz').gt.0 .or.
     .              index(infile(1),'esa').gt.0 ) then
                    num = 1
                    sig(j,ns) = 1.d-9
                else
                    num = 2
                    sig(j,ns) = values(2)
                endif
                if( sig(j,ns).eq.0.d0  ) then
                    sig(j,ns) = 200.0d-9
                end if
*               Start making a list of code in first file
                if( j.eq.1 ) then
                    found = .false.
                    do k = 1, num_code
                       if( uniq_code(k).eq.code(j,ns) ) then
                           found = .true.
                       end if
                    enddo
                    if( .not.found ) then
                        num_code = num_code + 1
                        uniq_code(num_code) = code(j,ns)
                    end if
                end if
 
            end if
         end do
         num_ent(j) = ns
         write(*,240) infile(j)(1:trimlen(infile(j))), num_ent(j)
 240     format('Input ',a15,' has ',i8,' entries')
      end do

****  Now loop over the entries matching times and types
      do j = 1, 2
         start(j) = 1
         end(j) = 0
      end do
      do j = 1,num_code
         do k = 1,3
            stats(k,j) = 0.d0
         end do
      end do 

      done = .false.
      do while (.not.done)

****     Get the next block of data in reference file
         end(1) = start(1)
         epc(1) = epoch(1,end(1))
         ref(1) = 0
         if( code(1,end(1)).eq.ref_code ) ref(1) = end(1)
         do while ( abs(epc(1)-epoch(1,end(1))).lt.1.d-5 .and.
     .              end(1).lt.num_ent(1) )
             end(1) = end(1) + 1
             if( code(1,end(1)).eq.ref_code .and.
     .           abs(epc(1)-epoch(1,end(1))).lt.1d-5 ) ref(1) = end(1)
         end do

         if( end(1).ge.num_ent(1) ) then
             done = .true.
         else
             end(1) = end(1) -1 
         end if

****     Now see if we can find matching epoch
         beg = end(2) + 1
         start(2) = 0
         ref(2) = 0
         do i = beg, num_ent(2)
            if( abs(epoch(2,i)-epc(1)).lt.1.d-5 .and. 
     .                                           start(2).eq.0 ) then
                start(2) = i
            end if
            if( abs(epoch(2,i)-epc(1)).lt.1.d-5 ) end(2) = i
            if( abs(epoch(2,i)-epc(1)).lt.1.d-5  .and. 
     .          code(2,i).eq.ref_code ) ref(2) = i
            if( epoch(2,i).gt.epc(1) ) exit
         end do

****     If start(2) is set then we have a matching epoch so
*        output
         if( start(2).ne.0 ) then

*            Get the clock correction to allign the clocks
             dclk = 0
             if( ref(1).ne.0 ) then 
                dclk = clk(1,ref(1))
             end if
             if( ref(2).ne.0 ) then
                dclk = dclk - clk(2,ref(2))
             end if
   
             call jd_to_yds(epc(1), date(1), date(2), date(3))
             doy = date(2) + date(3)/86400.d0 

*            Now loop over file 1 finding entries in file 2
             do i = start(1), end(1)
                do j = start(2),end(2)
                   if( code(1,i).eq.code(2,j) ) then
                       oclk = (clk(1,i) - clk(2,j) - dclk)*1.d9
                       if( sig(1,i).lt. 10.d-9 .and. 
     .                     sig(2,j).lt.10.d-9       ) then
                           write(200,310) doy, oclk, sig(2,j)*1e+9,  
     .                                   code(1,i), type(1,i)
 310                       format(f20.6,2F20.3, 1x, a4,1x,a2)

*                          Accumlate statistics
                           nc = -1
                           do k = 1, num_code
                              if( code(1,i).eq.uniq_code(k) ) then
                                 nc = k
                              end if
                           end do
                           if( nc.gt.0 ) then
                              stats(1,nc) = stats(1,nc) + 1
                              stats(2,nc) = stats(2,nc) + oclk
                              stats(3,nc) = stats(3,nc) + oclk**2
                           endif

                       end if
                   end if
                end do
             end do
          end if
          start(1) = end(1) + 1
          start(2) = end(2) + 1
      end do

****  Now write out the statistics of the fit
      write(*,410)
 410  format('Clock fit statistics',/,   
     .       '  # CODE  Num       Mean       RMS     StdDev',/,
     .       '                    (ns)       (ns)     (ns)')

      do j = 1, num_code
         if( stats(1,j).gt.0 ) then
             mean = stats(2,j)/stats(1,j)
             rms = sqrt(stats(3,j)/stats(1,j))
             std = sqrt(stats(3,j)/stats(1,j)-mean**2)
         else
             mean = 0.d0
             rms = 0.d0
             std = 0.d0
         endif
         write(*,450) j, uniq_code(j), int(stats(1,j)), mean,
     .                rms, std
 450     format(i3,1x,a4,1x,i4,1x,3F10.3)
      end do



****  Thats all
      close(200)
      end



              


 



