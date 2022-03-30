
      program tsjumps

      implicit none 

*     Routine to read a values file and analyze the time series of
*     individual series and detect and remove jumps.  Depending on
*     options passed a new values file can be written

      include 'tsjumps.h'

* MAIN PROGRAM VARIABLES

*  ierr -- iostat error variable

       integer*4 ierr, i

*  need_out -- Logical set if need to write the values file out
*      for a specific time series

       logical need_out

       character*8 access


****  Initialize the program variables

      call tsinit

****  Get the runstring to see options passed by used
      call get_runstring

****  If we passed the names of the org file; read it now so that we
*     have epochs 
      
      if( in_org ) then
          call read_org
      endif

****  Now loop over the values file, reading each of the time series
      open(100, file=values_in, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open', values_in, 1, 'TSJUMPS')

****  See if we will write an output values file
      if( out_vals ) then
          open(200,file=values_out, iostat=ierr, status='unknown')
          call report_error('IOSTAT',ierr,'creat', values_out, 1,
     .                      'TSJUMPS')
      endif

****  See if we writting a rename file and if we are to append.
      if( out_rn ) then
* MOD TAH 131111: Alloeq ~ substituion
          call subhome(rename_file)
          if( app_rn ) then 
              access = 'append'
              open(201,file=rename_file, iostat=ierr,status='unknown',
     .                 position=access)
          else
              open(201,file=rename_file, iostat=ierr,status='unknown')
          endif
          call report_error('IOSTAT',ierr,'open', rename_file, 1,
     .                      'TSJUMPS')
      endif  

**** Read the jumps file (eq file format)
        call read_jmp_file
        do i = 1, 10
          print*,'rn_codes, rn_times, rn_dpos, rn_types',rn_codes(1,i),
     .         rn_codes(2,i), rn_times(1,i),rn_times(2,i),rn_dpos(1,i),
     .         rn_dpos(2,i),rn_dpos(3,i),rn_types(i)
        enddo

****  Loop over values file
      do while ( ierr.eq. 0 )

*        Get the next time series
         call get_ts(100, ierr )

*        If no error then analyze
         if( ierr.eq.0 ) then
             call tsanal  

*            If we writing the time series, do so now
             need_out = out_vals
             if( out_new .and. num_jmp.eq.0 ) need_out = .false.

             if( need_out ) then
                  call write_vals( 200 )
             end if
             if( num_rn.gt.0 .and. comp_name(1:1).eq.'U' 
     .                       .and. out_rn ) then
                 call write_rn ( 201 )
             endif
         end if
      end do

****  Now close all the files
      close(100, iostat=ierr )
      close(200, iostat=ierr )
      close(201, iostat=ierr ) 

      end 

CTITLE TSINIT

      subroutine tsinit

      implicit none 

*     Routine to initialize the variables in tsjumps

      include 'tsjumps.h'

****  Set the default values

      values_in   = ' '
      values_out  = ' '
      org_in      = ' '
      rename_file = ' '
      jumps_file  = ' '

      out_vals = .false.
      out_rn   = .false.
      jmp_fl   = .false.
      in_org   = .false.
      app_rn   = .false.
      out_new  = .false.
      debug    = .false.
      rename_site = .false.

      scale_dres = 4.d0
      scale_sig  = 3.d0
      abs_dres   = 25.d0
      abs_slope  = 0.25d0
      max_sig    = 0.5d0

      num_org  = 0

****  Thats all 
      return
      end

CTITLE GET_RUNSTRING

      subroutine get_runstring

      implicit none 

*     Routine to decode the runstring and print help if needed

      include 'tsjumps.h'

* LOCAL VARIABLES
* ---------------

* na      --  Counter for runstring
* rcpar   -- Function to return lenght of runstring
* len_run -- Length of runstring  
* trimlen -- Length is string
* jerr    -- Error reading line

      integer*4 na, rcpar, len_run,  trimlen, jerr, i

* value   -- Value read from run-string 
      real*8 value
                                             
* runstring -- String returned from command line
      character*128 runstring

****  Write out header
      write(*,120)
 120  format(/,' Program TSJUMPS',/,' ---------------',/)

****  Loop over runstring
      na = 0
      len_run = rcpar(1, runstring )
      do while ( len_run.gt.0 )
         na = na + 1
         len_run = rcpar(na, runstring)

*        Now see what we should do: Check each argument types
*        -v[alues] file
         if( runstring(1:2).eq.'-v' ) then
*            Get the values file
             na = na + 1
             len_run = rcpar(na, values_in)

*        -p[rt/org] file
         else if( runstring(1:2).eq.'-p' ) then
*            Get the org/prt file
             na = na + 1
             len_run = rcpar(na, org_in)
             in_org = .true.

*        -o[ut] values
         else if( runstring(1:2).eq.'-o' ) then
*            Get the output values file
             na = na + 1
             len_run = rcpar(na, values_out)
             out_vals = .true.

*        -r[ename] file
         else if( runstring(1:2).eq.'-r' ) then
*            Get the output rename file
             na = na + 1
             len_run = rcpar(na, rename_file )
             out_rn = .true.

*        -a[ppend rename] -- Append to rename file
         else if( runstring(1:2).eq.'-a' ) then
             app_rn = .true.
         
*        -n[ew] values only
         else if( runstring(1:2).eq.'-n' ) then
             out_new = .true.

*        -s[cale] criteria for detecting breaks.  (Negative values will
*        use defualts
         else if( runstring(1:2).eq.'-s' ) then 
             do i = 1, 3
                na = na + 1
                len_run = rcpar(na,runstring)
                read(runstring,*,iostat=jerr) value
                if( len_run.gt.0 .and. jerr.eq.0 .and.value.gt.0) then
                    if( i.eq.1 ) scale_dres = value**2
                    if( i.eq.2 ) scale_sig = value
                end if
             end do

*        -m[ax] absolute values
         else if( runstring(1:2).eq.'-m' ) then 
             do i = 1, 3
                na = na + 1
                len_run = rcpar(na,runstring)
                read(runstring,*,iostat=jerr) value
                if( len_run.gt.0 .and. jerr.eq.0 .and.value.gt.0) then
                    if( i.eq.1 ) abs_dres = value**2
                    if( i.eq.2 ) abs_slope = value/1000.d0
                    if( i.eq.3 ) max_sig = value
                end if
             end do
         else if( runstring(1:2).eq.'-d' ) then
             debug = .true.

*        -f[irst] rename number. Sets the starting rename number
         else if( runstring(1:2).eq.'-f') then
             na = na + 1
             len_run = rcpar(na,runstring)
             read(runstring,*,iostat=jerr) value
             rn_off = nint(value)
             rename_site = .true.    

*        -j[umps] file
         else if( runstring(1:2).eq.'-j' ) then
*            Get the input jumps rename file
             na = na + 1
             len_run = rcpar(na, jumps_file )
             jmp_fl = .true.

         else
             if( len_run.gt.0 ) write(*,220) runstring(1:len_run)
 220         format('Unknown command line option: ',a)     
         endif
      end do

****  Now check status of options and print any warnings
      if( na.eq.0 ) then
          call proper_runstring('tsjumps.hlp','tsjumps',1)
      end if

****  See if rename and no prt/org file
      if( out_rn .and. .not. in_org ) then
         write(*,300) rename_file(1:max(1,trimlen(rename_file)))
 300     format('**WARNING** Rename file specified but no ',
     .          'prt/org file given',/,
     .          'Rename file ',a,' will not be created')
         out_rn = .false.
      end if

****  Write out options
      write(*,400) sqrt(scale_dres),  scale_sig,
     .             sqrt(abs_dres), abs_slope*1000.d0, max_sig
 400  format('OPTIONS:',/,'Scaled parameters: DRES ',f6.2,
     .       ' SIGNIFICANCE ',f6.2,/,
     .       'Absolutes: ABSDRES ',f6.2,' SLOPE ',f6.1,' mm/yr, SIGMA ',
     .       F6.3,' m',/)
 

****  Thats all
      return
      end

CTITLE READ_ORG

      subroutine read_org

      implicit none 

*     Routine to read a prt or org file and get the epochs
*     writing a rename file.

      include 'tsjumps.h'

* LOCAL VARIABLES
* ---------------

* ierr        -- IOSTAT errors
* date(5)     -- Date
* no          -- Number of entries found
* trimlen     -- Length of string

      integer*4 ierr,  date(5), no, trimlen

* sectag      -- Seconds tag
      real*8 sectag

* line        -- Line read from file
      character*256 line


****  Open the org/prt file and read
      open(100, file = org_in, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',org_in, 1, 
     .                  'read_org')

      no = 0
      sectag = 0.d0

      do while ( ierr.eq.0 ) 
         read(100,'(a)',iostat=ierr) line

*        See if these contain dates we need
         if( line(1:25).eq.' Solution commenced with:') then
             read(line,'(26x,i4,4(1x,i2))') date
             no = no + 1
             call ymdhms_to_jd(date,sectag,rn_ep(1,no))
         else if( line(1:25).eq.' Solution ended with    :') then
             read(line,'(26x,i4,4(1x,i2))') date
             call ymdhms_to_jd(date,sectag,rn_ep(2,no))
         else if( line(1:25).eq.' Solution refers to     :') then
             read(line,'(26x,i4,4(1x,i2))') date
             call ymdhms_to_jd(date,sectag,rn_ep(3,no))
         endif
      end do

*     Tell user
      num_org = no
      write(*,300) num_org, org_in(1:trimlen(org_in))
 300  format(' Found ',i4,' triplets of dates in ',a)

      close(100)

      return
      end

CTITLE GET_TS

      subroutine get_ts( unit, ierr )

      implicit none 

*     Routine to return the next block of time series.  When the
*     file is completely read ierr returns as -1

      include 'tsjumps.h'

* PASSED VARIABLES
* ----------------
* unit   -- Fortran unit number
* ierr   -- IOSTAT Error

      integer*4 unit, ierr

* LOCAL VARIABLES
* ---------------
* jerr  -- IOSTAT error
* nv    -- Counter for number of values
* i     -- Loop counter
* date(5)  -- Date of values estimate
* trimlen  -- Length of string

      integer*4 jerr, nv, date(5), i, trimlen

* sectag -- Seconds tag
* value, sigma -- TS value and sigma

      real*8 sectag, value, sigma

* line   -- Line read from file
      character*128 line 

* done  -- Logical to indicate that we have finished block read
      logical done


****  Read the next group of header lines
      nv = 0
      read(unit,'(a)', iostat=ierr) header(1)

*     If ERROR return
      if( ierr.ne.0 ) RETURN

****  Now get the next two lines
      read(unit,'(a)', iostat=ierr) header(2)
      read(unit,'(a)', iostat=ierr) header(3)

****  Now start reading the time series values
      done = .false.

      do while ( .not.done )
         read(100,'(a)',iostat=ierr) line

*        If line is not blank, then decode
         if( trimlen(line).gt.1 .and. line(1:1).eq.' ' ) then
             read(line,*,iostat=jerr) date, value, sigma
             if( jerr.eq.0 .and. nv.lt.max_ent) then
                nv = nv + 1
                call ymdhms_to_jd(date, sectag, ts_ep(nv))
                ts_val(nv) = value
                ts_sig(nv) = sigma
             end if
         else
             done = .true.
         end if
         if( ierr.ne.0 ) done = .true.
      end do

****  Now skip the next four lines with statistics
      do i = 1, 3
         read(unit,'(a)',iostat=jerr) line
      end do
      
      num_ent = nv 
      site_name = header(2)(1:8)
      comp_name = header(2)(13:13)
      
      if( debug ) write(*,*) ' '
      write(*,300) site_name, comp_name, num_ent
 300  format('For ',a8,1x,a1,1x,i5,' values found')

*     Now sort the time series into ascending order.
      call tssort

****  Thats all
      return
      end 

CTITLE TSANAL

      subroutine tsanal

      implicit none 

*     Routine to analyze the time series and detect jumps

      include 'tsjumps.h'


* LOCAL VARIABLES
* --------------- 

* np   -- Number of parameters being estimated
* iter -- Iteration counts 

      integer*4 np, iter, i, j, k 

* chi_tol -- Tolerance for chi**2: set to be mimimum of 10.d0
*     and the maxiumum of 1.0d0 and actual chi**2
* sig     -- Significance of jump

      real*8 chi_tol, sig

* jmp_removed -- Logical to indicate that jumps are removed
      logical jmp_removed

* done -- Logical to indicate when iteration has convergef

      logical done

****  Make sure we have values
      num_jmp = 0
      np      = 2
      if( num_ent-(num_jmp+2).lt.3 ) RETURN

****  Now iterate the analysis of the time series finding all the 
*     breaks.
      done    = .false.
      iter    = 0

      if ( .not.jmp_fl ) then
      do while ( .not.done )

*         Initilize the fitting of the data; clear normal equations
*         and stats.
          call init_est( np )

*         Now fit a linear regression and breaks to the time
*         series.
          call formsoln ( np )

*         Analyze the residuals to the time series to find and
*         breaks.  When there seems to be no more breaks, done returns
*         true.
          call find_jmp( np, done )

      end do 

****  Now do the final check to see which jumps are signficant.
      jmp_removed = .false.
      chi_tol = min(max(pos_chi,1.d0),4.d0)
      write(*,250) site_name, comp_name, pos_chi, 
     .             sol_vec(2),sqrt(norm_eq(2,2)), 
     .             (ts_jmp(j-2), sol_vec(j),sqrt(norm_eq(j,j)),j=3,np)
 250  format('Initial ',a,1x,a1,1x,'Chi**2/f ',F8.2,' Estimates ',
     .        F7.5,1x,F6.4, 1x, 200(i4,1x,F7.4,1x,F6.4,' : '))

      do i = 1, num_jmp
         sig = abs(sol_vec(i+2)/sqrt(norm_eq(i+2,i+2)))
         if( sig.lt. sqrt(chi_tol)*scale_sig ) then

*            This is not a significant jump, so remove it from the
*            list.
             ts_jmp(i) = -ts_jmp(i)
             jmp_removed = .true.
             if( debug )
     .       write(*,*) 'Removing Jump ',i,' at value ',ts_jmp(i),
     .                 ' Significance ', sig,' sigma'
         endif
      end do

****  See if any jumps removed
      if( jmp_removed ) then
         i = 0
         do while ( i.lt.num_jmp )
            i = i + 1
            if( ts_jmp(i).le.0 ) then
               do j = i+1, num_jmp
                  ts_jmp(j-1) = ts_jmp(j)
               end do
               num_jmp = num_jmp - 1
               i = i - 1
            end if
         end do
         if ( debug ) 
     .   write(*,*) 'Jumps reduced from ',np-2,  ' to ', num_jmp 
         np = num_jmp + 2
         call init_est( np )
         call formsoln ( np )   
      end if

****  Tell user final status:
      write(*,300) site_name, comp_name, pos_chi, 
     .             sol_vec(2),sqrt(norm_eq(2,2)), 
     .             (ts_jmp(j-2), sol_vec(j),sqrt(norm_eq(j,j)),j=3,np)
 300  format('Final   ',a,1x,a1,1x,'Chi**2/f ',F8.2,' Estimates ',
     .        F7.5,1x,F6.4, 1x, 200(i4,1x,F7.4,1x,F6.4,' : '))

****  Now save and merge the rename information
      if( out_rn ) call save_rn
      
      else 
         num_jmp = 0 
         if ( comp_name .eq. 'N' ) then
            k=1
         elseif ( comp_name .eq. 'E') then
            k=2
         elseif ( comp_name .eq. 'U') then
            k=3
         endif
         do i = 1,num_ent 
            do j = 1,num_renames  
c             print *,'rn_codes(1,j) site_name',rn_codes(1,j), site_name
               if ( rn_codes(1,j) .eq. site_name ) then
c                print*,'rn_times(1,j) ts_ep(i)',rn_times(1,j), ts_ep(i)
                  if ( ts_ep(i) .ge. rn_times(1,j) .and. 
     .               (.not. done_jmp(k,j)) ) then 
                     num_jmp = num_jmp + 1
                     ts_jmp(num_jmp) = i                    
                     print*,'comp_name',comp_name 
                     if ( comp_name .eq. 'N' ) then
                       done_jmp(k,j) = .true.
                       sol_vec(num_jmp+2) = rn_dpos(1,j)
                     elseif ( comp_name .eq. 'E') then
                       done_jmp(k,j) = .true.
                       sol_vec(num_jmp+2) = rn_dpos(2,j)
                     elseif ( comp_name .eq. 'U') then
                       done_jmp(k,j) = .true.
                       sol_vec(num_jmp+2) = rn_dpos(3,j)
                     endif
      print*,'num_jmp,comp_name,sol_vec(num_jmp+2)',num_jmp,comp_name,
     .        sol_vec(num_jmp+2) 
                  endif 
               endif
            enddo
         enddo 
      endif

****  Thats all
      return
      end

CTITLE INIT_EST

      subroutine init_est( np )

      implicit none 

*     Routine to initialize the estimates of the linear fit and
*     statistics

      include 'tsjumps.h'

 
* PASSED VARIABLES
* ----------------
* np  -- Number of parameters being estimated (2+num_jmp)

      integer*4 np

* LOCAL VARIABLES
* --------------- 
* i    -- Loop counter

       integer*4 i

****  Clear the normal equations
      call dwint( 0.d0, sol_vec, 1, np )
      do i = 1, np 
         call dwint( 0.d0, norm_eq(1,i), 1, np)
      end do

****  Clear the statistics accumulations
      call dwint( 0.d0, stats,    1, 4 )

****  Thats all
      return
      end

 
CTITLE FORMSOLN

      subroutine formsoln( np )

      implicit none 

*     Routine to accumulate normal equations and find the statistics
*     of the fit

      include 'tsjumps.h'

 
* PASSED VARIABLES
* ----------------
* np  -- Number of parameters being estimated (2+num_jmp)

      integer*4 np

* LOCAL VARIABLES
* ---------------
* i, j, k  -- Loop variables
* ipivot(max_jmp+2) -- Pivot elements in invert_vis
 
      integer*4 i,j,k, ipivot(max_jmp+2)

* wgh       -- Weight for each measurement 
* residual  -- Prefit residual (difference from first value).
* dstat     -- Change in chi**2 due to fit
* apart(max_jmp+2) -- Partials array: Note Jumps accumulate in value
* scale(max_jmp+2) -- Scaling used by inversion
* bvec(max_jmp+2)  -- Copy of sol_vec: used to compute change in chi**2

      real*8 wgh, residual, dstat, apart(max_jmp+2), 
     .       scale(max_jmp+2), bvec(max_jmp+2)


***** Loop over the data forming the normal equations
      do i = 1, num_ent

****     Form the partials array
         apart(1) = 1.d0
         apart(2) = (ts_ep(i)-ts_ep(1))/365.25d0

*        Now do the partials for the jumps.  The jumps are 
*        accumulative (ie each adding to the previous one)
         do j = 1, num_jmp
            apart(j+2) = 0.d0
            if( i.ge.ts_jmp(j) ) apart(j+2) = 1.d0
         end do

*        Now form up the contribution to the normal equations
         wgh = 1.d0/ts_sig(i)**2

         residual = ts_val(i)-ts_val(1)
         do j = 1, np
            sol_vec(j) = sol_vec(j) + residual*apart(j)*wgh
            do k = 1, np
               norm_eq(j,k) = norm_eq(j,k) + apart(j)*apart(k)*wgh
            end do
         end do

*****    Now accumulate the statistics
         stats(1) = stats(1) + residual*wgh
         stats(2) = stats(2) + wgh
         stats(3) = stats(3) + residual**2*wgh
         stats(4) = stats(4) + 1.d0

      end do

****  Now solve the system of equations
      call dwmov(sol_vec, 1, bvec, 1, np)
          
      call invert_vis(norm_eq, sol_vec, scale,ipivot,np,max_jmp+2,1)

****  Check the size of the rate:
      if( abs(sol_vec(2)).gt. abs_slope ) then 
          write(*,120) sol_vec(2)*1000.d0
 120      format('Rate too large (',F10.2,' mm/yr); setting to zero')
          sol_vec(2) = 0.d0
      end if
 
*     Now update the statistics to give values about fit
      call dwdot(dstat, sol_vec,1, bvec,1, np)

      pre_chi = (stats(3)-(stats(1)/stats(2))**2)/(num_ent-np)
      pos_chi = (stats(3)-dstat)/ (num_ent-np)
      
*     Now compute the Chi**2/f of the differnces between points
*     after removing the slope.
      do i = 1, num_ent-1
      end do

*     Now write out values (temporary)
      if( debug )
     .write(*,300) pos_chi, sol_vec(2),sqrt(norm_eq(2,2)), 
     .             (ts_jmp(j-2), sol_vec(j),sqrt(norm_eq(j,j)),j=3,np)
 300  format('Chi**2/f ',F8.2,' Estimates ',F7.5,1x,F6.4, 1x,
     .        200(i4,1x,F7.4,1x,F6.4,' : '))

****  Thats all
      return
      end

 
CTITLE FIND_JMP

      subroutine find_jmp( np, done )

      implicit none 

*     Routine to find jumps in time series by looking at the 
*     post-fit chi**2 and the magnitude of the largest jumps

      include 'tsjumps.h'

 
* PASSED VARIABLES
* ----------------
* np  -- Number of parameters being estimated (2+num_jmp)

      integer*4 np
* done -- Logical to indicate no more jumps

      logical done

* LOCAL VARIABLES
* ---------------

* i, j  -- Loop counters
* imax  -- Epoch with largest jump
* jmax  -- Epoch with largest residual 
* umax  -- Actual max jump or outlier epoch to use (either imax or jmax
*          depending on which is larger).
* ip    -- Index of previous point incase we skip a point due to
*          large sigma.

      integer*4 i,j, imax, jmax, umax, ip

* res, prev_res  -- Current residual and previous residual
* dres           -- Change in residual scaled by sigma
* max_dres       -- Maxiumum value of dres
* act_dres       -- Actual value of residual (not scaled by sigma)
* sres           -- Scaled residual
* max_sres       -- Maximum value of the residual itself
* next_res       -- Residual for next data point (used to to look for
*                   outliers).
* chi_tol        -- Chi**2 tolerance (between 1-4 depending on chi**2)

      real*8 res, prev_res, dres, max_dres, act_dres, sres, max_sres,
     .       next_res, chi_tol

* outlier  -- Set true is point seems to be outlier because next
*     timeseries value un-does jump
* brk_known -- Logical that checks to see if we already have a break
*     at a point with a large residual.

      logical outlier, brk_known

****  First loop at chi**2/f.  If this is less than 1 then say we
*     are done:
      if( pos_chi.lt. 1.d0 ) then
          done = .true.
      end if

****  Now scan residuals for largest jump
      max_dres = 0.d0
      max_sres = 0.d0

      if ( .not.done ) then
         prev_res = -sol_vec(1)
         ip = 1
         do i = 2, num_ent

*****       Ignored points with large sigmas (> max_sig)
            if( ts_sig(i).lt.max_sig ) then
               res = (ts_val(i)-ts_val(1)) - sol_vec(1) - 
     .                    sol_vec(2)*(ts_ep(i)-ts_ep(1))/365.25d0
               if( i.lt.num_ent ) then
                   next_res = (ts_val(i+1)-ts_val(1)) - sol_vec(1) - 
     .                    sol_vec(2)*(ts_ep(i+1)-ts_ep(1))/365.25d0
               end if

               do j = 1, num_jmp
                  if( i.ge.ts_jmp(j) ) then
                      res = res - sol_vec(j+2)
                  end if
                  if( i+1.ge.ts_jmp(j) ) then
                      next_res = next_res - sol_vec(j+2)
                  end if
 
               end do

               dres = (res-prev_res)**2/(ts_sig(i)**2+
     .                                    ts_sig(ip)**2)
               if( i.lt. num_ent ) then
                  sres = (res-(prev_res+next_res)/2)**2/
     .                (ts_sig(i)**2+(ts_sig(ip)**2+ts_sig(i+1)**2)/2)
               else
                  sres = 0.d0
               endif

C              if( debug ) 
C    .         write(*,220) 'Entry ',i,' DT ', ts_ep(i)-ts_ep(1),
C    .                    ' days, Residual ', res,
C    .                    ' m, DRES and ABSRES ', dres, sres
C220           format(a,i5,a,F10.2,a,F10.4,a,2F10.2)


               if( dres.gt.max_dres  ) then

*                  Make sure that we do not already have break at this 
*                  point
                   brk_known = .false.
                   do j = 1, num_jmp
                       if( i.eq.ts_jmp(j) ) brk_known = .true.
                   end do
                   if( .not.brk_known ) then 
                       max_dres = dres
                       act_dres = (res-prev_res)
                       imax = i
                   end if
               end if

*              Check the scaled residual
               if( sres.gt.max_sres  ) then

*                  Make sure that we do not already have break at this 
*                  point
                   brk_known = .false.
                   do j = 1, num_jmp
                       if( i.eq.ts_jmp(j) ) brk_known = .true.
                   end do
                   if( .not.brk_known ) then 
                       max_sres = sres
                       jmax = i
                   end if
               end if
 
               prev_res = res
               ip = i
            end if 

         end do
      end if

****  Now see how large the jump is:  If less than 4 sigma
*     then say we are done
      chi_tol =  min(max(pos_chi,1.d0),2.d0)
      if( debug )
     .write(*,*) 'MAX DRES ', sqrt(max_dres), ' MAX ABSRES ', 
     .            sqrt(max_sres),  ' Entries ', imax, jmax, chi_tol
      if( max_dres.lt.scale_dres*chi_tol .and. 
     .    max_dres.lt.abs_dres           .and.
     .    max_sres.lt.scale_dres*chi_tol .and. 
     .    max_sres.lt.abs_dres                 ) then
         done = .true.
      else
         if( max_dres.gt.max_sres ) then
             umax = imax
             res  = max_dres
             outlier = .false.
         else
             umax = jmax
             res  = max_sres
             outlier = .true.
         end if
         brk_known = .false.
         do j = 1, num_jmp
            if( ts_jmp(j).eq.umax ) brk_known = .true.
         end do
         
         if( .not. brk_known ) then   
            num_jmp = num_jmp + 1
            ts_jmp(num_jmp) = umax
            if ( debug )
     .      write(*,300) num_jmp, ts_jmp(num_jmp), res
 300        format('Found jump ',i4,' at Entry ',i4,' Size ', F10.2)
         end if

         if( outlier ) then
            brk_known = .false.
            do j = 1, num_jmp
               if( ts_jmp(j).eq.umax+1 ) brk_known = .true.
            end do
            
            if( .not. brk_known ) then 
                num_jmp = num_jmp + 1
                ts_jmp(num_jmp) = umax+1
                if( debug )
     .          write(*,320) num_jmp, ts_jmp(num_jmp)
 320            format('Found jump ',i4,' at Entry ',I4,
     .                 ' Added as outlier')
            end if
         end if

*        Update number of parameters estimated.  Make sure that we actually
*        added entries.
         if( np.eq.num_jmp+2 ) done = .true.
         np      = num_jmp + 2
       
      end if 

****  Now sort the jumps to keep them in time order
      do i = 1, num_jmp-1
         do j = 1, num_jmp -1
            if( ts_jmp(j).gt.ts_jmp(j+1) ) then
                jmax        = ts_jmp(j)
                ts_jmp(j)   = ts_jmp(j+1)
                ts_jmp(j+1) = jmax
            end if
         end do
      end do


****  Thats all
      return 
      end 

CTITLE TSSORT

      subroutine tssort

      implicit none 

*     Routine to sort the time series into ascending order.  (Needed
*     for the jump detection)

      include 'tsjumps.h'

* LOCAL VARIABLES
* i,j  -- Loop counters
      integer*4 i,j

* sort_needed -- true if sort is needed
      logical sort_needed

****  Use a bubble sort method.  First check to see if needed
      sort_needed = .false.
      do i = 2, num_ent
         if( ts_ep(i).lt.ts_ep(i-1) ) sort_needed = .true.
      end do

      if( sort_needed ) then
          do i = 1, num_ent - 1
              do j = 1, num_ent - 1

*                 See if we need to switch 
                  if( ts_ep(j).gt.ts_ep(j+1) ) then
                      call dwswp(ts_ep(j), 1, ts_ep(j+1), 1,1)
                      call dwswp(ts_val(j),1, ts_val(j+1),1,1)
                      call dwswp(ts_sig(j),1, ts_sig(j+1),1,1)
                  end if
              end do
          end do 
      end if

***** Thats all
      return
      end

CTITLE WRITE_VALS

      subroutine write_vals( unit )

      implicit none 

*     Routine to write out the values file with the jumps removed
*

      include 'tsjumps.h'

* PASSED VARIABLES
* unit  -- Unit number 

      integer*4 unit

* LOCAL VARIABLES 

* mean_jd  -- JD that means referrs to
* mean_year -- Corresponding as deciminl year
* duration  -- Duration in years
* mean_out  -- Mean length of timeseries
* offset    -- Offset in series as function of time
* sectag    -- Seconds tag
* sig_mean  -- Sigma of mean
* nrms      -- NRMS scatter of data
* wrms      -- WRMS scatter of data
* slope_out -- Output slope
* sig_slope -- Sigma of slope
                              
      real*8 mean_jd,  mean_year, duration, mean_out, offset, sectag,
     .       sig_mean, nrms, wrms, slope_out, sig_slope

* i,j  -- Loop counters
* date(5) -- Calender date
* trimlen -- Length of string

      integer*4 i,j, date(5), trimlen

***** Write out the header lines.  Use the third header to record
*     epochs of breaks
      write(unit,'(a,1x,a)' ) header(1)(1:max(1,trimlen(header(1)))),
     .                   'Breaks removed with TSJUMPS'
      write(unit,'(a)' ) header(2)(1:max(1,trimlen(header(2))))

      write(unit,120) num_jmp, (ts_jmp(i),i=1,num_jmp)
 120  format('* ',i4,' Jumps removed at entries ',128i6)

*     Now write the time-series
      do i = 1, num_ent
         call jd_to_ymdhms( ts_ep(i), date, sectag)
         offset = 0.d0
         do j = 1, num_jmp
            if( i.ge.ts_jmp(j) ) then
                offset = offset + sol_vec(j+2)
            end if
         end do

         ts_val(i) = ts_val(i) - offset
         write(unit, 150) date, ts_val(i), ts_sig(i), offset
 150     format(i5,4i3,f15.4,1x,f10.4,1x,2f10.4)
      end do

*     Now write the summary at the bottom
      mean_jd  = ts_ep(1) + (norm_eq(1,2)/norm_eq(2,2))*365.25d0
      call jd_to_decyrs(mean_jd, mean_year)
      duration = (ts_ep(num_ent)-ts_ep(1))/365.25d0 

      mean_out = ts_val(1) + sol_vec(1) + 
     .           sol_vec(2)*(ts_ep(1)-mean_jd)/365.25d0
      sig_mean = sqrt(norm_eq(1,1)*pos_chi)
      nrms     = sqrt(pos_chi)
      wrms     = sqrt((1.d0*num_ent)/stats(2))*nrms*1000.d0

      slope_out = sol_vec(2)*1000.d0
      sig_slope =  sqrt(norm_eq(2,2)*pos_chi)*1000.d0
      
      write(unit,200) mean_out, sig_mean, num_ent, wrms, nrms
 200  format(/,"Wmean ",f13.4," m +- ",f6.4," from ",i4," data.",
     .   " WRMS ",f5.1," mm, NRMS ",f5.2)
      write(unit,220) slope_out, sig_slope, wrms,
     .   nrms, duration, mean_year
 220  format("Slope ",f7.1," +- ",f6.1," mm/yr, WRMS ",f5.1,
     .   " mm, NRMS ",f5.2,", dur ",f4.1," <> ",f6.1," yr",/)

****  Thats all
      return
      end

CTITLE SAVE_RN

      subroutine save_rn

      implicit none 

*     Routine to save and merge the rename information about jumps
*

      include 'tsjumps.h'

* PASSED VARIABLES
* none


* LOCAL VARIABLES 
* chi_tol  -- Chi-square to use for writting sigmas
* jdep     -- Julian date of epoch

      real*8 chi_tol, jdep

* i,j,k    -- Loop counters
* imax, jmax -- Positions in arrays needed

      integer*4  i,j,k, imax, jmax

* found    -- Logical to indicate that value found
      logical found


****  If this is the N component then this is first and so we just save.
      chi_tol = min(max(pos_chi,1.d0),4.d0)

      if( comp_name(1:1).eq.'N' ) then
          do i = 1, num_jmp
             rn_jest(1,1,i) = sol_vec(i+2)
             rn_jest(2,1,i) = sqrt(norm_eq(i+2,i+2))*chi_tol

*            Clear the E and U values
             do j = 2,3
                rn_jest(1,j,i) = 0.d0
                rn_jest(2,j,i) = 0.d0
             end do

*            Now find this epoch in the list of org times.  If not present
*            then add to list.
             jdep = ts_ep(ts_jmp(i))
             rn_jmp(i) = 0

             do j = 1, num_org
                if( jdep.eq.rn_ep(3,j) ) then

*                   There is no reliable method for handling duplicate
*                   epochs, so just assign the epochs
                    rn_jmp(i) = j
                end if
             end do

****         See if we found a match
             if( rn_jmp(i).eq.0 ) then
                 write(*,*) 'WARNING -- No org/prt entry for JD ',
     .                      jdep, site_name, comp_name(1:1), i
                 num_org = num_org + 1
                 rn_ep(3,num_org) = jdep
                 rn_ep(1,num_org) = jdep - 0.5008d0
                 rn_ep(2,num_org) = jdep + 0.5008d0
                 rn_jmp(i) = num_org
             end if
          end do
          num_rn = num_jmp
      else

****      We have the E or U component. See if can find a value in the 
*         current list.
          do i = 1, num_jmp
             found = .false.
             imax  = 0
             do j = 1, num_rn
                jdep = rn_ep(3,rn_jmp(j))
                if( jdep.lt.ts_ep(ts_jmp(i)) ) imax = j 
                if( jdep.eq.ts_ep(ts_jmp(i)) ) then
*                   OK, match found so add jump estimate
                    found = .true. 
                    if( comp_name(1:1).eq.'E' ) then
                        rn_jest(1,2,j) = sol_vec(i+2)
                        rn_jest(2,2,j) = sqrt(norm_eq(i+2,i+2))*chi_tol
                    else
                        rn_jest(1,3,j) = sol_vec(i+2)
                        rn_jest(2,3,j) = sqrt(norm_eq(i+2,i+2))*chi_tol
                    end if
                end if
             end do

*            See if we found it:  If not we will need to add to list
             if( .not.found ) then

****             Scan the prt/org times to find value.  If not found then
*                add to list
                 jdep = ts_ep(ts_jmp(i))
                 jmax = 0

                 do j = 1, num_org
                    if( jdep.eq.rn_ep(3,j) ) then

*                       There is no reliable method for handling duplicate
*                       epochs, so just assign the epochs
                        jmax = j
                    end if
                 end do

****             See if we found a match
                 if( jmax.eq.0 ) then
                     write(*,*) 'WARNING -- No org/prt entry for JD ',
     .                          jdep, site_name, comp_name(1:1), i
                     num_org = num_org + 1
                     rn_ep(3,num_org) = jdep
                     rn_ep(1,num_org) = jdep - 0.5008d0
                     rn_ep(2,num_org) = jdep + 0.5008d0
                     jmax = num_org
                 end if

*                Now add values.  Insert immediately after the imax entry
*                which the last entry before the one we need
                 do j = num_rn, imax+1, -1
                    do k = 1,3
                       rn_jest(1,k,j+1) = rn_jest(1,k,j)
                       rn_jest(2,k,j+1) = rn_jest(2,k,j)
                    end do
                    rn_jmp(j+1) = rn_jmp(j)
                 end do

****             Now add current entry
                 imax = imax + 1
                 rn_jmp(imax) = jmax
                 rn_jest(1,1,imax) = 0.d0
                 rn_jest(2,1,imax) = 0.d0
                 if( comp_name(1:1).eq.'E' ) then
                     rn_jest(1,1,imax) = 0.d0
                     rn_jest(2,1,imax) = 0.d0 
                     rn_jest(1,2,imax) = sol_vec(i+2)
                     rn_jest(2,2,imax) = sqrt(norm_eq(i+2,i+2))*chi_tol 
                     rn_jest(1,3,imax) = 0.d0
                     rn_jest(2,3,imax) = 0.d0
                 else
                     rn_jest(1,1,imax) = 0.d0
                     rn_jest(2,1,imax) = 0.d0
                     rn_jest(1,2,imax) = 0.d0
                     rn_jest(2,2,imax) = 0.d0
                     rn_jest(1,3,imax) = sol_vec(i+2)
                     rn_jest(2,3,imax) = sqrt(norm_eq(i+2,i+2))*chi_tol 
                 end if

*                Save number of renames 
                 num_rn = num_rn + 1
              end if
          end do
      end if
 
****  Thats all
      return
      end 

CTITLE WRITE_RN

      subroutine write_rn( unit )

      implicit none 

*     Routine to write the rename entries into the rename file.
*

      include 'tsjumps.h'

* PASSED VARIABLES
* unit  -- Unit number for write

      integer*4 unit

* LOCAL VARIABLES

* date(5), start_date(5), end_date(5) -- Various dates
* i,j,k  -- Looop counters
      integer*4 date(5), start_date(5), end_date(5), i,j,k

* sectag  -- Seconds tag
* offs(3) -- NEU offset of site
* start_ep, end_ep -- Start and stop epoch for rename 

      real*8 sectag, start_ep, end_ep, offs(3)

* out_name -- Name for renamed site

      character*8 out_name 



****  Get the time the updates written
      if( num_rn.eq.0 ) RETURN

      call systime(date, sectag)
      write(unit,120) site_name, date
 120  format('* Rename updates for ',a,' by TSJUMPS ',i5,4(':',i2.2))

*     Now loop over the renames we have:  We have pointers to epochs
*     when jumps occurr.  What we want to do know is compute the position
*     changes to line up all of the values.  We would like the last points
*     to the coordinates adopted.
      do i = 1, num_rn

*        Get the total offset to this point 
         do j = 1,3
            offs(j) = 0.d0
            do k = i, num_rn
               offs(j) = offs(j) + rn_jest(1,j,k)
            end do
         end do

*        Now get the epoch range.  For first epoch start 1 minute before 
*        first experiment
         if ( i.eq.1 ) then
             start_ep =  rn_ep(1,1) - 0.0008d0 
         else
*            use the start epoch of the experiment before.
             start_ep = rn_ep(1,rn_jmp(i-1))
         end if
         end_ep = rn_ep(2,rn_jmp(i)-1) + 0.0008d0

*****    Now convert back to date
         call jd_to_ymdhms(start_ep, start_date, sectag)
         call jd_to_ymdhms(end_ep  , end_date,   sectag)

*****    If we are renaming then get new name:
         out_name = site_name
         if( rename_site ) then
            write(out_name(7:8),'(i2.2)') i+rn_off
         end if
         write(unit, 220) site_name, out_name, start_date, end_date,
     .                    offs
 220     format(1x,'RENAME ',a8,1x,a8,1x,i4,4(1x,i2),2x,i4,4(1x,i2),2x,
     .          3F12.4,' NEU ')
      end do
            

****  Now set the number of renames to zero for the next site
      num_rn = 0

****  Thats all
      return
      end
 
     
CTITLE READ_JMP_FILE
 
      subroutine read_jmp_file

      implicit none 
 
      include 'tsjumps.h'
 
* LOCAL VARIABLES
 
*   ierr    - IOSTAT error
*   trimlen - Length of string
*   iel     - Command number from list of earthquake commands
*   indx    - Position in string that we have decoded to.
 
      integer*4 ierr, trimlen, iel, indx, max_jmp_cmd
 
      character*8 jmp_commands

*   buffer  - Line read from input file
      character*256 buffer  

      num_renames = 0
 
****  See if the file name has been given
      if( trimlen(jumps_file).gt.0 ) then
* MOD TAH 131111: Allow ~ substitution
          call subhome( jumps_file )
          open(106,file=jumps_file, status='old', iostat=ierr)
          call report_error('IOSTAT',ierr,'open',jumps_file,0,
     .                    'read_jmp_file/blsum')
          if( ierr.ne.0 ) jumps_file = ' '
      end if
 
***   Now start reading the file and extracting commands
 
      do while ( ierr.eq.0 )
          read(106,'(a)',iostat=ierr) buffer
          if( ierr.eq.0 .and. buffer(1:1).eq.' ' .and.
     .        trimlen(buffer).gt.0 )  then
 
*             OK, Valid line for command.  Process:
              indx = 0                   
              jmp_commands = 'RENAME'
              max_jmp_cmd  = 1
              call get_cmd(buffer, jmp_commands, max_jmp_cmd, iel,
     .                        indx)
 
*             Process the command
              call proc_jmp_cmd(1, indx, buffer )
 
          end if
      end do
 
****  Close the input file
      if( ierr.eq.-1 ) close(106)
      return
      end
 
CTITLE PROC_EQ_CMD
 
      subroutine proc_jmp_cmd(iel, indx, buffer ) 

      implicit none 

      include 'tsjumps.h'

*     Routine to process the earthquake commands
 
* LOCAL VARIABLES
 
*   ierr    - IOSTAT error
*   trimlen - Length of string
*   iel     - Command number from list of earthquake commands
*   indx    - Position in string that we have decoded to.
 
*   jel     - Earthquake number
*   jndx    - Index when looking up code
*   date(5) - YMDHM of earthquake
*   jerr    - Error to test if hfile name passed in rename
*   indx_curr - Current value of indx incase no hfile name passed
 
      integer*4 ierr, trimlen, iel, indx, date(5), j
 
*   geod(3) - Geodetic postion of earthquake (co-lat, long,
*           - ellipsiodal height.
*   sectag  - Seconds tag
*   lat, long       - Decimal degrees for epicenter
*   rad_km, dep_km  - Radius of influence and depth in
*           - km (converted to meters for internal use)
 
      real*8 sectag
 
*   buffer  - Line read from input file
 
      character*(*) buffer
 
*   test_code   - TEst code for Earth quake to see if
*           - duplicate or new event.
*	cdum    - Dummy string for readline
*   rn_test     - Test string to see if hfile name passed in rename
 
      character*8 cdum
 
****  See if there was an error
      if( iel.le.0 ) then
          call report_error('Get_command',ierr,'decod',buffer,0,
     .            'proc_jmp_cmd')
          RETURN
      end if
 
****  Process the command:
      goto ( 600 ) iel
 
 600  CONTINUE
 
*             Make sure we do not have too many renames already.
              if( num_renames+1.gt.max_rn ) then
                  write(*,610) max_rn, (buffer(1:trimlen(buffer)))
 610              format('**WARNING** Too many renames specified.',
     .                ' Maximum allowed ',i5,/,
     .                ' Ignoring ',a)
                  RETURN
              end if
 
              num_renames = num_renames + 1
              call getword(buffer, rn_codes(1,num_renames),indx)
              call getword(buffer, rn_codes(2,num_renames),indx)
              call casefold(rn_codes(1,num_renames))
              call casefold(rn_codes(2,num_renames))

* MOD TAH 971112: Check to see if a name restriction has been included
*             in the line
c              indx_curr = indx
c              call getword(buffer, rn_test, indx )
*             See if number or string
c              call check_num( rn_test, jerr )
c              if( jerr.ne.0 ) then
c                  rn_hfiles(num_renames) = rn_test
c              else
c                  rn_hfiles(num_renames) = ' '
c                  indx = indx_curr
c              end if
 
*             Now start reading the optional parts of the command.
*             Start time for the rename (ymdhms)
              do j = 1,5
                  call read_line(buffer,indx,'I4',ierr,date(j),cdum)
                  if( ierr.ne.0 ) date(j) = 1
              end do
*             If there is an error, default the start to 1900
              if( ierr.ne.0 ) date(1) = 1900
 
*             Process the starting data
              sectag = 0.0d0
              call ymdhms_to_jd(date, sectag, rn_times(1,num_renames))
 
****          Get the end time:
              do j = 1,5
                  call read_line(buffer,indx,'I4',ierr,date(j),cdum)
                  if( ierr.ne.0 ) date(j) = 1
              end do
*             If there is an error, default the end to 2100
              if( ierr.ne.0 ) date(1) = 2100
 
*             Process the starting data
              sectag = 0.0d0
              call ymdhms_to_jd(date, sectag, rn_times(2,num_renames))
 
****          See if there is a displacement associated with this
*             rename:
              do j = 1,3
                  call read_line(buffer,indx,'R8',ierr,
     .                rn_dpos(j,num_renames),cdum)
                  if( ierr.ne.0 ) rn_dpos(j,num_renames) = 0.d0
              end do
 
****          Now see what frame this is (XYZ/NEU)
              call read_line(buffer,indx,'CH',ierr,sectag,
     .                rn_types(num_renames))
              call casefold(rn_types(num_renames))
              if( ierr.ne.0 ) rn_types(num_renames) = 'XYZ'
              RETURN
 
***** Thats all
          END

