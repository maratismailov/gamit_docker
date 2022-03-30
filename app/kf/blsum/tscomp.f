      program tscomp

      implicit none

*     Program to perform math on PBO times series.  First operations are
*     differences and mean values.
*
*     Usages:
*     tscomp <dir> <prod_id> <option> <pbo .pos files>
*     
*     where <dir>  -- directory to put the time series in.  
*           <prod_id> is product id with form: jpl.diffs_frame.  Characters
*                  5:9 will used for time series type (normally rapid or final)
*           <option> -- either DIFF or MEAN
*           <list of pbo .pos files>
*     
*     PROD_ID types valid in PBO:
*     Format <cen>.<series>_<frame>+<type>
*     <cen>     - Center 3 characters (bsl cwu pbo mit) or special product
*                 code e.g., aug for Augustine Volcano, aku for Akutan volcano
*     <series>  - "diffs" or "means" would be common choice.
*     <frame>   - Frame type: loose or frame 
*     <type>    - Optional type.  If not given series name is used.  Additional
*                 type put in the final series is supplemental run (suppl).
*     
*     Stsndard PROD_ID
*     pbo.rapid_frame
*     pbo.final_frame
*     pbo.final_frame+suppl

      integer*4 max_files
      parameter ( max_files = 4096 ) 

      include 'tsfit.h'
      include 'tscon.h'

      integer*4 len_run, nr, ierr, ns, rcpar, trimlen, i,j 
     .,      unc   ! Unit number for command file (0 is not given or error)
     .,      fl   ! File name length
     .,      indx ! String pointer
     .,      num_files  ! Number of files read
     .,      ent_start(max_files)   ! First entry for each file read
     .,      cent(max_files)        ! Entry numbers that align the times in the 
                  ! epoch being processed.
     .,      eent(max_files)        ! Last entry for each site
     .,      nm   ! Number of values available for mean calculation

      logical list  ! Set true is reading list of files
     .,       cont  ! Set true if we should continue loopsing
     .,       done  ! Set true no more files given
     .,       use_ref  ! Set true if a pbo ts files is read and
                    ! has a reference coordinate in it.
     .,       skip  ! Set true of next entry should be skipped while we
                    ! allign the epochs of the data.
     .,       mean_all   ! Set true but adding ALL to MEAN option and requiring
                    ! data set be available to take the mean value
     .,       mean_diff  ! Set true but adding DIFF to MEAN option.  This
                    ! is used to find the mean of the differences between
                    ! a set of differenced files that are assumed to be zero
                    ! mean and rate.
     .,       OK(max_files)   ! Set true when entry can be used in mean value

      logical edit       ! Set true if edits to be applied from editfile
     .,       app_jump   ! Set true if jumps are be applied as well.

      character*256 runname ! Name of file read from runstring
      character*256 editfile   ! Name of eq_formatted with edits and jumps
      character*16 cdum
      character*4 comp_op   ! Option DIFF or MEAN

      real*8 emjd  ! Evaluation MJD for stablizing coordinates
      real*8 curr_mjd   ! Current MJD of epoch being processed
      real*8 end_mjd    ! Last MJD of epochs being processed.

      character *80 line    ! DEBUG 

****  OK, Read the runsting (for the moment just by position)
      len_run = rcpar(1,tsdir)
      ln_tsdir = len_run
      if( len_run.eq.0 ) then
        call proper_runstring('tscomp.hlp','tscomp',1)
      endif
      len_run = rcpar(2,prod_id)

* MOD TAH 051129: Extract type from prod_id: Format pbo.final_frame+suppl
      ts_ref_type = prod_id(5:9)
      if( len_run.gt.15 ) then
          ts_ref_type = prod_id(17:21)
          prod_id = prod_id(1:15)
      endif

      ln_prod_id = trimlen(prod_id)

      ts_refresh = .true.
      edit = .false.
      app_jump = .false.

***** Get command file name
      nr = 3
      len_run = rcpar(nr,cdum)
      comp_op = 'DIFF'
      if( len_run.gt.0 ) then
          call casefold(cdum)
          if( cdum(1:1).eq.'M' ) then
             comp_op = 'MEAN'
             mean_all = .false.
             mean_diff = .false.
*            See if all needed from mean
             if( index(cdum,'AL').gt.0 ) then
                 mean_all = .true.
             endif
*            See if mean of differences
             if( index(cdum,'DIF').gt.0 ) then
                 mean_diff = .true.
             endif
*            See if edit/jump file
             if( index(cdum,'ED').gt.0 ) then
                 edit = .true.
                 if( index(cdum,'JU').gt.1 ) then
                     app_jump = .true.
                 endif
                 ! Get file name (tsfit rep_edits file)
                 nr = nr + 1
                 len_run = rcpar(nr,editfile)
                 open(101,file=editfile,iostat=ierr,status='old')
                 call report_error('IOSTAT',ierr,'open',editfile,1,
     .                             'TSCOMP')
             endif 
          endif
      endif 
      
***** OK Loop over input files
      num_ent = 0
      num_site = 0
      num_code = 0 
      len_run = 1
      first_mjd = 1.d20
      last_mjd = 0.d20
      save_ref_xyz = 0
      sigma_scale = 1.d0
      tsprog = 'tscomp'
      reference_frame = ''

****  See if file with list or if list passed
      done = .false.
      do while ( .not.done )
         nr = nr + 1
         len_run = rcpar(nr,in_file)
         if( len_run.gt.0 ) then 
             num_files = nr - 3
             ent_start(num_files) = num_ent + 1
             call read_in_pbo(num_files)
             eent(num_files) = num_ent
         else
             done = .true.
         endif 
      end do

***** OK say what is happening
      write(*,210) num_site, num_files, comp_op, mean_all, mean_diff
 210  format('TSCOMP: ',i5,' sites found in ',i5,' files.  Option ',a,
     .      ' MEANALL ',L1,' MEANDIFF ',L1)
      if( edit ) then
          write(*,220) trim(editfile), app_jump
 220      format('TCOMP: Edits from ',a,' applied.  APP_JUMP ',L1)
      endif

***** Now start looping (difference first)
      curr_mjd = in_mjd(1) 
      end_mjd  = in_mjd(1) 
      do j = 1, num_files
         cent(j) = ent_start(j)-1
         if( in_mjd(cent(j)+1).lt.curr_mjd ) 
     .      curr_mjd = in_mjd(cent(j)+1)     ! Needed for mean values
         if( in_mjd(eent(j)).gt. end_mjd ) 
     .      end_mjd = in_mjd(eent(j))
      end do

***   CODE for differences
      if ( comp_op(1:1).eq.'D' ) then 
         num_ts = 0
         do i = 1, ent_start(2)-1
            cent(1) = i
            cent(2) = cent(2)+1
            if( cent(2).gt. num_ent ) then
               exit    ! Exit do loop; we are done
            endif 
            if( abs(in_mjd(cent(1))-in_mjd(cent(2))).lt.0.001 ) then
               call diff_xyz(cent) 
            elseif( in_mjd(cent(2)).lt. in_mjd(cent(1))) then
*              need to advance cent(2)
               skip = .true.
               do while ( skip )
                  cent(2) = cent(2)+1
                  if( cent(2).gt. num_ent ) then
                      skip = .false.   
                      exit      ! Exit do loop; we are done
                  endif 
                  if( abs(in_mjd(cent(1))-in_mjd(cent(2))).lt.0.001 ) 
     .                                                          then
                      skip = .false.
                      call diff_xyz(cent)
                  elseif( in_mjd(cent(2)).gt. in_mjd(cent(1)) ) then
                      skip = .false.
                  end if
               end do
            else     ! Hold cent(2) at current value
               cent(2) = cent(2) - 1    ! Will be incremented next loop
            end if
         end do
      else
*        Generate the MEAN of the time series provided.
         num_ts = 0
*        Start at curr_mjd and loop until no data left
         do j = 1,num_files
             cent(j) = cent(j)+1
         enddo


         do while ( curr_mjd.le.end_mjd )
*            See what entries are available at this time.
             nm = 0
             do j = 1, num_files
                OK(j) = .false.
                if( cent(j).gt.0 .and. 
     .              abs(in_mjd(cent(j))-curr_mjd).lt.0.001 ) then
                    OK(j) = .true.
                    nm = nm + 1   ! Number of values in the mean
                endif
             end do
             if( nm.gt.0 ) then
                 if( (mean_all .and. nm .eq.num_files) .or. 
     .                .not.mean_all ) then
!                    call mean_xyz(cent, OK, num_files, mean_diff)
* MOD TAH 210309: Added passing num_files instead of nm
                     call mean_xyz(cent, OK, num_files, mean_diff)
                 endif
             end if
*            Now advance to next entry
             curr_mjd = end_mjd + 1
             do j = 1, num_files
* MOD TAH 210310: Changed logic here.  Advance counter cent(j)
*               for each site used but keep the others at their
*               current value (since they have not been used)
                if( OK(j) ) then
*                   Advance to next entry
                    cent(j) = cent(j) + 1
                endif
* MOD TAH 210310: Added test for cent(j).eq.0 so that curr_mjd
*               won't get reset.
                if( cent(j).gt.eent(j) .or. cent(j).eq. 0 ) then
                    cent(j) = 0 
                else
                    if( in_mjd(cent(j)).lt.curr_mjd ) 
     .                  curr_mjd = in_mjd(cent(j))
                endif
             end do
         end do
      endif
                
****  Now write out differenced time series.
*     Generate TS file name
      ns = 1
      ts_file = tsdir(1:ln_tsdir) // '/' // in_code(ns) // '.' // 
     .          prod_id(1:ln_prod_id) // '.pos'

      open(200,file=ts_file,iostat=ierr,status='unknown')
      ref_xyz = save_ref_xyz(:,ns)
      call write_ts(200, ns)

      end

CTITLE DIFF_XYZ

      subroutine diff_xyz( cent )

*     Subroutine to difference the XYZ coordinates of entries at
*     ent given in the cent array 
 
      implicit none

      include '../includes/const_param.h'
      
      include 'tsfit.h'
      include 'tscon.h'

* PASSED
      integer*4 cent(2)    ! Entry numbers for the values at the same 
                 ! epoch

* LOCAL
      integer*4 nt, i, j, k   ! Counters
      real*8  gmjd, pos_xyz_fin(3),xyz_std(6), unc_geod(3),llu_std(3),
     .        pos_neu_fin(3), neu_std(6)
      real*8 ucov_neu(3,3), ucov_xyz(3,3) ! "Unit" NEU convarinace matrix
                   ! transformed to XYZ to scale sigma
      real*8 tcov(3,3)  ! Scratch space.
      real*8 unc_llu(3)  ! Computed Lat/long/height
      real*8 loc_coord(3) 
      real*8 rot_mat(3,3) ! Rotation from XYZ to NEU (and transpose)

      real*8 in_cov(3,3,2)   ! Covariance matrices of the two inputs
 
      integer*4 zone  ! Zone for UTM coordinates
     
      character*1 hemi    ! Hemisphere for UTM coordinates


***   For the two epochs; form the difference.
      num_ts = num_ts + 1
      nt = num_ts          ! Short version
      ts_mjd(nt)   = in_mjd(cent(1))
      ts_xyz(:,nt) = (in_xyz(:,cent(1)) - in_xyz(:,cent(2))) +
     .                save_ref_xyz(:,1)

***   Form the covariance matrix
      do j = 1,2       ! Each input
         do k = 1,3    ! Variances
             in_cov(k,k,j) = in_xyz_std(k,cent(j))**2
         end do
         in_cov(1,2,j) = in_xyz_std(1,cent(j))*in_xyz_std(2,cent(j))*
     .                   in_xyz_std(4,cent(j))            
         in_cov(1,3,j) = in_xyz_std(1,cent(j))*in_xyz_std(3,cent(j))*
     .                   in_xyz_std(5,cent(j))            
         in_cov(2,3,j) = in_xyz_std(2,cent(j))*in_xyz_std(3,cent(j))*
     .                   in_xyz_std(6,cent(j))  
         in_cov(2,1,j) = in_cov(1,2,j) 
         in_cov(3,1,j) = in_cov(1,3,j) 
         in_cov(3,2,j) = in_cov(2,3,j) 
      end do
****  Now sum to get the combined covariances
      ucov_xyz = in_cov(:,:,1) + in_cov(:,:,2)

*     Create output
      do k = 1,3
         ts_xyz_std(k,nt) = sqrt(ucov_xyz(k,k))
      end do
      ts_xyz_std(4,nt) = ucov_xyz(1,2)/
     .                  (ts_xyz_std(1,nt)*ts_xyz_std(2,nt))
      ts_xyz_std(5,nt) = ucov_xyz(1,3)/
     .                  (ts_xyz_std(1,nt)*ts_xyz_std(3,nt))
      ts_xyz_std(6,nt) = ucov_xyz(2,3)/
     .                  (ts_xyz_std(2,nt)*ts_xyz_std(3,nt))

****  Now convert XYZ to lat/long/height
      pos_xyz_fin = ts_xyz(:,nt)
      call geod_to_geod(pos_xyz_fin, unc_geod, 
     .      'XYZ', 'GEOD','WGS84','WGS84',zone,hemi)

*     Convert co-lat, long in radians to degrees
      unc_llu(1) = (pi/2-unc_geod(1))*180/pi
      unc_llu(2) = unc_geod(2)*180/pi 
      unc_llu(3) = unc_geod(3)           

****  Now get scale to convert XYZ sigmas to NEU
!     We know the correlations so compute directly.
      call XYZ_to_GEOD(rot_mat, pos_xyz_fin, loc_coord )

      call var_comp(rot_mat,ucov_xyz, ucov_neu, tcov, 3,3,1)
      do i = 1,3
         neu_std(i) = sqrt(ucov_neu(i,i))
      end do
      neu_std(4) = ucov_neu(1,2)/(neu_std(1)*neu_std(2))
      neu_std(5) = ucov_neu(1,3)/(neu_std(1)*neu_std(3))
      neu_std(6) = ucov_neu(2,3)/(neu_std(2)*neu_std(3))

****  Convert mm NEU to lat/long sigmas                
      llu_std(1) = neu_std(1)/Earth_rad*180/pi*1.d9
      llu_std(2) = neu_std(2)/Earth_rad*180/pi/
     .          sin(unc_geod(1))*1.d9
      llu_std(3) = neu_std(3)

****  Now get NEU values
 
      call loc_to_geod(unc_geod, pos_neu_fin)

*     Save values
      do j = 1,3
           ts_neu(j,nt) = pos_neu_fin(j)
           ts_llu(j,nt) = unc_llu(j)
      end do
      do j = 1,6
         ts_neu_std(j,nt) = neu_std(j)
      end do

      ts_type(nt) = ts_ref_type

****  Thats all
      return
      end

 
 
CTITLE MEAN_XYZ

      subroutine mean_xyz( cent, OK, num_files, mean_diff  )

*     Subroutine to mean the XYZ coordinates of entries at
*     ent given in the cent array.  If site is available OK(j) set
*     tru
 
      implicit none

      include '../includes/const_param.h'
      
      include 'tsfit.h'
      include 'tscon.h'

* PASSED
      integer*4 num_files,         ! Number of files to average
     .          cent(num_files)    ! Entry numbers for the values e 
                                   ! at the same epoch
      logical OK(num_files)        ! Set true if epoch OK to average
      logical mean_diff    ! Set true to take mean of differences where
                           ! reference coordinate for each file is removed

* LOCAL
      integer*4 nt, nm, i, j, k   ! Counters
      real*8  gmjd, pos_xyz_fin(3),xyz_std(6), unc_geod(3),llu_std(3),
     .        pos_neu_fin(3), neu_std(6)
      real*8 ucov_neu(3,3), ucov_xyz(3,3) ! "Unit" NEU convarinace matrix
                   ! transformed to XYZ to scale sigma
      real*8 tcov(3,3)  ! Scratch space.
      real*8 scale(3), bsol(3)    ! Scale factor for invert_vis and solution
      real*8 dxyz(3)              ! Difference of XYZ from ref_location
      real*8 unc_llu(3)  ! Computed Lat/long/height
      real*8 loc_coord(3) 
      real*8 rot_mat(3,3) ! Rotation from XYZ to NEU (and transpose)

      real*8 in_cov(3,3)   ! Covariance matrices of the inputs
 
      integer*4 zone      ! Zone for UTM coordinates
      integer*4 pivot(3)  ! Pivot elements
      character*1 hemi    ! Hemisphere for UTM coordinates


***   For the two epochs; form the difference. 
      num_ts = num_ts + 1
      nt = num_ts          ! Short version

****  Now loop over files and take XYZ weighted average
      ts_mjd(nt) = 0
      ucov_xyz = 0
      bsol = 0
      nm = 0
      do j = 1, num_files
         if( OK(j) ) then       ! Increment the normal equations
            if( ts_mjd(nt).eq.0 ) then   ! Get epoch and apriori
               ts_mjd(nt)   = in_mjd(cent(j))
            end if
*           Form the covariance matrix for this site
            do k = 1,3    ! Variances
                in_cov(k,k) = in_xyz_std(k,cent(j))**2
            end do
            in_cov(1,2) = in_xyz_std(1,cent(j))*in_xyz_std(2,cent(j))*
     .                    in_xyz_std(4,cent(j))            
            in_cov(1,3) = in_xyz_std(1,cent(j))*in_xyz_std(3,cent(j))*
     .                    in_xyz_std(5,cent(j))            
            in_cov(2,3) = in_xyz_std(2,cent(j))*in_xyz_std(3,cent(j))*
     .                    in_xyz_std(6,cent(j))  
            in_cov(2,1) = in_cov(1,2) 
            in_cov(3,1) = in_cov(1,3) 
            in_cov(3,2) = in_cov(2,3) 
              
*           Get the difference XYZ values
* MOD TAH 160225: If mean_diff then remove reference value for this site
*           otherwize use a common value
            if( mean_diff ) then 
               dxyz = in_xyz(:,cent(j)) - save_ref_xyz(:,j)
            else   ! Orginal mean option
               dxyz = in_xyz(:,cent(j)) - save_ref_xyz(:,1)
            endif 
            call invert_vis(in_cov, dxyz, scale, pivot, 3,3,1)
            ucov_xyz = ucov_xyz + in_cov
            bsol = bsol + dxyz
            nm   = nm + 1
         end if
      end do

*     Now solve the problems
      call invert_vis(ucov_xyz, bsol, scale, pivot, 3,3,1 )
      ts_xyz(:,nt) = save_ref_xyz(:,1) + bsol 
      call check_ucov(ucov_xyz, ts_mjd(nt) ) 

*     Create output
      do k = 1,3
         ts_xyz_std(k,nt) = sqrt(ucov_xyz(k,k))
      end do
      ts_xyz_std(4,nt) = ucov_xyz(1,2)/
     .               (ts_xyz_std(1,nt)*ts_xyz_std(2,nt))
      ts_xyz_std(5,nt) = ucov_xyz(1,3)/
     .               (ts_xyz_std(1,nt)*ts_xyz_std(3,nt))
      ts_xyz_std(6,nt) = ucov_xyz(2,3)/
     .               (ts_xyz_std(2,nt)*ts_xyz_std(3,nt))

****  Now convert XYZ to lat/long/height
      pos_xyz_fin = ts_xyz(:,nt)
      call geod_to_geod(pos_xyz_fin, unc_geod, 
     .      'XYZ', 'GEOD','WGS84','WGS84',zone,hemi)

*     Convert co-lat, long in radians to degrees
      unc_llu(1) = (pi/2-unc_geod(1))*180/pi
      unc_llu(2) = unc_geod(2)*180/pi 
      unc_llu(3) = unc_geod(3)           

****  Now get scale to convert XYZ sigmas to NEU
!     We know the correlations so compute directly.
      call XYZ_to_GEOD(rot_mat, pos_xyz_fin, loc_coord )

      call var_comp(rot_mat,ucov_xyz, ucov_neu, tcov, 3,3,1)
      do i = 1,3
         neu_std(i) = sqrt(ucov_neu(i,i))
      end do
      neu_std(4) = ucov_neu(1,2)/(neu_std(1)*neu_std(2))
      neu_std(5) = ucov_neu(1,3)/(neu_std(1)*neu_std(3))
      neu_std(6) = ucov_neu(2,3)/(neu_std(2)*neu_std(3))

****  Convert mm NEU to lat/long sigmas                
      llu_std(1) = neu_std(1)/Earth_rad*180/pi*1.d9
      llu_std(2) = neu_std(2)/Earth_rad*180/pi/
     .       sin(unc_geod(1))*1.d9
      llu_std(3) = neu_std(3)

****  Now get NEU values
 
      call loc_to_geod(unc_geod, pos_neu_fin)

*     Save values
      do j = 1,3
        ts_neu(j,nt) = pos_neu_fin(j)
        ts_llu(j,nt) = unc_llu(j)
      end do
      do j = 1,6
         ts_neu_std(j,nt) = neu_std(j)
      end do

* MOD TAH 180206: Encode number of values into description
*     ts_type(nt) = ts_ref_type
      write(ts_type(nt),'(A3,I2.2)') ts_ref_type(1:3), nm

****  Thats all
      return
      end

CTITLE CHECK_UCOV
   
      subroutine check_ucov(ucov, ts_mjd) 

      implicit none

      real*8 ucov(3,3), ts_mjd

      integer*4 i,j, bad

      bad = 0 
      if( abs(ucov(1,2)/sqrt(ucov(1,1)*ucov(2,2))).gt.1.d0 ) then
           bad = bad + 1
      endif
      if( abs(ucov(1,3)/sqrt(ucov(1,1)*ucov(3,3))).gt.1.d0 ) then
           bad = bad + 1
      endif
      if( abs(ucov(2,3)/sqrt(ucov(2,2)*ucov(3,3))).gt.1.d0 ) then
           bad = bad + 1
      endif

      if( bad.gt.0 ) then
         write(*,210) ts_mjd, bad
 210     format('WARNING: At MJD ',F10.4,1x,i2,' correlations > 1')
         ucov(1,2) = 0 
         ucov(1,3) = 0 
         ucov(2,3) = 0    
      endif
      end
  




