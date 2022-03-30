      subroutine decode_data(field, bak_array, values, point)
c
c     Routine to decode the data string read from the data
c     file.
c
c Include files
c -------------
*                                   ! the parameter file
      include 'plot_param.h'
c
*                                   ! the common block
      include 'plot_com.h'
      include 'pltsl_com.h'
c
c Variables
c ---------
c field  -- the information on the fields of data
c     field(1) -- data type: 0 = time, 1 = normal
c     field(2) -- start column
c     field(3) -- number of columns in field.  For time fields
c        these are interpretted as year,month,day, hour, min, and
c        seconds.  For normal fields these are interpretted as
c        value and standard deviation.
c     field(4) -- the point type
c bak_array -- the ema file of data
c values   -- the returns values from the strings.  This will be
c     julian date for time fields, for normal fields these will be
c     the values.
c point -- the type of point to be plotted and edit flag
c
      integer*4 field(1), point(2)
 
c
      real*8 values(2)
 
c
      integer*4 bak_array(*)
 
c
c Local variables
c ---------------
c unw_flag  -- the unweight flag from the bak_array
c
 
      integer*4 unw_flag
 
*  i        - Loop counter
*  ib       - Baseline number
*  is       - Site number of current site dependent quantity
*  kel(2)   - baseline sites (used in base_entnum)
*  nstep    - Number of data types per step
*  tunw     - Test version of unw flag when site dependent data is
*           - plotted (all baselines to the site are tested)
 
      integer*4 i, ib, is, kel(2), nstep, tunw
 
c
c.... Intialize the values
      values(1) = 0.d0
      values(2) = 0.d0
c
c.... Start decoding string
*                                   ! read time arguments
      if( field(1).eq.0 ) then
c
         call get_ema_8(bak_array, values(1))
c
*                                     ! read normal data
      else
c
         call get_ema_4(bak_array(imar_val), values(1), field(2))
c
*                                  ! get the sigma as well
         if( field(3).eq.2 ) then
            call get_ema_4(bak_array(imar_val), values(2), field(2)+1)
         end if
c
      end if
c
c.... Get the type of point to be plotted
      if( p_field(1).ne.2 ) then
         point(1) = bak_array(ibak_sou)
c
c....    See if downweighted
*                                   ! change point for unweight flag
         if( p_field(2).gt.0 ) then
            unw_flag = bak_array(ibak_unw+p_field(2)-1)
            point(2) = unw_flag
*                                    ! convert point so that it will be lower
            if( unw_flag.ne.0 ) then
c                                   case
               point(1) = point(1) + 32
c
            end if
c
c....       See if data were truly there.  To do this we check to see if
c           if the KalObs record number is zero.
*                                                              ! no data
            if( bak_array(ibak_recs+p_field(2)-1).eq.0 ) then
               point(1) = -1
            end if
c
         end if
 
*                                     ! Get all unw flags for this site.
         if( p_field(2).lt.0 ) then
            is = -p_field(2)
            call num_step(data_type, nstep)
*                                     ! KEL contains the baseline
            kel(1) = is
            unw_flag = 32767
            do i = 1, num_site
*                                     ! Get the baseline (except for site
               if( i.ne.is ) then
*                                     ! with itself.
                  kel(2) = i
                  call base_entnum( kel, 1, nstep, ib)
*                                                            ! Baseline there
                  if( bak_array(ibak_recs+ib-1).ne.0 ) then
                      tunw = bak_array(ibak_unw+ib-1)
                      if( abs(tunw).lt.abs(unw_flag) ) unw_flag = tunw
                  end if
              end if
            end do
 
*                                          ! No observations at this epoch
            if( unw_flag.eq.32767 ) then
                unw_flag =  0
*                                ! Will cause point not to be plotted
                point(1) = -1
            end if
 
            point(2) = unw_flag
            if( unw_flag.ne.0 ) point(1) = point(1) + 32
         end if
c
*                  ! set point symbol to match unweight
      else
c
*                                   ! set point
         if( p_field(2).ne.0 ) then
            point(1) = bak_array(ibak_unw+p_field(2)-1)
            point(2) = 0
         end if
c
      end if
c
      return
      end
 
c.......................................................................
 
      subroutine get_ema_8(bak_array, value)
c
c     Routine to return real*8 value for bak_array
c
c Variables
c ---------
c bak_array -- the record from bak_array (only the first element is
c     real*8
c value -- the real*8 value to be returned
c
      real*8 bak_array(1), value
 
c
c
c.... Just save the value
      value = bak_array(1)
c
      return
      end
 
c.......................................................................
 
      subroutine get_ema_4(bak_array, value, field)
c
c     Routine to return real*4 values from bak array
c
c Variables
c ---------
c bak_array -- the bak array (may be only part of it)
c value  -- the value read from the bak_array
c field   -- the field information
c
      real*4 bak_array(1)
 
c
      real*8 value
 
c
      integer*4 field
 
c
c
c.... Get the value
      value = bak_array(field)
c
      return
      end
 
c.....................................................................
 
      subroutine read_file(ema_data, header_only)
c
c     Routine to read data from file into ema.  It will also set
c     up the dynamic mapping of ema for this data set.  The
c     site_names and souc_names are read from the header record
c     along with the codes for the markov processes which are included
c     in the file.
c
c     This routine will also change the units of the elements in
c     in bakfile to match the output units.
c
c Include files
c -------------
      include '../includes/kalman_param.h'
*                                 ! the parameter file
      include 'plot_param.h'
c
*                                 ! the common block
      include 'plot_com.h'
      include 'pltsl_com.h'
 
c Variables
c ---------
c ema_data -- an integer array which will used to store the contents
c     of the file and the data to be plotted.  The data to be plotted
c     will be extracted by read_data.
c
      integer*4 ema_data(max_plt_space)

* header_only -- Tells program to read only the header of the backfile
*                and not the data

      logical header_only
 
c
c
c Local variables
c ---------------
c ierr -- general error variable
c
*   dcb_buffers - The number of dcb buffers available for reading
*               - the BAKFILE.
      integer*4 ierr, dcb_buffers
 
c
c
c Variables used in reading file
c ------------------------------
c idcb  -- the dcb buffer used to read the file
c data_array -- the array of data from the file
c idata_array -- an integer array equivalenced to data_array
c cdata_array -- a character array equivalenced to data_array
c
c headr_len -- the length of the headr site/souc and markov codes
* num_headr -- Number of header records.
c
      character*8 cdata_array(1)
 
*   char_name       - Name for KalObs file
*   Fullname        - Full name of backfile

      character*128 char_name, fullname
c
      integer*4 idcb(dcb_size), idata_array(max_headr_size),
     .          num_headr
 
*   len_data_file   - Length of the KalObs file name
*   len_rec         - Length of record read
*   num_name_recs   - Number of records needed for the KalObs file
*                   - name
*   trimlen         - HP function for length of string

      integer*4 len_data_file, len_rec, num_name_recs, trimlen
c
      integer*4 headr_len, temprecl, lenr, i, j
 
c
      real*4 data_array(1)
 
c
      equivalence (idata_array,data_array)
      equivalence (cdata_array,data_array)
      equivalence (char_name  ,data_array)
c
c Some decoding information
c -------------------------
c islcat -- an interger array equivalanced to slcat_tem
c slcat_tem -- the charcater array for the slcat key read from BCK_FILE
c
      integer*4 islcat(2), FmpClose
 
c
      character*4 slcat_tem
 
c
      equivalence (slcat_tem,islcat)
c
c.... Set file read false
      file_read = .false.
c
c.... open the input file (if one has been specified)
*                                        ! no file name
      if( input_file(1:1).eq.' ' ) then
         write(termlu,'(" READ_FILE Error: No data file has",
     .     " been given")')
*                              ! open and read file
      else
*                                    ! convert name to upper case
* MOD TAH 891228: Get the record length and reopen file as type
*     2 with the correct record length.

         call FmpOpen( idcb, ierr, input_file, 'ro', 1)
         call readd(idcb, ierr, idata_array, 128, lenr, 1)
         temprecl = idata_array(7)
         call fullfilename(input_file, 2, 1, temprecl, fullname)
         ierr = FmpClose(idcb)
         call FmpOpen(idcb,ierr,fullname, 'RO',1)

         call report_error('FMPOPEN',ierr,'open',input_file,
     .                      0, 'READ_FILE')
c
*                                ! no use continuing
         if( ierr.lt.0 ) return
c
c....    Get the headr
         call readd(idcb,ierr,data_array,temprecl,len_rec,1)
c
         num_headr  = idata_array(1)
         num_epochs = idata_array(2)
         num_mar    = idata_array(3)
         num_site   = idata_array(4)
         num_souc   = idata_array(5)
         data_type  = idata_array(6)
         bak_recl   = idata_array(7)
         len_data_file  = idata_array(8)
c
c....    Get the run_time for solvk
         do i = 1,6
            run_time(i) = idata_array(9+i)
         end do
c
*****    Get the KalObs name from file (if it is there)
         num_name_recs = (len_data_file-1)/(4*bak_recl) + 1
         len_rec = min(33,bak_recl)
         do i = 1, num_name_recs
              call readd(idcb,ierr,idata_array((i-1)*bak_recl+1),
     .                   len_rec, lenr, 0)
         end do
 
         data_file = char_name(1:len_data_file)
 
c....    Get the rest of the header information (site names, souces
c        names and markov element iformation)
         headr_len = min(max_headr_size,bak_recl)    
c
         do i = 1,num_headr
            call readd(idcb,ierr,idata_array((i-1)*bak_recl+1),
     .         headr_len, lenr,0)
         end do
c
c....    Now get site and source names
         do i = 1, num_site
            site_names(i) = cdata_array(i)
         end do
c
         do i = 1, num_souc
            souc_names(i) = cdata_array(num_site+i)
         end do
c
c....    Get the markov types
         do i = 1,num_mar
            mar_elements(i) = idata_array(2*(num_site+num_souc)+i)
         end do
c
c....    Now list to the user what is going on
         call sum_file
c
****     If we are reading the data read now
         if( .not. header_only ) then

c....         Now set up the mapping for ema
              call plt_ema_map
c
c....         Now read the file data
              call read_markov(idcb, ema_data(ibak_array))
c
c....         Set the file read flag
              file_read = .true.
         end if
c
c....    Setup some basic header records
         write(headers(1),100) data_file(1:trimlen(data_file)),
     .                        (run_time(i), i = 1,5 )
 100     format(' Data file ',a,' Runtime: ',i4,'/',i2,'/',i2,1x,
     .          i2,':',i2)
          write(headers(2),120) (site_names(j),j=1,num_site)
 120      format(' SITES: ',16(:a8,1x))
 
          actual_num_headers = 2
 
*              ! file name available
      end if
c
c.... Thats all
      ierr = FmpClose( idcb )
      call report_error('FmpClose',ierr,'clos',input_file,0,'PLTSL')
 
      return
      end
 
c........................................................................
 
      subroutine read_markov(idcb,bak_array)
c
c     This routine will read the markov elements from the  bckfile
c
c Include files
c -------------
*                                 ! the parameter file
      include 'plot_param.h'
c
*                                 ! the plot_com common block
      include 'plot_com.h'
      include 'pltsl_com.h'
c
 
*         max_bak_recl - Largest bak recl length expected
 
      integer*4 max_bak_recl
 
      parameter ( max_bak_recl = 4*max_mar + 3*(max_site+1)*max_site)
 
*
*   scr_buffer(max_bak_recl)    - A buffer for reading the BakFile
*               - records into. (Needed since VREAD seems to require
*               - two MSEG pages)
 
      integer*4 scr_buffer(max_bak_recl)
 
c Variables
c ---------
c idcb  -- the dcb buffer
c bak_array -- the ema bak array storage
c len_rec -- Length of record read
 
      integer*4 idcb(1), bak_array(bak_recl,1), len_rec
 
c
c
c Local variables
c ---------------
c num -- the site/source number for the markov parameter
c code -- the type of parameter
c
*   i,j     - Loop counters

      integer*4 num, code, i,j, ierr
c
      do i = 1, num_epochs
c
         call readd(idcb,ierr,bak_array(1,i), bak_recl, len_rec,0 )
         call report_error('FmpRead',ierr,'read','BakFile',0,
     .                     'Read_markov')
      end do
c
c.... Now do the conversion of the units in bak_array entries
c     Loop over the markov elements first
      do i = 1, num_mar
c
c....    Decode the type of markov element
         call decode_mar(mar_elements(i), num, code)
c
c....    Now convert units. (Note that we need an even record length
c        in bak_file for this to work)
         if( conv_mar(code).ne.1.) then
c
c....       Convert the values
            call psmy(conv_mar(code), bak_array(imar_val+2*i-2,1),
     .         bak_recl, bak_array(imar_val+2*i-2,1),bak_recl,
     .         num_epochs)
c
c....       Convert the sigmas
            call psmy(conv_mar(code), bak_array(imar_val+2*i-1,1),
     .         bak_recl, bak_array(imar_val+2*i-1,1),bak_recl,
     .         num_epochs)
c
*               ! conversion was not one
         end if
*               ! looping over elevation angles
      end do
c
c.... Now convert elevation angles to deg
      do i = 1, num_site
         call psmy(57.295780,bak_array(ibelev+(i-1),1), bak_recl,
     .      bak_array(ibelev+(i-1),1), bak_recl, num_epochs)
      end do
c
      return
      end
 
c........................................................................
 
      subroutine decode_mar(element, num, code)
c
c     This routine will decode the markov element code which is stored
c     in packed form.  The site/source number is stored in the right-
c     hand 8 bits, and the code is in the lefthand 8 bits.  (See &plkbd
c     for listing of meaning of the codes)
c
c Variables
c ---------
c element -- the markov element in packed form
c num -- the site/source number
c code -- the markov element code number
c
      integer*4 element, num, code
 
c
c.... Depack the code, get the site/source number
      num = element/256
c
c.... Get the code
      code = element - num*256
c
c.... Thats all
      return
      end
 
c......................................................................
 
      subroutine plt_ema_map
c
c     This routine will step up the plotting package mapping of the
c     ema area.  The first (bak_recl*num_epochs) is used for the bak_file,
c     the remaining area is used for the x, y, and pt data
c
c Include files
c -------------
*                                   ! the plot parameter file
      include 'plot_param.h'
c
*                                   ! the plot control common block
      include 'plot_com.h'
      include 'pltsl_com.h'
c
c All the information we need for this routine is stored in plot_com
c
c
c.... Start setting up the address space
      ibak_array = 1
c
      if( .not.xdiff .and. .not.ydiff ) then
*                                                        ! add bak_array space
      ix_array   = ibak_array + bak_recl*num_epochs
c
*                                                      ! 2 real*4 elements for
      iy_array   = ix_array   + 2*num_epochs
c                                                       each epoch
c
*                                                      ! same amount as y_array
      ipt_array   = iy_array   + 2*num_epochs
      end if
c
c.... Make sure we have not gone outside ema space
*                                                         ! not enough room
      if( ipt_array+2*num_epochs.gt.max_plt_space ) then
         write(termlu,100) ipt_array+num_epochs, max_plt_space
  100    format(/" Out of ema space, space need is ",i4," words",/,
     .      " Only ",i4," words of ema space available.",
     .      " Increase max_plt_space in &plkpa (the parameter file)")
c
         stop ' Out of ema space in pltsl'
      end if
c
c.... Now set up the mapping for bak_array
c
      call bak_file_map(ibepoch, imar_val, ipost_res, ibelev, ibak_recs,
     .   ibak_sou, ibak_unw, num_mar, data_type, num_site, bak_recl)
c
c.... Thats all
      return
      end
 
c......................................................................
 
      subroutine sum_file
c
c     This routine summarizes the contents of the bak_file.  It lists
c     the sites sources and the markov elements which are available
c
c Include files
c -------------
*                                  ! the parameter file
      include 'plot_param.h'
c
*                                  ! the common file
      include 'plot_com.h'
      include 'pltsl_com.h'
c
c All the information we need is in the common file
c
c Functions
c ---------
c trimlen -- HP string length utility
c kbit  -- SOLVE bit checking function
c
      integer*4 trimlen
 
c
      logical kbit
 
c
c.... Start by giving the summary of the solvk solution
 
*                                             ! No file yet so do not
      if( trimlen(input_file).eq.0 ) RETURN
*                                             ! summaries
 
      write(termlu,100) input_file(1:trimlen(input_file))
 100  format(/" Summary of data in ",a," ", $)
c
      call out_run_time(termlu,' SOLVK ',run_time)
c
c.... Give SLCAT key and data types used
      write(termlu,150) data_file(1:trimlen(data_file))
 150  format(" KalObs file ",a,". ",$)
c
      if( kbit(data_type,1) ) then
         write(termlu,'("Group Delay "$)')
      end if
c
      if( kbit(data_type,2) ) then
         write(termlu,'("Phase Delay "$)')
      end if
 
      if( kbit(data_type,3) ) then
         write(termlu,'("SB Delay "$)')
      end if
 
      if( kbit(data_type,4) ) then
         write(termlu,'("Phase delay rate "$)')
      end if
c
      write(termlu,'("data types used")')
c
c.... Output more about volumn of data
      write(termlu,200) num_epochs, num_site, num_souc
 200  format(" There are ",i4," epochs of data using ",i3," sites",
     .   " and ",i3," sources")
c
c.... Output the sites
      call out_names(termlu,'sites', site_names, num_site)
c
c.... Output the soucres
      call out_names(termlu,'sources ', souc_names, num_souc)
c
c.... Output the makov elements which are avaliable
      call out_markov(termlu)
c
c.... Thats all
      return
      end
 
c........................................................................
 
      subroutine out_markov(iout)
c
c     Routine to output the markov element descriptions
c
c Include files
c -------------
*                                 ! the parameter file
      include 'plot_param.h'
c
*                                 ! the common block
      include 'plot_com.h'
      include 'pltsl_com.h'
c
c Variables
c ---------
c iout -- the output lu number
c
      integer*4 iout
 
c
c Local variables
c ---------------
c num -- the site/souce number for the markov element
c code -- the markov element code (see &plkbd)
c
      integer*4 num, code, i
 
c
c.... Loop over the markov codes and output descripion
      write(iout,100) num_mar
 100  format(" There are ",i3," markov elements in the plot file:")
c
      do i = 1, num_mar
c
c....    decode the markov element code
         call decode_mar(mar_elements(i), num, code)
c
c....    If num is not zero then either a site/source goes with the
c        the markov element
*                              ! see if site or source
         if( num.ne.0 ) then
c
*                                                  ! not a soucre
            if( code.ne.10 .and. code.ne.11 ) then
               write(iout,200) plot_types(code), site_names(num)
*                                                  ! then a source
            else
               write(iout,200) plot_types(code), souc_names(num)
            end if
c
*                              ! no site/source name for this code
         else
            write(iout,200) plot_types(code)
         end if
      end do
c
c.... The output format
 200  format(1x,a,1x,a)
c
c.... Thats all
      return
      end
 
c......................................................................
 
      subroutine get_field(field, axis_label)
c
c     Subroutine to decode the field information for the plot
c
c Include files
c -------------
*                                 ! the parameter file
      include 'plot_param.h'
c
*                                 ! the common block
      include 'plot_com.h'
      include 'pltsl_com.h'
c
c Variables
c ---------
c field  -- the field array:
c     field(1) = type of data (0 is time, 1 is normal)
c     field(2) = the index in bak_array to find the quanitity to be
c        plotted
c     field(3) = number of colmumns in time field [at the moment it
c        is hardwired to 5 (yr,mth,day, ihr,min)]
c axis_label -- the default label to be put on the axis
c
      character*(*) axis_label
 
c
      integer*4 field(1)
 
c
c local variables
c ---------------
c iel -- an index returned from get_cmd which is the position in
c     the list of commands
c jel -- a second index used to get site/source
c kel -- a double entry index for baselines
c test_mar -- the mar_elemnt code (use to check against the available
c     markov elements)
c nstep -- the number of data types in this solution
c
      integer*4 iel, jel, test_mar, kel(2), nstep, i
 
c
c Functions
c ---------
c trimlen -- HP string length utility
c kbit    -- SOLVE bit checking routine
c duse    -- True if a delay or rate data type used
 
      integer*4 trimlen, ierr
 
c
      logical kbit, duse
 
c
c.... Clear the axis label
      axis_label = ' '
c
c.... Set the field type to -1 so that we can later trap non valid fields
      field(1) = -1
c
c.... Start decoding; see what type of plot
      call get_cmd(buffer(9:), plot_types, num_types, iel)
      call command_error(buffer(9:), iel)
c
c.... See if valid plot type
*                         ! valid type
      if( iel.gt.0 ) then
c
c....    Start constructing the axis label
         axis_label = '"' // plot_types(iel)
c
         jel = 0
         kel(1) = 1
         kel(2) = 1
c
c....    Check some specific fields; see if site dependent markov
*                                             ! get the site
         if( iel.ge.1 .and. iel.le.10 ) then
            call get_cmd(buffer(17:), site_names, num_site, jel)
            call command_error(buffer(17:), jel)
c
c....       Add site name to label
            if( jel.gt.0 ) then
               axis_label(trimlen(axis_label)+2:) = site_names(jel)
            end if
c
         end if
c
c....    see if source dependent
*                                              ! get the source
         if( iel.ge.11 .and. iel.le.12 ) then
            call get_cmd(buffer(17:), souc_names, num_souc, jel)
            call command_error(buffer(17:), jel)
c
c....       Add source name to lable
            if( jel.gt.0 ) then
               axis_label(trimlen(axis_label)+2:) = souc_names(jel)
            end if
c
         end if
c
c....    check tides (these may be site dependent or independent)
*                                              ! get the site
         if( (iel.ge.18 .and. iel.le.28) .or.
     .        iel.eq.47 ) then
            call get_cmd(buffer(17:) , site_names, num_site, jel)
            call command_error(buffer(17:), jel)
c
c....       Add site names to axis label
            if( jel.gt. 0 ) then
               axis_label(trimlen(axis_label)+2:) = site_names(jel)
            end if
c
*                                ! string was empty (no site given,
            if( jel.lt.0 ) then
c                                  so assume global
               jel = 0
            end if
         end if
c
c....    see if elevation angle
*                               ! get site
         if( iel.eq.49 ) then
            call get_cmd(buffer(17:), site_names, num_site, jel)
            call command_error(buffer(17:), jel)
c
c....       Add site name to label
            if( jel.gt.0 ) then
               axis_label(trimlen(axis_label)+2:) = site_names(jel)
            end if
         end if
c
c....    See if residual
*                                            ! get pair of sites
         if( iel.ge.50 .and. iel.le.51) then
            call get_cmd(buffer(17:), site_names, num_site, kel(1))
            call command_error(buffer(17:), kel(1))
            call get_cmd(buffer(25:), site_names, num_site, kel(2))
            call command_error(buffer(25:), kel(2))
c
c...        Add baseline to axis_label
            if( kel(1).gt.0 .and. kel(2).gt.0 ) then
               axis_label(trimlen(axis_label)+2:) = site_names(kel(1))
     .            // '-' // site_names(kel(2))
            end if
c
         end if
c
c
c....    Before check whether we can plot the selected quantities, add
c        the units to the label
c
         if( iel.ne.0 ) then
            axis_label(trimlen(axis_label)+2:) =
     .         unit_label(iel) // '"'
         else
            axis_label(trimlen(axis_label)+1:) = '"'
         end if
c
c.....   See if we found all the needed information, check markov
*                              ! these are markov
         if( iel.le.47 ) then
            test_mar = jel*256 + iel
c
c....       See if we have this markov element
            field(2) = -1
            do i = 1, num_mar
*                                                       ! yes we have
               if( test_mar.eq. mar_elements(i) ) then
*                                          ! allow space for sigmas
                  field(2) = 2*(i-1) + 1
               end if
            end do
c
c....       See if we found a match
*                                    ! no we didnot
            if( field(2).lt.1 ) then
               write(termlu,100) buffer(9:trimlen(buffer))
  100          format(" Failed to find markov field ",a)
c
*                                    ! set up rest of field
            else
c
*                                  ! normal data
               field(1) = 1
*                                  ! get sigma as well
               field(3) = 2
*                                  ! Save the site number for getting the
               p_field(2) = -jel
*                                  ! edit flag.
c
            end if
            return
c
*                              ! this is not a markov field
         else
c
c....       See if time
*                                 ! yes it is
            if( iel.eq.48 ) then
c
c....          Clear the axis label
               axis_label = ' '
c
*                                 ! set time field
               field(1) = 0
*                                 ! set of field
               field(2) = 1
*                                 ! five elements in field for readin scales
               field(3) = 5
               p_field(2) = 0
c
               return
            end if
c
c....       See if elevation angle
*                                 ! yes it is
            if( iel.eq.49 ) then
c
c....          Check we found site
*                                  ! yes we found site
               if( jel.gt.0 ) then
*                                  ! normal data
                  field(1) = 1
c                           | # r*4s past | |# r*4s in |
c                           | markov start| |elev field|
                  field(2) = ibelev-2 +    jel - 1
                  field(3) = 0
                  p_field(2) = -jel
*                                  ! site not found
               else
                  write(termlu,200) buffer(17:trimlen(buffer))
 200              format(" Could not find site ",a)
               end if
c
               return
            end if
 
c....       See if postfit residual
*                                                 ! postfit residual
            if( iel.ge.50 .and. iel.le.51 ) then
c
c....          see if we have baseline
*                                                       ! yes we do
               if( kel(1).gt.0 .and. kel(2).gt.0 ) then
c
c....             Compute the baseline number
                  call num_step(data_type, nstep)
*                                       ! get index if data avaiable
                  if( iel.eq.50 ) then
*                                                   ! get index
                     if( duse(data_type,1) ) then
                        call base_entnum(kel, 1, nstep, jel)
                        field(1) = 1
                     else
                       write(termlu,'(1x,a," data type not",
     .                 " avaiable")') plot_types(iel)
                     end if
                  end if
c
*                                        ! get index  if data available
                  if( iel.eq.51 ) then
*                                                   ! get index
                     if( kbit(data_type,4) ) then
                        call base_entnum(kel, 2, nstep, jel)
                        field(1) = 1
                     else
                       write(termlu,'(1x,a," data type not",
     .                 " available",i4)') plot_types(iel)
                     end if
                  end if
c
c....             Now compute the entry
                  field(2) = ipost_res-2 + 2*(jel-1) 
                  field(3) = 2
                  p_field(2) = jel
c
*                                    ! no match on sites
               else
c
                  write(termlu,300) buffer(17:trimlen(buffer))
  300             format(" Could not find baseline ",a)
c
               end if
c
               return
            end if
c
*                           ! either markov or not plot types
         end if
c
*                           ! we did not find plot command
      else
c
         call command_error(buffer(9:),ierr)
c
      end if
c
      return
      end
 
c......................................................................
 
      subroutine psmy(scalar, v1, inc1, v2, inc2, num)
c
c     Subroutine to multiple v1 by scalar and save in v2
c
c Variables
c ---------
c scalar -- real*4 multiplication factor
c v1  -- ema vector to be mulitplied
c inc1 -- increment for stepping through v1
c v2  -- ema vector for results
c inc2 -- increment for stepping through v2
c num -- the number of times we should loop
c
      integer*4 inc1, inc2, num, i
 
c
      real*4 scalar, v1(inc1,1), v2(inc2,1)
c
c.... Loop over the vector
      do i = 1, num
         v2(1,i) = v1(1,i)*scalar
      end do
c
      return
      end
 
c......................................................................
 
      subroutine help
c
c
c     This routine will list available commands, and give help information
c     from the help_file on specific commands
c
c Include files
c -------------
*                         ! the pltsl parameter file
      include 'plot_param.h'
c
*                         ! the pltsl common block
      include 'plot_com.h'
      include 'pltsl_com.h'
c
c Variables
c ---------
c iel -- the command number from get_cmd.
c ierr -- the error return from reading files
c eof  -- end of file indicator
c found -- set true when command found
c finished -- set true when when help has finished output from file
c line -- the line read from the help file
c i    -- loop counter
c len  -- length of 'line' or 1 if line is null
c
      integer*4 iel, ierr, i, len
c
      logical eof, found, finished
c
      character*79 line

* FullHelpName - Full name of help file

      character*128 FullHelpName
 
c
c Functions
c ---------
c trimlen -- HP string length utility
c
      integer*4 trimlen
 
c
c.... See if help is requested for a specific command
      call get_cmd(buffer(9:), commands, num_commands, iel)
c
c.... If iel is greater then 0 then give specific help information
      if( iel.gt.0 ) then
c
          call gen_help_name(pltsl_help_file, fullhelpname)
          open(300,file=fullhelpname,iostat=ierr, status='old')
          if( ierr.ne.0 ) then
              call report_error('IOSTAT',ierr,'open',fullhelpname,
     .            0,'HELP')
c
c....         Set iel = 0 so that standard help massage will be output
              iel = 0
c
*                 ! file openened OK, find help message
          else
c
              eof = .false.
              finished = .false.
              found = .false.
c
c....         Read file until we reach EOF or we have finished outputting
c             the help meassge
              do while ( .not.eof .and. .not.finished )
c
                  read(300,'(a)',iostat=ierr, err=100) line
  100             continue
*                                       ! EOF
                  if( ierr.eq.-1 ) then
                      eof = .true.
*                                       ! see if other error
                  else
                      if( ierr.ne.0 ) then
                          call report_error('IOSTAT',ierr,'read',
     .                        help_file,0,'HELP')
*                                   ! force standard massage
                          iel = 0
                      end if
                  end if
c
c....             See if help command found
*                                       ! check
                  if( ierr.eq.0 ) then
*                                           ! see if command line
                      if( .not.found ) then
                          if( line(3:10).eq.commands(iel) ) then
*                                              ! we found the command
                              found = .true.
                          end if
*                                          ! see if next command found ( end
                      else
*                                          ! of help for current command)
*                                                       ! end
                          if( line(1:2).eq.'**' ) then
                              finished = .true.
*                                          ! write help line
                          else
                              len = max(1,trimlen(line))
                              write(termlu,'(a)') line(1:len)
                         end if
                      end if
*                                  ! no file reading error
                   end if
*                                  ! looping over file
              end do
              close (300)
*                                  ! file opened OK
          end if
*                                  ! iel > 0 (may have been set no to zero
      end if
*                                  ! if file error
c
c.... List the availabke commands to the termlu if iel<=0
      if( iel.le.0 ) then
          write(termlu,300)
 300      format(/" HELP for PLTSL:",/,
     .            " The runstring for PLTSL is:",//,
     .            " CI> PLTSL,<control file name or lu>, <plot lu>,",
     .            " <file to be plotted>",//,
     .            " All runstring parameters are optional")
 
          write(termlu,310)
 310      format(/" All commands must be preceeded by at least one",
     .            " blank character.",/,
     .            " Commands may be abbreviated")
 
          write(termlu,'(" The commands available in plot are:")')
          write(termlu,320) (commands(i),i=1,num_commands)
 320      format(5(1x,a,3x))
 
          write(termlu,'(" Some field types which are available:")')
          write(termlu,320) (plot_types(i),i=1,num_types)
c
          write(termlu,330)
 330      format(" Further help may be obtained using",/,
     .           "?  HELP <command>"/)
c
          call sum_file
c
      end if
 
c
      return
      end
 
c......................................................................
 
CTITLE PLOT_FINAL_SETUP
      subroutine plot_final_setup
 
 
*     Routine to do the final setup for the plot programs
 
      include 'plot_param.h'
      include 'plot_com.h'
 
c.... Tell the user we are here
      call report(
     .     ' PLTSL: (preceed all commands by at least one blank)')
      call report('        (use HELP for instructions, END to quit)')
 
      return
      end
 
