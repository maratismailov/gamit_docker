 
      program velcom

      implicit none 
 
*     Program to compare and difference from a reference velocity
*     field a number of velocity fields extracted from GLOBK output
*     of from a getrel run.
 
*     The runstring of the program is:
*     % velcom <ext> <average> <ref> [Input fields .... ]
*     where <ext> is the next extent to the added to the file names
*                 to get the name of the file for the differences.
*           <average> is the name of average field and the rms of the
*                 input fields
*           <ref> if the name of the reference field to be removed from
*                 all fields
*           [Input fields .... ] is the list of input fields.
 
      include 'velcom.h'
 
* LOCAL VARIABLES
 
*   i       - Loop counter for looping over fields
 
 
      integer*4 i
 
***** First read the runstring string and get all of the files.  If not
*     enough information is given then stop program and print help file.
 
      write(*,100)
 100  format(/' VELCOM: Program to compare velocity fields',/)
 
      call get_velcom_runstring
 
***** Now loop over all of the files and read the inputs.  The first
*     field is the reference field.
 
      do i = 1, num_files
 
          call read_vel_fields( i )
 
      end do
 
****  Now loop over all the fields and compute the residuals, output
*     the residual field and compute the average field
 
      call clear_av_field
      do i = 2, num_files
          call gen_res_field(i)
      end do
 
****  Now finish up the statistics on the average field and output this
*     field
 
      call gen_av_field
 
****  Thats all
      end
 
CTITLE  GET_VELCOM_RUNSTRING
 
      subroutine get_velcom_runstring

      implicit none 
 
*     This routine will read the runstring of the program.  If the
*     runstring is not complete then the help file is listed.
 
      include 'velcom.h'
 
* LOCAL VARIABLES
 
*   i       - Loop counter
*   rcpar   - Read runstring routine
*   len_run - Length of runstring
 
      integer*4 i, rcpar, len_run
 
****  Try to get the extent for the difference files
 
      len_run = rcpar(1, extent )
*                             ! Print Help and stop
      if( len_run.le.0 ) then
          call proper_runstring('velcom.hlp', 'velcom/extent', 1)
      end if
 
*     Try to get name of file for average field
      len_run = rcpar(2, average_file)
*                             ! Print Help and stop
      if( len_run.le.0 ) then
          call proper_runstring('velcom.hlp', 'velcom/average file',
     .                          1)
      end if
 
*     Now start looping over the input files
 
      i = 0
      write(*,100)
 100  format(' Input files ',/,
     .      ' ------------'  )
      do while ( i.lt.max_files .and. len_run.gt.0 )
          i = i + 1
          len_run = rcpar(2+i, file_names(i) )
          if( len_run.gt.0 ) then
              write(*,110) i, file_names(i)(1:len_run)
 110          format(I3,'. ',a)
          end if
      end do
 
      num_files = i - 1
 
*     See if we have enough data
      if( num_files.lt.2 ) then
          write(*,150)
 150      format(/' VELCOM: Not enough input files in runstring')
          call proper_runstring('velcom.hlp', 'velcom/files', 1)
      end if
 
****  Thats all
      return
      end
 
CTITLE READ_VEL_FIELDS
 
      subroutine read_vel_fields( nf )

      implicit none 
 
*     Routine to read velocity field file.  The first is the reference
*     field and the name of the sites are taken from this file.  Sites
*     the other files will only be considered if they are not in the
*     reference file.
 
      include 'velcom.h'
 
* PASSED VARIBALES
 
*       nf      - File number to be read.  If this value is 1
*               - then site list is made, otherwize the site
*               - number is determined.
 
      integer*4 nf
 
* LOCAL VARIABLES
 
*   iel         - Site number determined from get_cmd
*   ierr, jerr  - IOSTAT errors on file read and string decode
*   terr        - total error in decoding read_line.
*   indx        - Position in string
*   trimlen     - Length of string (used portion)
*   nt          - Temporary (short) name for the number of
*               - site entries for a particular file.
 
      integer*4 iel, ierr, jerr, terr, indx, trimlen, nt
 
*   tlong, tlat - Read latitude and longitude
*   tEvel, tNvel    - read East and North velocity (mm/yr)
*   tEadj, tNadj    - Read East and North adjustment (not used)
*   tEsig, tNsig    - read East and North sigma (mm/yr)
*   tNErho      - Read East North correlation
*   tUvel, tUadj, tUsig - Read Up velocity, adjustent (not
*               - used) and sigma
 
      real*4 tlong, tlat, tEvel, tNvel, tEadj, tNadj, tEsig, tNsig,
     .    tNErho, tUvel, tUadj, tUsig
 
*   tsite_name  - Read site name
*   cd          - Dummy charcater string for read_line
 
      character*8 tsite_name, cd
 
*   line        - Line read from file.
 
      character*256 line
 
****  Initialize the number of sites of this is first file
 
      if( nf.eq.1 ) then
          tot_sites = 0
      end if
 
*     Initialize number of sites in file
      num_sites(nf) = 0
 
****  Try to open file.
 
      open(100, file=file_names(nf), iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',file_names(nf), 0,
     .                'read_vel_fields')
 
 
*     Skip this file if it is not the first. If reference field then
*     kill job
      if( ierr.ne.0 .and. nf.gt.1 ) RETURN
      if( ierr.ne.0 .and. nf.eq.1 ) then
          stop ' VELCOM: Could not open reference velocity field'
      end if
 
*     Start reading the velocity field
      read(100,'(a)', iostat=ierr) header(nf)
      if ( trimlen(header(nf)).eq.0 ) then
          header(nf) = ' UNKOWN SOLUTION TYPE: File ' //
     .                file_names(nf)
      else if( header(nf)(1:1).eq.' ' ) then
          header(nf) = ' No Header: File '  //
     .                file_names(nf)
          rewind(100)
      end if
 
*     Now loop over the file.
      do while ( ierr.eq.0 )
          read(100,'(a)', iostat=ierr ) line
 
*         Process if no error and the first character is blank
          if( ierr.eq.0  .and. line(1:1).eq.' ' .and.
     .        trimlen(line).gt.0                    ) then
 
*             Try to read long and lat first.  If problem then skip
*             line
              indx = 0
              terr = 0
*             Now each of the values from the line
              call read_line(line, indx, 'R4', jerr, tlong, cd)
              terr = terr + jerr
              call read_line(line, indx, 'R4', jerr, tlat, cd)
              terr = terr + jerr
              call read_line(line, indx, 'R4', jerr, tEvel, cd)
              terr = terr + jerr
              call read_line(line, indx, 'R4', jerr, tNvel, cd)
              terr = terr + jerr
              call read_line(line, indx, 'R4', jerr, tEadj, cd)
              terr = terr + jerr
              call read_line(line, indx, 'R4', jerr, tNadj, cd)
              terr = terr + jerr
              call read_line(line, indx, 'R4', jerr, tEsig, cd)
              terr = terr + jerr
              call read_line(line, indx, 'R4', jerr, tNsig, cd)
              terr = terr + jerr
              call read_line(line, indx, 'R4', jerr, tNErho, cd)
              terr = terr + jerr
              call read_line(line, indx, 'R4', jerr, tUvel, cd)
              terr = terr + jerr
              call read_line(line, indx, 'R4', jerr, tUadj, cd)
              terr = terr + jerr
              call read_line(line, indx, 'R4', jerr, tUsig, cd)
              terr = terr + jerr
*              Get the site name (tUadj is a dummy here: will not be
*             changed by call.
              call read_line(line, indx, 'CH', jerr, tUadj, tsite_name)
              call casefold( tsite_name )
              terr = terr + jerr
 
*             Now see if any error.  If there are error assume that
*             this is comment line, and ignore.
              if ( terr.eq.0 ) then
 
*                 Get the site number.  If this file 1 then simply
*                 add site to list
                  if( nf.eq.1 ) then
                      tot_sites = tot_sites + 1
                      site_names(tot_sites) = tsite_name
*                     Save the site position
                      long(tot_sites) = tlong
                      lat(tot_sites)  = tlat
                      iel = tot_sites
                  else
                      indx = 1
                      call get_cmd(tsite_name, site_names,
     .                                num_sites, iel, indx )
                  end if
 
*                 Now save the values if the site name was found
                  if( iel.gt.0 ) then
                      num_sites(nf) = num_sites(nf) + 1
                      nt = num_sites(nf)
                      site(nt,nf)  = iel
                      Nvel(nt,nf)  = tNvel
                      Evel(nt,nf)  = tEvel
                      Nsig(nt,nf)  = tNsig
                      Esig(nt,nf)  = tEsig
                      NErho(nt,nf) = tNErho
                      Uvel(nt,nf)  = tUvel
                      Usig(nt,nf)  = tUsig
*                                 ! All values have been saved
                  end if
*                                 ! There was no error decoding line
              end if
*                                 ! No error reading file
          end if
*                                 ! Looping over the input file
      end do
 
***** Now just summarize what has happened.
      write(*,300) num_sites(nf), nf,
     .            file_names(nf)(1:trimlen(file_names(nf)))
 300  format(' VELCOM: ',I4,' Sites found in file ',i3,1x,a)
 
      close ( 100 )
****  Thats all
      return
      end
 
CTITLE CLEAR_AV_FIELD
 
      subroutine clear_av_field

      implicit none 
 
*     Routine to clear the statistics variables for the average field
*     computation.
 
      include 'velcom.h'
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
 
      integer*4 i,j
 
*     Loop over all summation values and clear
      do i = 1, tot_sites
          do j = 1,8
              sum_stat(j,i) = 0.d0
          end do
      end do
 
***** Thats all
      return
      end
 
CTITLE GEN_RES_FIELD
 
      subroutine gen_res_field (nf)

      implicit none 
 
*     Routine to compute and output the differences in the velocity
*     fields between file NF and the reference field (NF=1)
*     The routine also accumulated the statistics of the difference
*     which will be output in the average_file.
 
      include 'velcom.h'
 
* PASSED VARIABLES
 
*   nf      - File number being processed.
 
      integer*4 nf
 
* LOCAL VARIABLES
 
*   iel         - Site number determined from get_cmd
*   ierr, jerr  - IOSTAT errors on file write.
*   trimlen     - Length of string (used portion)
*   nt          - Temporary (short) name for the number of
*               - site entries for a particular file.
 
      integer*4  ierr, trimlen, nt, i
 
*   rEres, rNres    - Residual East and North velocity (mm/yr)
*   rEtot, rNtot    - total East and North adjustment (not used)
*   rEsig, rNsig    - residual East and North sigma (mm/yr) computed
*               - as the difference between the reference and
*               - current field.  (If less than 0.1 value set
*               - to 0.1)
*   rNErho      - residual East North correlation
*   rUres, rUtot, rUsig - residual Up velocity,total and sigma
*   covnf(2,2), covrf(2,2)  - Covariance matrices for the
*               - Current file and reference file
 
 
      real*4 rEres, rNres, rEtot, rNtot, rEsig, rNsig, rNErho,
     .    rUres, rUtot, rUsig, covnf(2,2), covrf(2,2)
 
*   outfile     - Name of the output file
 
 
      character*256 outfile
 
*     First see if we have any data
      if( num_sites(nf).eq.0 .or. nf.eq.1 ) RETURN
 
*     Generate the output file name (remove old extent and add new one)
      call gen_res_name ( file_names(nf), extent, outfile )
 
      open(200, file=outfile, iostat=ierr, status='unknown')
      write(*,120) nf, outfile(1:trimlen(outfile))
 120  format(' Generating residuals to file ',i3,'. Output ',a)
 
*     Write the headers for the output
      write(200,150) file_names(nf)(1:trimlen(file_names(nf))),
     .               file_names(1)(1:trimlen(file_names(1))),
     .              header(nf)(1:trimlen(header(nf)))
 150  format('* VELCOM: Residuals for file ',a,' from ',a,/,a,/,
     .    '*  Long.  ',7x,'Lat.',8x,'E & N Res',6x,'E & N Total',6x,
     .    'E & N +-',4x,'RHO',7x,'H res',1x,'H Total',2x,'+-', 2x,
     .    'Site',/,'*  (deg)',7x,
     .    ' (deg)  ',7x,'(mm/yr)',8x,'(mm/yr)',8x,'(mm/yr)',20x,
     .    '(mm/yr)')
 
****  Now loop over of then entries
 
      do i = 1, num_sites(nf)
 
*         Get the site number in the reference field
          nt = site(i,nf)
 
*         Now compute the residual values and save the totals
          rNres = Nvel(i,nf) - Nvel(nt,1)
          rEres = Evel(i,nf) - Evel(nt,1)
          rUres = Uvel(i,nf) - Uvel(nt,1)
          rNtot = Nvel(i,nf)
          rEtot = Evel(i,nf)
          rUtot = Uvel(i,nf)
 
*         Now compute the covariance matrix of the difference
          covnf(1,1) = Nsig(i,nf)**2
          covnf(2,2) = Esig(i,nf)**2
          covnf(1,2) = NErho(i,nf)*Nsig(i,nf)*Esig(i,nf)
 
*         Now for the reference field
          covrf(1,1) = Nsig(nt,1)**2
          covrf(2,2) = Esig(nt,1)**2
          covrf(1,2) = NErho(nt,1)*Nsig(nt,1)*Esig(nt,1)
 
*****     Now difference the values.  (Don't worry about the sign since
*         we don't know which field is the most accurate)
* Current code seems to add the sigma.  Changed height to do this.
          covnf(1,1) = abs(covnf(1,1)+covrf(1,1))
          covnf(2,2) = abs(covnf(2,2)+covrf(2,2))
          covnf(1,2) =     covnf(1,2)+covrf(1,2)
 
*         Now compute the sigmas
* MOD TAH 920904: Lowered limit of small values for new format
          rNsig = sqrt(covnf(1,1))
          if( rNsig.lt.0.005 ) rNsig = 0.005
          rEsig = sqrt(covnf(2,2))
          if( rEsig.lt.0.005 ) rEsig = 0.005
          rNErho = covnf(1,2)/(rNsig*rEsig)
          rUsig  = sqrt(abs(Usig(i,nf)+Usig(nt,1)))
          if( rUsig.lt.0.005 ) rUsig = 0.005
 
****      Now write out the writes
          write(200,200) long(nt), lat(nt), rEres, rNres,
     .                    rEtot, rNtot, rEsig, rNsig, rNErho,
     .                    rUres, rUtot, rUsig, site_names(nt)
* MOD TAH 930211 Changed format to add one more digit.
 200      format(2(1x,f10.5),1x,6(1x,f7.2),1x,f6.3,2x,
     .           3(1x,f6.2), 1x,a8,a1)
 
****      Now accumulate the statistics for the average velocity field
          call accum_stat( sum_stat(1,nt),  rNres, rEres, rUres)
*             ! Looping over sites in this file
      end do
 
****  Thats all
      close ( 200 )
      return
      end
 
CTITLE ACCUM_STAT
 
      subroutine accum_stat( sum_stat, rNres, rEres, rUres)

      implicit none 
 
*     Routine to accumulate the statistics for the the average solution.
 
* PASSED VARIABLES
 
 
*   rNres, rEres, rUres       ! Differences in N, E and U
      real*4  rNres, rEres, rUres
 
 
*   sum_stat(8) - Accumulation variable for statistcs (see
*               - velcom.h for contents)
 
      real*8 sum_stat(8)
 
****  Now sum the values into the array
      sum_stat(1) = sum_stat(1) + 1
      sum_stat(2) = sum_stat(2) + rNres
      sum_stat(3) = sum_stat(3) + rNres**2
      sum_stat(4) = sum_stat(4) + rEres
      sum_stat(5) = sum_stat(5) + rEres**2
      sum_stat(6) = sum_stat(6) + rNres*rEres
      sum_stat(7) = sum_stat(7) + rUres
      sum_stat(8) = sum_stat(8) + rUres**2
 
***** Thats all
      return
      end
 
CTITLE GEN_RES_NAME
 
      subroutine gen_res_name ( file_name, extent, outfile )

      implicit none 
 
*     Routine to generate the name of the file for the residual
*     field.  Any characters after the last '.' in the name are
*     replaced with the 'extent' string.  If there is not . then
*     the extent is just added to the name
 
* PASSED VARIABLES
 
*   file_name   - Original file name
*   extent      - New extent to be used
*   outfile     - Output file name
 
      character*(*) file_name, extent, outfile
 
* LOCAL VARIABLES
 
*   trimlen     - Length of string
*   i           - Loop counter
 
      integer*4 trimlen, i
 
****  Get the length of the file names
      i = trimlen(file_name)
      do while (i.gt.2 .and. file_name(i:i).ne.'.' .and.
     .                       file_name(i:1).ne.'/' )
          i = i - 1
      end do
 
****  Now see what we found
      if( file_name(i:i).eq.'.' ) then
          if( extent(1:1).eq.'.' ) then
              outfile = file_name(:i) // extent(2:)
          else
              outfile = file_name(:i) // extent
          end if
*                 ! Just concatinate on the end
      else
          outfile = file_name(:trimlen(file_name)) // extent
      end if
 
***** Thats all
      return
      end
 
CTITLE GEN_AV_FIELD
 
      subroutine gen_av_field

      implicit none 
 
*     Routine to generate the average velocity field with error bars
*     given by the rms correlation between the residuals to the
*     reference field.
 
      include 'velcom.h'
 
*   ierr    - IOSTAT error
*   trimlen     - Length of string (used portion)
*   i       - Loop counter
*   nt      - Number of values in determination of average
 
      integer*4 ierr, trimlen, i, nt
 
*   aEres, aNres    - Average Residual East and North velocity  
*   aEtot, aNtot    ! Average total East and North adjustment
*               - from the rms of the values
*   aNErho      - average residual East North correlation
*   aUres, aUtot, aUsig - Average Up velocity,total and sigma
 
 
      real*4 aEres, aNres,aEtot, aNtot, aEsig, aNsig,  aNErho, 
     .       aUres, aUtot, aUsig
 
****  First open the average field and write the header for it.
      open(200, file= average_file , iostat=ierr, status='unknown')
      call report_error('IOSTAT',ierr,'open',average_file, 1,
     .                'gen_av_field' )
 
****  Write out the header
      write(200,100) file_names(1)(1:trimlen(file_names(1))),
     .              header(1)(1:trimlen(header(1)))
 100  format('* VELCOM: Average field relative to file ',a,/,
     .    '*',a,/,
     .    '*  Long. ',7x,'Lat.',7x,'E & N Avg',6x,'E & N Res.',3x,
     .    'E & N +-',5x,'RHO',7x,'H Avg',3x,'H Res',2x,'+-', 2x,
     .    'Site',3x,'#',/,'*  (deg)',5x,
     .    ' (deg)  ',7x,'(mm/yr)',6x,'(mm/yr)',3x,'(mm/yr)',12x,
     .    '(mm/yr)')
 
****  Now start looping over the sites
 
      do i = 1, tot_sites
 
*         Compute the means and rms values if we have enough values
          if( sum_stat(1,i).ge.2 ) then
              nt = sum_stat(1,i)
              aNres = sum_stat(2,i)/nt
              aNsig = sqrt(abs((sum_stat(3,i) - aNres**2*nt)/(nt-1)))
              if( aNsig.lt.0.2 ) aNsig = 0.2
              aNtot = Nvel(i,1) + aNres
              aEres = sum_stat(4,i)/nt
              aEsig = sqrt(abs((sum_stat(5,i) - aEres**2*nt)/(nt-1)))
              if( aEsig.lt.0.2 ) aEsig = 0.2
              aEtot = Evel(i,1) + aEres
              aUres = sum_stat(7,i)/nt
              aUsig = sqrt(abs((sum_stat(8,i) - aUres**2*nt)/(nt-1)))
              if( aUsig.lt.0.2 ) aUsig = 0.2
              aUtot = Uvel(i,1) + aUres
*             Compute NE correlation
              aNErho = (sum_stat(6,i) - aNres*aEres*nt)/(nt-1)/
     .                                  (aNsig*aEsig)
 
*             Now write out the results
              write(200,200) long(i), lat(i), aEtot, aNtot,
     .                    aEres, aNres, aEsig, aNsig, aNErho,
     .                    aUres, aUtot, aUsig, site_names(i),
     .                    int(sum_stat(1,i))
c200          format(1x,f7.3,1x,f7.3,1x,6(f7.1,1x),f6.3,1x,3(f7.1,1x),
c    .                1x,a,1x,I3)
* MOD TAH 021005: Increase format
 200          format(2(1x,f10.5),1x,6(1x,f7.2),1x,f6.3,2x,
     .               3(1x,f6.2), 1x,a8,a1)
          end if
      end do
 
***** Thats all
      return
      end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
