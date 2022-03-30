 
      program genrelref

      implicit none 

*     Program to generate a relative veloctity field by subracting
*     the velocity of the first site in the list from all others in
*     the list
 
*     The runstring of the program is:
*     % genrelref <input> <output>
*     where <input> is the input file and
*          <output> is output file name.  (If first character is
*                  '.' then extent will be changed
 
* LOCAL VARIABLES
 
*   i       - Loop counter
*   rcpar   - Read runstring routine
*   len_run - Length of runstring
*   iel         - Site number determined from get_cmd
*   ierr, jerr  - IOSTAT errors on file read and string decode
*   terr        - total error in decoding read_line.
*   indx        - Position in string
*   trimlen     - Length of string (used portion)
*   nt          - Temporary (short) name for the number of
*               - site entries for a particular file.
 
      integer*4  rcpar, len_run, ierr, jerr, terr, indx,
     .    trimlen
 
*   tlong, tlat - Read latitude and longitude
*   tEvel, tNvel    - read East and North velocity (mm/yr)
*   tEadj, tNadj    - Read East and North adjustment (not used)
*   tEsig, tNsig    - read East and North sigma (mm/yr)
*   tNErho      - Read East North correlation
*   tUvel, tUadj, tUsig - Read Up velocity, adjustent (not
*               - used) and sigma
*   REvel, RNvel    - Reference E and North velocity
*   RUvel       - Reference Up velocity
 
      real*4 tlong, tlat, tEvel, tNvel, tEadj, tNadj, tEsig, tNsig,
     .    tNErho, tUvel, tUadj, tUsig, REvel, RNvel, RUvel
 
*   ref_found   - Indicates that reference site has been found.
 
      logical ref_found
 
*   tsite_name  - Read site name
*   Rsite_name  - Name of reference site name
*   cd          - Dummy charcater string for read_line
 
      character*8 tsite_name, Rsite_name, cd
 
*   line        - Line read from file.
*   input       - Name of input file
*   output      - Name of output file
*   outentry        - Output name entered (may be extent only)
*   header      - Header line of input file
 
      character*256 line, input, output, outentry, header

      rEvel = 0.d0
      rNvel = 0.d0
      rUvel = 0.d0
 
****  Try to get the extent for the difference files
 
      len_run = rcpar(1, input )
*                             ! Print Help and stop
      if( len_run.le.0 ) then
          call proper_runstring('genrelref.hlp', 'genrelref/input',1)
      end if
 
*     Try to get name of file for outentry field
      len_run = rcpar(2, outentry)
*                             ! Print Help and stop
      if( len_run.le.0 ) then
          call proper_runstring('genrelref.hlp',
     .            'genrelref/outentry',  1)
      end if
 
 
****  Try to open file.
 
      open(100, file=input, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',input, 1,
     .                'genrelref/input')
 
*     Generate name of output file
      if( outentry(1:1).eq.'.' ) then
          call gen_res_name ( input, outentry, output  )
      else
          output  = outentry
      end if
 
      open(200, file=output, iostat=ierr, status='unknown' )
      call report_error('IOSTAT',ierr,'open',output, 1,
     .                'genrelref/output')
 
*     Start reading the velocity field
      read(100,'(a)', iostat=ierr) header
      if ( trimlen(header).eq.0 ) then
          header = ' UNKOWN SOLUTION TYPE: File ' //
     .                input
      end if
 
      ref_found = .false.
 
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
C              terr = terr + jerr
              if( jerr.ne.0 ) tsite_name = ' '
 
*             Now see if any error.  If there are error assume that
*             this is comment line, and ignore.
              if ( terr.eq.0 ) then
 
*                 We have found valid line.  See if this is reference
*                 site
                  if( .not. ref_found ) then
                      rsite_name = tsite_name
                      rEvel = tEvel
                      rNvel = tNvel
                      rUvel = tUvel
 
****                  Write the headers for the output file
                      write(200,200) rsite_name,
     .                    header(1:trimlen(header))
 200                  format('* GENRELREF: Reference relative',
     .                    ' velocity field for ',a8,/,a)
                      write(200,210)
 210                  format('   Long.       Lat. ',7x,'E & N Rate ',3x,
     .                       ' E & N Adj. ',2x,' E & N +-',1x,
     .                    ' RHO ',5x,' H Rate  H adj.   +-',1x,'SITE',/,
     .                   2x,' (deg)      (deg) ',2x,3(6x,'(mm/yr)'),15x,
     .                    '(mm/yr)' )
                      ref_found = .true.
                  end if
 
*                 Now write out the values.
                  write(200,220) tlong, tlat,
     .                    (tEvel-rEvel), (tNvel-rNvel),
     .                    tEadj, tNadj, tEsig, tNsig, tNErho,
     .                    (tUvel-rUvel), tUadj, tUsig, tsite_name
* MOD TAH 930215 Changed format to add one more digit.
 220              format(2(1x,f10.5),1x,6(1x,f7.2),1x,f6.3,2x,
     .                   3(1x,f7.2), 1x,a8,a1)
*
 
 
*                                 ! There was no error decoding line
              end if
*                                 ! No error reading file
          end if
*                                 ! Looping over the input file
      end do
 
****  Thats all
      close(100)
      close(200)
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
     .                      file_name(i:1).ne.'/' )
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
 
