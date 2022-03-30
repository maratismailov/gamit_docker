CTITLE DECODE_PRT_OPT

      subroutine decode_prt_opt( line, indx, opt )

      implicit none 

*     Routine to decode the print options for globk and glorg

* num_opt  - Number of options supported

      integer*4 num_opt
      parameter ( num_opt = 32 )

* PASSED VARIABLES
*  line  - Line read from command file or runstring
*  indx  - Current position in line
*  opt   - Bit mapped options values

      integer*4 indx, opt
      character*(*) line

* LOCAL VARIABLES

* trimlen  - Length of string
* is       - OPtion number    
* ierr     - IOSTAT error.
* indx_save - Saved value of index

      integer*4 trimlen, is, ierr, indx_save

* opt_names(num_opt)  - OPtions for output

      character*4 opt_names(num_opt) 

      character*8 cdum
* MOD TAH 030109: Added two options:
* RNRP  -- Rename Report printed in org/prt file, and equate lines
*          reported in <org root>.eqs
* FIXA  -- Fixes different apriori values when equates are made 
*          except for positions that differ by more than 1.0 m
* MOD TAH 030223: Added NAPP option.  This option only applies
*          the rotations and does not cnage covariance matrix.
* MOD TAH 030730: Introduced GEOD (geodetic coordinates, WGS84) and
*          UTM (WGS84) outputs
* MOD TAH 070823: Introduced MIDP which will set the globk output
*          epoch to the midpoint of the data and not the end.
*          (-M option in glsave will do the same).
* MOD TAH 110223: Added NEUC (North/east/up correlation matrix
*          output (bit 31 used).
* MOD TAH 130310: Forced PBOP option on as standard (needs prt_opt/org_opt call)
* MOD TAH 210509: Added BALN option to Balance numbers of stations when
*          multiple files are combined on same day (used for IGS
*          combination mostly).

      data opt_names / 'CORR', 'BLEN', 'BRAT', 'CMDS',
     .                 'VSUM', 'NUS6', 'NUS7', 'NUS8',
     .                 'NUS9', 'NU10', 'RAFX', 'MOUT',
     .                 'COVA', 'PSUM', 'GDLF', 'DBUG',
     .                 'ERAS', 'NOPR', 'SDET', 'RNRP',
     .                 'FIXA', 'PLST', 'NAPP', 'RPRN',
     .                 'GEOD', 'UTM ', 'SMAR', 'PBOP',
     .                 'MIDP', 'RENM', 'NEUC', 'BALN' /

*     Remove delimiters that may be present from runstring 
*     passing
      call casefold(line)
      call sub_char(line,':',' ')
      call sub_char(line,'=',' ')

****  First try to read an integer value
      opt = 0
      indx_save = indx 
      call read_line(line, indx, 'I4', ierr, opt, cdum ) 
      if( ierr.ne.0 ) then
          opt = 0
          indx = indx_save
      end if

      call sbit(opt,28,1)   ! Force PBO output on MOD TAH 130310.

*     See if there options.
      do while ( indx.lt.trimlen(line))
         call get_cmd(line, opt_names, num_opt, is, indx)
* MOD TAH 030610: If if the RPRN option set (bit 24), change it
*        to bit 20 
         if( is.eq.24 ) is = 20
         if( is.gt.0 ) call sbit(opt, is, 1)
      end do

****  Thats all
      return
      end

