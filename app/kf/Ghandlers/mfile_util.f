 
CTITLE OPEN_MF
 
      subroutine open_mf( mf_unit, mf_name, ierr)

      implicit none 
 
*     Routine to open the mfile and return the open status
*
 
* PASSED VARIABLES
 
*   mf_unit - Unit number for the mfile.
*   ierr    - IOSTAT Error
 
      integer*4 mf_unit, ierr
 
*   mf_name - Name of the mfile
 
 
      character*(*) mf_name
 
****  Try to open the mfile
      open( mf_unit, file=mf_name, iostat=ierr, status = 'old',
     .     access='sequential', form='unformatted')
      call report_error('IOSTAT',ierr,'open',mf_name,0,'open_mf')
 
****  Thats all
      return
      end
 
CTITLE READ_MF
 
      subroutine read_mf( mf_unit, ierr )

      implicit none 
 
*     This routine will read the whole (single session) mfile and
*     save the results in the mf common block
 
      include '../includes/kalman_param.h'
      include '../includes/mfile_def.h'
 
* PASSED VARIABLES
 
*   mf_unit - Unit number for mfile
*   ierr    - IOSTAT error
 
      integer*4 mf_unit, ierr
 
* LOCAL VARIABLES
 
*   i   - Loop counter
*   jerr    - Dimension checking error flag.
*   cf      - Cfile number counter
 
*   iflag   - 1,2 or 3 depending on record
*   nwds    - Number of words in record
*   nversn  - Version number for mfile

      integer*4 i, jerr, cf, iflag, nwds, nversn
 
*   message     - String for warning messages.
 
      character*128 message
 
***** Start by read the first "parameter" record
 
      read(mf_unit,iostat=ierr)
     .    iflag, nwds, nversn,
     .    mf_ndy, mf_nepch, mf_mtpart,
     .    (mf_alabel(i), i = 1, min(mf_mtpart,mf_maxprm)),
     .    (mf_idms(i), i = 1, min(mf_mtpart,mf_maxprm)),
     .    (mf_islot_all(i), i = 1, min(mf_mtpart,mf_maxprm)),
     .    (mf_aprval(i), i = 1, min(mf_mtpart,mf_maxprm)),
     .    (mf_adjust(i), i = 1, min(mf_mtpart,mf_maxprm)),
     .    mf_nsat, (mf_isat(i), i = 1, min(mf_nsat,mf_maxsat)),
     .    mf_nsite, (mf_sitet(i), i = 1, min(mf_nsite, mf_maxcfl)),
     .    mf_nrfile, (mf_rfname(i), i = 1, min(mf_nrfile,mf_maxcfl)),
     .    mf_ntfile, (mf_tfname(i), i = 1, min(mf_ntfile,mf_maxtfl)),
     .    mf_norb

* MOD TAH 190605: Save the mfile version number
      mf_nversn = nversn

* MOD TAH 050211: Casefold xfile names to be consistent with older m-files
      do i = 1, min(mf_nrfile,mf_maxcfl)
         call casefold(mf_rfname(i))
      end do
 
      call report_error('IOSTAT',ierr,'read','Block 1 of mfile',
     .                0,'read_mfile')
 
*     Check that the flags and versions are OK.  Also check the
*     maximum values
 
      if( iflag.ne.1 ) then
          write(*,120) iflag
 120      format('Incorrect FLAG read from mfile: Value read ',i8)
          ierr = -1
      end if
 
*     Now check the version
      if( nversn.lt. 921 ) then
          write(message,130) nversn/100.d0
 130      format('MFile version ',f5.2,' too old.  Cannot read')
          call report_stat('warning','autcln','read_mf','Mfile',
     .                    message,0)
          ierr = -1
          RETURN
      end if
 
***** Check the number of sessions
      if( mf_ndy.gt.1 ) then
*         Too many sessions
          write(message,140) mf_ndy
 140      format('MFile: has ',i3,' sessions; only single session ',
     .            'supported')
          call report_stat('warning','autcln','read_mf','Mfile',
     .                    message,0)
          ierr = -1
          RETURN
      end if
 
 
***** Now check sizes:
      call check_msize( mf_mtpart, mf_maxprm, 'Number of parameters',
     .                jerr )
      ierr = min(ierr,jerr)
      call check_msize( mf_nsat, mf_maxcfl, 'Number of sites',
     .                jerr )
      ierr = min(ierr,jerr)
      call check_msize( mf_nsat, mf_maxcfl, 'Number of satellies',
     .                jerr )
      ierr = min(ierr,jerr)
 
***** Check to see if errors
      if( ierr.ne.0 ) call report_stat('FATAL','AUTCLN','read_mfile',
     .                'M-file','Dimensions too large ',ierr)     
 
****  Now continue reading by getting the session information.
 
      read(mf_unit,iostat=ierr) iflag, nwds, mf_idy,
     .    (mf_mdy(i), i = 1,3),
     .    (mf_hms(i), i = 1,3),
     .    mf_nepch_sess,
     .    mf_inter,
     .    mf_skip,
     .    mf_ncfile, (mf_cfname(i),i=1,min(mf_ncfile,mf_maxcfl)),
     .    (mf_stawgh(i), i = 1, min(mf_ncfile,mf_maxcfl)),
     .    mf_nsat_sess,
     .    (mf_satwgh(i), i = 1, min(mf_nsat_sess,mf_maxsat))
 
      call report_error('IOSTAT',ierr,'read','Block 2 of mfile',
     .    0, 'read_mf')

      if( iflag.ne.2 ) then
          write(*,*) 'Flag not correct in mfile record 2 ',
     .               iflag,cf
          ierr = -1
      end if
 
      if( ierr.ne.0 ) call report_stat('FATAL','AUTCLN','read_mfile',
     .                'M-file','Format not correct',ierr)     
 
***** Now loop over the sites reading the individual entries
      do cf = 1, mf_ncfile

*        Check the version of the mfile
         if( nversn.lt.940 ) then       
            read(mf_unit,iostat=ierr) iflag, nwds,
     .          (mf_mdy_cf(i), i = 1,3),
     .          (mf_hms_cf(i), i = 1,3),
     .           mf_kpart,
     .          (mf_islot_cf(i), i = 1,min(mf_kpart, mf_maxprm_cf)),
     .           mf_elvcut(cf), mf_zenmod,
     .           mf_numzen,(mf_idtzen(i), 
     .                        i = 1,min(mf_numzen,mf_maxatm))
*           Set the values which are used by m-file version 980
            mf_gradmod = 'CON'
            mf_numgrad = 1
            mf_idtgrad(1) = 1   
         else
            read(mf_unit,iostat=ierr) iflag, nwds,
     .          (mf_mdy_cf(i), i = 1,3),
     .          (mf_hms_cf(i), i = 1,3),
     .           mf_kpart,
     .          (mf_islot_cf(i), i = 1,min(mf_kpart, mf_maxprm_cf)),
     .           mf_elvcut(cf), mf_zenmod,
     .           mf_numzen,(mf_idtzen(i), 
     .                        i = 1,min(mf_numzen,mf_maxatm)),
     .           mf_gradmod, mf_numgrad, 
     .          (mf_idtgrad(i), i = 1,min(mf_numgrad,mf_maxgrad))
          endif
            
          if( iflag.ne.3 ) then
             write(*,*) 'Flag not correct in mfile record 3 ',
     .                  iflag,cf
             ierr = -1
             RETURN
          end if
 
          call report_error('IOSTAT',ierr,'read','Block 3 of mfile',
     .        0, 'read_mf')
      end do
      call check_msize( mf_numzen, mf_maxatm, 
     .                 'Number of zeniths delays', jerr )
      ierr = min(ierr, jerr)
      
      call check_msize( mf_numgrad, mf_maxgrad, 
     .                 'Number of gradients', jerr )
      ierr = min(ierr, jerr)
 
      if( ierr.ne.0 ) call report_stat('FATAL','AUTCLN','read_mfile',
     .                'M-file','Dimensions too large ',ierr)     
****  Thats all
      return
      end
 
CTITLE CHECK_MSIZE

      subroutine check_msize( mf_num, mf_max, type, jerr )

      implicit none 

*     Routine to check that we do not violate size limits

* PASSED VARIABLES

*  mf_num  - size in mfile
*  mf_max  - Max allowed
*  jerr    - error return if dimension is too small

      integer*4 mf_num, mf_max, jerr

*  type    - Type of parameter checked

      character*(*) type

* LOCAL VARIABLES

*  trimlen - length of string

      integer*4 trimlen

****  See if OK
      jerr = 0
      if( mf_num.gt.mf_max ) then
          write(*,120) type(1:trimlen(type)), mf_num, mf_max
 120      format('MFILE Dimension too large for ',a,' File has ',
     .           i6,' Max allowed ',i6)
          jerr = -1
      end if

      return
      end

