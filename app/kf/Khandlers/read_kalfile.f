ctitle
 
      subroutine read_kalfile(kal_array, phy_rec, nrecs, kalfile_open)
 

      implicit none
 
c
c     Routine to actually read the kalfil.  The routine will
c     check to see is the kalfil is open.  If the number of records to be
c     read (nrecs) is positive, the file will be read forward starting at
c     at record phy_rec. If nrecs is negative, the file be read backwards
c     starting at phy_rec.
 
*                                  ! the solvk parameter file
      include '../includes/kalman_param.h'
c
*                                  ! the forsl control common file
      include '../includes/forsl_comm.h'
*                                  ! the solvk control common block
      include '../includes/solvk_cntl.h'
*                                  ! the kalfil common block desription
      include '../includes/kalfile_def.h'
c
c Variables
c ---------
c kal_array -- the array in ema in to which the kalfile is read
c phy_rec -- the current physical record to be read first
c nrecs  -- the number of records to be read
c kalfile_open -- indicates if kalfile is open
c
      integer*4 kal_array(krec_len, max_krec)
 
c
c
      integer*4 phy_rec, nrecs
 
c
 
      logical kalfile_open
 
c
c Local variables
c ---------------
c ikaldc -- a copy of the dcb buffer for kalfil
c kal_name -- the name of kalfile
c cart -- the cartridge for kalfile
c nbuf -- the number of 128 word buffers available for reading kalfil
c ema_rec -- the record number in ema.  Always starts from 1 when
c     a new read is carried out
c st_rec  -- the starting record of the file to be read
c en_rec  -- the ending record number for the read
c
c ierr  -- error flag
c
 
      integer*4 ikaldc(144), nbuf, ema_rec, st_rec,
     .    en_rec, ierr
 
*   i,j,k   - Loop counters
*   len     - length of the kalfile name
 
      integer*4 i, len

*   FullKalname - Full name of the kalfile, with the file type
*     and size of record added.

      character*132 fullKalname
 
c
c ..................
c Scratch common area
c ..................
c scr_data -- a scratch common area used in rading kalfil
c ident -- an identifier saying that get_epoch has used the scratch
c     read.
c
      common ident, scr_data
 
c
      integer*4 ident, scr_data(scr_forsl)
 
c
c
c.... see if kalfil is open
*                                    ! open kalfil
      if( .not. kalfile_open ) then
c
c
c....    compute number of buffers available
         nbuf = (scr_forsl-17)/128

*        Generate the full Kalfile name
         call FullFileName( kal_file, 2, 1, krec_len, FullKalName)
*        Forces FmpOpen to use next unit number.
         scr_data(1) = 0
         call FmpOpen(scr_data, ierr, FullKalName,'rwo', nbuf)
c
c....    Check error
         call report_error('FMGR',ierr,'open',FullKalName,1,
     .      'READ_KALFIL')
c
c....    set open indicator
         kalfile_open = .true.
c
c....    Move the dcb into local memory
         call wmove(scr_data, ikaldc, 144)
c
      end if
c
c.... Move the dcb back to the scratch common area
      call wmove(ikaldc, scr_data, 144)
      ident = 3
c
c.... Now read the data from kalfil
      ema_rec = 0
c
c.... See if we are reading forwards or backwards
*                              ! read forward
      if( nrecs.gt.0 ) then
 
         st_rec = phy_rec
         en_rec = phy_rec + nrecs - 1
*                              ! read backwards
      else
         st_rec = phy_rec + nrecs + 1
         en_rec = phy_rec
      end if
c
c
*                            ! read records (Note: file is physically read
c                              forward.
      do i = st_rec, en_rec
c
         ema_rec = ema_rec + 1
         call readd(scr_data,ierr,kal_array(1,ema_rec),krec_len, len,i)
c
      end do
c
c.... save the dcb
      call wmove(scr_data, ikaldc, 144)
c
c.... Move the first ema record into main memory
*                             ! move first record
      if( nrecs.gt.0 ) then
         call wmove(kal_array(1,1), ksite, krec_len)
*                             ! move last record
      else
         call wmove(kal_array(1,ema_rec), ksite, krec_len)
      end if
c
c.... Thats all
      return
      end
 
c......................................................................
