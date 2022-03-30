CTITLE    ..................................................................
 
      subroutine label_ax0(idate,indx, lf, ntic, idep)
c
c     Routine to write time axis label
c
c Variables
c ---------
c idate -- the date to be written
c indx -- the type of incrementing we are doing
c space_cntrl -- gives the x and y coordinates for a back space
c     and line feed
c ntic -- the tic number
c idep -- the depth of the label of the axis
c
      integer*4 idate(6), indx, ntic, idep
 
c
*      lf(2)    - Indicates how to move in x and y to get a line feed
      real*4 lf(2)
 
c
c Local variables
c ---------------
c label -- the label to be be written
c months -- lists the months of the years
c
      character*15 label
 
*            lower_label    - The second label to be written on the
*                           - axis.  This string is checked for each
*                           - tic and only output if it changes
      character*15 lower_label
 
c
      character*3 month(12)
 
c
      data month / 'Jan','Feb','Mar','Apr','May','Jun','Jul',
     .             'Aug','Sep','Oct','Nov','Dec' /
c
c.... Write into the label
      label = ' '
*                            ! just output year
      if( indx.eq.1 ) then
*                   ! only one line of labels
         idep = 1
         call swh( 4, 1)
 
         write(label,'(i4)') idate(1)
         call save_label(label)
         return
      end if
c
*                             ! output month and year
      if( indx.eq.2 ) then
*                   ! two lines of labels
         idep = 2
         call swh( 4, 2 )
 
         write(label,'(a)') month(idate(2))
c
         call save_label(label)
c
c....    Add year if first tic of Jan of a year
         label = ' '
         write(label,'(i4)') idate(1)
c
         if( label.ne.lower_label .or. ntic.eq.1 ) then
            call sr2mv(lf(1), lf(2) )
            call save_label(label)
            lower_label = label
         end if
c
      end if
c
c.... Day labels
*                            ! output day, month year
      if( indx.eq.3 ) then
*                            ! two lines of labels
         idep = 2
         call swh( 3, 2)
         write(label,'(i2)') idate(3)
 
         call save_label(label)
c
c....    add month and year
         label = ' '
         write(label,'(a,", ",i4)') month(idate(2)),idate(1)
 
         if( lower_label.ne.label .or. ntic.eq.1 ) then
 
*                                          ! Save full width
            if( ntic.ne.1 ) call swh( 9,2)
c
            call sr2mv(lf(1), lf(2) )
            call save_label(label)
            lower_label = label
 
         end if
c
      end if
c
c.... hour min labels
*                            ! output hour min, day, month year
      if( indx.eq.4 ) then
*                            ! two lines of labels
         idep = 2
*                            ! Set with for just hours and minites part.
         call swh( 5, 2)
c
c MOD JLD To print leading '0' in minutes
         write(label,'(i2,":",i2.2)') idate(4), idate(5)
C        call write_label(label)
         call save_label(label)
c
c....    add day month and year
         label = ' '
c
         write(label,'(i2,1x,a,", ",i4)') idate(3),
     .       month(idate(2)),idate(1)
 
         if( lower_label.ne.label .or. ntic.eq.1 ) then
 
*                                            ! If large label inside region
            if( ntic.ne.1 ) call swh( 12, 2)
*                                            ! axis then set full width
c
c.....      Move to correct position
C           call jr2mv(2.*space_cntrl(1,1),2.*space_cntrl(1,2)) ! back 2 spac
C           call jr2mv(space_cntrl(2,1),space_cntrl(2,2))       ! line feed
 
            call sr2mv(lf(1), lf(2) )
            call save_label(label)
            lower_label = label
         end if
c
      end if
c
c.... hour min labels (spaced by minutes)
*                            ! output hour min, day, month year
      if( indx.eq.5 ) then
*                            ! two lines of labels
         idep = 2
         call swh( 5, 2)
c
c MOD JLD To print leading '0' in minutes
         write(label,'(i2,":",i2.2)') idate(4), idate(5)
C        call write_label(label)
         call save_label(label)
c
c....    add day month and year
         label = ' '
c
         write(label,'(i2,1x,a,", ",i4)') idate(3),
     .         month(idate(2)),idate(1)
 
         if( lower_label.ne.label .or. ntic.eq.1 ) then
 
*                                           ! Save full width
            if( ntic.ne.1 ) call swh(12, 2)
c
c.....      Move to correct position
            call sr2mv(lf(1), lf(2) )
            call save_label(label)
            lower_label = label
         end if
c
      end if

c.... hour min sec labels (spaced by sec)
*                            ! output hour min, day, month year
      if( indx.eq.6 ) then
*                            ! two lines of labels
         idep = 2
         call swh( 5, 2)
c
c MOD JLD To print leading '0' in minutes
         write(label,'(i2,":",i2.2,":",I2.2)') idate(4), idate(5), 
     .            idate(6)
C        call write_label(label)
         call save_label(label)
c
c....    add day month and year
         label = ' '
c
         write(label,'(i2,1x,a,", ",i4)') idate(3),
     .         month(idate(2)),idate(1)
 
         if( lower_label.ne.label .or. ntic.eq.1 ) then
 
*                                           ! Save full width
            if( ntic.ne.1 ) call swh(12, 2)
c
c.....      Move to correct position
            call sr2mv(lf(1), lf(2) )
            call save_label(label)
            lower_label = label
         end if
c
      end if

      return
      end
 
