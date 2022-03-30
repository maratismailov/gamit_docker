CTITLE ASSIGN_CONRT
 
      subroutine assign_conrt(itype, icon)

      implicit none
 
c
c     routine to determine the type of control line and to set the
c     icon variable to match.
c
c Variables
c ---------
c itype  -- The number of the string in con_names which was found
c           get_cmd
c icon   -- control variable for what to do with the rest of the
c        buffer and the next lines in the file.
c        icon = 1 for observation files to be read next from subsequent
c           records.
c        icon = 2 for station data to be read next from subsequent
c           records.
c        icon = 3 for GPS sat. data to be read next from subsequent
c           records.
c        icon = 4 for miscellaneous data to be read next.  For this
c           catagory the rest of the string is read to get the
c           desired information.
c
*                                        ! the MAKEXK parameter file
      include 'trackRT_cmds.h'
c
      integer*4 itype,  icon
 

c
c.... set icon based on the value of itype.
      icon = -1
      if( itype.le.file_end .and. itype.ge.file_start ) icon = 1
      if( itype.le.site_end .and. itype.ge.site_start ) icon = 2
c      if( itype.le.souc_end  .and. itype.ge.souc_start  ) icon = 3
      if( itype.le.mis_end  .and. itype.ge.mis_start  ) icon = 4
c
      return
      end
 
