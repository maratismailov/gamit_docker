CTITLE jend
 
      subroutine jend     
c
c     Routine to cleanly exit from GKS

      include 'g1000.h'
c
c
c trimlen - Length of string
c fmprename - renames a file
c ierr      - IOSTAT error on rename.

      integer*4 trimlen, fmprename, ierr 
c
c.... Flush the graphics 1000 buffers
      call jmcur
c
*     close seesion
      call clsgks

*     Rename the meta file created by NCAR to the desired name
      if( trimlen( gmetaf ) .gt.0 .and.
     .    gmetaf(1:1).ne.char(0)        ) then
          ierr = fmprename( 'gmeta', gmetaf, '' )
          call report_error('FmpRename',ierr,'renam',gmetaf, 0,
     .                      'PLOT/Main')
      end if

c
c.... Thats all
      return
      end
 
