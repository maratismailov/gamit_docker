c    Program (for shell scripts) to call gamit/lib subroutine to 
c    determine whether the L-file is old-style (spherical) or GLOBK apr.  
c    R. King 15 Feb 2003   

c    Input: l-file name, on command line

c    Output: 'old' or 'apr'   

      implicit none

      integer*4 lflag,iarg,iclarg

      character*256 lfname 
      character*3 ltype
                  
      lflag = -1

      iarg = iclarg(1,lfname)
      if( iarg.le.0 ) then 
         write(*,*) 'Missing filename argument for check_lfile'
         stop
      endif

      call crd_file_type(lfname,lflag) 

      if( lflag.eq.0 ) then
         ltype = 'old'
      elseif( lflag.eq.1 ) then
         ltype = 'apr'
      else
         write(*,*) 'Invalid lfile type in check_lfile'
      endif

      write(*,'(a3)') ltype
      end


