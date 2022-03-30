CTITLE RUN_GLINIT
 
      subroutine run_glinit

      implicit none 
 
*     Routine for use in globc: Simply call glinit
*     The SUBR passed tell the routine that it is a subroutine
*     of globc rather than a main program.
 
      call glinit('SUBR')
 
      return
      end
 
CTITLE RUN_GLFOR
 
      subroutine run_glfor

      implicit none 
 
*     Routine for use in globc: Simply call glfor
*     The SUBR passed tell the routine that it is a subroutine
*     of globc rather than a main program.

      call glfor('SUBR')
 
      return
      end
 
CTITLE RUN_GLBAK
 
      subroutine run_glbak

      implicit none 
 
*     Routine for use in globc: Simply call glbak
*     The SUBR passed tell the routine that it is a subroutine
*     of globc rather than a main program.
 
      call glbak('SUBR')
 
      return
      end
 
CTITLE RUN_GLOUT
 
      subroutine run_glout( out_file, opts)

      implicit none 
 
*     Routine for use in globc: Simply call glout
*     The SUBR passed tell the routine that it is a subroutine
*     of globc rather than a main program.

      integer*4 opts
      character*(*) out_file

      call glout('SUBR', out_file, opts)
 
      return
      end
      
CTITLE RUN_glorg
 
      subroutine run_glorg( cmd_file, out_file, opts)

      implicit none 
 
*     Routine for use in globc: Simply call glorg
*     The SUBR passed tell the routine that it is a subroutine
*     of globc rather than a main program.

      integer*4 opts
      character*(*) out_file, cmd_file

      call glorg('SUBR', cmd_file, out_file, opts)
 
      return
      end
 
CTITLE RUN_GLSAVE
 
      subroutine run_glsave

      implicit none 
 
*     Routine for use in globc: Simply call glsave
*     The SUBR passed tell the routine that it is a subroutine
*     of globc rather than a main program.
 
      call glsave('SUBR','NO')
 
      return
      end
 
 
