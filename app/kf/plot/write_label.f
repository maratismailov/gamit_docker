CTITLE    ...............................................................
 
      subroutine write_label(label)
c
c     Routine to write all labels on the plot.  It will reset the
c     scales and view of the plot so that the aspect ratio of the
c     world coordinates martch that of the device.  In this way the
c     character size will not depend on orientation.
c
c     The routine assumes that the current position for the label
c     is set correctly and that the orientation has also been set.
c
c Include files
c -------------
*                         ! the solvk parameter file
      include 'plot_param.h'
c
*                         ! the solvk common block
      include 'plot_com.h'
c
c Variables
c ---------
c label  -- the label to be output to the plot
c
      character*(*) label
c
c
c Local variables
c ---------------
c new_label -- Label modified so that lower case letters can be
*     handled.
c width -- width of character in new world cooridinates
c height -- height in new world coordinates
c gap   -- character gap in fraction of width
* nlen  -- Length of new label as it is being made.
* inquote -- Logical indicating that we are in a quote command - not used
* inupper -- Logical indicating that we are in upper case  - not used
c
      character*150 new_label
c      character*1 one_char - not used
c
      integer*4 ilen, nlen
 
c      logical inquote, inupper  - not used
c
      real*4 width, height, gap
 
c
c Functions
c ---------
c trimlen -- HP utility
c
      integer*4 trimlen
 
c
c.... Firstly check length of label
      ilen = trimlen(label)
*                               ! no label to be output
      if( ilen.le.0 ) return
c
      width  = charsz_x
      height = charsz_y
      gap    = 0.05
c
c.... Now set the character size
      call jcsiz(width, height, gap)
 
*     Select the font we want to use
      call null_terminate( font_name )
* MOD TAH 131018: Return font size in pixels
      call jfont( font_name, font_size )

*     Procress the string (mainly for NCAR graphics)
      if( Trimlen( label ).le.0 ) RETURN

      call jtexp( label, new_label, font_type )
 
      nlen = trimlen(new_label)
*                                 ! use high quality text
      call jtexh(nlen,new_label)
c
      return
      end
 
