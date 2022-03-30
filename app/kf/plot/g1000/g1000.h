*
*     Include file for the g1000 emulation program
 
*   jcenter     - Defines the centering value to be used
*               - for text strings
*   jor_deg     - Set the orientation of text
*   jsze        - Character size for text
*   jdash       - Dash pattern to use (if <=1 then solid line
*               - used)
*   jwidth      - Sets the width of the line to be used.
*   jfnt        - Current font number set by user. 1 - Cartographic,
*               - 2 - Roman, 3 - greek.  Most important here is
*               - type 1 which doubles the font size (because this
*               - is a small font.
 
      integer*4 jcenter, jor_deg, jsze, jdash, jwidth, jfnt

*   jfntset     - Indicates that the font is set.
*   jmode       - TRue if we are mapping mode (set with jmapm call)

      logical   jfntset, jmode
 
*   gvl, gvr, gvb, gvt  - Virtual coordinates of plotter area
*               - to use.
*   gwl, gwr, gwb, gwt  - World coordinates of the edges of the
*               - area to use.
*   gheight     - Heigth justificaation (value between -0.5 and 0.5)
 
      real*4 gvl, gvr, gvb, gvt, gwl, gwr, gwb, gwt, gheight

*   gmetaf      - Name of the metafile to use

      character*256 gmetaf
 
*--------- Common declaration -------------
 
      common / g1000_comm / jcenter, jor_deg, jsze, jdash, jwidth,
     .    jfnt, jmode, gvl, gvr, gvb, gvt, gwl, gwr, gwb, gwt, gheight,
     .    jfntset, gmetaf
 
*--------------------------------------------------------------------
 
