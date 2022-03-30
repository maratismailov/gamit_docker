#include "x1000.h"

void
jend()

{
   /*  terminate graphics package */
      XUnloadFont(display, font_info->fid);
      XFreeGC(display, gc);
      XCloseDisplay(display);
      fflush(stdout);
      Curx=0;
      Cury=0;
}
