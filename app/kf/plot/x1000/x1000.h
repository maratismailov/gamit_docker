#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <stdio.h>
/* #include <math.h> */
#include "x1000.icon"
#define WINFRAC 2      /* x1k window fraction of largest window (1/WINFRAC)*/
#define MAXSTRLEN 256
  /* maximum length of a label */

Display     *display;
int          screen;
GC           gc;
XGCValues    gcv;
XFontStruct *font_info;
XEvent       report;
XSizeHints   size_hints;
Window       win;
Pixmap       icon_pixmap ;
Pixmap       clip_mask ;
XSetWindowAttributes setwinattr;
unsigned long valuemask;

/* Definitions for defining the plotting area and the size of this 
   area */

float xvl, xvr, xvb, xvt  ; 
      /* Virtual corrdinates of the cornors of the
         area for plot (left, right, bottom, top */
float xwl, xwr, xwb, xwt; 
      /* World coordinates of the cornors of the
         virtual area */
int Window_width; /* current width of window in pixels */
int Window_height; /* current height of window in pixels */
int curr_charsize[2];  /* current character size for actual font used 
                          width and height (pixels) */
int usr_charsize    ;  /* User character size in pixels */
int Curx;           /* current x coord of cursor (viewport scale) */
int Cury;           /* current y coord of cursor (viewport scale) */

/* Current font style and info */

float jcenter;   /* Defines the centering to be used 0-1       
                    for left through right justification */
float jheight;  /* Justification in height for strings 
                   (0-1.0 for bottom to top justifcation */
int jor_tan;     /* Defines orientation of baseline of the text
                    strings given as tan of angle */
int jzse;        /* Current character size as specified by user
                    which is different from current font size 
                   (pixels) */
int jdash;      /* Current dash pattern (init to solid line)*/
int jwidth;      /* Current width of lines (init to 1 pixel)*/

int jclip ;     /* Set to one to enable clipping, 0 for no clipping*/

char font[MAXSTRLEN] ; /* Name of current font to be used */



