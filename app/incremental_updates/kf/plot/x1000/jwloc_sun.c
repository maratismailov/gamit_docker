#include "x1000.h"
#include <X11/cursorfont.h>
#include <X11/keysym.h>

void jwloc_( unit, num, rtchar, xv, yv )

float *xv, *yv ;  /* Virtual coordinates of point sensed */
int *unit, *num ; /* Unit number and number of characters (both
                    ignored and never changed */
char *rtchar    ; /* Return character.  For mouse bottons LMR are
                     returned fro Left, Middle and Right buttons */

{
int cursor_shape = XC_crosshair ;
Cursor cursor;

int ix, iy ; /* Location of cursor in pixel coorinates */
int button_press ; /* Button pressed */
void check_window() ; /* MOD TAH 210523: Added declaration */

/* Create and display the cursor  */
/* MOD TAH 210523: Removed void; this is function call */
check_window(); 

cursor = XCreateFontCursor(display, cursor_shape);
XDefineCursor(display, win, cursor);
XFlush(display);
XSelectInput(display,win, ButtonPressMask );
while (1) {
   XNextEvent(display, &report);
   switch (report.type) {
      case ButtonPress :  
          ix = report.xbutton.x ; iy = report.xbutton.y ; /* Position of currsor */
          button_press = report.xbutton.button ;
          /* See which button was pressed */
          if( button_press == Button1 ) *rtchar = 'L' ;
          if( button_press == Button2 ) *rtchar = 'M' ;
          if( button_press == Button3 ) *rtchar = 'R' ;
          *xv = (float) ix / (float) Window_width; 
          *yv = (float) (Window_height-iy)/(float) Window_height;
          XBell(display, 50 );
          XUndefineCursor(display,win);
          XFlush(display);
          XFreeCursor(display, cursor);
      /*    XSelectInput(display, win, ExposureMask | ButtonPressMask);*/
          return;
     default: ;
      }
   }
}
