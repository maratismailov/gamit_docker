#include "x1000.h"

void check_window() 
{
/* This routine will check if a window has been resized.  If it
   has then the window will be updated.  A call to jclr should be
   made to clear the screen after resizing */

XEvent report;

XFlush(display);
XSelectInput( display, win, StructureNotifyMask ); 

/* XNextEvent ( display, &report ) ; */
XCheckMaskEvent( display, StructureNotifyMask, &report );

switch ( report.type ) {
    case ConfigureNotify:   /* Window resized */
        if( (report.xconfigure.send_event == 0 ||
             report.xconfigure.send_event == 1) &&
             (int) report.xconfigure.width < 2000 &&
             (int) report.xconfigure.height < 2000 ) {
        Window_width = (int) report.xconfigure.width ;
        Window_height = (int) report.xconfigure.height ;
       /* printf("Window width and height %d %d %d\n", 
                Window_width, Window_height,
                report.xconfigure.send_event); */ }
    }     /* End switch */
}

