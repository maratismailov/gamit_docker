#include "x1000.h"

void
jclr()

/* Routine to clear the screen */

{
XClearWindow( display, win ) ;
XFlush( display) ;
}


