#include "x1000.h"

void
DrawChar( ix, iy, type )

int ix, iy ;   /* Screen coordinates of mark to be drawn */
int type   ;   /* ASCII Code of character to be output */
{

/* Local variables */
char mark  ;   /* The character to be written to the window */
XCharStruct mark_size ; /* Dimensions of character, to be used
                   to center the character */
int  mx, my        ; /* Start point of character so it is centered  */
int  mwidth, mheight ; /* Mark width and height */

/* Move the point so that the center of the character will be
   al the mark */
mx = ix - curr_charsize[0]/2 ;
my = iy + curr_charsize[1]/2 ;

/* Get the character */
mark = (char) type ;
XDrawString( display, win, gc, mx, my, &mark, 1 ) ;

}
          
