#include <stdlib.h>
#include "x1000.h"
  
void
jtexh_( ilen, text )

int *ilen ; /* Length of string to be output */
char *text ; /* Text of string to be output */

{
/* Routine to write out text using current justification,
   orientation.  The cursor is left positioned at the last
   point of the text (ready to write next character) */

int dcx, dcy ; /* Space in pixel for each character.  X and
                  Y corrdinate step */
int ocx, ocy ; /* Offset from start position for the current
                  string with current justification */
int scx, scy ; /* Start position for text */

int pcx, pcy ; /* Position of current character */

int i        ; /* Loop counter for stepping through string */

char oc      ; /* Character currently being output */

/* In this version we use just the width and height of the
   font to compute justification (this is to handle sloped
   text. */

/* Compute the spacing of the charcaters */
if( abs(jor_tan) <=  1.0 )
    {
     dcx = curr_charsize[0] ; dcy = curr_charsize[1]* jor_tan ;
     ocx = ((*ilen)*dcx+curr_charsize[0])*jcenter ; 
     ocy = ((*ilen)*dcy+curr_charsize[1])*jheight ;
    }
else
    {
     dcx = curr_charsize[0]/jor_tan ; dcy = curr_charsize[1] ;
     ocx = ((*ilen)*dcx+curr_charsize[0])*jheight ;
     ocy = ((*ilen)*dcy+curr_charsize[1])*jcenter ;
    } ;

/* Move the cursor the correct position and start writing */

scx = Curx - ocx ; scy = Cury + ocy ;

/* Now write the text */
for ( i=0; i<*ilen ; ++i )
     {
      pcx = scx + i*dcx ; pcy = scy - i*dcy ;
      oc = *(text+i) ;
      XDrawString( display, win, gc, pcx, pcy, &oc, 1 );
     } ;

/* Now move cursor to start of next character */
Curx = pcx + dcx ; Cury = pcy - dcy ;

}

