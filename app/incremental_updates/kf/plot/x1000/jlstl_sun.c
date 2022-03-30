#include "x1000.h"
 
void
jlstl_( line_type )

/* Routine to set line style.  Mod 10 sets the dash type,
   part greater then 10 set the width. */

int *line_type ; /* Line type passed from plotting program */

{
int dash_style ; /* Dash style to use */
int line_width ; /* Width of line to use */
int line_style ; /* Style of line */
int cap_style = CapRound ;
int join_style = JoinRound ;

static const char dotted[2] = {1,3};
static const char dot_dash[4] = {1,3,2,4} ;
static const char short_dash[2] = {4,4} ;
static const char long_dash[2] = {4,7} ;

int dash_offset = 0 ;

/* Increment pixel width at half the rate of NCAR so that the
   appearance of the plots will be more similar */
line_width = *line_type / 20 + 1 ;
dash_style = *line_type - (line_width-1)*20 ;


/* Set the dash style */
switch ( dash_style ) {
      case 1:              /* Solid line */
        /*  gcv.line_style = LineSolid ; */
          XSetLineAttributes(display, gc, line_width, LineSolid ,
                             cap_style, join_style ) ;
          break ;
      case 2:             /* dotted line */
          XSetDashes( display, gc, dash_offset, dotted, 2);
          XSetLineAttributes(display, gc, line_width, LineOnOffDash,
                             cap_style, join_style ) ;
          break ;
      case 3:             /* dot dash */
          XSetDashes( display, gc, dash_offset, dot_dash, 4);
          XSetLineAttributes(display, gc, line_width, LineOnOffDash,
                             cap_style, join_style ) ;
          break ;
      case 4:             /* short dash */
          XSetDashes( display, gc, dash_offset, short_dash, 2);
          XSetLineAttributes(display, gc, line_width, LineOnOffDash,
                             cap_style, join_style ) ;
          break ;
      case 5:             /* long dash */
          XSetDashes( display, gc, dash_offset, long_dash, 2);
          XSetLineAttributes(display, gc, line_width, LineOnOffDash,
                             cap_style, join_style ) ;
          break;
      default:            /*Solid line*/
          XSetLineAttributes(display, gc, line_width, LineSolid ,
                             cap_style, join_style ) ;
	}
}



