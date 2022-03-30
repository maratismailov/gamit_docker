#include "x1000.h"

void
gsclip( clip )

int *clip ;  /* Clipping to set.  If =1 then clipping will
                be applied to current viewport; if 0 then
                clipping will be turned off */
{
int width_clip ;   /*Width of clip region (+1 to make we see lines) */
int height_clip ;  /*height of clip region (+1 to make we see lines) */
int xorg_clip, yorg_clip ; /* Orgin of clip region (-1 to make we are in) */
void check_window90 ; /* MOD TAH 210523: Added declaration */

XRectangle clip_field ;    /* Rectangle defined by viewport */
XRectangle *clip_test ; 
Pixmap clip_mask ;         /* Pixmap definition of mask */

/* see which clipping is to be applied */
/* MOD TAH 210523: Removed void; this is function call */
check_window();

if( *clip ) 
     /* OK we want to do clipping, First create the 
        Pixmap that we need */
    { 
      width_clip = (int) ((float) Window_width * (xvr-xvl) + 1) ;
      height_clip = (int) ((float) Window_height * (xvt-xvb) + 2) ;
      xorg_clip = (int) ((float) Window_width * xvl -0 );
      yorg_clip = (int) ((float) Window_height* (1.0-xvt) -0 );

/*     Set the field for the rectangle we have  */
      clip_field.x = 0 ; clip_field.y = 0 ;
      clip_field.width = (short) width_clip ;
      clip_field.height = (short) height_clip ; 
      clip_test = &clip_field ;

     /* Create the Pixmap for the width and height 
      clip_mask = XCreatePixmap(display, win, width_clip, height_clip, 1); */
      XSetClipRectangles( display, gc, xorg_clip, yorg_clip,
                                      clip_test, 1, Unsorted ) ;

     /* Now set the mask 
      XSetClipMask( display, gc, clip_mask );
      XSetClipOrigin ( display, gc, xorg_clip, yorg_clip) ; */
    }
else       /* Set the clip mask to zero */
    {
      XSetClipMask( display, gc, None ) ;
    } 

/* Save the clip mask incase we need it */
jclip = *clip ;

}



