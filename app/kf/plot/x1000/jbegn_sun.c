
#include <stdlib.h>
#include "x1000.h"

void jbegn_(display_name, x0, y0, wdt0, hgt0)

char *display_name; /* Name of the Xhost for the display; if NULL use
                       the name from DISPLAY Environment variable */
int *x0, *y0, *wdt0, *hgt0 ;  /* optional passed values for x.y, origin
                       and width and height.  If x0 is zero defaults used */

{
/* set up the graphics package */

   char    *window_name = "X1000";
   char    *icon_name = "x1000";
   unsigned int width,height,winwidth,winheight;
   int border_width = 0;
   int x, y;
   Pixmap icon_pixmap;
   XSizeHints size_hints;
   char   **argv;
   int      argc;
   char  *default_font = "6x10";
   short csize[2];
   int nullterm(), lend ; /* Function to null terminate a string and
                             then length of the dsplay name */

   argv = (char **) malloc (1 * sizeof(char));
   argv[0] = (char *) malloc(5 * sizeof(char));
   strcpy(argv[0],"X1000");
   argc = 1;

/* connect to X server */
   lend = nullterm(display_name) ;
   
   if ((display=XOpenDisplay(display_name)) == NULL)
      {
      (void) fprintf(stderr, "x1000:  cannot connect to X server %s \n",
		    XDisplayName(display_name));
      exit(-1);
      }

   screen = DefaultScreen(display);

   width  = DisplayWidth(display, screen);
   height = DisplayHeight(display,screen);

   winheight = 3*height/5;
   if( winheight > 560 ) winheight = 560 ;
   winwidth  = winheight ;
   x=winheight/10;
   y=winheight/10;


/* Now see if user has passed values, if so use these */
   if( *x0 != 0 ) {
       x = *x0 ; y = *y0; 
       winwidth = *wdt0 ; winheight = *hgt0 ;
   }

   Window_width=winwidth;
   Window_height=winheight;
       
   printf ("Screen Size %d %d  Origin %d %d\n",Window_width, Window_height, x, y );

   win = XCreateSimpleWindow(display, RootWindow(display,screen),
	     x, y, winwidth, winheight,
 	     border_width,  BlackPixel(display,screen),  
        WhitePixel(display,screen)); 


   icon_pixmap = XCreateBitmapFromData(display, win, x1000_bits,
	   x1000_width, x1000_height);

   size_hints.flags = USPosition | PSize | PMinSize;
   size_hints.x = x;
   size_hints.y = y;
   size_hints.width = winwidth;
   size_hints.height = winheight;
   size_hints.min_width = 0;
   size_hints.min_height = 0;

   XSetStandardProperties(display, win, window_name, icon_name,
	 icon_pixmap, argv, argc, &size_hints);

   XSelectInput(display, win, ExposureMask | ButtonPressMask | 
                StructureNotifyMask | KeyPressMask );

   gc = XCreateGC(display,win, 0, &gcv ); 

   XSetForeground(display, gc, BlackPixel(display,screen));
   XSetBackground(display, gc, WhitePixel(display,screen));  


/* initialize current cursor position */
/* does not necessarily correspond to position of mouse cursor ! */
   Curx=0;
   Cury=0;

/* always keep backing store so that never lose contents of window */
   valuemask = CWBackingStore | CWBitGravity;
   setwinattr.backing_store=Always;
   setwinattr.bit_gravity = CenterGravity ;
   XChangeWindowAttributes(display,win,valuemask,&setwinattr);

   XMapWindow(display, win);
   XNextEvent(display, &report);

/* Set up the default font selection */ 
  
   void jfont_( char default_font ) ;

/* Set no clipping default */
   jclip = 0 ;
   XSetClipMask(display, gc, None) ;

}

