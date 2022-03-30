/* X wrappers for CVIEW -- System 5 (HP, DEC, etc) version
   Written by Andrea Donnellan 7/24/1990  
   E-mail andrea@seismo.gps.caltech.edu
   Caltech */

/* Color option introduced by M. Burc Oral 2/19/1992 using a modified
 * get_color_ro.c from a package from 1989 O'Reilly and Associates, Inc.
 * if you do not like the color scheme, it can be changed in get_color_ro.c. 
 * E-mail oral@erl.mit.edu
 * MIT- ERL  */

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <stdio.h>
#include <stdlib.h>

#include "cview_icon"

Display     *display;
int          screen;
GC           gc;
XFontStruct *font_info;
XEvent       report;
XSizeHints   size_hints;
Window       win;
unsigned long foreground_pixel, background_pixel, border_pixel;
unsigned int seconds=0;
int get_colors();

int gsetup (ibmx0,ibmy0,bmwdth,bmhght,lcolor,csize,ipos)
short *ibmx0,*ibmy0,*bmwdth,*bmhght,*csize,*ipos;
unsigned *lcolor;
{
/* set up the graphics package
     returns:
        bmwdth bit map width
        bmhght bit map height
        csize  width of a character in screen units          
     all other variables are ignored. */

   char    *window_name = "CVIEW";
   char    *icon_name = "CVIEW";
   char    *display_name = NULL;
   unsigned int width,height,winwidth,winheight;
   int border_width = 0;
   int x=0, y=0;
   Pixmap icon_pixmap;
   XSizeHints size_hints;
   char   **argv;
   int      argc;
   char    *fontname = "9x15";

   argv = (char **) malloc (1 * sizeof(char));
   argv[0] = (char *) malloc(5 * sizeof(char));
   strcpy(argv[0],"cview");
   argc = 1;

/* connect to X server */
   if ((display=XOpenDisplay(display_name)) == NULL)
      {
      (void) fprintf(stderr, "cview:  cannot connect to X server %s \n",
		    XDisplayName(display_name));
      exit(-1);
      }

   screen = DefaultScreen(display);

   width  = DisplayWidth(display, screen);
   height = DisplayHeight(display,screen);

/* MOD TAH 031218/080925: Limit maximim window to 1280x1024 */
   if ( width > 1280 && height > 1024) {
      width = 1280 ; height = 1024;
   }
   winwidth = width-60; 
/* MOD TAH/RWK 090508: Increase height allowance from 60 to 120 for high aspect ratio screens with ribbons
   winheight = height-120; 

/* If the screen size is <= 1024 pixels wide, use a smaller font size */
   if ( width <= 1025 ) {
      fontname = "7x13";
   }

   *bmwdth = (short) winwidth;
   *bmhght = (short) winheight;
   get_colors();
   win = XCreateSimpleWindow(display, RootWindow(display,screen),
	  x, y, winwidth, winheight,
	  border_width, border_pixel,
	   background_pixel);
/* 	  border_width, BlackPixel(display,screen),  */
/* 	  WhitePixel(display,screen)); */


   icon_pixmap = XCreateBitmapFromData(display, win, cview_icon_bits,
	   cview_icon_width, cview_icon_height);

   size_hints.flags = PSize | PMinSize;
   size_hints.width = winwidth;
   size_hints.height = winheight;
   size_hints.min_width = 350;
   size_hints.min_height = 250;

   XSetStandardProperties(display, win, window_name, icon_name,
	 icon_pixmap, argv, argc, &size_hints);

   XSelectInput(display, win, ExposureMask | ButtonPressMask);

   if ((font_info = XLoadQueryFont(display,fontname)) == NULL)
      {
      (void) fprintf(stderr, "CVIEW:  Cannot open %s font \n",fontname);
      exit(-1);
      }

   gc = XCreateGC(display,win, 0, NULL);
   XSetForeground(display, gc, foreground_pixel);
   XSetBackground(display, gc, background_pixel);
/* XSetForeground(display, gc, BlackPixel(display,screen));  */
/* XSetBackground(display, gc, WhitePixel(display,screen));  */

   XSetFont(display, gc, font_info->fid);

   csize[0] = (short) XTextWidth(font_info, " ", 1);
   csize[1] = (short) (font_info->ascent + font_info->descent);
   XMapWindow(display, win);

   XNextEvent(display, &report);
/*   XSelectInput(display, win, ButtonPressMask);*/

}

void gmsgc (iwin,msg,size)
short *iwin,*size;
char  *msg;
{
/*  prints a message in centered in the designated window */

   int length, x, y;
   int igerr;
   int x1,y1;
   unsigned int width,height;

   x1 = (int) iwin[0]; 
   y1 = (int) iwin[1];
   width = (unsigned int) iwin[2];
   height = (unsigned int) iwin[3];

   length = strlen(msg);

   x = (int)(iwin[0] + iwin[2]/2 - size[0]/2);
   y = (int)(iwin[1] + iwin[3]/2 + size[1]/2);

   XSetForeground(display, gc, background_pixel); 
   XFillRectangle(display, win, gc, x1, y1, width, height);
   XSetForeground(display, gc, foreground_pixel); 
/*   XSetForeground(display, gc, WhitePixel(display,screen));  */
/*   XFillRectangle(display, win, gc, x1, y1, width, height);  */
/*   XSetForeground(display, gc, BlackPixel(display,screen));  */
  

   XDrawString(display, win, gc, x, y, msg, length);
}

        
void gtsizec (msg,size,len)
char *msg;
short *size;
int   *len;
{
/* return x and y extents of text */

   while ((msg[*len] == ' ') && (*len > 0))
      --*len;
   msg[*len+1]='\0';

   size[0] = (short) XTextWidth(font_info, msg, *len+1);
   size[1] = (short) (font_info->ascent + font_info->descent);
}                           

                          
void gmvcur (ipos,igerr)
short *ipos;
int   *igerr;
{
/* move the cursor to position x=ipos(1),y=ipos(2)
   Not used.   */
}

void gclip(iwin,lon)
short *iwin;
unsigned int *lon;
{
/* turn graphics clipping on or off
   lon is 1 to turn it on
   the iwin is defined as x,y,dx,dy
   I don't think this function is necessary
   and it is therefore not used.*/

}

void gmouse (key,pos,wait)
char  *key;
short *pos;
unsigned *wait;
{
/* Get input from mouse
    input:
      wait = 1 wait for next event
    output:
      pos  x,y of cursor
      key  mouse button */

   Window  root, child;
   int     root_x, root_y, win_x, win_y;
   unsigned int button;
   int     cwait;

   cwait = (int)*wait;

   if (!cwait)
      {
      cwait = XCheckMaskEvent(display, ButtonPressMask, &report);
      if (cwait)
	 {
         switch (report.xbutton.button)
	    {
	    case 1:
	       seconds += 500000;
	       cwait = 0;
	       break;
	    case 2:
	       seconds = 0;
	       break;
	    case 3:
	       if (seconds > 0)
	          seconds -= 500000;
	       cwait = 0;
	       break;
	    default:
	       break;
	    }
	 }
      sleep(seconds);
      }
   else
      {
      XNextEvent(display, &report);
      seconds = 0;
      }

   pos[0] = (short) report.xbutton.x;
   pos[1] = (short) report.xbutton.y;
   if (!cwait)
      *key = '.';
   else
      {
      sprintf(key,"%1d",(int)(report.xbutton.button+1));
      if (report.type == Expose)
	 {
	 pos[0] = 0;
	 pos[1] = 0;
	 }
      }
}

         
void gline (ix1,iy1,ix2,iy2,igerr)
short *ix1,*iy1,*ix2,*iy2;
int    *igerr;
{
/*  draw a line from (ix1,iy1) to (ix2,iy2) */

   XDrawLine(display, win, gc, (int)*ix1, (int)*iy1, 
			       (int)*ix2, (int)*iy2);
}


void gbox (iwin,igerr)
short *iwin;
int   *igerr;
{
/*  fill in window iwin */

   int x1,y1;
   unsigned int width,height;

   x1 = (int) iwin[0]; 
   y1 = (int) iwin[1];
   width = (unsigned int) iwin[2];
   height = (unsigned int) iwin[3];
   XDrawRectangle(display, win, gc, x1, y1, width, height);
   XFillRectangle(display, win, gc, x1, y1, width, height);
}  

void grect (ix,iy, iradi, igerr)
short *ix, *iy, *iradi;
int    igerr;
{
/* draw an open rectangle of height and width iradi centred at ix,iy */

   int x,y;

   x = (int)*ix-(unsigned int)*iradi;
   y = (int)*iy-(unsigned int)*iradi;

   XDrawRectangle(display, win, gc, x, y, (unsigned int)*iradi*2,
	     (unsigned int)*iradi*2);
}  


void gcirc (ix,iy, iradi, igerr)
short *ix, *iy, *iradi;
int    igerr;
{
/* draw an open circle of radius iradi at ix,iy */

   int x,y;

   x = (int)*ix-(unsigned int)*iradi;
   y = (int)*iy-(unsigned int)*iradi;

   XDrawArc(display, win, gc, x, y, (unsigned int)*iradi*2,
	     (unsigned int)*iradi*2, 0, 360*64);
}  

void gdot (ix,iy, iradi, igerr)
short *ix, *iy, *iradi;
int    igerr;
{
/* draw an filled circle of radius iradi at ix,iy */
                     
   int x,y;

   x = (int)*ix-(unsigned int)*iradi;
   y = (int)*iy-(unsigned int)*iradi;

   XDrawArc(display, win, gc, x, y, (unsigned int)*iradi*2,
	     (unsigned int)*iradi*2, 0, 360*64);
   XFillArc(display, win, gc, x, y, (unsigned int)*iradi*2,
	     (unsigned int)*iradi*2, 0, 360*64);
}

void gplus (ix,iy, iradi, igerr)
short *ix,*iy,*iradi;
int   *igerr;
{
/* draw a plus sign of radius iradi at ptcent(1),ptcent(2) */

   XDrawLine(display, win, gc, (int)(*ix+*iradi), (int)*iy, 
			       (int)(*ix-*iradi), (int)*iy);
   XDrawLine(display, win, gc, (int)*ix, (int)(*iy+*iradi), 
			       (int)*ix, (int)(*iy-*iradi));
}
                       
void gend (igerr)
int *igerr;
{
/*  terminate graphics package */

   XUnloadFont(display, font_info->fid);
   XFreeGC(display, gc);
   XCloseDisplay(display);
} 
