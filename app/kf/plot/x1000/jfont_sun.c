#include "x1000.h"
   
void
jfont_( name, font_size)
char *name ; /* Name of the font to be set */
int font_size[2] ; /* Font size in pixels */

{
int nullterm() ;    /* Function to return length of string
                       and null terminate it */

   if( nullterm(name)<=2 ) return  ; /* Null terminate the name string */
   if ((font_info = XLoadQueryFont(display,name)) == NULL)
      {
      (void) fprintf(stderr, "X1000: Cannot open %s font \n",name);
      return;
      }

   XSetFont(display, gc, font_info->fid);
   curr_charsize[0] = (int) XTextWidth(font_info, "_", 1);
   curr_charsize[1] = (int) (font_info->ascent + font_info->descent);
   font_size[0] = curr_charsize[0];
   font_size[1] = curr_charsize[1];

   /* printf("Character size pixels Width %d Height %d =  %d + %d Font %s\n",
        curr_charsize[0], curr_charsize[1],  
        (int) font_info->ascent , (int) font_info->descent, name); */
   strcpy(&font[0],name) ;

}


