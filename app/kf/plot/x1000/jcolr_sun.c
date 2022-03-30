#include "x1000.h"

void jcolr_( pen )

/* Emulation of the color routines */

int *pen;   /* Pen number which set color.  The first color
               is used to set the background color. */
{
XColor setcol_col;
Colormap colormap;
#define MAX_COLORS 35


char colors[][MAX_COLORS] = {
"#FFFFFF",  // White
"#000000",  // Black
"#FF0000",  // Red
"#008000",  // Green
"#0000FF",  // Blue
"#FF00FF",  // Magenta
"#00FFFF",  // Cyan
"#FFFF00",  // Yellow
"#808080",  // Gray
"#C0C0C0",  // Silver
"#8B0000",  // Dark Red
"#FF4500",  // Orange Red
"#D2691E",  // Chocolate
"#FF8C00",  // Dark Orange
"#FFA500",  // Orange
"#FFD700",  // Gold
"#808000",  // Olive
"#006400",  // Dark Green
"#00FF00",  // Lime
"#98EE90",  // Light Green
"#66CFAA",  // Aquamarine
"#008080",  // Teal
"#2F4F4F",  // Dark Slate Gray
"#AFEEEE",  // Pale Turquoise
"#5F9EA0",  // Cadet Blue
"#00BFFF",  // Deep Sky Blue
"#4682B4",  // Steel Blue
"#778899",  // Light Slate Gray
"#000080",  // Navy
"#9400D3",  // Dark Violet
"#800080",  // Purple
"#DDA0DD",  // Plum
"#C71585",  // Medium Violet Red
"#DC143C",  // Crimson
"#FFB6C1"   // LightPink
}; 
char setcol[8];

colormap = DefaultColormap(display, 0);

strcpy(setcol,colors[*pen]);
printf("Setting pen %d, colors %s\n",*pen,setcol);

XParseColor(display, colormap, setcol, &setcol_col);
XAllocColor(display, colormap, &setcol_col);
XSetForeground(display, gc, setcol_col.pixel);
}


