#include <fcntl.h>
#include <termios.h>
#include <unistd.h>

/* INKEY

   Rewritten 910204 J.L. Davis
   Updated to use revised system calls and termios structure.

   See termio(7), tcattribute(3C), and /usr/include/sys/termio.h

   Routine to read a single key from the keyboard using non-canonical
   input.

   Returned arguments:   -1 = Error
                          0 = Not terminal type, character not read
                          1 = Character read successfully, returned in ikey

   NOTE: Character is not echoed to screen.
	Modified to force "inkey_" as the name of the function.
	Lada L. Dimitrova          June 30 1998 
                                                                            */

/* LLD>> */
/*
#ifdef _NEEDED
int inkey_(ikey)
#else
int inkey(ikey)
#endif

  IRIX:	C functions used in FORTRAN 77 programs need name ending in _
 	in the FORTRAN 77 code they are refered by the name without _
		ex. if FORTRAN 77 needs a function inkey, 
		the C function should be named inkey_
*/

int inkey_(ikey)
/* LLD << */

    char *ikey; 

{
  void perror();
  struct termios tbuf, tbufsave;
  int  i;
  char c;
  static int stdio_fd = 0;

  *ikey = '\0';

  /* Check that this is a terminal */

  if(!isatty(stdio_fd))
    return(0);

  /* Get current attributes */

  if (tcgetattr(stdio_fd,&tbuf) == -1) 
  {
    perror("tcgetattr");
    return(-1);
  }

  /* Save old values */

  tbufsave = tbuf;

  /* Set TIME and MIN */

  tbuf.c_cc[VMIN]  = 0;
  tbuf.c_cc[VTIME] = 0;

  /* Set Non-canonical input processing, no echo. */

  tbuf.c_lflag = IEXTEN + ECHOK + ECHOE + ISIG;

  if (tcsetattr(stdio_fd,TCSAFLUSH,&tbuf) == -1)
  {
    perror("tcsetattr");
    return (-1);
  }

  /* Loop until a non-null is read */

  c = '\0';
  while (c == '\0')
    read (stdio_fd, &c,1);

  /* Restore old values */

  tbuf = tbufsave;
  if (tcsetattr(stdio_fd,TCSAFLUSH,&tbuf) == -1)
  {
    perror("tcsetattr");
    return (-1);
  }

  /* Save value and return */

  *ikey = c;

  return(1);
}
