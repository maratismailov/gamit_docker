/*
  inkey - return a single character from file descriptor 0 if its tty
          ikey returns 1: one character read in *ikey
                       0: no characters available or not tty
                      -1: error
  NOTE: calling this routine will destroy type-ahead
        so before doing another read you better:

                     short inkey,ikey;
                     while(inkey(&ikey)>0) ;
        or:
                     INTEGER*2 INKEY,IKEY
                     DO WHILE(INKEY(IKEY).GT.0)
                     ENDDO
*/
 
#include <fcntl.h>
#include <sgtty.h>
#include <unistd.h>

int inkey_(ikey)
    int *ikey; 
{
    void perror();
    struct sgttyb tbuf, tbufsave;
    int iret, arg;
    char c='\0';
    static int tty_fd=0;

    *ikey=0;
    if(!isatty(tty_fd)) return(0);

    if(ioctl(tty_fd,TIOCGETP,&tbuf) == -1) {
      perror("ioctl");
      return(-1);
    }

    tbufsave = tbuf;		/* save old terminal parameters  */
    tbuf.sg_flags = 218   ;    /* turn off record processing */

    if(ioctl(tty_fd,TIOCSETN,&tbuf) == -1) {
	perror("ioctl 2");
	return(-1);
    }   
    if( ioctl(tty_fd,FIONREAD,&arg) == -1) {
	perror("ioctl 2");
	return(-1);
    }

    if( arg>0 ) { 
    if((iret=read(tty_fd,&c,1))==-1) {
      perror("read");
      return(-1);
    };}
			
    if(ioctl(tty_fd,TIOCSETN,&tbufsave) == -1) {
	perror("ioctl 3");
	return(-1);
    }   

    *ikey=c;
    return(arg);
}
