/* pickfn is a function callable from FORTRAN to pick up
   a user choice of wildcard matched files in current
   directory (see test fortran code below)

   Has been tested on Solaris 5.5.1, HP-UX 10.20, Linux 2.0.34.
   The Sun/LINUX/DEC version differsn form the HP/IBM only in the
   underscore at the end of the function name.

   by pfang@ucsd.edu July 9, 1998, and January 3, 1999, `
   with mods by S. McClusky July 20, 1998  */

#include <dirent.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <fnmatch.h>
#include <ctype.h>
#include <stdlib.h>
#define	MAXSTR 64
/* MOD TAH 080512: Increased number of allowed files to 512 */
#define MAXFIL 512
  
void
pickfn_ (retval_p, len, ch_ptr, n_ptr, ch_len)
char *retval_p, *ch_ptr;
int len, *n_ptr, ch_len;

{
	DIR *dirp;
	struct dirent *direntp;
	char aline[MAXSTR], matched[MAXFIL][MAXSTR];
	int i, n, off_stdin;
                 
/*   off_stdin = ftell(stdin);  */
/*   printf("Unit 5 offset: %d \n",off_stdin); */
	ch_ptr[*n_ptr] = (char) NULL; 
/*   printf("Ch_ptr %s len %d ch_len %d\n",ch_ptr,*n_ptr,ch_len); */
        retval_p[len] = (char) NULL ; 
/*   printf("RetVal %s len %d\n", retval_p, len); */
/* removing trailing blank characters */
	while (ch_ptr[--*n_ptr] == ' ') {
		ch_ptr[*n_ptr] = (char) NULL;
	}
	printf("\nAvailable files matching: %s are\n",ch_ptr);
	n = 1;
	dirp = opendir(".");
	while ((direntp = readdir(dirp)) != NULL) {
		if ((fnmatch(ch_ptr,direntp->d_name,0)) == 0) {
			strcpy(matched[n],direntp->d_name);
			printf("%4d %s\n", n++, direntp->d_name);
			if (n > MAXFIL-1) {
				printf("WARNING pickfn: # of files exceeds the limit %d\n", MAXFIL);
				break;
			}
		}
	}
	closedir(dirp);
	printf("\nEnter a file name or pick a number (0 for none): ");
	scanf("%s",aline);
	if (isdigit(aline[0])) {
		i = atoi(aline);
		if (i > 0 && i < n) {
			strcpy(retval_p,matched[i]);
		} else {
			strcpy(retval_p,"0");
		}
	} else {
		strcpy(retval_p,aline);
	}
	*n_ptr = strlen(retval_p);
/*   if (fseek(stdin,off_stdin,0) < 0) {
     printf("Error in PICKFN.c fseek for stdin %d \n",off_stdin);
     exit(0);            
   } */

}

