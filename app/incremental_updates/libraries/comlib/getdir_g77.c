/*

   getdir is C-routine to replace the slow and system-dependent getdir.f in 
   gamit/lib.  The Sun/DEC/LINUX version differs from the HP/IBM version only in the
   use of the underscore following the subroutine name, to designate a C-routine. 

   get returns a set of files matching wildcard in current directory
   note: MAXLEN must match the string length in Fortran
         MAXFIL must be greater or equal to argument maxfil
   (see test fortran code below)

   has been tested on Solaris 5.5.1, HP-UX 10.20, Linux 2.0.34

   by pfang@ucsd.edu July 13, 1998 

	program testget
c program to test the C subroutine getdir

	integer maxfil, nret
	character*80 wildcard,returned(500)
	data maxfil/500/

	print *, 'Enter wildcard :'
	read(*,'(a)') wildcard

	call getdir(wildcard,maxfil,returned,nret)

	do i = 1, nret
		print *,returned(i)
	enddo
	end

*/

#include <dirent.h> 
#include <sys/types.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <fnmatch.h>
#define MAXLEN 80
#define MAXFIL 500

void getdir_ (wild, maxfil, ret, nret)
char wild[MAXLEN], ret[MAXFIL][MAXLEN];
int *maxfil, *nret;
{
	DIR *dirp;
	struct dirent *direntp;
	char pth[MAXLEN];
	int i,j,len,off_stdin;

/* It appears that on LINUX use of ftell/fseek on the stdin file causes problems. Simon McClusky 980916
   off_stdin = ftell(stdin);  
   printf("Unit 5 offset: %d \n",off_stdin); */

	strtok(wild," ");

	len = strlen(wild);
	j = 0;
	for(i = 0; i < len; i++) {
		if (wild[i] == '/')  j = i;
	}
	if (j == 0) {
		strcpy(pth,".");
	} else {
		j++;
		strncpy(pth,wild,j);
		pth[j] = (char) NULL;
		for(i = 0; i <= len-j; i++) {
			wild[i] = wild[j+i];
		}
	}

	i = 0;
	dirp = opendir(pth); 
	while ((direntp = readdir(dirp)) != NULL) {
		if ((fnmatch(wild,direntp->d_name,0)) == 0) {
			strcpy(ret[i++],direntp->d_name);
			if (i == MAXFIL || i == *maxfil) {
				printf("WARNING getdir: # of files exceeds the limit %d with wildcard %s\n", MAXFIL, wild);
				break;
			}
		}
	}
	closedir(dirp);
	*nret = i;

/* It appears that on LINUX use of ftell/fseek on the stdin file causes problems. Simon McClusky 980916
   if (fseek(stdin,off_stdin,0) < 0) {
     printf("Error in GETDIR.c fseek for stdin %d \n",off_stdin);
     exit(0);            
   } */

}
