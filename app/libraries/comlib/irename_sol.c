/* 
   irename is a wrapper written in C to force calls of 'rename' to
   use the C-library routine rather than the Fortran-library routine,
   which doesn't exist in Fortran 90

   The Sun/DEC/LINUX version differs from the HP/IBM version only in the
   use of the underscore following the subroutine name, to designate a C-routine.

   Written by P. Fang; installed by R. King  980720

*/

#include <stdio.h>
#include <string.h>
void
irename_ (old, new, istat)
char *old, *new;
int *istat;
{
	strtok(old," ");
	strtok(new," ");
	*istat = rename(old, new);
}

