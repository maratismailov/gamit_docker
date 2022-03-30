/*
   it is also possible to get hostname etc if wish so.
   buf length depends on <sys/utsname.h> which varies
   from system to system, 2048 is safer.
 
   by peng fang, pfang@pgga.ucsd.edu  March 21, 1993; November 1, 1995.  
  
   Some MODS TAH 201026: Pretty dangerous routine that needs re-doing.
   Further MOD TAH 210112: Used the structure from utsname rather than
      the int buf array which latest versions of c-compilers don't
      support anymore (clang version 12.0.0) 

   Sun OS version with underscore after routine name.
   HP Version: No _ in name. 
*/
#include <stdio.h>
#include <string.h>
#include <sys/utsname.h>

void 
getmac (retval_p, len)
char *retval_p;
int len;
{
        char buf[10]; // Reduced size amount needed TAH 21012
        struct utsname unameData;
        int istat, i;
        /* istat = uname(buf); // Old form which is now outdated */
        istat = uname( &unameData);

        /* strcpy(buf,str1); */
        /* Use strncpy to make sure only length of buf is copied 
           TAH 210112 */
        strncpy(buf,unameData.sysname,10);
        /* Changed to 5 from 4 and to 6 from 5 below to handle Darwin */ 
        for (i=0; i <= 5; i++) {
             if( buf[1] == 0 ) {buf[i]=32;};
             *retval_p++ = buf[i];
        }
        for (i=6; i <= 9; i++) {
             /* set null to 32 for space in Fortran TAH 201026*/ 
             *retval_p++ = 32;
        }
}
