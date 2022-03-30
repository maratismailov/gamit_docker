/* collection of useful c-commands for FONDA */

/* change a character to upper case  */
#include <ctype.h>
short upperc(c)
short c;
{
   return(islower(c)? c+'A'-'a' : c);
}

/* change a character to lower case  */
#include <ctype.h>
short lowerc(c)
short c;
{
   return(isupper(c)? c+'a'-'A' : c);
}

/* sort a string and output argument number */
#include <strings.h>
int argsrt(s)
char s[];
{ 
   int i,j,ia;
   j = 0; ia = 0;
   for (i = 0; s[i] != '\n'; i ++) {
      if (strcmp(s[i]," ") == 0 || strcmp(s[i],"\t") == 0) 
         j = 0;
      else
         j += 1;
      if (j == 1) ia += 1;
   }
   return(ia);
}      

/* get a specified argument from a string */
#include <stdio.h>
#include <strings.h>
int argext(s,id,arg)
char s[], arg[];
int id;
{
   int anum, i, j, ia;
   anum = argsrt(s);
   if (anum < id) errexit(" No enough arguments. ", NULL);
   j = 0; ia = 0;
   for (i = 0; s[i] != '\n'; i ++) {
      if (strcmp(s[i]," ") == 0 || strcmp(s[i],"\t") == 0) 
         j = 0;
      else
         j += 1;
      if (j == 1) ia += 1;
      if (ia == id && j > 0) strcpy(arg[j-1],s[i]);
      if (ia > id) break;
   }
   return(0);
}

/* exit when error happens */
#include <stdio.h> 
errexit(s1,s2)
char *s1,*s2;
{
   printf(s2==NULL?"%s\n":"%s%s\n",s1,s2);
   exit(-1);
}
   
