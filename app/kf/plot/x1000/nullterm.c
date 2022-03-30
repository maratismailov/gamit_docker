
int nullterm( name )

/* Function to put a null at the first blank or \0 in string and
   return the length */

char *name ; /* Name of the character string to be null terminated */
{

int i       ; /* Loop counter for stepping through string */

for ( i=0 ; i <=256 ; ++i )
    { 
      if ( *(name+i) == '\0' ) 
         { return(i);
         }
      if ( *(name+i)==' '  ) 
       { *(name+i) = '\0' ; 
         return(i);
       }
    }
return(i) ;
}
