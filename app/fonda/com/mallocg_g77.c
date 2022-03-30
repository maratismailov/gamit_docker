/* 
   mallocg.c

     Created by Simon McClusky

     Tested w/ gcc-2.7.2.1f/SunOS 4.1.3 
            w/ gcc-2.7.2.1f/Linux 2.0.18

     Update on 97/01/14

*/
/*-------------------------------------------------------------malloc*/
int mallocg_(int *size) 
{ return (int)malloc(*size); }
