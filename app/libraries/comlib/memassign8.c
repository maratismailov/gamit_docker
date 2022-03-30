/* 
   memassign8.c

     Created by Thomas Herring

     Two arguments are passed:
     mem_i8  -- Number of I*4 works to be assigned
     vma_start -- First element in the array that the memory is to
                be assigned to.  The element number in this routine
                where the memory is available is returned.

     Update on 04/17/19

*/
/*-------------------------------------------------------------malloc*/
#include <stdlib.h>
#include <stdio.h>
#ifdef ADD64BIT
#define F_PTR_SIZE long
#else
#define F_PTR_SIZE long long
#endif

F_PTR_SIZE memassign8_(long *numwords, int *sizeword, long *vma_data) {
/* Get the meomory address as 64-but integer for large memory
  machines */
int *memad ; /* Address where memory is available */
long size  ; /* Number of bytes to be assigned */
F_PTR_SIZE vma_start ; /* Location in vma_data where memory is
                   available */

/* Compute the memory needed in bytes and get the address
printf("Values passed %d %d %d\n",*mem_i8, *vma_data, vma_data); */

/* size = *mem_i8*4; */
size = (long) *numwords  * (long)(*sizeword) * 4 ;
memad = (int *) malloc(size);
/* printf("Trying to allocate memory %lu bytes of memory %d %d \n",size, *numwords, *sizeword); */
if( memad == 0 ) {
   printf("Unable to allocate memory %lu bytes\n",size);
   vma_start = 0;
   return vma_start;
   }
/* Compute where we need to start in array */
vma_start = ((int *) memad - (int *) (*vma_data)) + 1; 
/* printf("Adress %lu Offset %lu \n",(long)memad, vma_start); */
 
return vma_start; }

/* Version for HP 
   (remove _ from name of function */

/* 
   memassign8.c

     Created by Thomas Herring

     Two arguments are passed:
     mem_i8  -- Number of I*4 works to be assigned
     vma_start -- First element in the array that the memory is to
                be assigned to.  The element number in this routine
                where the memory is available is returned.

     Update on 04/17/19

*/
/*-------------------------------------------------------------malloc*/
#include <stdlib.h>
#ifdef ADD64BIT
#define F_PTR_SIZE long
#else
#define F_PTR_SIZE long long
#endif

F_PTR_SIZE memassign8(long *numwords, int *sizeword, long *vma_data) {
/* Get the meomory address as 64-but integer for large memory
  machines */
int *memad ; /* Address where memory is available */
long size  ; /* Number of bytes to be assigned */
F_PTR_SIZE vma_start ; /* Location in vma_data where memory is
                   available */

/* Compute the memory needed in bytes and get the address
printf("Values passed %d %d %d\n",*mem_i8, *vma_data, vma_data); */

/* size = *mem_i8*4; */
size = (long) *numwords  * (long)(*sizeword) * 4 ;
memad = (int *) malloc(size);
/* printf("Trying to allocate memory %lu bytes of memory\n",size); */
if( memad == 0 ) {
   printf("Unable to allocate memory %lu bytes\n",size);
   vma_start = 0;
   return vma_start;
   }
/* Compute where we need to start in array */
vma_start = ((int *) memad - (int *) (*vma_data)) + 1; 
/* printf("Adress %lu Offset %lu \n",(long)memad, vma_start); */
 
return  vma_start; }
