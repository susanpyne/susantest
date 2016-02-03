#include <stdio.h>
int main()
{

 int i ;
 int nums[10] ={10,20,30,40,50,60,70,80,90,100} ;
 int *ptr = nums ; /* sets pointer to 1st element of array*/
 
 printf("At address %p is value %d\n", ptr, *ptr) ;

 ptr++ ;
 printf("At address %p is value %d\n", ptr, *ptr) ;
 ptr++ ;
 printf("At address %p is value %d\n", ptr, *ptr) ;
 ptr-=2 ;
 
 for (i=0; i<10; i++)
 {printf("Element %d contains value %d\n", i, *ptr) ;
  ptr++ ;
 }
  
 

return 0;
}
