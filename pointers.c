#include <stdio.h>
int main()
{

 int num = 8 ;
 int *ptr=&num ;
 printf("Variable contains %d\n", num) ;
 printf("Pointer contains address 0x%p\n", ptr) ;
 printf("Pointer points to value %d\n", *ptr) ;
 *ptr =12 ;
 printf("Variable contains %d\n", num) ;
 printf("Pointer contains address 0x%p\n", ptr) ;
 printf("Pointer points to value %d\n", *ptr) ;

return 0;
}
