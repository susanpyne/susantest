#include <stdio.h>
void twice (int *ptr) ;
void thrice (int *ptr) ;
int main()
{

 int num = 5 ;
 int *ptr =&num ;/* value of num */
 printf("ptr stores address %p\n", ptr) ;
 printf("*ptr holds value %d\n", *ptr) ;
 printf("Value of num is %d\n", num) ;
 
 twice(ptr);
 printf("After twice num has value %d\n", num) ;
 thrice(ptr) ;
 printf("After thrice num has value %d\n", num) ;
 

return 0;
}
void twice(int *number)
{
  *number=(*number*2) ;
}
void thrice (int *number)
{
  *number=(*number*3) ;
}
