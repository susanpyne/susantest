#include <stdio.h>

int main()
{
 int x = 1;
 int y = 1;
 int *px = &x;
 int *py = &y;
 
 printf("Value of x = %d, address of x = %d\n", x, &x) ;
 x++ ;
 printf("px points to address %d which holds value %d\n", px, *px ) ;
 *px =10 ;
 printf("px points to address %d which holds value %d\n", px, *px ) ;

return 0;
}

