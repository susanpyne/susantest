#include <stdio.h>
void thing() ;/* void because use arguments to return values */

int main()
{
 int i=0 ;
 int j=0 ;
 printf ("\nStarting main: i = %d, j = %d", i,j) ;

 thing(i,&j) ; /* first argument is a value, second is a pointer/*

 printf("\nEnd of main: i= %d, j = %d\n", i,j) ;/* i has not been incremented but j has becasue argument as ref to address*/
 
return 0;
}

void thing(int x, int *y) 
{
 printf("\nStarting thing: x = %d, y = %d", x,*y) ;
 x ++ ; /* increments local variable within function */
 (*y) ++ ; /* increments value held in address pointed to */
 printf("\nEnd of thing: x = %d, y = %d", x, *y) ; 

}
