#include <stdio.h>
void first(char str[]) ;
int square(int x );
int cube(int y) ;

int main() 
{
 int num, a ;
 char msg[10]="Hello" ;
 printf("Enter a ") ;
 scanf("%d", &a) ;
 printf("a= %d\n",a) ;
 first(msg) ;
 num = square(a) ;
 printf("square = %d\n", num) ;
 printf("cube = %d\n", cube(a)) ;
   
 return 0;
}
void  first(char str[])
{
 printf("%s\n", str) ;
}
int square( int x) 
{
 
 return (x*x) ;
}

int cube(int y)
{ 
  
  return (y*y*y);
}



