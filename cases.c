#include <stdio.h>
int main()
{
  int a,b,c,d ;
  printf("Enter a ") ;
  scanf("%d",&a) ;
  printf("Enter b ") ;  
  scanf("%d",&b) ;
  printf("Enter c ") ;
  scanf("%d",&c) ;
  printf("Enter d ") ;  
  scanf("%d",&d) ;
  printf("a= %d\n",a) ;
  printf("b= %d\n", b) ;
  printf("c= %d\n",c) ;
  printf("d= %d\n", d) ;



  switch( a)
{
 case 1 : printf("a=one\n") ;
 case 2 : printf("a=two\n") ;
}
  
  return 0;
}

