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



  if ( a>b )
  {
   if (c>d)
   {printf("a>b and c>d\n") ; }
   else
   {printf("a>b and c<=d\n") ; }
  }
  else
  {
   if (c>d)
   {printf("a<=b and c>d\n") ; }
   else 
   {printf("a<=b and c<=d\n") ; }
  }
  
  return 0;
}

