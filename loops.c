#include <stdio.h>
int main()
{
  int  i, j , arr[3]= {10,20,30};
  
  for ( i=1 ; i <4 ;i++)
  {
    for ( j=1; j<3 ; j++)
     {
       printf("i= %d, j= %d\n", i,j) ;
     }
  }
  i=0 ;
  while (i<3)
   { 
    printf("arr[%d] = %d\n", i, arr[i]) ;
    i++ ;
   }
   
  
  return 0;
}

