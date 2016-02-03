#include <stdio.h>

typedef struct /* structure consisting of char array */
{
  char str[5] ;
} Arrtype ;

typedef struct /* structure consisting of char pointer */
{
  char *str ; /* * means this variable is a pointer*/
} Ptrtype ;

Arrtype arr = {'B', 'a', 'd', ' ', '\0'} ; /* intialise arr which has type Arrtype*/

Ptrtype ptr = {"Good"} ; /* initialise pointer which has type Ptrtype*/

int main()
{
  printf("\nArray string is a %s\n", arr.str ) ;
   
 return 0 ;
}




