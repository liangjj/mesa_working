/*deck @(#)second.c	1.1 5/18/95 */
#include <time.h>
double second()
{

  clock_t junk;
  junk=clock();
  return ( (((double)junk)/((double) CLOCKS_PER_SEC)) );
}

