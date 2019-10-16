/* Program to compute the sizes of various data type and determine needed
 * FEAP configuration parameters
*/

#include <stdio.h>

int main(int argc, char **argv)
{
  char c;
  int  i;
  float f;
  double d;
  long l;
  long long ll;
 
  printf("\n **** FEAP Memory Prober: Start ****\n");
  printf("\nYour machine's data-types use the following number of bits\n");
  printf("\n char %ld, int %ld, long %ld, long long %ld, float %ld, \
double %ld, pointers %ld\n",
  8*sizeof(c),
  8*sizeof(i),
  8*sizeof(l),
  8*sizeof(ll),
  8*sizeof(f),
  8*sizeof(d),
  8*sizeof(&c));

  i = 8*sizeof(&c);

  if(i  == 32 ) printf("\n Select the integer4 include files in \
the top level makefile.in\n");
  else if( i == 64 ) printf("\n Select the integer8 include files \
in the top level makefile.in\n");
  else printf("\n Your pointer size (%d bit pointers) in not \
supported by FEAP\n",i);

  i = sizeof(d)/sizeof(i);

  if ( i == 1 || i == 2) printf("\n Set ipr to %d in main/feapXX.f\n",i);
  else printf("\n The real to integer ratio of %d is not supported \
by FEAP\n",i);
  printf("\n **** FEAP Memory Prober: End ****\n");
}

