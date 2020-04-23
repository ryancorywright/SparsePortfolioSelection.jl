#include <stdio.h>
#include <stdlib.h>
/*--------------------------------------------------------------------------*/
// COMPILE with -lm option (math library)
/*--------------------------------------------------------------------------*/
// Parameter file must contain:
// n
// average lower bound
// average upper bound
// maximum bound variation
// min expected return
// max expected return
// seed

/*--------------------------------------------------------------------------*/
// parameters
int n ;
double avgminbd,avgmaxbd,varbd ;
double minret, maxret ;
unsigned int seed ;
// generated values
double *minbd , *maxbd ;
double *retvalue ;
double rho ;
/*--------------------------------------------------------------------------*/

void read_param( char *probname );
void gen_data( );
void write_data( char *probname );

/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
  if ( argc < 2 )
    {
      fprintf(stderr,"Error: specify test name\n");
      fprintf(stderr,"Usage: genfin <test name>\n");
      exit(10);
    }
  read_param( argv[1] ); 
  gen_data( );
  write_data( argv[1] );
}
/*--------------------------------------------------------------------------*/
void read_param( char *probname )
{
  char filename[100];
  sprintf(filename,"%s.param",probname);
  FILE *file = fopen(filename,"r");
  fscanf(file,"%d",&n);
  fscanf(file,"%lf",&avgminbd);
  fscanf(file,"%lf",&avgmaxbd);
  fscanf(file,"%lf",&varbd);
  fscanf(file,"%lf",&minret);
  fscanf(file,"%lf",&maxret);
  fscanf(file,"%d",&seed);
  fclose(file);
}
/*--------------------------------------------------------------------------*/
void  gen_data( )
{
  int i = 0 ;
  srand( seed );
  minbd = ( double *) calloc( n , sizeof(double));
  maxbd = ( double *) calloc( n , sizeof(double));
  for (; i < n ; i ++ )
    {
      int x = rand();
      double tmp = avgminbd - varbd/2 + varbd * (((double) x)/RAND_MAX);
      if ( tmp < 0 ) tmp = 0 ;
      if ( tmp > 1 ) tmp = 1 ;
      minbd[ i ] = tmp ;
    }
  for ( i = 0 ; i < n ; i ++ )
    {
      int x = rand();
      double tmp = avgmaxbd - varbd/2 + varbd * (((double) x)/RAND_MAX);
      if ( tmp < 0 ) tmp = 0 ;
      if ( tmp > 1 ) tmp = 1 ;
      maxbd[ i ] = tmp ;
      if ( maxbd[i] < minbd[i] ) 
	{
	  maxbd[i] = minbd[i];
	  minbd[i] = tmp ;
	}
    }

  retvalue = ( double *) calloc( n , sizeof(double));
  for ( i = 0 ; i < n ; i ++ )
    {
      int x = rand();
      retvalue[ i ] = minret + (maxret-minret) * (((double) x)/RAND_MAX);
    }
  rho = minret + (maxret-minret) * (((double) rand())/RAND_MAX); 
  fprintf(stderr,"Proposed seed %d\n",rand());
}
/*--------------------------------------------------------------------------*/
void write_data( char *probname )
{
  char filename[100];
  sprintf(filename,"%s.bds",probname);
  FILE *file = fopen(filename,"w");
  int i = 0 ;
  for (; i < n ; i ++ )
    fprintf(file,"%10.8f %10.8f\n", minbd[i], maxbd[i]);
  fclose(file);
  
  sprintf(filename,"%s.txt",probname);
  file = fopen(filename,"w");
  fprintf( file , "%d\n", n );
  for ( i = 0; i < n ; i ++ )
    fprintf(file,"%10.8f 0.0\n", retvalue[i]);  
  fclose(file);

  sprintf(filename,"%s.rho",probname);
  file = fopen(filename,"w");
  fprintf( file , "%10.8f\n", rho );
  fclose(file);
}
/*--------------------------------------------------------------------------*/


