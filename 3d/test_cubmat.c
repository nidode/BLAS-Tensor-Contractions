// ============================================================================
// CONTRIBUTORS
//
// Edoardo Di Napoli
//   Juelich Supercomputing Centre, Institute for Advanced Simulation
//   Forschungszentrum Juelich
//   Wilhelm-Johnen strasse, 52425 Juelich, Germany
//   e.di.napoli@fz-juelich.de
// Diego Fabregat, Paolo Bientinesi
//   AICES, RWTH-Aachen University, 52056-Aachen, Germany
//   {fabregat,pauldj}@aices.rwth-aachen.de
// G. Quintana-Orti
//   Depto. de Ingenieria y Ciencia de Computadores 
//   Universitat Jaume I, 12.080 Castellon, Spain
//   gquintan@icc.uji.es
//
// ============================================================================
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "FLAME.h"
#include "time_cubmat.h"


#define case1
#define case2
#define case3
#define case4
#define case5


// ============================================================================
// Declaration of local prototypes.

int test_read_info( int * ndim, int ** mvalues, int ** nvalues, int ** kvalues,
        int ** lvalues, int * nexecs, int * matrix_type, 
        int * print_data, int * check_result, int * seed );


// ============================================================================
int main( int argc, char *argv[] ) {
  FLA_Obj  C, C_copy;
  int      ndim, * mvalues, * nvalues, * kvalues, * lvalues, idim, m, n, k, l,
           info, nexecs, matrix_type, print_data, check_result, seed[4], 
           variant, nit;
  double   dtime, gflops, resid, max_resid;

  // Initialize FLAME.
  FLA_Init( );

  // Checking the number of arguments.
  if ( argc != 1 ) {
    printf( "Usage:  %s\n\n", argv[0] );
    exit( -1 );
  }

  // Reading the input arguments.
  info = test_read_info( & ndim, & mvalues, & nvalues, & kvalues, & lvalues,
                         & nexecs, & matrix_type, & print_data, & check_result,
                         seed );

  if( info != 0 ) {
    printf( "Error in test_read_info\n" );
    exit( -1 );
  }

  printf( "%%                        m    n    k   " );
  printf( " l    time      gflops  resid \n" );
  printf( "%%--------------------------------------" );
  printf( "--------------------------------------\n" );

  max_resid = 0.0;
  nit = 1;
  for ( idim = 0; idim < ndim; idim++ ) {
    m = mvalues[ idim ];
    n = nvalues[ idim ];
    k = kvalues[ idim ];
    l = lvalues[ idim ];

    // Create new objects.
    FLA_Obj_create( FLA_DOUBLE, m, n, 0, 0, & C );
    FLA_Obj_create( FLA_DOUBLE, m, n, 0, 0, & C_copy );

#ifdef case1
    // ====================================================================
    // Timing case 1.
    // ====================================================================
    variant = 1;
    time_cubmat( variant, nexecs, matrix_type, print_data, check_result, 
                 m, n, k, l, C, C_copy,
                 & dtime, & gflops, & resid );
    if( resid > max_resid ) max_resid = resid;
    printf( " case1( %3d, 1:6 ) = [ %4d %4d %4d %4d  ",
            nit, m, n, k, l );
    printf( " %7.2le  %6.2lf ];\n", 
            dtime, gflops, resid );
    fflush( stdout );
#endif

#ifdef case2
    // ====================================================================
    // Timing case 2.
    // ====================================================================
    variant = 2;
    time_cubmat( variant, nexecs, matrix_type, print_data, check_result, 
                 m, n, k, l, C, C_copy,
                 & dtime, & gflops, & resid );
    if( resid > max_resid ) max_resid = resid;
    printf( " case2( %3d, 1:7 ) = [ %4d %4d %4d %4d  ",
            nit, m, n, k, l );
    printf( " %7.2le  %6.2lf %8.1le ];\n", 
            dtime, gflops, resid );
    fflush( stdout );
#endif

#ifdef case3
    // ====================================================================
    // Timing case 3.
    // ====================================================================
    variant = 3;
    time_cubmat( variant, nexecs, matrix_type, print_data, check_result, 
                 m, n, k, l, C, C_copy,
                 & dtime, & gflops, & resid );
    if( resid > max_resid ) max_resid = resid;
    printf( " case3( %3d, 1:7 ) = [ %4d %4d %4d %4d  ",
            nit, m, n, k, l );
    printf( " %7.2le  %6.2lf %8.1le ];\n", 
            dtime, gflops, resid );
    fflush( stdout );
#endif

#ifdef case4
    // ====================================================================
    // Timing case 4.
    // ====================================================================
    variant = 4;
    time_cubmat( variant, nexecs, matrix_type, print_data, check_result, 
                 m, n, k, l, C, C_copy,
                 & dtime, & gflops, & resid );
    if( resid > max_resid ) max_resid = resid;
    printf( " case4( %3d, 1:7 ) = [ %4d %4d %4d %4d  ",
            nit, m, n, k, l );
    printf( " %7.2le  %6.2lf %8.1le ];\n", 
            dtime, gflops, resid );
    fflush( stdout );
#endif

#ifdef case5
    // ====================================================================
    // Timing case 5.
    // ====================================================================
    variant = 5;
    time_cubmat( variant, nexecs, matrix_type, print_data, check_result, 
                 m, n, k, l, C, C_copy,
                 & dtime, & gflops, & resid );
    if( resid > max_resid ) max_resid = resid;
    printf( " case5( %3d, 1:7 ) = [ %4d %4d %4d %4d  ",
            nit, m, n, k, l );
    printf( " %7.2le  %6.2lf %8.1le ];\n", 
            dtime, gflops, resid );
    fflush( stdout );
#endif

    // Remove objects.
    FLA_Obj_free( & C );
    FLA_Obj_free( & C_copy );

    printf( "%%--------------------------------------" );
    printf( "----------------------------------------\n" );

    nit++;
  }

  // Remove vectors: mvalues, nvalues, kvalues.
  free( mvalues );
  free( nvalues );
  free( kvalues );

  // Finalize FLAME.
  printf( "\n" );
  printf( "%% Maximum residual: %10.1le\n", max_resid );
  printf( "%%End of Program\n" );
  FLA_Finalize( );

  return 0;
}


// ============================================================================
int test_read_info( int * ndim, int ** mvalues, int ** nvalues, int ** kvalues,
        int ** lvalues, int * nexecs, int * matrix_type, 
        int * print_data,  int * check_result, int * seed ) {
  FILE       * fp;
  const int  MAX_LINE_LENGTH = 2048;
  char       myLine[ MAX_LINE_LENGTH ];
  int        i;

  if ( ( fp = fopen( "test.in", "r" ) ) == NULL )
    return -1;

  // Read the number of dimensions to test.
  fscanf( fp, "%d", ndim );
  fgets( myLine, MAX_LINE_LENGTH, fp );
  if( * ndim == 0 ) {
    fprintf( stderr, "+++ ERROR in test_read_info: " );
    fprintf( stderr, "Number of dimensions is zero.\n" );
    return (-1);
  }

  // Create the vector to hold the dimensions "m".
  *mvalues = (int *) malloc( *ndim * sizeof( int ) );

  // Read the dimensions.
  for( i = 0; i < * ndim; i++ )
    fscanf( fp, "%d", (*mvalues+i) );
  fgets( myLine, MAX_LINE_LENGTH, fp );

  // Write the dimensions.
  printf( "%% Test %d dimensions:\n", * ndim );
  printf( "%%   m: " );
  for( i = 0; i < * ndim; i++ )
    printf( "%4d ", *(*mvalues+i) );
  printf( "\n" );


  // Create the vector to hold the dimensions "n".
  *nvalues = (int *) malloc( *ndim * sizeof( int ) );

  // Read the dimensions.
  for( i = 0; i < * ndim; i++ )
    fscanf( fp, "%d", (*nvalues+i) );
  fgets( myLine, MAX_LINE_LENGTH, fp );

  // Write the dimensions.
  printf( "%% Test %d dimensions:\n", * ndim );
  printf( "%%   n: " );
  for( i = 0; i < * ndim; i++ )
    printf( "%4d ", *(*nvalues+i) );
  printf( "\n" );


  // Create the vector to hold the dimensions "k".
  *kvalues = (int *) malloc( *ndim * sizeof( int ) );

  // Read the dimensions.
  for( i = 0; i < * ndim; i++ )
    fscanf( fp, "%d", (*kvalues+i) );
  fgets( myLine, MAX_LINE_LENGTH, fp );

  // Write the dimensions.
  printf( "%% Test %d dimensions:\n", * ndim );
  printf( "%%   k: " );
  for( i = 0; i < * ndim; i++ )
    printf( "%4d ", *(*kvalues+i) );
  printf( "\n" );

  
  // Create the vector to hold the dimensions "l".
  *lvalues = (int *) malloc( *ndim * sizeof( int ) );

  // Read the dimensions.
  for( i = 0; i < * ndim; i++ )
    fscanf( fp, "%d", (*lvalues+i) );
  fgets( myLine, MAX_LINE_LENGTH, fp );

  // Write the dimensions.
  printf( "%% Test %d dimensions:\n", * ndim );
  printf( "%%   l: " );
  for( i = 0; i < * ndim; i++ )
    printf( "%4d ", *(*lvalues+i) );
  printf( "\n" );


  // Read the number of executions.
  fscanf (fp, "%d", nexecs );
  fgets( myLine, MAX_LINE_LENGTH, fp );
  printf( "%% Number of executions:             %d\n", *nexecs );

  // Read the matrix type.
  fscanf (fp, "%d", matrix_type );
  fgets( myLine, MAX_LINE_LENGTH, fp );
  printf( "%% Matrix type:                      %d\n", *matrix_type );

  // Read whether print data.
  fscanf( fp, "%d", print_data );
  fgets( myLine, MAX_LINE_LENGTH, fp );
  printf( "%% Print_data (0=no;1=yes):          %d\n", *print_data );

  // Read whether print data.
  fscanf( fp, "%d", check_result );
  fgets( myLine, MAX_LINE_LENGTH, fp );
  printf( "%% check_result (0=no;1=yes):        %d\n", *check_result );


  // Reading seeds. Last value must be odd.
  fscanf( fp, "%d", seed );
  fscanf( fp, "%d", seed+1 );
  fscanf( fp, "%d", seed+2 );
  fscanf( fp, "%d", seed+3 );
  fgets( myLine, MAX_LINE_LENGTH, fp );
  if( *( seed+3 ) % 2 == 0 ) {
    printf( "ERROR: Last value of seeds must be odd.\n" );
    return (-1);
  }
  printf( "%% Seeds:                            %d %d %d %d\n",
          *seed, *(seed+1), *(seed+2), *(seed+3) );

  fclose( fp );
  return 0;
}


