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
#include "MyFLA_Utils.h"
#include "time_cubmat.h"


#define case1ref
//// #define case1a
#define case1b
#define case2ref
//// #define case2a
#define case2b
//// #define case2e
//// #define case2f
#define case3ref
//// #define case3a
#define case3b
#define case4ref
//// #define case4a
#define case4b


// ============================================================================
// Declaration of local prototypes.

int test_read_info( int * ndim, int ** avalues, int ** bvalues, int ** cvalues,
        int ** dvalues, int ** ivalues, int ** jvalues,
        int * nexecs, int * matrix_type, 
        int * print_data, int * check_result, int * seed );

void matrix_generate( int matrix_type, int seed, FLA_Obj A );


// ============================================================================
int main( int argc, char *argv[] ) {
  int      ndim, * avalues, * bvalues, * cvalues, * dvalues, * ivalues,
           * jvalues, idim, size_a, size_b, size_c, size_d, size_i, size_j,
           info, nexecs, matrix_type, print_data, check_result, seed[4], 
           variant, nit, size_cb_A, size_cb_B, size_cb_C, ref_executed;
  double   dtime, gflops, resid, max_resid;
  FLA_Obj  cb_A, cb_B, cb_C, cb_C_ref_result;

  // Initialize FLAME.
  FLA_Init( );

  // Checking the number of arguments.
  if ( argc != 1 ) {
    printf( "Usage:  %s\n\n", argv[0] );
    exit( -1 );
  }

  // Reading the input arguments.
  info = test_read_info( & ndim, & avalues, & bvalues, & cvalues, & dvalues,
             & ivalues, & jvalues, & nexecs, & matrix_type, & print_data, 
             & check_result, seed );

  if( info != 0 ) {
    printf( "Error in test_read_info\n" );
    exit( -1 );
  }

  printf( "%%                        a   b   c   d " );
  printf( "   i    j    time    gflops  resid \n" );
  printf( "%%--------------------------------------" );
  printf( "--------------------------------------\n" );

  max_resid = 0.0;
  nit = 1;
  for ( idim = 0; idim < ndim; idim++ ) {
    size_a = avalues[ idim ];
    size_b = bvalues[ idim ];
    size_c = cvalues[ idim ];
    size_d = dvalues[ idim ];
    size_i = ivalues[ idim ];
    size_j = jvalues[ idim ];


    // --------------------------------------------------------------------
    // Initializations for case 1.
    // --------------------------------------------------------------------

    ref_executed = 0;

    // Allocate space for cubes.
    size_cb_A = size_i * size_a * size_j * size_b;
    FLA_Obj_create( FLA_DOUBLE, size_cb_A, 1, 0, 0, & cb_A );
  
    size_cb_B = size_i * size_c * size_j * size_d;
    FLA_Obj_create( FLA_DOUBLE, size_cb_B, 1, 0, 0, & cb_B );
  
    size_cb_C = size_a * size_b * size_c * size_d;
    FLA_Obj_create( FLA_DOUBLE, size_cb_C, 1, 0, 0, & cb_C );
    FLA_Obj_create( FLA_DOUBLE, size_cb_C, 1, 0, 0, & cb_C_ref_result );
  
    // Generate data.
    matrix_generate( matrix_type, 1,             cb_A );
    matrix_generate( matrix_type, size_cb_A + 1, cb_B );
    MyFLA_Obj_set_to_one( cb_C );
    MyFLA_Obj_set_to_one( cb_C_ref_result );

#ifdef case1ref
    // ====================================================================
    // Timing case 1. Reference.
    // ====================================================================
    variant = 100;
    ref_executed = 1;
    time_cubmat( variant, nexecs, print_data, check_result, 
                 size_a, size_b, size_c, size_d, size_i, size_j,
                 cb_A, cb_B, cb_C, cb_C_ref_result,
                 & dtime, & gflops, & resid );
    if( resid > max_resid ) max_resid = resid;
    printf( "case1ref( %2d, 1:9 ) = [ %3d %3d %3d %3d %4d %4d ",
            nit, size_a, size_b, size_c, size_d, size_i, size_j );
    printf( " %7.2le %6.2lf %8.1le ];\n",
            dtime, gflops, resid );
    fflush( stdout );
#endif

#ifdef case1a
    // ====================================================================
    // Timing var1a.
    // ====================================================================
    variant = 101;
    if( ref_executed == 0 ) check_result = 0;
    time_cubmat( variant, nexecs, print_data, check_result, 
                 size_a, size_b, size_c, size_d, size_i, size_j,
                 cb_A, cb_B, cb_C, cb_C_ref_result,
                 & dtime, & gflops, & resid );
    if( resid > max_resid ) max_resid = resid;
    printf( "case1a( %2d, 1:9 )   = [ %3d %3d %3d %3d %4d %4d ",
            nit, size_a, size_b, size_c, size_d, size_i, size_j );
    printf( " %7.2le %6.2lf %8.1le ];\n", 
            dtime, gflops, resid );
    fflush( stdout );
#endif

#ifdef case1b
    // ====================================================================
    // Timing var1b.
    // ====================================================================
    variant = 102;
    if( ref_executed == 0 ) check_result = 0;
    time_cubmat( variant, nexecs, print_data, check_result, 
                 size_a, size_b, size_c, size_d, size_i, size_j,
                 cb_A, cb_B, cb_C, cb_C_ref_result,
                 & dtime, & gflops, & resid );
    if( resid > max_resid ) max_resid = resid;
    printf( "case1b( %2d, 1:9 )   = [ %3d %3d %3d %3d %4d %4d ",
            nit, size_a, size_b, size_c, size_d, size_i, size_j );
    printf( " %7.2le %6.2lf %8.1le ];\n", 
            dtime, gflops, resid );
    fflush( stdout );
#endif

    // --------------------------------------------------------------------
    // Final clearings for case 1.
    // --------------------------------------------------------------------
    printf( "\n" );

    // Remove space of cubes.
    FLA_Obj_free( & cb_A );
    FLA_Obj_free( & cb_B );
    FLA_Obj_free( & cb_C );
    FLA_Obj_free( & cb_C_ref_result );



    // --------------------------------------------------------------------
    // Initializations for case 2.
    // --------------------------------------------------------------------

    ref_executed = 0;

    // Allocate space for cubes.
    size_cb_A = size_i * size_a * size_j * size_b;
    FLA_Obj_create( FLA_DOUBLE, size_cb_A, 1, 0, 0, & cb_A );
  
    size_cb_B = size_j * size_c * size_i * size_d;
    FLA_Obj_create( FLA_DOUBLE, size_cb_B, 1, 0, 0, & cb_B );
  
    size_cb_C = size_a * size_b * size_c * size_d;
    FLA_Obj_create( FLA_DOUBLE, size_cb_C, 1, 0, 0, & cb_C );
    FLA_Obj_create( FLA_DOUBLE, size_cb_C, 1, 0, 0, & cb_C_ref_result );
  
    // Generate data.
    matrix_generate( matrix_type, 1,             cb_A );
    matrix_generate( matrix_type, size_cb_A + 1, cb_B );
    MyFLA_Obj_set_to_one( cb_C );
    MyFLA_Obj_set_to_one( cb_C_ref_result );

#ifdef case2ref
    // ====================================================================
    // Timing case 2. Reference.
    // ====================================================================
    variant = 200;
    ref_executed = 1;
    time_cubmat( variant, nexecs, print_data, check_result, 
                 size_a, size_b, size_c, size_d, size_i, size_j,
                 cb_A, cb_B, cb_C, cb_C_ref_result,
                 & dtime, & gflops, & resid );
    if( resid > max_resid ) max_resid = resid;
    printf( "case2ref( %2d, 1:9 ) = [ %3d %3d %3d %3d %4d %4d ",
            nit, size_a, size_b, size_c, size_d, size_i, size_j );
    printf( " %7.2le %6.2lf %8.1le ];\n",
            dtime, gflops, resid );
    fflush( stdout );
#endif

#ifdef case2a
    // ====================================================================
    // Timing var2a.
    // ====================================================================
    variant = 201;
    if( ref_executed == 0 ) check_result = 0;
    time_cubmat( variant, nexecs, print_data, check_result, 
                 size_a, size_b, size_c, size_d, size_i, size_j,
                 cb_A, cb_B, cb_C, cb_C_ref_result,
                 & dtime, & gflops, & resid );
    if( resid > max_resid ) max_resid = resid;
    printf( "case2a( %2d, 1:9 )   = [ %3d %3d %3d %3d %4d %4d ",
            nit, size_a, size_b, size_c, size_d, size_i, size_j );
    printf( " %7.2le %6.2lf %8.1le ];\n", 
            dtime, gflops, resid );
    fflush( stdout );
#endif

#ifdef case2b
    // ====================================================================
    // Timing var2b.
    // ====================================================================
    variant = 202;
    if( ref_executed == 0 ) check_result = 0;
    time_cubmat( variant, nexecs, print_data, check_result, 
                 size_a, size_b, size_c, size_d, size_i, size_j,
                 cb_A, cb_B, cb_C, cb_C_ref_result,
                 & dtime, & gflops, & resid );
    if( resid > max_resid ) max_resid = resid;
    printf( "case2b( %2d, 1:9 )   = [ %3d %3d %3d %3d %4d %4d ",
            nit, size_a, size_b, size_c, size_d, size_i, size_j );
    printf( " %7.2le %6.2lf %8.1le ];\n", 
            dtime, gflops, resid );
    fflush( stdout );
#endif

    // --------------------------------------------------------------------
    // Final clearings for case 2.
    // --------------------------------------------------------------------
    printf( "\n" );

    // Remove space of cubes.
    FLA_Obj_free( & cb_A );
    FLA_Obj_free( & cb_B );
    FLA_Obj_free( & cb_C );
    FLA_Obj_free( & cb_C_ref_result );




    // --------------------------------------------------------------------
    // Initializations for case 3.
    // --------------------------------------------------------------------

    ref_executed = 0;

    // Allocate space for cubes.
    size_cb_A = size_i * size_a * size_j * size_b;
    FLA_Obj_create( FLA_DOUBLE, size_cb_A, 1, 0, 0, & cb_A );
  
    size_cb_B = size_j * size_c * size_i * size_d;
    FLA_Obj_create( FLA_DOUBLE, size_cb_B, 1, 0, 0, & cb_B );
  
    size_cb_C = size_a * size_b * size_c * size_d;
    FLA_Obj_create( FLA_DOUBLE, size_cb_C, 1, 0, 0, & cb_C );
    FLA_Obj_create( FLA_DOUBLE, size_cb_C, 1, 0, 0, & cb_C_ref_result );
  
    // Generate data.
    matrix_generate( matrix_type, 1,             cb_A );
    matrix_generate( matrix_type, size_cb_A + 1, cb_B );
    MyFLA_Obj_set_to_one( cb_C );
    MyFLA_Obj_set_to_one( cb_C_ref_result );

#ifdef case3ref
    // ====================================================================
    // Timing case 3. Reference.
    // ====================================================================
    variant = 300;
    ref_executed = 1;
    time_cubmat( variant, nexecs, print_data, check_result, 
                 size_a, size_b, size_c, size_d, size_i, size_j,
                 cb_A, cb_B, cb_C, cb_C_ref_result,
                 & dtime, & gflops, & resid );
    if( resid > max_resid ) max_resid = resid;
    printf( "case3ref( %2d, 1:9 ) = [ %3d %3d %3d %3d %4d %4d ",
            nit, size_a, size_b, size_c, size_d, size_i, size_j );
    printf( " %7.2le %6.2lf %8.1le ];\n",
            dtime, gflops, resid );
    fflush( stdout );
#endif

#ifdef case3a
    // ====================================================================
    // Timing var3a.
    // ====================================================================
    variant = 301;
    if( ref_executed == 0 ) check_result = 0;
    time_cubmat( variant, nexecs, print_data, check_result, 
                 size_a, size_b, size_c, size_d, size_i, size_j,
                 cb_A, cb_B, cb_C, cb_C_ref_result,
                 & dtime, & gflops, & resid );
    if( resid > max_resid ) max_resid = resid;
    printf( "case3a( %2d, 1:9 )   = [ %3d %3d %3d %3d %4d %4d ",
            nit, size_a, size_b, size_c, size_d, size_i, size_j );
    printf( " %7.2le %6.2lf %8.1le ];\n", 
            dtime, gflops, resid );
    fflush( stdout );
#endif

#ifdef case3b
    // ====================================================================
    // Timing var3b.
    // ====================================================================
    variant = 302;
    if( ref_executed == 0 ) check_result = 0;
    time_cubmat( variant, nexecs, print_data, check_result, 
                 size_a, size_b, size_c, size_d, size_i, size_j,
                 cb_A, cb_B, cb_C, cb_C_ref_result,
                 & dtime, & gflops, & resid );
    if( resid > max_resid ) max_resid = resid;
    printf( "case3b( %2d, 1:9 )   = [ %3d %3d %3d %3d %4d %4d ",
            nit, size_a, size_b, size_c, size_d, size_i, size_j );
    printf( " %7.2le %6.2lf %8.1le ];\n", 
            dtime, gflops, resid );
    fflush( stdout );
#endif

    // --------------------------------------------------------------------
    // Final clearings for case 3.
    // --------------------------------------------------------------------
    printf( "\n" );

    // Remove space of cubes.
    FLA_Obj_free( & cb_A );
    FLA_Obj_free( & cb_B );
    FLA_Obj_free( & cb_C );
    FLA_Obj_free( & cb_C_ref_result );



    // --------------------------------------------------------------------
    // Initializations for case 4.
    // --------------------------------------------------------------------

    ref_executed = 0;

    // Allocate space for cubes.
    size_cb_A = size_i * size_a * size_j * size_b;
    FLA_Obj_create( FLA_DOUBLE, size_cb_A, 1, 0, 0, & cb_A );
  
    size_cb_B = size_j * size_c * size_i * size_d;
    FLA_Obj_create( FLA_DOUBLE, size_cb_B, 1, 0, 0, & cb_B );
  
    size_cb_C = size_a * size_b * size_c * size_d;
    FLA_Obj_create( FLA_DOUBLE, size_cb_C, 1, 0, 0, & cb_C );
    FLA_Obj_create( FLA_DOUBLE, size_cb_C, 1, 0, 0, & cb_C_ref_result );
  
    // Generate data.
    matrix_generate( matrix_type, 1,             cb_A );
    matrix_generate( matrix_type, size_cb_A + 1, cb_B );
    MyFLA_Obj_set_to_one( cb_C );
    MyFLA_Obj_set_to_one( cb_C_ref_result );

#ifdef case4ref
    // ====================================================================
    // Timing case 4. Reference.
    // ====================================================================
    variant = 400;
    ref_executed = 1;
    time_cubmat( variant, nexecs, print_data, check_result, 
                 size_a, size_b, size_c, size_d, size_i, size_j,
                 cb_A, cb_B, cb_C, cb_C_ref_result,
                 & dtime, & gflops, & resid );
    if( resid > max_resid ) max_resid = resid;
    printf( "case4ref( %2d, 1:9 ) = [ %3d %3d %3d %3d %4d %4d ",
            nit, size_a, size_b, size_c, size_d, size_i, size_j );
    printf( " %7.2le %6.2lf %8.1le ];\n",
            dtime, gflops, resid );
    fflush( stdout );
#endif

#ifdef case4a
    // ====================================================================
    // Timing var4a.
    // ====================================================================
    variant = 401;
    if( ref_executed == 0 ) check_result = 0;
    time_cubmat( variant, nexecs, print_data, check_result, 
                 size_a, size_b, size_c, size_d, size_i, size_j,
                 cb_A, cb_B, cb_C, cb_C_ref_result,
                 & dtime, & gflops, & resid );
    if( resid > max_resid ) max_resid = resid;
    printf( "case4a( %2d, 1:9 )   = [ %3d %3d %3d %3d %4d %4d ",
            nit, size_a, size_b, size_c, size_d, size_i, size_j );
    printf( " %7.2le %6.2lf %8.1le ];\n", 
            dtime, gflops, resid );
    fflush( stdout );
#endif

#ifdef case4b
    // ====================================================================
    // Timing var4b.
    // ====================================================================
    variant = 402;
    if( ref_executed == 0 ) check_result = 0;
    time_cubmat( variant, nexecs, print_data, check_result, 
                 size_a, size_b, size_c, size_d, size_i, size_j,
                 cb_A, cb_B, cb_C, cb_C_ref_result,
                 & dtime, & gflops, & resid );
    if( resid > max_resid ) max_resid = resid;
    printf( "case4b( %2d, 1:9 )   = [ %3d %3d %3d %3d %4d %4d ",
            nit, size_a, size_b, size_c, size_d, size_i, size_j );
    printf( " %7.2le %6.2lf %8.1le ];\n", 
            dtime, gflops, resid );
    fflush( stdout );
#endif

    // --------------------------------------------------------------------
    // Final clearings for case 4.
    // --------------------------------------------------------------------
    printf( "\n" );

    // Remove space of cubes.
    FLA_Obj_free( & cb_A );
    FLA_Obj_free( & cb_B );
    FLA_Obj_free( & cb_C );
    FLA_Obj_free( & cb_C_ref_result );





    printf( "%%--------------------------------------" );
    printf( "----------------------------------------\n" );
    nit++;
  }

  // Remove vectors: mvalues, nvalues, kvalues.
  free( avalues );
  free( bvalues );
  free( cvalues );
  free( dvalues );
  free( ivalues );
  free( jvalues );

  // Finalize FLAME.
  printf( "\n" );
  printf( "%% Maximum residual: %10.1le\n", max_resid );
  printf( "%%End of Program\n" );
  FLA_Finalize( );

  return 0;
}


// ============================================================================
int test_read_info( int * ndim, int ** avalues, int ** bvalues, int ** cvalues,
        int ** dvalues, int ** ivalues, int ** jvalues,
        int * nexecs, int * matrix_type, 
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
  *avalues = (int *) malloc( *ndim * sizeof( int ) );

  // Read the dimensions.
  for( i = 0; i < * ndim; i++ )
    fscanf( fp, "%d", (*avalues+i) );
  fgets( myLine, MAX_LINE_LENGTH, fp );

  // Write the dimensions.
  printf( "%% Test %d dimensions:\n", * ndim );
  printf( "%%   a: " );
  for( i = 0; i < * ndim; i++ )
    printf( "%4d ", *(*avalues+i) );
  printf( "\n" );

  // Create the vector to hold the dimensions "m".
  *bvalues = (int *) malloc( *ndim * sizeof( int ) );

  // Read the dimensions.
  for( i = 0; i < * ndim; i++ )
    fscanf( fp, "%d", (*bvalues+i) );
  fgets( myLine, MAX_LINE_LENGTH, fp );

  // Write the dimensions.
  printf( "%% Test %d dimensions:\n", * ndim );
  printf( "%%   b: " );
  for( i = 0; i < * ndim; i++ )
    printf( "%4d ", *(*bvalues+i) );
  printf( "\n" );

  // Create the vector to hold the dimensions "m".
  *cvalues = (int *) malloc( *ndim * sizeof( int ) );

  // Read the dimensions.
  for( i = 0; i < * ndim; i++ )
    fscanf( fp, "%d", (*cvalues+i) );
  fgets( myLine, MAX_LINE_LENGTH, fp );

  // Write the dimensions.
  printf( "%% Test %d dimensions:\n", * ndim );
  printf( "%%   c: " );
  for( i = 0; i < * ndim; i++ )
    printf( "%4d ", *(*cvalues+i) );
  printf( "\n" );

  // Create the vector to hold the dimensions "m".
  *dvalues = (int *) malloc( *ndim * sizeof( int ) );

  // Read the dimensions.
  for( i = 0; i < * ndim; i++ )
    fscanf( fp, "%d", (*dvalues+i) );
  fgets( myLine, MAX_LINE_LENGTH, fp );

  // Write the dimensions.
  printf( "%% Test %d dimensions:\n", * ndim );
  printf( "%%   d: " );
  for( i = 0; i < * ndim; i++ )
    printf( "%4d ", *(*dvalues+i) );
  printf( "\n" );

  // Create the vector to hold the dimensions "m".
  *ivalues = (int *) malloc( *ndim * sizeof( int ) );

  // Read the dimensions.
  for( i = 0; i < * ndim; i++ )
    fscanf( fp, "%d", (*ivalues+i) );
  fgets( myLine, MAX_LINE_LENGTH, fp );

  // Write the dimensions.
  printf( "%% Test %d dimensions:\n", * ndim );
  printf( "%%   i: " );
  for( i = 0; i < * ndim; i++ )
    printf( "%4d ", *(*ivalues+i) );
  printf( "\n" );

  // Create the vector to hold the dimensions "m".
  *jvalues = (int *) malloc( *ndim * sizeof( int ) );

  // Read the dimensions.
  for( i = 0; i < * ndim; i++ )
    fscanf( fp, "%d", (*jvalues+i) );
  fgets( myLine, MAX_LINE_LENGTH, fp );

  // Write the dimensions.
  printf( "%% Test %d dimensions:\n", * ndim );
  printf( "%%   j: " );
  for( i = 0; i < * ndim; i++ )
    printf( "%4d ", *(*jvalues+i) );
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
  printf( "%%\n" );

  fclose( fp );
  return 0;
}


// ============================================================================
void matrix_generate( int matrix_type, int intseed, FLA_Obj A ) {
  int     ldim_A, m_A, n_A, i, j;
  double  * buff_A, num, scale;

  buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
  ldim_A  = FLA_Obj_col_stride( A );
  m_A     = FLA_Obj_length( A );
  n_A     = FLA_Obj_width( A );

  if( matrix_type == 1 ) {
    //
    // matrix_type: 1
    // ---------------
    //
    // Initialize matrix.
    num = ( double ) intseed;
    for ( j = 0; j < n_A; j++ ) {
      for ( i = 0; i < m_A; i++ ) {
        buff_A[ i + j * ldim_A ] = num;
        num++;
      }
    }

  } else if( matrix_type == 2 ) {
    //
    // matrix_type: 2
    // ---------------
    //
    // Initialize matrix.
    num = ( double ) intseed;
    for ( j = 0; j < n_A; j++ ) {
      for ( i = 0; i < m_A; i++ ) {
        buff_A[ i + j * ldim_A ] = num;
        num++;
      }
    }
    // Scale down matrix.
    scale = 1.0 + ( double ) m_A * ( double ) n_A;
    if( scale != 0.0 ) {
      scale = 1.0 / scale;
      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_A[ i + j * ldim_A ] = buff_A[ i + j * ldim_A ] * scale;
        }
      }
    }

  } else {
    //
    // Unknown matrix_type
    // --------------------
    //
    fprintf( stderr, "ERROR in matrix_generate: Unknown matrix type: %d\n", 
             matrix_type );

  }
}


