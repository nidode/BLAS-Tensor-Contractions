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
#include "FLAME.h"
#include "MyFLA_Utils.h"
#include "case1.h"
#include "case2.h"
#include "case3.h"
#include "case4.h"
#include "time_cubmat.h"


// ============================================================================
// Declaration of local prototypes.

double compute_resid( FLA_Obj A_orig, FLA_Obj A_new );


// =============================================================================
void time_cubmat( int variant, int nexecs, int print_data, int check_result,
         int size_a, int size_b, int size_c, int size_d, int size_i, int size_j,
         FLA_Obj cb_A, FLA_Obj cb_B, FLA_Obj cb_C, FLA_Obj cb_C_ref_result,
         double * dtime, double * gflops, double * resid ) {

  int      irep;
  double   dtime_1, dtime_2;


  // Loop for performing "nexecs" executions.
  for ( irep = 0 ; irep < nexecs; irep++ ) {
    dtime_1 = FLA_Clock();

    switch( variant ){

    case 100:
      compute_case1ref( size_a, size_b, size_c, size_d, size_i, size_j,
                        cb_A, cb_B, cb_C, print_data );
      break;

    case 101:
      compute_case1a( size_a, size_b, size_c, size_d, size_i, size_j,
                      cb_A, cb_B, cb_C, print_data );
      break;

    case 102:
      compute_case1b( size_a, size_b, size_c, size_d, size_i, size_j,
                      cb_A, cb_B, cb_C, print_data );
      break;

    case 200:
      compute_case2ref( size_a, size_b, size_c, size_d, size_i, size_j,
                        cb_A, cb_B, cb_C, print_data );
      break;

    case 201:
      compute_case2a( size_a, size_b, size_c, size_d, size_i, size_j,
                      cb_A, cb_B, cb_C, print_data );
      break;

    case 202:
      compute_case2b( size_a, size_b, size_c, size_d, size_i, size_j,
                      cb_A, cb_B, cb_C, print_data );
      break;

    case 300:
      compute_case3ref( size_a, size_b, size_c, size_d, size_i, size_j,
                        cb_A, cb_B, cb_C, print_data );
      break;

    case 301:
      compute_case3a( size_a, size_b, size_c, size_d, size_i, size_j,
                      cb_A, cb_B, cb_C, print_data );
      break;

    case 302:
      compute_case3b( size_a, size_b, size_c, size_d, size_i, size_j,
                      cb_A, cb_B, cb_C, print_data );
      break;

    case 400:
      compute_case4ref( size_a, size_b, size_c, size_d, size_i, size_j,
                        cb_A, cb_B, cb_C, print_data );
      break;

    case 401:
      compute_case4a( size_a, size_b, size_c, size_d, size_i, size_j,
                      cb_A, cb_B, cb_C, print_data );
      break;

    case 402:
      compute_case4b( size_a, size_b, size_c, size_d, size_i, size_j,
                      cb_A, cb_B, cb_C, print_data );
      break;

    default:
      fprintf( stderr, "ERROR in time_cubmat: Variant not implemented.\n" );
    }

    dtime_2 = FLA_Clock();
    if ( irep == 0 )
      * dtime = ( dtime_2 > dtime_1 ? dtime_2 - dtime_1 : 0.0 );
    else
      * dtime += ( dtime_2 > dtime_1 ? dtime_2 - dtime_1 : 0.0 );
  }

  // Check residuals.
  if( check_result == 0 ) {
    * resid = -1.0;
  } else {
    // Compute residuals.
    switch( variant ){
    case 100:
    case 200:
    case 300:
    case 400:
      FLA_Copy( cb_C, cb_C_ref_result );
      * resid = -1.0;
      break;
    case 101:
    case 102:
    case 103:
    case 104:
    case 201:
    case 202:
    case 203:
    case 204:
    case 205:
    case 206:
    case 301:
    case 302:
    case 401:
    case 402:
      * resid = compute_resid( cb_C, cb_C_ref_result );
      break;
    default:
      printf( "ERROR: variant not implemented.\n" );
    }
  }

  // Compute timings.
  * dtime = * dtime / nexecs;

  // Compute gigaflops.
  * gflops = ( 2.0 * ( ( double ) size_a ) * ( ( double ) size_b ) * 
                     ( ( double ) size_c ) * ( ( double ) size_d ) * 
                     ( ( double ) size_i ) * ( ( double ) size_j ) ) /
             ( * dtime * 1.0e+9 );
}


// ============================================================================
double compute_resid( FLA_Obj A_orig, FLA_Obj A_new ) {
// Returns residual of: norm( A_orig - A_new, 1 ).
  FLA_Obj Diff, alpha;
  double  nrm_diff, nrm_a, resid;

  // Create auxiliary objects.
  FLA_Obj_create( FLA_Obj_datatype( A_orig ), 1, 1, 0, 0, & alpha );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A_orig, & Diff );

  // Diff := A_orig - A_new.
  FLA_Copy( A_orig, Diff );
  FLA_Axpy( FLA_MINUS_ONE, A_new, Diff );

  // nrm_diff := norm( Diff, 1 ).
  MyFLA_Nrm1( Diff, alpha );
  nrm_diff = *( ( double * ) FLA_Obj_buffer_at_view( alpha ) );

  // nrm_a := norm( A_orig, 1 ).
  MyFLA_Nrm1( A_orig, alpha );
  nrm_a = *( ( double * ) FLA_Obj_buffer_at_view( alpha ) );

  // Free auxiliary objects.
  FLA_Obj_free( & Diff );
  FLA_Obj_free( & alpha );

  // Return relative residual.
  if( nrm_a != 0.0 ) {
    resid = nrm_diff / nrm_a;
  } else {
    resid = nrm_diff;
  }
  return resid;
}


