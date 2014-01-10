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
#include "case5.h"
#include "time_cubmat.h"


// ============================================================================
// Declaration of local prototypes.

double compute_resid( FLA_Obj A_orig, FLA_Obj A_new );

void matrix_generate( int matrix_type, int seed, FLA_Obj A );


// ============================================================================
void time_cubmat( int variant, int nexecs, int matrix_type, 
         int print_data, int check_result,
         int m, int n, int k, int l, FLA_Obj C, FLA_Obj C_copy,
         double * dtime, double * gflops, double * resid ) {
  FLA_Obj  cb_A, cb_B;
  int      irep, datatype;
  double   dtime_1, dtime_2, d_m, d_n, d_k, d_l;

  // Some initializations.
  datatype = FLA_Obj_datatype( C );

  // Allocate space for cubes and matrices.
  FLA_Obj_create( datatype, m * k * l, 1,  0, 0,  & cb_A );
  FLA_Obj_create( datatype, n * k * l, 1,  0, 0,  & cb_B );

  // Generate data.
  matrix_generate( matrix_type,             1, cb_A );
  matrix_generate( matrix_type, m * k * l + 1, cb_B );
  MyFLA_Obj_set_to_zero( C );

  // Loop for performing "nexecs" executions.
  for ( irep = 0 ; irep < nexecs; irep++ ) {
    dtime_1 = FLA_Clock();

    switch( variant ){

    case 1:
      // Time case 1.
      compute_case1( m, n, k, l, cb_A, cb_B, C, print_data );
      break;

    case 2:
      // Time case 2.
      compute_case2( m, n, k, l, cb_A, cb_B, C, print_data );
      break;

    case 3:
      // Time case 3.
      compute_case3( m, n, k, l, cb_A, cb_B, C, print_data );
      break;

    case 4:
      // Time case 4.
      compute_case4( m, n, k, l, cb_A, cb_B, C, print_data );
      break;

    case 5:
      // Time case 5.
      compute_case5( m, n, k, l, cb_A, cb_B, C, print_data );
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
    // Compute residuals for dense cholesky factorizations.
    // Residual checking only works if all variants are computed.
    switch( variant ){
    case 1:
      FLA_Copy( C, C_copy );
      * resid = -1.0;
      break;
    case 2:
      * resid = compute_resid( C, C_copy );
      break;
    case 3:
      FLA_Copy( C, C_copy );
      * resid = -1.0;
      break;
    case 4:
      * resid = compute_resid( C, C_copy );
      break;
    case 5:
      * resid = compute_resid( C, C_copy );
      break;
    case 6:
      * resid = compute_resid( C, C_copy );
      break;
    default:
      * resid = - 1.0;
      printf( "ERROR: variant not implemented.\n" );
    }
  }

  // Remove objects.
  FLA_Obj_free( & cb_A );
  FLA_Obj_free( & cb_B );

  // Compute timings.
  * dtime = * dtime / nexecs;

  // Compute gigaflops.
  d_m = ( double ) m;
  d_n = ( double ) n;
  d_k = ( double ) k;
  d_l = ( double ) l;
  switch( variant ){
  case 1:
    * gflops = ( 2.0 * d_m * d_n * d_k * d_l )/( * dtime * 1.0e+9 );
    break;
  case 2:
    * gflops = ( 2.0 * d_m * d_n * d_k * d_l )/( * dtime * 1.0e+9 );
    break;
  case 3:
    * gflops = ( 2.0 * d_m * d_n * d_k * d_l )/( * dtime * 1.0e+9 );
    break;
  case 4:
    * gflops = ( 2.0 * d_m * d_n * d_k * d_l )/( * dtime * 1.0e+9 );
    break;
  case 5:
    * gflops = ( 2.0 * d_m * d_n * d_k * d_l )/( * dtime * 1.0e+9 );
    break;
  case 6:
    * gflops = ( 2.0 * d_m * d_n * d_k * d_l )/( * dtime * 1.0e+9 );
    break;
  default:
    * gflops = - 1.0;
    printf( "ERROR: variant not implemented.\n" );
  }

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

