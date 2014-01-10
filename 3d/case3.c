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
#include "extract_flat_matrices.h"
#include "case3.h"


// ============================================================================
void compute_case3( int m, int n, int k, int l, 
         FLA_Obj cb_A, FLA_Obj cb_B, FLA_Obj C, int print_data ) {
  FLA_Obj  slice_A, slice_B;
  int      datatype, h;
  double   * buff_cb_A, * buff_cb_B;

  // Some initializations.
  datatype  = FLA_Obj_datatype( cb_A );
  buff_cb_A = ( double * ) FLA_Obj_buffer_at_view( cb_A );
  buff_cb_B = ( double * ) FLA_Obj_buffer_at_view( cb_B );

  // Prepare temporal slices.
  FLA_Obj_create_without_buffer( datatype, l, m, & slice_A );
  FLA_Obj_create( datatype, l, n, 0, 0, & slice_B );

  // Initialize matrix C for the result.
  MyFLA_Obj_set_to_zero( C );

  // Show data.
  if( print_data == 1 ) {
    FLA_Obj_show( " Ci = [ ", C, "%le", " ];" );
    FLA_Obj_show( " cb_A = [ ", cb_A, "%le", " ];" );
    FLA_Obj_show( " cb_B = [ ", cb_B, "%le", " ];" );
  }

  // Perform computation.
  for( h = 0; h < k; h++ ) {
    FLA_Obj_attach_buffer( buff_cb_A + l * m * h, 1, l, & slice_A );
    extract_flat_matrix_from_cubic_with_plane_nk( k, l, n, cb_B, h, slice_B );
    FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_ONE, slice_A, slice_B, FLA_ONE, C );
  }
  
  // Remove temporal slices.
  FLA_Obj_free_without_buffer( & slice_A );
  FLA_Obj_free( & slice_B );

  // Show data.
  if( print_data == 1 ) {
    FLA_Obj_show( " Cf = [ ", C, "%le", " ];" );
  }
}

