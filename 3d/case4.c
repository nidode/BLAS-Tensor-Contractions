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
#include "case4.h"


// ============================================================================
void compute_case4( int m, int n, int k, int l, 
         FLA_Obj cb_A, FLA_Obj cb_B, FLA_Obj C, int print_data ) {
  int      datatype, h, j, i_one = 1;
  double   * buff_cb_A, * buff_cb_B, d_one = 1.0;

  // Some initializations.
  datatype  = FLA_Obj_datatype( cb_A );
  buff_cb_A = ( double * ) FLA_Obj_buffer_at_view( cb_A );
  buff_cb_B = ( double * ) FLA_Obj_buffer_at_view( cb_B );

  // Initialize matrix C for the result.
  MyFLA_Obj_set_to_zero( C );

  // Show data.
  if( print_data == 1 ) {
    FLA_Obj_show( " Ci = [ ", C, "%le", " ];" );
    FLA_Obj_show( " cb_A = [ ", cb_A, "%le", " ];" );
    FLA_Obj_show( " cb_B = [ ", cb_B, "%le", " ];" );
  }

  // Perform computation.
  for( j = 0; j < n; j++ ) {
    for( h = 0; h < k; h++ ) {
      dgemv_( "Transpose", & l, & m,
              & d_one, ( double * ) buff_cb_A + l * m * h, & l,
                       ( double * ) buff_cb_B + k * l * j + h, & k,
              & d_one, ( double * ) FLA_Obj_buffer_at_view( C ) + m * j, 
                       & i_one );
    }
  }
  
  // Show data.
  if( print_data == 1 ) {
    FLA_Obj_show( " Cf = [ ", C, "%le", " ];" );
  }
}


