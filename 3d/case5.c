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
#include "case5.h"


// ============================================================================
void compute_case5( int m, int n, int k, int l, 
         FLA_Obj cb_A, FLA_Obj cb_B, FLA_Obj C, int print_data ) {
  int      datatype, i, j, kl;
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
  kl = k * l;
  for( i = 0; i < k; i++ ) {
    for( j = 0; j < l; j++ ) {
      dger_( & m, & n,
             & d_one, ( double * ) buff_cb_A + l * m * i + j,
                      & l,
                      ( double * ) buff_cb_B + i + k * j,
                      & kl,
             ( double * ) FLA_Obj_buffer_at_view( C ), 
             & m );
    }
  }
  
  // Show data.
  if( print_data == 1 ) {
    FLA_Obj_show( " Cf = [ ", C, "%le", " ];" );
  }
}

