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
#include "FLAME.h"
#include "extract_flat_matrices.h"


// ============================================================================
void extract_flat_matrix_from_cubic_with_plane_nk(
         int m, int n, int k, FLA_Obj cubmat, int idx, FLA_Obj A ) {
  double  * buff_cubmat, * buff_A;
  int     ldim_A, i, j;

  buff_cubmat = ( double * ) FLA_Obj_buffer_at_view( cubmat );
  buff_A      = ( double * ) FLA_Obj_buffer_at_view( A );
  ldim_A      = FLA_Obj_col_stride( A );

  for ( j = 0 ; j < k; j++ ) {
    for ( i = 0 ; i < n; i++ ) {
      buff_A[ i + j * ldim_A ] = 
        buff_cubmat[ idx + i * m + j * m * n ];
    }
  }
}

